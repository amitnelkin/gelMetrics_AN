function static_control(workRoot, datasets, outName, versions, measure)
% STATIC_CONTROL  Pixel-level negative control for the extraction (RESULTS.md
% section 23): does the extraction itself (per-frame or steady multithresh
% thresholds) create apparent granule dynamics?
% For every real track (>= 24 points): a static granule = median of its first
% 31 real 15x15 crops (offset removed), carried along the track (the granule
% stays centred, so tracking and motion play no role), bleached with the movie's
% own background trend (cubic fit to the mean of pixels > 10 px from any track),
% plus fresh Gaussian camera + photon noise per pixel and frame with variance
% a*(I - 200) + c (calibrated per movie as in noise_from_background.m). Nothing
% about the granule changes. Each synthetic crop sequence is measured with the
% collect_intensity.m lines (multithresh(4), spot = top two classes, cell =
% middle two, spotSum = sum(spot - cell mean)), thresholds per frame (thrWin 0,
% published) or their running median over 31 frames (steady), and the traces
% get exactly the real-data treatment: burst detection (analysis.m,
% window 13 and 17, SP 12) on the traces and on burst-free surrogates (shared
% bleaching shape + white noise, as detection_grid.m), and the structure
% function D(1 frame), D(30 frames) of y = spotSum / cubicTrend - 1 (as
% noise_methods_real.m). Seeded. Tracks whose crop touches the image edge are skipped.
% Writes analysis_outputs/steady/static_control.csv (or outName).
% Optional (section 39): measure = 'sum' (default, spot sum + cell mean as in vivo) or
% 'mean' (spot mean, offset removed, as the in vitro traces); D1/D30 then use the spot mean.
if nargin < 2, datasets = {'4PP7', '4PP7-FUS'}; end
if nargin < 3, outName = 'static_control.csv'; end
if nargin < 4, versions = [0 31]; end
if nargin < 5, measure = 'sum'; end   % thresholds: 0 per frame, >0 running median, -1 one fixed mask (section 33)
vname = @(w) [repmat('vfix', 1, double(w < 0)), repmat(sprintf('v%d', w), 1, double(w >= 0))];
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'checks'));
off = 200; rng(2023);
tVec = (0:360) * 10/60; rows = cell(0, 15);
for d = 1:numel(datasets)
    ds = datasets{d}; wd = fullfile(workRoot, ds, 'work');
    tifs = dir(fullfile(wd, '*.tif')); [~, o] = sort({tifs.name}); tifs = tifs(o);
    X = struct(); for w = versions, X.(vname(w)) = {}; end
    Bk = X; mv = [];
    for f = 1:numel(tifs)
        [~, base] = fileparts(tifs(f).name);
        mov = generate_movArray([wd filesep], tifs(f).name);
        st = single(cat(3, mov.data)); clear mov
        xy = load(fullfile(wd, [base '_xy.mat']));
        near = false(size(st, 1), size(st, 2));
        for t = 1:numel(xy.x)
            r = round(xy.y{t} + 1); c = round(xy.x{t} + 1);
            ok = r >= 1 & r <= size(st,1) & c >= 1 & c <= size(st,2);
            near(sub2ind(size(near), r(ok), c(ok))) = true;
        end
        near = imdilate(near, strel('disk', 10));
        m = mean(st, 3) - off; v = var(diff(st, 1, 3), 0, 3) / 2; mm = m(~near); vv = v(~near);
        edges = unique(quantile(mm, linspace(0, 1, 31))); [~, ~, bin] = histcounts(mm, edges);
        pf = polyfit(accumarray(bin(bin > 0), mm(bin > 0), [], @median), accumarray(bin(bin > 0), vv(bin > 0), [], @median), 1);
        nf = size(st, 3); fm = zeros(nf, 1);
        for k = 1:nf, fr = st(:, :, k); fm(k) = mean(fr(~near)) - off; end
        tt = (1:nf)'; bl = polyval(polyfit(tt, fm, 3), tt); bl = bl / bl(1);
        % one static granule per track, carried along the track: template = median
        % of its first 31 crops; every frame = template bleached + fresh noise
        for t = 1:numel(xy.x)
            fr0 = xy.frame{t} + 1; xs = xy.x{t} + 1; ys = xy.y{t} + 1;
            if numel(xs) < 24, continue, end
            idx = fr0 - fr0(1) + 1; n = idx(end); C = zeros(15, 15, numel(xs), 'single'); okc = true;
            for q = 1:numel(xs)
                bb = round([ys(q)-7, ys(q)+7, xs(q)-7, xs(q)+7]);
                if bb(1) < 1 || bb(3) < 1 || bb(2) > size(st,1) || bb(4) > size(st,2), okc = false; break, end
                C(:, :, q) = st(bb(1):bb(2), bb(3):bb(4), fr0(q));
            end
            nT = min(31, numel(xs)); T = median(C(:, :, 1:nT), 3) - off; bl0 = mean(bl(fr0(1:nT)));
            for w = versions
                fn = vname(w);
                if ~okc, X.(fn){end+1} = []; Bk.(fn){end+1} = []; continue, end
                Y = zeros(15, 15, numel(xs), 'single');
                for q = 1:numel(xs)
                    I = T * bl(fr0(q)) / bl0;
                    Y(:, :, q) = max(off + I + randn(15, 15, 'single') .* sqrt(max(pf(1) * I + pf(2), 1)), 0);
                end
                TH = nan(n, 4);
                if w < 0   % one fixed mask from the granule's average crop
                    TH(idx, :) = repmat(double(multithresh(mean(Y, 3), 4)), numel(xs), 1);
                else
                    for q = 1:numel(xs), TH(idx(q), :) = double(multithresh(Y(:, :, q), 4)); end
                end
                if w > 0, TH = movmedian(TH, w, 1, 'omitnan'); end
                ss = nan(n, 1); cm = nan(n, 1);
                for q = 1:numel(xs)
                    sub = double(Y(:, :, q)); th = TH(idx(q), :);
                    sp = sub >= th(3); ce = sub >= th(1) & sub < th(3);
                    cm(idx(q)) = mean(sub(ce)); ss(idx(q)) = sum(sub(sp) - cm(idx(q))); cm(idx(q)) = cm(idx(q)) - off;
                    if strcmp(measure, 'mean'), ss(idx(q)) = mean(sub(sp)) - off; end
                end
                X.(fn){end+1} = fillmissing(ss, 'linear'); Bk.(fn){end+1} = fillmissing(cm, 'linear');
            end
            mv(end+1) = f; %#ok<AGROW>
        end
        clear st
        fprintf('%s %s done\n', ds, base);
    end
    for fn = arrayfun(vname, versions, 'UniformOutput', false)
        Ssum = X.(fn{1})(:); Bc = Bk.(fn{1})(:);
        keep = cellfun(@numel, Ssum) >= 180 & ~cellfun(@(c) any(isnan(c)), Ssum); Ssum = Ssum(keep); Bc = Bc(keep); mk = mv(keep);
        S = cellfun(@(s, b) s(:) + b(:), Ssum, Bc, 'UniformOutput', false);
        if strcmp(measure, 'mean'), S = cellfun(@(s) s(:), Ssum, 'UniformOutput', false); end
        Lm = min(cellfun(@numel, S)); shp = zeros(numel(S), Lm);
        for i = 1:numel(S)
            s = S{i}(:); x = (1:numel(s))'; trd = polyval(polyfit(x, s, 3), x); shp(i,:) = trd(1:Lm)' / mean(trd(1:Lm));
        end
        shape = median(shp, 1); Ssur = S;
        for i = 1:numel(S)
            s = S{i}(:); n = numel(s); sh = [shape, repmat(shape(end), 1, max(n - Lm, 0))]'; sh = sh(1:n);
            sgm = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
            Ssur{i} = mean(s) / mean(sh) * sh + sgm * randn(n, 1);
        end
        res = zeros(numel(S), 8);
        for wq = 1:2
            win = [13 17]; win = win(wq);
            for kind = 1:2
                Z = S; if kind == 2, Z = Ssur; end
                ws = warning('off', 'all'); [~, ~, ~, ~, ~, ~, seg] = analysis(Z, Bc, tVec, win, 12, 0, 0); warning(ws);
                seg = seg(seg(:,2) < 0 & seg(:,6) > 0, :);
                res(:, (wq - 1) * 4 + (kind - 1) * 2 + 1) = accumarray(seg(:,1), 1, [numel(S) 1]);
                res(:, (wq - 1) * 4 + (kind - 1) * 2 + 2) = accumarray(seg(:,1), abs(seg(:,5)) ./ seg(:,6), [numel(S) 1]);
            end
        end
        for i = 1:numel(S)
            g = Ssum{i}(:)'; t = 1:numel(g); trd = polyval(polyfit(t, g, 3), t); D1 = NaN; D30 = NaN;
            if all(trd > 0)
                y = g ./ trd - 1; D1 = mean((y(2:end) - y(1:end-1)).^2); D30 = mean((y(31:end) - y(1:end-30)).^2);
            end
            rows(end+1,:) = [{ds, fn{1}, i, mk(i), numel(g)}, num2cell(res(i,:)), {D1, D30}]; %#ok<AGROW>
        end
    end
end
R = cell2table(rows, 'VariableNames', {'dataset','version','trace','movie','lengthFrames', ...
    'nReal13','fracReal13','nSurr13','fracSurr13','nReal17','fracReal17','nSurr17','fracSurr17','D1','D30'});
writetable(R, fullfile(repoDir, 'analysis_outputs', 'steady', outName));
disp('STATIC_DONE');
end
