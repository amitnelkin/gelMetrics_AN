function noise_from_background(datasets, outSub)
% NOISE_FROM_BACKGROUND  Independent check of the per-trace noise level used for
% the burst-free surrogates (RESULTS.md section 21; reviewer point R2 6).
%
% The surrogates use the noise measured on each trace itself (1.4826*MAD of the
% frame-to-frame differences / sqrt(2)). That could include fast real granule
% changes. Here the expected camera + photon noise is predicted WITHOUT using the
% granule trace, from pixels away from all granules:
%  1. Calibration per movie: for every pixel farther than 10 px from any
%     tracked point, its mean (offset 200 subtracted) and its frame-to-frame
%     variance var(diff)/2 (removes slow drift and bleaching). Pixels are binned
%     by mean (30 quantile bins), and the bin medians are fitted with a line:
%     variance = a * intensity + c (EMCCD: photon noise grows linearly with
%     intensity, c = read/dark noise).
%  2. For each kept track (same crop and multithresh masks as
%     collect_intensity.m), the predicted noise variance per frame of
%       spot sum  = sum_spot v(I) + N^2 * sum_cell v(I) / M^2
%       spot mean = sum_spot v(I) / N^2
%       cell mean = sum_cell v(I) / M^2
%     (N spot, M cell pixels; v from the calibration line). Predicted sigma =
%     sqrt(median over frames).
%  3. Measured sigma of the same three series: 1.4826*MAD(diff)/sqrt(2), using
%     only consecutive frames that were both detected.
% Ratio measured / predicted ~ 1: the trace's frame-to-frame jitter is fully
% explained by camera + photon noise. > 1: something extra (mask jitter from the
% per-frame threshold, or fast real changes).
% Movies: ..\zigzag_scratch\reextract\<dataset>\work (fluorescence channel only).
% Writes analysis_outputs/noise_background/{calibration.csv, traces.csv}.
% Optional (section 39): outSub = output folder under analysis_outputs, e.g.
% noise_from_background({'vitro_1000ng_10ug', ...}, 'noise_background_vitro').
if nargin < 1, datasets = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}; end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir);
if nargin < 2, outSub = 'noise_background'; end
outDir = fullfile(repoDir, 'analysis_outputs', outSub);
if ~exist(outDir, 'dir'), mkdir(outDir); end
scratch = fullfile(fileparts(repoDir), 'zigzag_scratch', 'reextract');
off = 200; hw = 7; cal = cell(0, 8); tr = cell(0, 14);
for d = 1:numel(datasets)
    ds = datasets{d}; wd = fullfile(scratch, ds, 'work');
    tifs = dir(fullfile(wd, '*.tif')); [~, o] = sort({tifs.name}); tifs = tifs(o);
    k = 0;
    for f = 1:numel(tifs)
        [~, base] = fileparts(tifs(f).name);
        mov = generate_movArray([wd filesep], tifs(f).name);
        st = single(cat(3, mov.data)); clear mov
        xy = load(fullfile(wd, [base '_xy.mat']));
        % calibration mask: far from every tracked point
        near = false(size(st, 1), size(st, 2));
        for t = 1:numel(xy.x)
            r = round(xy.y{t} + 1); c = round(xy.x{t} + 1);
            ok = r >= 1 & r <= size(st,1) & c >= 1 & c <= size(st,2);
            near(sub2ind(size(near), r(ok), c(ok))) = true;
        end
        near = imdilate(near, strel('disk', 10));
        m = mean(st, 3) - off; v = var(diff(st, 1, 3), 0, 3) / 2;
        mm = m(~near); vv = v(~near);
        edges = unique(quantile(mm, linspace(0, 1, 31)));
        [~, ~, bin] = histcounts(mm, edges);
        bm = accumarray(bin(bin > 0), mm(bin > 0), [], @median);
        bv = accumarray(bin(bin > 0), vv(bin > 0), [], @median);
        pf = polyfit(bm, bv, 1); a = pf(1); c0 = pf(2);
        R2 = 1 - sum((bv - polyval(pf, bm)).^2) / sum((bv - mean(bv)).^2);
        vfun = @(I) max(a * (double(I) - off) + c0, 1); p99 = prctile(mm, 99);
        cal(end+1,:) = {ds, base, a, c0, R2, prctile(mm, 1), prctile(mm, 99), nnz(~near)}; %#ok<AGROW>
        for t = 1:numel(xy.x)
            fr = xy.frame{t} + 1; xs = xy.x{t} + 1; ys = xy.y{t} + 1;
            if numel(xs) < 24, continue, end
            k = k + 1;
            idx = fr - fr(1) + 1; n = idx(end);
            S = nan(n,1); Mn = S; Cm = S; vS = S; vM = S; vC = S; Np = S; above = S;
            for q = 1:numel(xs)
                fra = st(:, :, fr(q));
                b = round([max(ys(q)-hw,1), min(ys(q)+hw,size(fra,1)), max(xs(q)-hw,1), min(xs(q)+hw,size(fra,2))]);
                sub = fra(b(1):b(2), b(3):b(4));
                th = multithresh(sub, 4);
                sp = sub >= th(3); ce = sub >= th(1) & sub < th(3);
                Is = double(sub(sp)); Ic = double(sub(ce)); N = numel(Is); M = numel(Ic);
                cm = mean(Ic);
                S(idx(q)) = sum(Is - cm); Mn(idx(q)) = mean(Is) - off; Cm(idx(q)) = cm - off;
                sv = sum(vfun(Is)); cv = sum(vfun(Ic));
                vS(idx(q)) = sv + N^2 * cv / M^2; vM(idx(q)) = sv / N^2; vC(idx(q)) = cv / M^2;
                Np(idx(q)) = N; above(idx(q)) = mean(Is - off > p99);
            end
            tr(end+1,:) = {ds, base, k, n, median(S, 'omitnan'), median(Np, 'omitnan'), ...
                rsd(S), sqrt(median(vS, 'omitnan')), rsd(Mn), sqrt(median(vM, 'omitnan')), ...
                rsd(Cm), sqrt(median(vC, 'omitnan')), median(Mn, 'omitnan'), mean(above, 'omitnan')}; %#ok<AGROW>
        end
        clear st
        fprintf('%s %s done (a %.1f, c %.0f, R2 %.3f)\n', ds, base, a, c0, R2);
    end
end
writetable(cell2table(cal, 'VariableNames', {'dataset','movie','gain_a','readvar_c','R2','meanP1','meanP99','nPixels'}), ...
    fullfile(outDir, 'calibration.csv'));
writetable(cell2table(tr, 'VariableNames', {'dataset','movie','trace','lengthFrames','spotSum','nSpotPx', ...
    'sdSum_meas','sdSum_pred','sdMean_meas','sdMean_pred','sdCell_meas','sdCell_pred','spotMean','fracSpotPxAboveCal'}), ...
    fullfile(outDir, 'traces.csv'));
disp('NOISE_BG_DONE');
end

function s = rsd(x)
% robust frame-to-frame noise, consecutive detected frames only
d = diff(x); d = d(~isnan(d));
s = 1.4826 * median(abs(d - median(d))) / sqrt(2);
end
