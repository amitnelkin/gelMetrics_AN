function crop_sum_diagnosis(workRoot)
% CROP_SUM_DIAGNOSIS  Why does the crop-sum "fraction shed" not show the
% Qbeta vs PP7 difference that the mean and spot-sum measures show
% (RESULTS.md section 15)? Reviewer point: measure robustness.
% 1. Radial profiles: mean intensity in 1-px rings 0..20 px around tracked
%    granules (one movie per dataset, tracks >= 180 frames, every 20th frame),
%    minus the median of the 18-20 px ring.
% 2. Variant integrated intensities (collect_intensity opts.sumDiag): box sums
%    at half-widths 7/10/14/20 px (border-median background), a granule-only
%    sum (r <= 6 px minus the 8-11 px ring median), and the cell light just
%    around the granule (4-8 px ring minus the 8-11 px ring median).
% 3. Fraction shed per negative burst with each variant (analysis.m, SP 12,
%    signal = sum + background as in shed_fraction.m).
% 4. Where does the shed light go? For each negative burst on the spot-sum
%    trace: change of the spot sum vs change of the 4-8 px cell ring, over the
%    same frames.
% Writes analysis_outputs/crop_diagnosis/radial.csv, fractions.csv, nearcell.csv.
checksDir = fileparts(mfilename('fullpath')); repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir); addpath(fullfile(repoDir, 'analysis'));
outDir = fullfile(repoDir, 'analysis_outputs', 'crop_diagnosis');
if ~exist(outDir, 'dir'), mkdir(outDir); end
names = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'};
% 1. radial profiles
rad = {};
for d = 1:numel(names)
    w = fullfile(workRoot, names{d}, 'work'); tf = dir(fullfile(w, '*.tif'));
    [~, base] = fileparts(tf(1).name);
    M = generate_movArray([w filesep], tf(1).name); xy = load(fullfile(w, [base '_xy.mat']));
    prof = [];
    for t = 1:numel(xy.x)
        if numel(xy.x{t}) < 180, continue, end
        for n = 1:20:numel(xy.x{t})
            f = xy.frame{t}(n) + 1; x = xy.x{t}(n) + 1; y = xy.y{t}(n) + 1;
            I = double(M(f).data);
            r = max(round(y)-21,1):min(round(y)+21,size(I,1)); c = max(round(x)-21,1):min(round(x)+21,size(I,2));
            W = I(r,c); [cc, rr] = meshgrid(c, r); dd = hypot(cc - x, rr - y);
            bg = median(W(dd >= 18 & dd <= 20));
            p = arrayfun(@(k) mean(W(dd >= k & dd < k+1)) - bg, 0:20);
            prof = [prof; p]; %#ok<AGROW>
        end
    end
    mp = median(prof, 1, 'omitnan');
    for k = 0:20, rad(end+1,:) = {names{d}, k, mp(k+1), size(prof,1)}; end %#ok<AGROW>
    fprintf('%s radial: %d samples, half-max radius %.1f px\n', names{d}, size(prof,1), find(mp < mp(1)/2, 1) - 1);
end
writetable(cell2table(rad, 'VariableNames', {'dataset','radius','medianIntensity','nSamples'}), fullfile(outDir, 'radial.csv'));
% 2-4. variant sums, fractions, near-cell response
tVec = (0:360) * 10/60; fr = {}; nc = {};
for d = 1:numel(names)
    [I, B, ~, S] = collect_intensity_tests(fullfile(workRoot, names{d}, 'work'), struct('offset', 200, 'sumDiag', true));
    keep = [false, cellfun(@numel, I(2:end)) >= 180];
    Bk = B(keep)';
    V = {'spotsum', S.spotSum; 'box7', S.cropSumW{1}; 'box10', S.cropSumW{2}; 'box14', S.cropSumW{3}; ...
         'box20', S.cropSumW{4}; 'ringsum', S.ringSum};
    for v = 1:size(V,1)
        Sv = V{v,2}(keep)'; Sv = cellfun(@(s, b) s(:) + b(:), Sv, Bk, 'UniformOutput', false);
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(Sv, Bk, tVec, 13, 12, 0, 0);
        warning(ws);
        seg = seg(seg(:,2) < 0 & seg(:,6) > 0, :);
        for q = 1:size(seg,1)
            fr(end+1,:) = {names{d}, V{v,1}, seg(q,1), abs(seg(q,5)) / seg(q,6), abs(seg(q,5)), seg(q,6)}; %#ok<AGROW>
        end
        if strcmp(V{v,1}, 'spotsum')
            % near-cell response over the same frames (raw traces, 13-frame smoothing)
            Sk = S.spotSum(keep); Nk = S.nearCell(keep);
            for q = 1:size(seg,1)
                i = seg(q,1); a = max(1, round(seg(q,3))); b = min(numel(Sk{i}), round(seg(q,4)));
                s1 = smoothdata(Sk{i}, 'movmean', 13); n1 = smoothdata(Nk{i}, 'movmean', 13);
                nc(end+1,:) = {names{d}, i, s1(b) - s1(a), n1(b) - n1(a), n1(a)}; %#ok<AGROW>
            end
        end
    end
    fprintf('%s variants done\n', names{d});
end
writetable(cell2table(fr, 'VariableNames', {'dataset','measure','trace','fraction','absAmp','levelBefore'}), fullfile(outDir, 'fractions.csv'));
writetable(cell2table(nc, 'VariableNames', {'dataset','trace','dSpotSum','dNearCellMean','nearCellBefore'}), fullfile(outDir, 'nearcell.csv'));
disp('CROPDIAG_DONE');
end
