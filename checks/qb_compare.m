% QB_COMPARE  Are the Naor .nd2 movies (5Qb\250ms, 10Qb) the source of the
% stored Qbeta data, and does Qbeta behave the same after re-extraction?
% Compares the stored MS data (in_vivo_experiments-old) with the re-extraction
% (in_vivo_reextracted/after_R2_offset: red channel only, multithresh spot,
% offset 200; same pipeline as PP7). Reviewer points: general data validity,
% the Qbeta (lambda 4) vs PP7 (lambda 1) claim.
% Per dataset and version:
%  - traces: number, >= 180 frames, median level, length, alternation flag
%    (as checks/check_interleave.py: alt > 0.85 or |log(odd/even mean)| > 0.2);
%  - bursts at SP 10..15 (paper SP 14): counts per trace, amplitudes (csv);
%  - fits at SP 14: grid/WS (paper weights) both signs; likelihood fit (MLE,
%    100-resample granule bootstrap) negative bursts; sigma from ampVar.
% Writes analysis_outputs/qb_compare/traces.csv, bursts.csv, fits.csv, amplitudes.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'qb_compare');
if ~exist(outDir, 'dir'), mkdir(outDir); end
V = {'stored_MS', 'in_vivo_experiments-old'; 'R2_reextracted', 'in_vivo_reextracted/after_R2_offset'};
names = {'5Qb', '10Qb'}; tVec = (0:360) * 10/60; win = 13; trunc = 15;
tr = {}; bu = {}; fi = {}; am = cell(0, 6);
for v = 1:size(V,1)
    for d = 1:numel(names)
        [S, B] = gm_loadDataset(repoDir, V{v,2}, names{d});
        L = cellfun(@numel, S);
        flag = 0; n60 = 0;
        for i = 1:numel(S)
            x = double(S{i}(:));
            if numel(x) < 60, continue, end
            n60 = n60 + 1;
            s = sign(diff(x)); s = s(s ~= 0);
            alt = mean(s(2:end) ~= s(1:end-1)); r = mean(x(1:2:end)) / mean(x(2:2:end));
            flag = flag + (alt > 0.85 || abs(log(r)) > 0.2);
        end
        tr(end+1,:) = {V{v,1}, names{d}, numel(S), sum(L >= 180), median(L), ...
            median(cellfun(@(x) median(double(x)), S)), median(cellfun(@(x) median(double(x)), B)), n60, flag}; %#ok<SAGROW>
        for SP = 10:15
            T = gm_amplitudes_by_trace(S, B, tVec, win, SP);
            T = T([T.len] >= 180);
            nN = cellfun(@numel, {T.ampNeg}); nP = cellfun(@numel, {T.ampPos});
            bu(end+1,:) = {V{v,1}, names{d}, SP, numel(T), sum(nN), sum(nP), mean(nN + nP), ...
                median(abs(vertcat(T.ampNeg))), median(abs(vertcat(T.ampPos))), gm_sigma_from_ampvar(vertcat(T.ampVar))}; %#ok<SAGROW>
            for i = 1:numel(T)
                for sg = {'Neg', 'Pos'}
                    a = abs(T(i).(['amp' sg{1}]));
                    for k = 1:numel(a), am(end+1,:) = {V{v,1}, names{d}, SP, sg{1}, i, a(k)}; end %#ok<SAGROW>
                end
            end
            if SP == 14
                for sg = {'Neg', 'Pos'}
                    o = gm_fit_dataset(T, struct('sign', sg{1}, 'truncPercent', trunc, ...
                        'nBoot', 100 * strcmp(sg{1}, 'Neg'), 'seed', d, 'grid', true));
                    m = o.mle;
                    fi(end+1,:) = {V{v,1}, names{d}, sg{1}, o.nTraces, o.nBursts, m.n, o.sigma, o.grid.lambda, o.grid.k0, ...
                        m.lambda, o.boot.lamCI(1), o.boot.lamCI(2), m.k0, o.boot.k0CI(1), o.boot.k0CI(2), m.lamK0}; %#ok<SAGROW>
                end
            end
        end
        fprintf('%s %s done\n', V{v,1}, names{d});
    end
end
writetable(cell2table(tr, 'VariableNames', {'data','dataset','nTraces','nTraces180','medianLength', ...
    'medianSpotLevel','medianBackgroundLevel','nTraces60','nAlternationFlagged'}), fullfile(outDir, 'traces.csv'));
writetable(cell2table(bu, 'VariableNames', {'data','dataset','SP','nTraces180','nNeg','nPos','burstsPerTrace', ...
    'medianAbsNeg','medianAbsPos','sigma'}), fullfile(outDir, 'bursts.csv'));
writetable(cell2table(fi, 'VariableNames', {'data','dataset','sign','nTraces180','nBursts','nFit','sigma', ...
    'gridLam','gridK0','mleLam','mleLamLoBoot','mleLamHiBoot','mleK0','mleK0LoBoot','mleK0HiBoot','mleLamK0'}), fullfile(outDir, 'fits.csv'));
writetable(cell2table(am, 'VariableNames', {'data','dataset','SP','sign','trace','amp'}), fullfile(outDir, 'amplitudes.csv'));
disp('QB_COMPARE_DONE');
