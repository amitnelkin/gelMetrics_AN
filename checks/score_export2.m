function score_export2(tag, opts)
% SCORE_EXPORT2  Burst tables for the robustness analyses of the chosen scores
% (fraction shed per burst; absolute burst size; turnover), checks/robustness.py.
% Reviewer points R2 6 (noise floor, one SP), R2 4a (sensitivity), R2 5
% (nested bootstrap, time halves), R2 2b (smoothing window).
% Same detection as the runners (analysis.m). Datasets: in vivo R2 re-extraction
% (mean and spot-sum measures; movie id per trace from extraction_info) and in
% vitro stored data (mean; no movie information).
% opts.SPs (10:15), opts.win (13), opts.minLen (180),
% opts.surrogate:
%   'none'    real traces
%   'shuffle' trace -> 61-frame moving-mean trend + its residuals in random
%             order. NOT a valid burst-free control: the trend still contains the
%             real bursts (smeared), so it is kept only as a record.
%   'poly'    burst-free control: trace -> its cubic trend (the slow bleaching
%             curve, as analysis.m fits the background) + white Gaussian noise
%             at the trace's own noise level, sigma = 1.4826*MAD(diff(s))/sqrt(2)
%             (robust to occasional steps); seeded. Every detection is false.
%   The lag-1 autocorrelation of the real first differences (white noise:
%   -0.5) is written per trace to <tag>_traces.csv (column acDiff) to check the
%   white-noise assumption.
% Writes analysis_outputs/robustness/<tag>_bursts.csv, _traces.csv, _sigma.csv.
if nargin < 2, opts = struct(); end
d0 = struct('SPs', 10:15, 'win', 13, 'minLen', 180, 'surrogate', 'none');
fn = fieldnames(d0);
for i = 1:numel(fn), if ~isfield(opts, fn{i}), opts.(fn{i}) = d0.(fn{i}); end, end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'robustness');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset'; SU = 'in_vivo_reextracted/after_R2_sums';
info = [readtable(fullfile(repoDir, R2, 'extraction_info.csv')); readtable(fullfile(repoDir, R2, 'extraction_info_5Qb_10Qb.csv'))];
jobs = {};
for d = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    jobs(end+1,:) = {d{1}, 'mean', R2, d{1}, false}; %#ok<AGROW>
    jobs(end+1,:) = {d{1}, 'spotsum', SU, [d{1} '_spotsum'], true}; %#ok<AGROW>
end
for d = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'}
    jobs(end+1,:) = {d{1}, 'mean', 'in_vitro_experiments', d{1}, false}; %#ok<AGROW>
end
tVec = (0:360) * 10/60; bu = cell(0, 8); tr = cell(0, 6); sg = cell(0, 5);
for j = 1:size(jobs,1)
    [S, B] = gm_loadDataset(repoDir, jobs{j,3}, jobs{j,4});
    I = info(strcmp(info.dataset, jobs{j,1}), :);
    if isempty(I), mov = zeros(numel(S), 1); else, mov = repelem((1:height(I))', I.nKept); end
    keep = cellfun(@numel, S) >= opts.minLen; S = S(keep); B = B(keep); mov = mov(keep);
    if jobs{j,5}, S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false); end
    acd = nan(numel(S), 1);
    for i = 1:numel(S)
        dd = diff(S{i}(:)); if numel(dd) > 3, cc = corrcoef(dd(1:end-1), dd(2:end)); acd(i) = cc(1,2); end
    end
    if strcmp(opts.surrogate, 'poly')
        rng(11 + j);
        for i = 1:numel(S)
            s = S{i}(:); x = (1:numel(s))'; trd = polyval(polyfit(x, s, 3), x);
            sgm = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
            S{i} = trd + sgm * randn(size(s));
        end
    end
    if strcmp(opts.surrogate, 'shuffle')
        rng(7 + j);
        for i = 1:numel(S)
            s = S{i}(:); trd = smoothdata(s, 'movmean', 61); r = s - trd;
            S{i} = trd + r(randperm(numel(r)));
        end
    end
    for i = 1:numel(S), tr(end+1,:) = {jobs{j,1}, jobs{j,2}, i, numel(S{i}), mov(i), acd(i)}; end %#ok<AGROW>
    for SP = opts.SPs
        ws = warning('off', 'all');
        [~, ~, ~, ~, av, ~, seg] = analysis(S, B, tVec, opts.win, SP, 0, 0);
        warning(ws);
        sg(end+1,:) = {jobs{j,1}, jobs{j,2}, SP, gm_sigma_from_ampvar(av(2:end)), size(seg,1)}; %#ok<AGROW>
        for q = 1:size(seg,1)
            bu(end+1,:) = {jobs{j,1}, jobs{j,2}, SP, seg(q,1), seg(q,2), seg(q,5), seg(q,6), seg(q,3)}; %#ok<AGROW>
        end
    end
    fprintf('%s %s done\n', jobs{j,1}, jobs{j,2});
end
writetable(cell2table(bu, 'VariableNames', {'dataset','measure','SP','trace','sign','amp','levelBefore','startFrame'}), ...
    fullfile(outDir, [tag '_bursts.csv']));
writetable(cell2table(tr, 'VariableNames', {'dataset','measure','trace','lengthFrames','movie','acDiff'}), fullfile(outDir, [tag '_traces.csv']));
writetable(cell2table(sg, 'VariableNames', {'dataset','measure','SP','sigma','nBursts'}), fullfile(outDir, [tag '_sigma.csv']));
disp(['EXPORT2_' tag '_DONE']);
end
