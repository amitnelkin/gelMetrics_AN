function score_export(what)
% SCORE_EXPORT  Burst tables for comparing candidate summary scores
% (checks/score_compare.py): every classified burst (both signs) with its
% trace, amplitude, granule level before it and duration, plus per trace its
% length, and per dataset the quiescent noise sigma (ampVar). Same detection as
% the runners (analysis.m, window 13). Reviewer points R2 4a/4b/7.
% what = 'data': in vivo (R2 re-extraction: mean and spot-sum measures) and in
%        vitro (stored data, mean measure); SP 10..15.
% what = 'sim':  seeded simulations with known lambda (0.5,1,2,3,5) x k0 (60,120)
%        x noise SD (40,80) x seeds 1:2, 150 tracks, SP 12; measure = the
%        simulated signal (granule level known from the simulation).
% Writes analysis_outputs/scores/<what>_bursts.csv, <what>_traces.csv, <what>_sigma.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'scores');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tVec = (0:360) * 10/60; bu = cell(0, 9); tr = cell(0, 5); sg = cell(0, 5);
if strcmp(what, 'data')
    R2 = 'in_vivo_reextracted/after_R2_offset'; SU = 'in_vivo_reextracted/after_R2_sums';
    jobs = {};
    for d = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
        jobs(end+1,:) = {d{1}, 'mean', R2, d{1}, false}; %#ok<AGROW>
        jobs(end+1,:) = {d{1}, 'spotsum', SU, [d{1} '_spotsum'], true}; %#ok<AGROW>
    end
    for d = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'}
        jobs(end+1,:) = {d{1}, 'mean', 'in_vitro_experiments', d{1}, false}; %#ok<AGROW>
    end
    SPs = 10:15;
else
    jobs = {};
    for lam = [0.5 1 2 3 5], for k0 = [60 120], for nz = [40 80], for s = 1:2
        jobs(end+1,:) = {sprintf('sim_l%g_k%d_n%d_s%d', lam, k0, nz, s), 'sim', lam, k0, [nz s]}; %#ok<AGROW>
    end, end, end, end
    SPs = 12;
end
for j = 1:size(jobs,1)
    if strcmp(what, 'data')
        [S, B] = gm_loadDataset(repoDir, jobs{j,3}, jobs{j,4});
        keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep);
        if jobs{j,5}, S = cellfun(@(s, b) s(:) + b(:), S, B, 'UniformOutput', false); end
    else
        p = jobs{j,5}; rng(5000 + 100*p(2) + round(10*jobs{j,3}) + jobs{j,4} + p(1));
        [Sm, Bm] = generateSim(150, 361, 3, 0.0026, p(1), jobs{j,3}, jobs{j,4});
        S = num2cell(Sm, 2); B = num2cell(Bm, 2);
    end
    for i = 1:numel(S), tr(end+1,:) = {jobs{j,1}, jobs{j,2}, i, numel(S{i}), 0}; end %#ok<AGROW>
    for SP = SPs
        ws = warning('off', 'all');
        [~, ~, ~, ~, av, ~, seg] = analysis(S, B, tVec, 13, SP, 0, 0);
        warning(ws);
        sg(end+1,:) = {jobs{j,1}, jobs{j,2}, SP, gm_sigma_from_ampvar(av(2:end)), size(seg,1)}; %#ok<AGROW>
        for q = 1:size(seg,1)
            bu(end+1,:) = {jobs{j,1}, jobs{j,2}, SP, seg(q,1), seg(q,2), seg(q,5), seg(q,6), seg(q,4) - seg(q,3), 0}; %#ok<AGROW>
        end
    end
    fprintf('%s %s done\n', jobs{j,1}, jobs{j,2});
end
writetable(cell2table(bu(:,1:8), 'VariableNames', {'dataset','measure','SP','trace','sign','amp','levelBefore','durationFrames'}), ...
    fullfile(outDir, [what '_bursts.csv']));
writetable(cell2table(tr(:,1:4), 'VariableNames', {'dataset','measure','trace','lengthFrames'}), fullfile(outDir, [what '_traces.csv']));
writetable(cell2table(sg, 'VariableNames', {'dataset','measure','SP','sigma','nBursts'}), fullfile(outDir, [what '_sigma.csv']));
disp(['SCORE_EXPORT_' upper(what) '_DONE']);
end
