function export_amps_shape()
% EXPORT_AMPS_SHAPE  Negative-burst amplitudes in manuscript units (spot mean
% measure, as the manuscript's lambda/k0), window 13, SP 12, traces >= 180
% frames, for the shape-only lambda/k0 fit (checks/lambda_k0_shape.py,
% RESULTS.md section 26):
%   in vivo  R2 re-extraction (in_vivo_reextracted/after_R2_offset), with movie id;
%   in vitro stored data (in_vitro_experiments);
%   simulations generateSim type 3, k0 = 60, lambda 0.5 1 2 3 4 5 x noise 12 / 40
%   / 120 (SNR k0/noise 5, 1.5, 0.5), 200 tracks x 361 frames, expConst 0.0026,
%   seeded; the same code path (gm_amplitudes_by_trace -> analysis.m).
% Writes analysis_outputs/lambda_k0_shape/amplitudes.csv
% (kind, dataset, lambdaTrue, noise, trace, movie, amp).
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'lambda_k0_shape'); if ~exist(outDir, 'dir'), mkdir(outDir); end
R2 = 'in_vivo_reextracted/after_R2_offset';
info = [readtable(fullfile(repoDir, R2, 'extraction_info.csv')); readtable(fullfile(repoDir, R2, 'extraction_info_5Qb_10Qb.csv'))];
tVec = (0:360) * 10/60; C = cell(0, 7);
jobs = {'in vivo', R2, {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}; ...
        'in vitro', 'in_vitro_experiments', {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'}};
for j = 1:size(jobs, 1)
    for d = jobs{j,3}
        [S, B] = gm_loadDataset(repoDir, jobs{j,2}, d{1});
        I = info(strcmp(info.dataset, d{1}), :);
        if isempty(I), mov = zeros(numel(S), 1); else, mov = repelem((1:height(I))', I.nKept); end
        C = addAmps(C, jobs{j,1}, d{1}, NaN, NaN, S, B, mov, tVec);
        fprintf('%s done\n', d{1});
    end
end
lams = [0.5 1 2 3 4 5]; noises = [12 40 120]; k = 0;
for lam = lams
    for nz = noises
        k = k + 1; rng(700 + k);
        [S, B] = generateSim(200, 361, 3, 0.0026, nz, lam, 60);
        C = addAmps(C, 'simulated', sprintf('lambda %g, noise %d', lam, nz), lam, nz, num2cell(S, 2), num2cell(B, 2), zeros(200, 1), tVec);
    end
end
writetable(cell2table(C, 'VariableNames', {'kind','dataset','lambdaTrue','noise','trace','movie','amp'}), fullfile(outDir, 'amplitudes.csv'));
disp('AMPS_DONE');
end

function C = addAmps(C, kind, name, lam, nz, S, B, mov, tVec)
ws = warning('off', 'all');
T = gm_amplitudes_by_trace(S, B, tVec, 13, 12);
warning(ws);
for i = 1:numel(T)
    if T(i).len < 180, continue, end
    a = abs(T(i).ampNeg);
    for q = 1:numel(a)
        C(end+1,:) = {kind, name, lam, nz, i, mov(i), a(q)}; %#ok<AGROW>
    end
end
end
