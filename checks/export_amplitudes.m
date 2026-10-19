% EXPORT_AMPLITUDES  Writes the burst amplitudes that checks/mle_real_data.m
% fits (same datasets, settings and code path: gm_amplitudes_by_trace, traces
% >= 180 frames), for the histogram + fitted-model figure
% (checks/plot_mle_histograms.py). Reviewer points R2 4a/4b (goodness of fit).
% Writes analysis_outputs/mle_real/amplitudes.csv (dataset, sign, trace, amp).
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'mle_real');
R3 = 'in_vivo_reextracted/after_R3_aperture';
D = {'protein_RNA_1_1','in_vitro_experiments',12; 'protein_RNA_10_1','in_vitro_experiments',12;
     'protein_RNA_100_1','in_vitro_experiments',12; '4PP7',R3,10; '4PP7-FUS',R3,10;
     '30PP7',R3,12; '30PP7-FUS',R3,12; '5Qb','in_vivo_experiments',14; '10Qb','in_vivo_experiments',14};
tVec = (0:360) * 10/60; win = 13;
ds = {}; sg = {}; tr = []; amp = [];
for d = 1:size(D,1)
    [S, B] = gm_loadDataset(repoDir, D{d,2}, D{d,1});
    T = gm_amplitudes_by_trace(S, B, tVec, win, D{d,3});
    T = T([T.len] >= 180);
    for i = 1:numel(T)
        for s = {'Neg','Pos'}
            a = abs(T(i).(['amp' s{1}])); a = a(isfinite(a));
            ds = [ds; repmat(D(d,1), numel(a), 1)]; sg = [sg; repmat(s, numel(a), 1)]; %#ok<AGROW>
            tr = [tr; repmat(i, numel(a), 1)]; amp = [amp; a(:)]; %#ok<AGROW>
        end
    end
end
writetable(table(ds, sg, tr, amp, 'VariableNames', {'dataset','sign','trace','amp'}), fullfile(outDir, 'amplitudes.csv'));
disp('EXPORT_DONE');
