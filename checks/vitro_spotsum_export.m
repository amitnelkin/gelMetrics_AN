function vitro_spotsum_export()
% VITRO_SPOTSUM_EXPORT  In vitro data for the manuscript figures with the spot sum (the
% user's decision, RESULTS.md section 53): GelMetricsReport on gm_paper_dataset (spot sum),
% paper seeds, nCopies 1.
% Writes analysis_outputs/vitro_spotsum/
%   bursts.csv   dataset, kind (real / surr), f (burst size, fraction of the granule)
%   sp_grid.csv  dataset, SP, realNeg, falseNeg (bursts in granules / copies, window 13)
%   summary.csv  dataset, label, every GelMetricsReport measure with _lo/_hi, halves p
repo = fileparts(fileparts(mfilename('fullpath'))); addpath(repo, fullfile(repo, 'analysis'), fullfile(repo, 'legacy'));
out = fullfile(repo, 'analysis_outputs', 'vitro_spotsum'); if ~exist(out, 'dir'), mkdir(out); end
V = {'vitro_1_1', 'protein_RNA_1_1'; 'vitro_10_1', 'protein_RNA_10_1'; 'vitro_100_1', 'protein_RNA_100_1'};
F = {'brightness', 'brightnessPooled', 'frameNoise', 'meanBurst', 'medianBurst', 'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'noiseShare', 'net', 'turnover1', 'turnover2', 'net1', 'net2'};
bu = {}; sg = {}; su = {};
for k = 1:size(V, 1)
    D = gm_paper_dataset(V{k, 1});
    R = GelMetricsReport(D, struct('seed', D.paperSeed, 'nCopies', 1, 'verbose', false));
    r = vertcat(R.items.r0); c = vertcat(R.items.c0);
    bu = [bu; repmat(V(k, 2), numel(r), 1), repmat({'real'}, numel(r), 1), num2cell(r); ...
          repmat(V(k, 2), numel(c), 1), repmat({'surr'}, numel(c), 1), num2cell(c)]; %#ok<AGROW>
    row = {V{k, 2}, R.granules, R.movies};
    for q = 1:numel(F), row = [row, {R.(F{q}), R.([F{q} '_lo']), R.([F{q} '_hi'])}]; end %#ok<AGROW>
    su(end+1, :) = [row, {R.pHalvesTurnover, R.pHalvesNet}]; %#ok<AGROW>
    for SP = 10:20
        Rs = GelMetricsReport(D, struct('seed', D.paperSeed, 'nCopies', 1, 'SP', SP, 'nBoot', 1, 'verbose', false));
        sg(end+1, :) = {V{k, 2}, SP, numel(vertcat(Rs.items.r0)), numel(vertcat(Rs.items.c0))}; %#ok<AGROW>
    end
    fprintf('%s done\n', V{k, 1});
end
writetable(cell2table(bu, 'VariableNames', {'dataset', 'kind', 'f'}), fullfile(out, 'bursts.csv'));
writetable(cell2table(sg, 'VariableNames', {'dataset', 'SP', 'realNeg', 'falseNeg'}), fullfile(out, 'sp_grid.csv'));
names = {'dataset', 'granules', 'movies'};
for q = 1:numel(F), names = [names, {F{q}, [F{q} '_lo'], [F{q} '_hi']}]; end %#ok<AGROW>
writetable(cell2table(su, 'VariableNames', [names, {'pHalvesTurnover', 'pHalvesNet'}]), fullfile(out, 'summary.csv'));
disp('EXPORT_DONE');
end
