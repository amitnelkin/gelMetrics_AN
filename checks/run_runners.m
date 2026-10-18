function T = run_runners(stage, doExp, doSim)
% RUN_RUNNERS  Runs ExperimentalAnalysis.m on every dataset (paper settings
% from the comment on its line 12) and SimulationRunWithAnalysis.m, and
% writes analysis_outputs/<stage>/lambda_k0_table.csv.
if nargin < 2, doExp = true; end
if nargin < 3, doSim = true; end
checksDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(checksDir);
outDir = fullfile(repoDir, 'analysis_outputs', stage);
if ~exist(outDir, 'dir'), mkdir(outDir); end

% name, folder, SP, GOF  (ExperimentalAnalysis.m line 12)
datasets = { ...
    'protein_RNA_1_1',   'in_vitro_experiments', 12, 15;
    'protein_RNA_10_1',  'in_vitro_experiments', 12, 10;
    'protein_RNA_100_1', 'in_vitro_experiments', 12, 10;
    '4PP7',              'in_vivo_experiments',  10, 15;
    '4PP7-FUS',          'in_vivo_experiments',  10, 15;
    '30PP7',             'in_vivo_experiments',  12, 10;
    '30PP7-FUS',         'in_vivo_experiments',  12, 10;
    '5Qb',               'in_vivo_experiments',  14, 15;
    '10Qb',              'in_vivo_experiments',  14, 15};

rows = {};
if doExp
    for d = 1:size(datasets, 1)
        r = run_experimental(datasets{d,1}, datasets{d,2}, datasets{d,3}, datasets{d,4}, outDir);
        rows = [rows; r]; %#ok<AGROW>
    end
end
if doSim
    rows = [rows; run_simulation(outDir)];
end
T = cell2table(rows, 'VariableNames', {'dataset','SP','GOF','sign','nAmpAll','nFit', ...
    'lam1','k0_1','lamK0_1','score1','lam2','k0_2','score2','k0AtGridEdge','nTraces','nTraces180'});
T.stage = repmat({stage}, height(T), 1);
disp(T);
writetable(T, fullfile(outDir, 'lambda_k0_table.csv'));
end
