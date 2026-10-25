function rows = run_simulation(outDir)
% RUN_SIMULATION  Runs the repo's SimulationRunWithAnalysis.m without
% editing it: script text copied to a temp file with close/clear removed,
% rng(1) set first (the same seed the bug-fix commit adds to the file),
% plotting removed (plotMax = 0, plot lines stripped, modifiedPoissonFit
% via gm_noplot; figure rendering hangs in headless MATLAB here), and
% modifiedPoissonFit( routed through gm_logged_mpfit for logging.
% Saves sim_run.mat in outDir and returns one table row per sign.

checksDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir);
global GM_FITLOG GM_RUNOUT
GM_FITLOG = {}; GM_RUNOUT = [];

src = strrep(fileread(fullfile(repoDir, 'legacy', 'SimulationRunWithAnalysis_WS.m')), sprintf('\r'), '');
src = regexprep(src, '^\s*(close all|clearvars|clear|clc)\s*;[^\n]*$', '', 'lineanchors');
src = gm_strip_plot_lines(src);
src = regexprep(src, 'plotMax\s*=\s*\d+\s*;', 'plotMax = 0;', 'once');
src = strrep(src, 'modifiedPoissonFit(', 'gm_logged_mpfit(');
src = ['rng(1);' newline src newline 'global GM_RUNOUT; GM_RUNOUT = struct(''nPos'', numel(ampPos), ''nNeg'', numel(ampNeg), ''lambda_sim'', lambda_sim, ''SP'', stringencyParameter, ''numTracks'', numTracks, ''expConst'', expConst, ''ampPos'', ampPos, ''ampNeg'', ampNeg, ''ampVar'', ampVar, ''eventsPerSignal'', eventsPerSignal);' newline];
noplotDir = gm_noplot();

tmpScript = [tempname '.m'];
fid = fopen(tmpScript, 'w'); fwrite(fid, src); fclose(fid);
cleanup = onCleanup(@() localCleanup(tmpScript, noplotDir));

fprintf('\n######## SimulationRunWithAnalysis ########\n');
run(tmpScript);

if ~exist(outDir, 'dir'), mkdir(outDir); end
R = GM_RUNOUT; R.name = 'SIM'; R.fits = GM_FITLOG; %#ok<NASGU>
save(fullfile(outDir, 'SIM_run.mat'), '-struct', 'R');

signsTxt = {'Positive amplitudes', 'Negative amplitudes'};
nAll = [GM_RUNOUT.nPos, GM_RUNOUT.nNeg];
rows = {};
for s = 1:2
    e = GM_FITLOG{find(cellfun(@(x) strcmp(x.title, signsTxt{s}), GM_FITLOG), 1, 'last')};
    rows(end+1, :) = {sprintf('SIM(lam=%g,k0=60)', GM_RUNOUT.lambda_sim), GM_RUNOUT.SP, 15, signsTxt{s}(1:3), ...
        nAll(s), e.nFit, e.lamOrder(1), e.k0Order(1), e.lamOrder(1)*e.k0Order(1), e.score(1), ...
        e.lamOrder(2), e.k0Order(2), e.score(2), ...
        e.k0Order(1) == e.K0min || e.k0Order(1) == e.K0max, GM_RUNOUT.numTracks, GM_RUNOUT.numTracks}; %#ok<AGROW>
end
end

function localCleanup(tmpScript, noplotDir)
if exist(tmpScript, 'file'), delete(tmpScript); end
rmpath(noplotDir);
rmdir(noplotDir, 's');
end
