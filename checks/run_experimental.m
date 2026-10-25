function rows = run_experimental(name, folder, SP, GOF, outDir)
% RUN_EXPERIMENTAL  Runs the repo's ExperimentalAnalysis.m on one dataset
% without editing it. The script text is copied to a temp file with:
%   - close all / clearvars / clc removed,
%   - the dataset load inserted (the script has all loads commented out),
%   - stringencyParameter set to SP,
%   - the truncation percentile set to GOF (through truncPercent after the
%     bug-fix commit; before it, through gm_noplot's copy of
%     modifiedPoissonFit.m with "percent = 15" replaced, which is how the
%     paper figures with GOF = 10 were made),
%   - plotting lines removed and modifiedPoissonFit run without its plotting
%     section (gm_noplot): figure rendering hangs in headless MATLAB here,
%   - modifiedPoissonFit( routed through gm_logged_mpfit (logging only).
% Saves <name>_run.mat (amplitudes + fit score matrices) in outDir for
% plotting, and returns one table row per sign.

checksDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(checksDir);
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(checksDir);
global GM_FITLOG GM_RUNOUT
GM_FITLOG = {}; GM_RUNOUT = [];

src = strrep(fileread(fullfile(repoDir, 'legacy', 'ExperimentalAnalysis_WS.m')), sprintf('\r'), '');
src = regexprep(src, '^\s*(close all|clearvars|clear|clc)\s*;[^\n]*$', '', 'lineanchors');
src = gm_strip_plot_lines(src);
src = regexprep(src, 'stringencyParameter\s*=\s*\d+\s*;', sprintf('stringencyParameter = %d;', SP), 'once');
loadStr = sprintf(['signals = struct2cell(load(''%s''));\n' ...
                   'backgrounds = struct2cell(load(''%s''));\n'], ...
    fullfile(repoDir, folder, [name '_data.mat']), fullfile(repoDir, folder, [name '_background.mat']));
src = strrep(src, '% Signal plotting controls', [loadStr '% Signal plotting controls']);
src = strrep(src, 'modifiedPoissonFit(', 'gm_logged_mpfit(');
src = [src newline 'global GM_RUNOUT; GM_RUNOUT = struct(''nPos'', numel(ampPos), ''nNeg'', numel(ampNeg), ''nTraces'', numel(signals), ''nTraces180'', sum(cellfun(@numel, signals) >= 180), ''ampPos'', ampPos, ''ampNeg'', ampNeg, ''ampVar'', ampVar, ''timeBetweenPos'', timeBetweenPos, ''timeBetweenNeg'', timeBetweenNeg, ''eventsPerSignal'', eventsPerSignal);' newline];

if contains(src, 'truncPercent = 15;')            % after the bug-fix commit
    src = strrep(src, 'truncPercent = 15;', sprintf('truncPercent = %g;', GOF));
    noplotDir = gm_noplot();
else                                              % original code
    noplotDir = gm_noplot(GOF);
end

tmpScript = [tempname '.m'];
fid = fopen(tmpScript, 'w'); fwrite(fid, src); fclose(fid);
cleanup = onCleanup(@() localCleanup(tmpScript, noplotDir));

fprintf('\n######## %s (SP=%d, GOF=%g) ########\n', name, SP, GOF);
run(tmpScript);   % executes in this function's workspace

if ~exist(outDir, 'dir'), mkdir(outDir); end
R = GM_RUNOUT; R.name = name; R.SP = SP; R.GOF = GOF; R.fits = GM_FITLOG; %#ok<NASGU>
save(fullfile(outDir, [name '_run.mat']), '-struct', 'R');

signsTxt = {'Positive amplitudes', 'Negative amplitudes'};
nAll = [GM_RUNOUT.nPos, GM_RUNOUT.nNeg];
rows = {};
for s = 1:2
    e = GM_FITLOG{find(cellfun(@(x) strcmp(x.title, signsTxt{s}), GM_FITLOG), 1, 'last')};
    rows(end+1, :) = {name, SP, GOF, signsTxt{s}(1:3), nAll(s), e.nFit, ...
        e.lamOrder(1), e.k0Order(1), e.lamOrder(1)*e.k0Order(1), e.score(1), ...
        e.lamOrder(2), e.k0Order(2), e.score(2), ...
        e.k0Order(1) == e.K0min || e.k0Order(1) == e.K0max, GM_RUNOUT.nTraces, GM_RUNOUT.nTraces180}; %#ok<AGROW>
end
end

function localCleanup(tmpScript, noplotDir)
if exist(tmpScript, 'file'), delete(tmpScript); end
rmpath(noplotDir);
rmdir(noplotDir, 's');
end
