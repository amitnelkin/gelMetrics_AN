function d = gm_noplot(percentOverride)
% GM_NOPLOT  Puts on top of the path a temporary copy of the repo's current
% modifiedPoissonFit.m with ONLY its plotting section removed (everything
% from "% Plot a subset of the results" to the end), so it can run in
% headless MATLAB, where figure rendering hangs on this machine. The scoring
% code is untouched; figures are redrawn from the saved fit data in Python
% (checks/plot_runner_outputs.py).
%
% percentOverride (optional): for the ORIGINAL file only, replaces the
% hardcoded "percent = 15;" (how the paper figures with GOF = 10 were made).
% Returns the temp folder; remove it with rmpath(d); rmdir(d,'s').
repoDir = fileparts(fileparts(mfilename('fullpath')));
src = fileread(fullfile(repoDir, 'legacy', 'modifiedPoissonFit.m'));
i0 = strfind(src, '% Plot a subset of the results');
assert(numel(i0) == 1, 'gm_noplot: plotting section marker not found');
src = [src(1:i0-1) sprintf('pvec = []; %% plotting removed by gm_noplot\r\nend\r\n')];
if nargin >= 1 && ~isempty(percentOverride)
    assert(contains(src, 'percent = 15;'), 'gm_noplot: "percent = 15;" not found');
    src = regexprep(src, 'percent = 15;', sprintf('percent = %g;', percentOverride), 'once');
end
d = tempname; mkdir(d);
fid = fopen(fullfile(d, 'modifiedPoissonFit.m'), 'w'); fwrite(fid, src); fclose(fid);
addpath(d, '-begin');
rehash path;
end
