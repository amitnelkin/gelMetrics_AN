function reextract_pp7(workRoot, stages, datasets, outRoot)
% REEXTRACT_PP7  Re-extracts the PP7 datasets from the .nd2 movies with the
% fixed extraction (zigzag cause: analysis_outputs/zigzag_diagnosis; general
% data validity, R2 5). Never touches in_vivo_experiments/: writes
% in_vivo_reextracted/<stage>/<dataset>_data.mat and _background.mat.
%
% workRoot/<dataset>/work must hold <movie>.tif (plain bfconvert output) and
% <movie>.csv (Mosaic Traj_<movie>.tif.csv from the R1 BatchParticleTracker.ijm,
% renamed as the user did). read_tables.m is run there, then
% collect_intensity.m once per stage:
%   after_R1_redchannel  red channel only + frame/xy mapping (offset 0, thresholded spot)
%   after_R2_offset      + camera offset 200
%   (after_R3_aperture, a 5 px aperture spot, was reverted: not a stage any more)
% Also writes in_vivo_reextracted/<stage>/extraction_info.csv.
% outRoot (optional, default 'in_vivo_reextracted'): output folder in the repo, e.g.
% 'in_vitro_reextracted' for the titration movies (RESULTS.md section 38). The work
% folder may hold the fluorescence-only TIFF written by the tracking macro instead of
% the bfconvert output: the frames read are identical (checked on 4PP7 005, 30PP7 001).

if nargin < 2 || isempty(stages)
    stages = {'after_R1_redchannel', struct('offset',0,  'aperture',0);
              'after_R2_offset',     struct('offset',200,'aperture',0);
              };   % R3 (5 px aperture) was reverted at the user's request: multithresh spot kept
end
checksDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(checksDir);
addpath(repoDir);
if nargin < 3, datasets = {'4PP7','4PP7-FUS','30PP7','30PP7-FUS'}; end  % also used for 5Qb/10Qb
if nargin < 4, outRoot = 'in_vivo_reextracted'; end

here = pwd;
for d = 1:numel(datasets)
    runReadTables(fullfile(workRoot, datasets{d}, 'work'));
end
cd(here);

for s = 1:size(stages,1)
    outDir = fullfile(repoDir, outRoot, stages{s,1});
    if ~exist(outDir,'dir'), mkdir(outDir); end
    allInfo = table();
    for d = 1:numel(datasets)
        work = fullfile(workRoot, datasets{d}, 'work');
        fprintf('%s  %s\n', stages{s,1}, datasets{d});
        [intensity, background, info] = collect_intensity(work, stages{s,2});
        save(fullfile(outDir, [datasets{d} '_data.mat']), 'intensity');
        save(fullfile(outDir, [datasets{d} '_background.mat']), 'background');
        I = struct2table(info, 'AsArray', true);
        I.dataset = repmat(datasets(d), height(I), 1);
        allInfo = [allInfo; I]; %#ok<AGROW>
    end
    infoName = 'extraction_info.csv';
    if nargin >= 3, infoName = ['extraction_info_' strjoin(datasets, '_') '.csv']; end
    writetable(allInfo, fullfile(outDir, infoName));
end
end

function runReadTables(work)
% read_tables.m is a script that calls clear: run it in its own workspace
cd(work);
read_tables;                           % writes <movie>_xy.mat next to <movie>.csv
end
