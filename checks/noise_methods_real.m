function noise_methods_real(r2Folder, suFolder, outFile, doVitro)
% NOISE_METHODS_REAL  Detection-free fluctuation measure for the real data
% (RESULTS.md section 19; reviewer points R2 6, R2 4a, "comparison with other
% methods"). Granule intensity g per trace: the spot sum (in vivo, R2 sums) or
% the spot mean minus the cell level (mean measure; in vivo R2 and in vitro
% stored data). Relative intensity y = g / cubicTrend(g) - 1 (traces with a
% non-positive trend are skipped). Structure function D(tau) = mean (y(t+tau) -
% y(t))^2 for tau = 1 2 3 6 12 18 30 45 60 frames, for the real trace and for
% the mean of 3 burst-free surrogates (cubic trend + white noise at the trace's
% own level, seeded), so D_excess = D_real - D_surr is the fluctuation beyond noise.
% Movie id per trace as in score_export2.m. Traces >= 180 frames.
% Writes analysis_outputs/noise_methods/real_sf.csv.
% Optional (section 23): r2Folder, suFolder (mean and spot-sum data), outFile
% (path under analysis_outputs), doVitro (default true).
if nargin < 1, r2Folder = 'in_vivo_reextracted/after_R2_offset'; end
if nargin < 2, suFolder = 'in_vivo_reextracted/after_R2_sums'; end
if nargin < 3, outFile = fullfile('noise_methods', 'real_sf.csv'); end
if nargin < 4, doVitro = true; end
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'noise_methods');
if ~exist(outDir, 'dir'), mkdir(outDir); end
R2i = 'in_vivo_reextracted/after_R2_offset'; R2 = r2Folder; SU = suFolder;
info = [readtable(fullfile(repoDir, R2i, 'extraction_info.csv')); readtable(fullfile(repoDir, R2i, 'extraction_info_5Qb_10Qb.csv'))];
jobs = {};
for d = {'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb'}
    jobs(end+1,:) = {d{1}, 'mean', R2, d{1}, false}; %#ok<AGROW>
    jobs(end+1,:) = {d{1}, 'spotsum', SU, [d{1} '_spotsum'], true}; %#ok<AGROW>
end
for d = {'protein_RNA_1_1', 'protein_RNA_10_1', 'protein_RNA_100_1'}
    if doVitro, jobs(end+1,:) = {d{1}, 'mean', 'in_vitro_experiments', d{1}, false}; end %#ok<AGROW>
end
taus = [1 2 3 6 12 18 30 45 60]; rows = cell(0, 8);
for j = 1:size(jobs,1)
    [S, B] = gm_loadDataset(repoDir, jobs{j,3}, jobs{j,4});
    I = info(strcmp(info.dataset, jobs{j,1}), :);
    if isempty(I), mov = zeros(numel(S), 1); else, mov = repelem((1:height(I))', I.nKept); end
    keep = cellfun(@numel, S) >= 180; S = S(keep); B = B(keep); mov = mov(keep);
    rng(31 + j);
    for i = 1:numel(S)
        if jobs{j,5}, g = S{i}(:)'; else, g = S{i}(:)' - B{i}(:)'; end
        t = 1:numel(g); trd = polyval(polyfit(t, g, 3), t);
        if any(trd <= 0), continue, end
        Dr = sf(g ./ trd - 1, taus);
        sgm = 1.4826 * median(abs(diff(g) - median(diff(g)))) / sqrt(2);
        Ds = zeros(1, numel(taus));
        for k = 1:3
            gs = trd + sgm * randn(size(g));
            ts = polyval(polyfit(t, gs, 3), t);
            Ds = Ds + sf(gs ./ ts - 1, taus) / 3;
        end
        rows(end+1,:) = {jobs{j,1}, jobs{j,2}, i, mov(i), numel(g), mat2str(taus), mat2str(Dr, 6), mat2str(Ds, 6)}; %#ok<AGROW>
    end
    fprintf('%s %s done\n', jobs{j,1}, jobs{j,2});
end
writetable(cell2table(rows, 'VariableNames', {'dataset','measure','trace','movie','lengthFrames','taus','Dreal','Dsurr'}), ...
    fullfile(repoDir, 'analysis_outputs', outFile));
disp('NOISE_REAL_SF_DONE');
end

function D = sf(y, taus)
D = arrayfun(@(k) mean((y(1+k:end) - y(1:end-k)).^2), taus);
end
