function vitro_same_granules()
% VITRO_SAME_GRANULES  The stored in vitro granules, re-read with the corrected
% pipeline (RESULTS.md section 43; the user's request: same granules and traces as
% the original data, but with the offset and the corrected extraction).
% 1. Track of every stored trace: each stored point (spot, cell) is matched exactly
%    to a point of the original-code run on the re-tracked movies
%    (checks/vitro_orig_verify.m, ..\zigzag_scratch\orig_vitro), which gives its
%    movie, Mosaic frame and x, y. Unmatched points inside a trace get frame and x, y
%    by linear interpolation between matched neighbours; unmatched points at the
%    ends are dropped.
% 2. These tracks (one per stored trace, same points as stored) are written as
%    <movie>_xy.mat into ..\zigzag_scratch\same_granules\<folder> and read with
%    collect_intensity (R2: page = frame + 1, x + 1, y + 1, offset 200, per-frame
%    multithresh), exactly as the re-extracted data.
% Writes in_vitro_reextracted/same_granules_R2/<same_*>_data.mat, _background.mat
% (traces in movie order), extraction_info_same.csv (movie ids) and
% same_granules_map.csv (stored trace -> movie, points matched, length).
% Section 53: also <same_*>_spotsum_data.mat (variable intensity: integrated granule
% light above the cell mean, collect_intensity 4th output) and _spotsum_background.mat
% (the same cell mean as _background.mat), as in_vivo_reextracted/after_R2_sums. Existing
% _data.mat / _background.mat files are not rewritten: the re-run is checked to
% reproduce them exactly.
Z = 'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\zigzag_scratch';
repo = fileparts(fileparts(mfilename('fullpath'))); addpath(repo);
outDir = fullfile(repo, 'in_vitro_reextracted', 'same_granules_R2'); if ~exist(outDir, 'dir'), mkdir(outDir); end
J = {'protein_RNA_1_1', 'vitro_1000ngul_0p1mgml', 'same_1_1'; 'protein_RNA_10_1', 'vitro_1000ng_10ug', 'same_10_1'; ...
     'protein_RNA_100_1', 'vitro_1000ngul_1mgml', 'same_100_1'};
info = cell(0, 7); mp = cell(0, 8);
for j = 1:size(J, 1)
    w = fullfile(Z, 'orig_vitro', J{j,2}); sw = fullfile(Z, 'same_granules', J{j,2});
    d = load(fullfile(w, 'data.mat')); b = load(fullfile(w, 'background.mat'));
    D = d.intensity(2:end); Bg = b.background(2:end);
    xs = dir(fullfile(w, '*_xy.mat')); [~, o] = sort({xs.name}); xs = xs(o);
    P = zeros(0, 6); k = 0;                % spot, cell, movie, Mosaic frame, x, y
    for f = 1:numel(xs)
        xy = load(fullfile(w, xs(f).name));
        for t = 1:numel(xy.x)
            fr = xy.frame{t}(2:end); x = xy.x{t}(2:end); y = xy.y{t}(2:end);   % as the original code
            if numel(x) < 24, continue, end
            k = k + 1;
            P = [P; double(D{k}(:)), double(Bg{k}(:)), repmat(f, numel(x), 1), fr(:), x(:), y(:)]; %#ok<AGROW>
        end
    end
    assert(k == numel(D));
    S = struct2cell(load(fullfile(repo, 'in_vitro_experiments-old', [J{j,1} '_data.mat']))); S = S{1}(2:end);
    B = struct2cell(load(fullfile(repo, 'in_vitro_experiments-old', [J{j,1} '_background.mat']))); B = B{1}(2:end);
    tr = struct('movie', {}, 'frame', {}, 'x', {}, 'y', {}, 'stored', {});
    for i = 1:numel(S)
        [ok, loc] = ismember([double(S{i}(:)), double(B{i}(:))], P(:, 1:2), 'rows');
        if ~any(ok), mp(end+1,:) = {J{j,1}, i, NaN, numel(S{i}), 0, 0, 0, 'no match'}; continue, end %#ok<AGROW>
        mv = mode(P(loc(ok), 3)); ok(ok) = P(loc(ok), 3) == mv;
        n = (1:numel(S{i}))'; F = nan(size(n)); X = F; Y = F;
        F(ok) = P(loc(ok), 4); X(ok) = P(loc(ok), 5); Y(ok) = P(loc(ok), 6);
        a = find(ok, 1); e = find(ok, 1, 'last'); in = (a:e)';
        F = round(interp1(n(ok), F(ok), in)); X = interp1(n(ok), X(ok), in); Y = interp1(n(ok), Y(ok), in);
        note = '';
        if any(diff(F) <= 0), note = 'frames not increasing'; end
        tr(end+1) = struct('movie', mv, 'frame', F, 'x', X, 'y', Y, 'stored', i); %#ok<AGROW>
        mp(end+1,:) = {J{j,1}, i, mv, numel(S{i}), sum(ok), numel(in), F(end) - F(1) + 1, note}; %#ok<AGROW>
    end
    % one _xy.mat per movie, stored order inside each movie
    order = [];
    for f = 1:numel(xs)
        q = find([tr.movie] == f & arrayfun(@(r) numel(r.x) >= 24 && all(diff(r.frame) > 0), tr));
        x = arrayfun(@(r) r.x, tr(q), 'UniformOutput', false); y = arrayfun(@(r) r.y, tr(q), 'UniformOutput', false);
        frame = arrayfun(@(r) r.frame, tr(q), 'UniformOutput', false);
        save(fullfile(sw, xs(f).name), 'x', 'y', 'frame');
        order = [order, [tr(q).stored]]; %#ok<AGROW>
        info(end+1,:) = {xs(f).name(1:3), 1, 1, 361, numel(q), numel(q), J{j,3}}; %#ok<AGROW>
    end
    [intensity, background, ~, sums] = collect_intensity(sw, 200);
    assert(numel(intensity) - 1 == numel(order));
    fD = fullfile(outDir, [J{j,3} '_data.mat']); fB = fullfile(outDir, [J{j,3} '_background.mat']);
    if isfile(fD) && isfile(fB)
        o1 = load(fD); o2 = load(fB);
        same = isequal(o1.intensity, intensity) && isequal(o2.background, background);
        fprintf('%s: spot mean and background identical to the saved files: %d\n', J{j,3}, same);
        assert(same, 'the re-run does not reproduce %s', fD);
    else
        save(fD, 'intensity'); save(fB, 'background');
    end
    intensity = sums.spotSum; save(fullfile(outDir, [J{j,3} '_spotsum_data.mat']), 'intensity');
    save(fullfile(outDir, [J{j,3} '_spotsum_background.mat']), 'background');
    writematrix(order(:), fullfile(outDir, [J{j,3} '_stored_index.csv']));
    fprintf('%s: %d stored traces, %d tracks rebuilt, %d read\n', J{j,1}, numel(S), numel(tr), numel(order));
end
writetable(cell2table(info, 'VariableNames', {'movie','channel','nChannels','nFrames','nTracks','nKept','dataset'}), fullfile(outDir, 'extraction_info_same.csv'));
writetable(cell2table(mp, 'VariableNames', {'stored','trace','movie','storedPoints','pointsMatched','pointsKept','frameSpan','note'}), fullfile(outDir, 'same_granules_map.csv'));
disp('SAME_DONE');
end
