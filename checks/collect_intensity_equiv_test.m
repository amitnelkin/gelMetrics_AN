% Check (RESULTS.md section 64): collect_intensity.m (simplified) and checks/collect_intensity_tests.m (default options) give identical output on synthetic movies.
% Run with the repository and checks/ on the path.
% synthetic movies -> old (collect_intensity_tests, default options) vs new collect_intensity
D = fullfile(tempdir, 'ci_equiv'); if exist(D, 'dir'), rmdir(D, 's'); end; mkdir(D);
rng(7);
for m = 1:2
    H = 60 + 20 * m; W = 90; nF = 50;                % non-square images
    tracks = {[20 30], [70 15], [3 40], [45 H-2]};    % one granule near the left edge, one near the bottom
    fname = fullfile(D, sprintf('mov%d.tif', m));
    x = {}; y = {}; frame = {};
    for k = 1:numel(tracks)
        if k == 2, fr = 5:nF-1; elseif k == 3, fr = [0:10, 13:30, 34:45]; else, fr = 0:nF-1; end   % gaps in track 3
        frame{k} = fr; x{k} = tracks{k}(1) - 1 + 0.3 * randn(1, numel(fr)); y{k} = tracks{k}(2) - 1 + 0.3 * randn(1, numel(fr));
    end
    frame{end+1} = 0:10; x{end+1} = 50 + zeros(1, 11); y{end+1} = 20 + zeros(1, 11);   % short track (< 24 points, dropped)
    for t = 1:nF
        [cc, rr] = meshgrid(1:W, 1:H);
        img = 200 + 1000 + 50 * randn(H, W);
        for k = 1:numel(tracks)
            j = find(frame{k} == t - 1, 1);
            if ~isempty(j), img = img + 3000 * exp(-((cc - x{k}(j) - 1).^2 + (rr - y{k}(j) - 1).^2) / 4); end
        end
        if t == 1, imwrite(uint16(img), fname); else, imwrite(uint16(img), fname, 'WriteMode', 'append'); end
    end
    save(fullfile(D, sprintf('mov%d_xy.mat', m)), 'x', 'y', 'frame');
end
for off = [200 0]
    [I1, B1, inf1, S1] = collect_intensity_tests(D, struct('offset', off));
    [I2, B2, inf2, S2] = collect_intensity(D, off);
    fprintf('offset %d: tracks %d / %d | intensity %d background %d spotSum %d spotArea %d info %d\n', off, numel(I1) - 1, numel(I2) - 1, ...
        isequaln(I1, I2), isequaln(B1, B2), isequaln(S1.spotSum, S2.spotSum), isequaln(S1.spotArea, S2.spotArea), isequaln(inf1, inf2));
end
fprintf('track lengths: %s\n', mat2str(cellfun(@numel, I2(2:end))));
