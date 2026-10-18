% Step 1b: run the ORIGINAL read_tables.m and Collect_intensity_data.m
% (git show f7eded0:<file>, unchanged) on 4PP7 movies 005 (1-channel) and
% 007 (2-channel), then record which TIFF page every trace point is read from.
S = 'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\zigzag_scratch';
addpath(fullfile(S,'shim'));       % uigetdir -> run folder
addpath(fullfile(S,'orig_code'));  % original generate_movArray.m
cd(fullfile(S,'run'));
read_tables;              % original, copied into run folder (it calls clear)
Collect_intensity_data;   % original, copied into run folder (clearvars)

S = 'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\zigzag_scratch';
R = fullfile(S,'run'); O = fullfile(S,'out'); if ~exist(O,'dir'), mkdir(O); end
D = load(fullfile(R,'data.mat')); B = load(fullfile(R,'background.mat'));
movies = {'005','007'};

% --- page statistics: which pages are red, which DIC
pg = [];
for m = 1:2
    fn = fullfile(R,[movies{m} '.tif']); info = imfinfo(fn);
    for k = [1:6, numel(info)-1, numel(info)]
        I = double(imread(fn,k));
        pg = [pg; m, numel(info), k, min(I(:)), median(I(:)), max(I(:))]; %#ok<AGROW>
    end
end
writematrix(pg, fullfile(O,'page_stats.csv'));   % movie, nPages, page, min, median, max

% --- per-point table: trace id, movie, Mosaic frame, TIFF page read, intensity
rows = []; tr = 0;
for m = 1:2
    xy = load(fullfile(R,[movies{m} '_xy.mat']));
    fn = fullfile(R,[movies{m} '.tif']); nP = numel(imfinfo(fn));
    for t = 1:numel(xy.x)
        f = xy.frame{t}; x = xy.x{t}; y = xy.y{t};
        f = f(2:end); x = x(2:end); y = y(2:end);
        if numel(x) >= 24
            tr = tr + 1;
            sig = D.intensity{tr+1}; bg = B.background{tr+1};
            for n = 1:numel(f)
                rows = [rows; tr, m, t, f(n), f(n), x(n), y(n), sig(n), bg(n), nP]; %#ok<AGROW>
            end
        end
    end
end
% columns: trace, movie(1=005,2=007), mosaicTraj, mosaicFrame, tiffPageRead, x, y, spot, cell, nPages
writematrix(rows, fullfile(O,'points.csv'));

% --- where is the spot? sample the 2-channel movie at each point, on
% page f+1 (Mosaic's plane if it walked raw planes) and on the red/DIC pages of
% time point f (pages 2f+1, 2f+2), as max - median in a 7x7 box.
fn = fullfile(R,'007.tif'); nP = numel(imfinfo(fn));
cache = containers.Map('KeyType','double','ValueType','any');
getp = @(k) double(imread(fn,k));
sel = rows(rows(:,2)==2,:);
pk = nan(size(sel,1),3);
for i = 1:size(sel,1)
    f = sel(i,4); x = round(sel(i,6)); y = round(sel(i,7));
    ks = [f+1, 2*f+1, 2*f+2];
    for j = 1:3
        k = ks(j); if k > nP, continue; end
        if ~isKey(cache,k), cache(k) = getp(k); end
        I = cache(k);
        r = max(y-3,1):min(y+3,size(I,1)); c = max(x-3,1):min(x+3,size(I,2));
        b = I(r,c); pk(i,j) = max(b(:)) - median(I(:));
    end
end
writematrix([sel(:,[1 4]) pk], fullfile(O,'peak_check_007.csv')); % trace, frame, page f+1, page 2f+1 (red t=f?), page 2f+2
disp('DIAG_DONE');
