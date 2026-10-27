% Copy of the committed burstMeasures.m (trend from trenddecomp); test only (RESULTS.md section 60).
function M = burstMeasures_td(spotMean, background, spotSum, movie, nCopies, seed)
% BURSTMEASURES  Noise-controlled burst measures of one condition (MS Table I).
%   M = burstMeasures(spotMean, background, spotSum, movie, nCopies, seed)
%
% Inputs (cell arrays with one trace per granule, as saved by Collect_intensity_data.m):
%   spotMean    mean intensity of the granule pixels (camera offset already removed)
%   background  mean intensity of the surrounding pixels (cell or dilute phase)
%   spotSum     integrated granule intensity above the background; [] = use
%               spotMean - background instead
%   movie       movie (field of view) of every trace; [] = every granule is its own replicate
%   nCopies     burst-free copies per granule (the paper used 1; 20 gives a stable noise floor)
%   seed        random seed of the copies (the paper's seeds: 4PP7 18, 4PP7-FUS 19,
%               30PP7 20, 30PP7-FUS 21, 5Qb 22, 10Qb 23; in vitro 1:1 18, 10:1 19, 100:1 20)
%
% Output M: the per-granule data, every measure (see summarizeGranules.m) and its
% 95% confidence interval (fields name_lo, name_hi), and the bootstrap values (M.boot).
%
% Steps:
%   1. keep traces of at least 30 min
%   2. make burst-free copies of every trace: the same photobleaching trend and the
%      same noise level, but no bursts
%   3. detect bursts (analysis.m) in the granules and in their copies
%   4. per granule: brightness, noise, net change, hours measured
%   5. all measures, and their 95% CI from resampling movies, then granules within movies

% ---- settings (the values used in the paper) ----
window = 13;          % moving-average window (time points)
SP = 12;              % stringency parameter (time points)
frameInterval = 10;   % seconds between frames
minFrames = 180;      % 30 min
nBoot = 2000;         % bootstrap samples for the confidence intervals

framesPerHour = 3600 / frameInterval;
tVec = (0:400) * frameInterval / 60;           % time in minutes (long enough for any trace)

% ---- 1. keep traces of at least 30 min ----
nAll = numel(spotMean);
if isempty(movie), movie = 1:nAll; end
keep = [];
for i = 1:nAll
    if numel(spotMean{i}) >= minFrames, keep(end+1) = i; end %#ok<AGROW>
end
n = numel(keep);
S = cell(n, 1); B = cell(n, 1); G = cell(n, 1); X = cell(n, 1);
for k = 1:n
    i = keep(k);
    S{k} = double(spotMean{i}(:));                  % spot mean
    B{k} = double(background{i}(:));                % background
    if isempty(spotSum)
        G{k} = S{k} - B{k};                         % granule above background (from the spot mean)
    else
        G{k} = double(spotSum{i}(:));               % granule above background (spot sum)
    end
    % analysis.m expects "granule + background": it divides by the background to remove
    % photobleaching. So we add the background here, and analysis.m takes it away again.
    X{k} = G{k} + B{k};
end
M.movie = movie(keep); M.movie = M.movie(:);
M.nGranules = n; M.nCopies = nCopies;

% ---- 2. burst-free copies ----
% The slow trend of every trace (a cubic fit), normalised to its mean; the median over
% all traces is the bleaching shape of the condition. Each copy = that shape scaled to
% the trace + white noise with the trace's own frame-to-frame noise.
rng(seed);
L = zeros(n, 1);
for k = 1:n, L(k) = numel(X{k}); end
Lmin = min(L);
shapes = zeros(n, Lmin); sigma = zeros(n, 1);
for k = 1:n
    x = (1:L(k))';
    trend = trenddecomp(X{k});   % TEST ONLY: trenddecomp instead of the cubic
    shapes(k, :) = trend(1:Lmin)' / mean(trend(1:Lmin));
    d = diff(X{k});
    sigma(k) = 1.4826 * median(abs(d - median(d))) / sqrt(2);   % robust noise of one frame
end
shape = median(shapes, 1);
C = cell(n * nCopies, 1); BC = cell(n * nCopies, 1); copyOwner = zeros(n * nCopies, 1);
q = 0;
for c = 1:nCopies
    for k = 1:n
        sh = [shape, repmat(shape(end), 1, L(k) - Lmin)]';      % extend the shape to the trace length
        q = q + 1;
        C{q} = mean(X{k}) / mean(sh) * sh + sigma(k) * randn(L(k), 1);
        BC{q} = B{k};
        copyOwner(q) = k;                                       % which granule this copy belongs to
    end
end

% ---- 3. burst detection (the original analysis.m) ----
% each row of seg: [trace, sign, first frame, last frame, amplitude, level at start, level at end]
oldWarnings = warning('off', 'all');            % analysis.m warns about poorly fitted segments
[~, ~, ~, ~, ~, ~, segAU] = analysis(S, B, tVec, window, SP, 0, 0);    % intensity units
[~, ~, ~, ~, ~, ~, segReal] = analysis(X, B, tVec, window, SP, 0, 0);  % granules
[~, ~, ~, ~, ~, ~, segCopy] = analysis(C, BC, tVec, window, SP, 0, 0); % copies
warning(oldWarnings);
% negative bursts only; burst size as a fraction of the granule = amplitude / level before it
segAU = segAU(segAU(:, 2) < 0, :);
segReal = segReal(segReal(:, 2) < 0 & segReal(:, 6) > 0, :);
segCopy = segCopy(segCopy(:, 2) < 0 & segCopy(:, 6) > 0, :);
% the bursts of every granule, kept in one list per granule
M.ampAU = cell(n, 1); M.fracReal = cell(n, 1); M.startReal = cell(n, 1); M.fracCopy = cell(n, 1); M.startCopy = cell(n, 1);
copyGranule = copyOwner(segCopy(:, 1));
for k = 1:n
    M.ampAU{k} = abs(segAU(segAU(:, 1) == k, 5));                                   % a.u.
    M.fracReal{k} = abs(segReal(segReal(:, 1) == k, 5)) ./ segReal(segReal(:, 1) == k, 6);   % fraction of the granule
    M.startReal{k} = segReal(segReal(:, 1) == k, 3);                                % first frame of each burst
    M.fracCopy{k} = abs(segCopy(copyGranule == k, 5)) ./ segCopy(copyGranule == k, 6);       % the granule's copies
    M.startCopy{k} = segCopy(copyGranule == k, 3);
end

% ---- 4. per granule ----
M.level = zeros(n, 1); M.noise = zeros(n, 1);
M.hours = zeros(n, 1); M.hours1 = zeros(n, 1); M.hours2 = zeros(n, 1);
M.granuleNet = zeros(n, 1); M.granuleNet1 = zeros(n, 1); M.granuleNet2 = zeros(n, 1);
for k = 1:n
    M.level(k) = median(S{k});                                   % brightness
    if isempty(spotSum), sig = S{k}; else, sig = G{k}; end
    d = diff(sig);
    M.noise(k) = 1.4826 * median(abs(d - median(d))) / sqrt(2) / median(sig);   % frame noise, fraction
    M.hours(k) = L(k) / framesPerHour;                           % time measured
    M.hours1(k) = min(L(k), 180) / framesPerHour;                % first 30 min
    M.hours2(k) = max(L(k) - 180, 0) / framesPerHour;            % second 30 min
    ratio = G{k} ./ B{k};                                        % granule relative to its background
    M.granuleNet(k) = netChange(ratio, framesPerHour);
    M.granuleNet1(k) = netChange(ratio(1:min(180, end)), framesPerHour);
    if L(k) - 180 >= 60
        M.granuleNet2(k) = netChange(ratio(181:end), framesPerHour);
    else
        M.granuleNet2(k) = NaN;
    end
end

% ---- 5. measures and their 95% confidence intervals ----
V = summarizeGranules(M, 1:n);
names = fieldnames(V);
for j = 1:numel(names), M.(names{j}) = V.(names{j}); end
movieList = unique(M.movie); nMovies = numel(movieList);
M.nMovies = nMovies; M.nBursts = numel(vertcat(M.ampAU{:}));
for j = 1:numel(names), M.boot.(names{j}) = zeros(nBoot, 1); end
for b = 1:nBoot
    % resample movies, then granules within each chosen movie
    picks = [];
    chosenMovies = randi(nMovies, nMovies, 1);
    for m = 1:nMovies
        granulesOfMovie = find(M.movie == movieList(chosenMovies(m)));
        chosen = granulesOfMovie(randi(numel(granulesOfMovie), numel(granulesOfMovie), 1));
        picks = [picks; chosen]; %#ok<AGROW>
    end
    Vb = summarizeGranules(M, picks);
    for j = 1:numel(names), M.boot.(names{j})(b) = Vb.(names{j}); end
end
for j = 1:numel(names)
    ci = prctile(M.boot.(names{j}), [2.5 97.5]);
    M.([names{j} '_lo']) = ci(1); M.([names{j} '_hi']) = ci(2);
end
% does the condition change between the first and the second 30 min? (two-sided bootstrap p)
d = M.boot.turnover2 - M.boot.turnover1;
M.pHalvesTurnover = min(1, 2 * min(mean(d <= 0), mean(d >= 0)));
d = M.boot.net2 - M.boot.net1;
M.pHalvesNet = min(1, 2 * min(mean(d <= 0), mean(d >= 0)));
end


function r = netChange(y, framesPerHour)
% Net change beyond photobleaching, %/h: the robust (Theil-Sen) slope of
% ln(granule / background), i.e. the median of the slopes between all pairs of frames.
% NaN if more than 10% of the frames are <= 0.
ok = isfinite(y) & y > 0;
if mean(ok) < 0.9, r = NaN; return; end
t = find(ok); v = log(y(ok));
slopes = [];
for i = 1:numel(t) - 1
    slopes = [slopes; (v(i+1:end) - v(i)) ./ (t(i+1:end) - t(i))]; %#ok<AGROW>
end
r = median(slopes) * framesPerHour * 100;
end
