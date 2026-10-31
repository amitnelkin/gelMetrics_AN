function M = burstMeasures(spotMean, background, spotSum, replicate, nCopies, seed, window, sp, trackLength)
%   Noise-controlled burst measures of one condition
%
%   spotMean    Mean intensity of the granule pixels with camera offset removed (signals)
%   background  Mean intensity of the surrounding pixels cell or dilute phase (backgrounds)
%   spotSum     Integrated granule intensity above the background; when the input is [], defaults to (spotMean - background)
%   replicate   Replicate of every trace (e.g. its movie); when the input is [], every granule is its own replicate
%   nCopies     Burst-free copies per granule
%   seed        Random seed of the copies
%   window      Moving average window span for signal denoising
%   sp          Stringency parameter for burst events detection

% ---- settings
nBoot = 2000;                         % bootstrap samples for the confidence intervals
tVec = linspace(0, 60, trackLength);
showHistograms = false;               % true: plot the amplitude histograms of the granules and of the copies


% ----  burst-free copies ----
% The slow trend of every trace (a cubic fit), normalised to its mean; the median over
% all traces is the bleaching shape of the condition. Each copy = that shape scaled to
% the trace + white noise with the trace's own frame-to-frame noise.
rng(seed);
S = []; % spot mean: the average brightness of one granule pixel.
B = []; % background: the average brightness of one pixel of the cell (or dilute phase) around the granule.
X = []; % S + B (to be handed seamlessly through analysis.m)
G = []; % spot sum: all the granule's light above the background, added over its pixels.
M = struct("replicate", []);

if isempty(replicate)
    replicate = 1:numel(spotMean);
end

for k = 1:numel(spotMean)
    if length(spotMean{k}) < 180     % Discard signals with less than 180 frames (30 min at 10 s/frame) -> USER SHOULD CONSIDER ALTERING THIS SETTING
        continue
    end
    
    S{end+1} = double(spotMean{k}(:));                  % spot mean
    B{end+1} = double(background{k}(:));                % background

    if isempty(spotSum)
        G{end+1} = S{end} - B{end};                 % granule above background (from the spot mean)
    else
        G{end+1} = double(spotSum{k}(:));               % granule above background (spot sum)
    end

    % analysis.m expects "granule + background": it divides by the background to remove
    % photobleaching effect. So we add the background here, and analysis.m takes it away again.
    X{end+1} = G{end} + B{end};

    M.replicate(end+1, 1) = replicate(k);
end

S = S(:); 
B = B(:); 
G = G(:); 
X = X(:);

n = numel(X);                   % Number of traces larger with minimum 180 frames
if n < 10
    warning('burstMeasures:notEnoughTraces', 'Only %d traces of at least 180 frames.', n); 
end
M.nGranules = n; 
M.nCopies = nCopies;
L = cellfun(@numel,S);          % Number of frames for qualified movies
minMovLength = min(L); 
for j = 1:n
    x = (1:L(j))';              
    trend = polyval(polyfit(x, X{j}, 3), x);                              % Fit a cubic polynom to catch the overall trend (stemming from photobleaching)
    shapes(j, :) = trend(1:minMovLength)' / mean(trend(1:minMovLength));  % Divide by mean to convert to unitless shape
    d = diff(X{j});                                                       % Change between frames
    sigma(j) = 1.4826 * median(abs(d - median(d))) / sqrt(2);             % Noise of one frame - from the change between frames
end
shape = median(shapes, 1); % Median bleaching shape. Taking the median minimizes the effect of bursts.
q = 0;
for c = 1:nCopies
    for k = 1:n
        sh = [shape, repmat(shape(end), 1, L(k) - minMovLength)]';        % Extend the shape to the trace length. Longer traces repeat the last value.
        q = q + 1;
        C{q} = mean(X{k}) / mean(sh) * sh + sigma(k) * randn(L(k), 1);    % Scale to the granule average brightnes + noise
        BC{q} = B{k};                                                     % The copy uses real background to go through analysis.m seamlessly
        copyOwner(q) = k;                                                 % which granule this copy belongs to
    end
end

C = C(:); 
BC = BC(:);

% ---- burst detection----
[~,~,~,~,~,~, segAU] = analysis(S, B, tVec, window, sp, 0, 0);    % Intensity units (BURST = lambda * k0).

[~,~,ampPos,ampNeg,ampVar,~, segReal] = analysis(X, B, tVec, window, sp, 0, 0);  % Granules - using spotsum (follow the percentage out of total light no matter if the granule has changed dimensions).
if showHistograms 
    hists(ampPos, ampNeg, ampVar, 'Real Granules'); 
end

[~,~,ampPos,ampNeg,ampVar,~, segCopy] = analysis(C, BC, tVec, window, sp, 0, 0); % Copies - using spotsum.
if showHistograms 
    hists(ampPos, ampNeg, ampVar, 'Synthetic Granules'); 
end

% negative bursts; burst size as a fraction of the granule = amplitude / level before it
segAU = segAU(segAU(:, 2) < 0 & segAU(:, 6) > 0, :);
segReal = segReal(segReal(:, 2) < 0 & segReal(:, 6) > 0, :);
segCopy = segCopy(segCopy(:, 2) < 0 & segCopy(:, 6) > 0, :);

% the bursts of every granule, kept in one list per granule
copyGranule = copyOwner(segCopy(:, 1));
M.ampAU = cell(n, 1); 
M.fracReal = cell(n, 1);
M.startReal = cell(n, 1);
M.fracCopy = cell(n, 1); 
M.startCopy = cell(n, 1);
for k = 1:n
    M.ampAU{k} = abs(segAU(segAU(:, 1) == k, 5));                                            % a.u.
    M.fracReal{k} = abs(segReal(segReal(:, 1) == k, 5)) ./ segReal(segReal(:, 1) == k, 6);   % amplitude / level at first (fraction)
    M.startReal{k} = segReal(segReal(:, 1) == k, 3);                                         % first frame of burst
    M.fracCopy{k} = abs(segCopy(copyGranule == k, 5)) ./ segCopy(copyGranule == k, 6);       % amplitude / level at first (fraction)
    M.startCopy{k} = segCopy(copyGranule == k, 3);                                           % first frame of burst
end

% ---- per granule ----
M.level = zeros(n, 1); 
M.noise = zeros(n, 1);
M.hours = zeros(n, 1); 
M.hours1 = zeros(n, 1); 
M.hours2 = zeros(n, 1);
M.granuleNet = zeros(n, 1); 
M.granuleNet1 = zeros(n, 1); 
M.granuleNet2 = zeros(n, 1);
for k = 1:n
    M.level(k) = median(S{k});                                                  % brightness
    if isempty(spotSum)
        sig = S{k}; 
    else 
        sig = G{k}; 
    end
    d = diff(sig);
    M.noise(k) = 1.4826 * median(abs(d - median(d))) / sqrt(2) / median(sig);   % frame noise, fraction
    M.hours(k) = L(k) / 360;                                                    % time measured (in our experiments we took an image every 10s over 60 min)
    M.hours1(k) = min(L(k), 180) / 360;                                         % first 30 min
    M.hours2(k) = max(L(k) - 180, 0) / 360;                                     % second 30 min
    ratio = G{k} ./ B{k};                                                       % granule relative to its background
    M.granuleNet(k) = netChange(ratio, 360);
    M.granuleNet1(k) = netChange(ratio(1:min(180, end)), 360);
    if L(k) - 180 >= 60
        M.granuleNet2(k) = netChange(ratio(181:end), 360);
    else
        M.granuleNet2(k) = NaN;
    end
end

% ---- Output parameters and their 95% confidence intervals ----
V = summarizeGranules(M, 1:n);
names = fieldnames(V);
for j = 1:numel(names)
    M.(names{j}) = V.(names{j});
end

replicateList = unique(M.replicate); 
nReplicates = numel(replicateList);
M.nReplicates = nReplicates; 
M.nBursts = numel(vertcat(M.ampAU{:}));

for b = 1:nBoot
    % resample replicates, then granules within each chosen replicate
    picks = [];
    chosenReplicates = randi(nReplicates, nReplicates, 1);
    for m = 1:nReplicates
        granulesOfReplicate = find(M.replicate == replicateList(chosenReplicates(m)));
        chosen = granulesOfReplicate(randi(numel(granulesOfReplicate), numel(granulesOfReplicate), 1));
        picks = [picks; chosen]; %#ok<AGROW>
    end
    Vb = summarizeGranules(M, picks);
    for j = 1:numel(names)
        M.boot.(names{j})(b) = Vb.(names{j}); 
    end
end
for j = 1:numel(names)
    ci = prctile(M.boot.(names{j}), [2.5 97.5]); % 95 percent confidence intervals
    M.([names{j} '_lo']) = ci(1); 
    M.([names{j} '_hi']) = ci(2);
end

% does the condition change between the first and the second 30 min? (two-sided bootstrap p)
d = M.boot.turnover2 - M.boot.turnover1;
M.pHalvesTurnover = min(1, 2 * min(mean(d <= 0), mean(d >= 0)));
d = M.boot.net2 - M.boot.net1;
M.pHalvesNet = min(1, 2 * min(mean(d <= 0), mean(d >= 0)));
end

function [] = hists(ampPos,ampNeg,ampVar,title_str)
    % Remove the -Inf placeholder used to initialize the variable
    ampPos = ampPos(2:end);
    ampNeg = ampNeg(2:end);
    ampVar = ampVar(2:end);

    % Plot histograms and statistical analysis results
    figure;
    yyaxis left
    histogram(ampVar(abs(ampVar) < 250), 25);
    ylabel('Counts');
    yyaxis right
    var_gauss = pdf('Normal', -250:0.5:250, 0, std(ampVar(abs(ampVar) < 100)));
    plot(-250:0.5:250, var_gauss);
    ylabel('PDF');
    xlabel('Intensity [A.U]');
    title(sprintf('Quiescent Amplitudes - %s',title_str));


    % Plot histogram of all amplitudes with 3 Gaussians overlaid corresponding
    % to Negative, Zero, Positive
    figure; 
    histogram(ampPos, 'FaceColor', '#009E73', 'BinWidth', 5, 'Normalization', 'probability'); hold on;
    histogram(ampVar, 'FaceColor', '#0072B2', 'BinWidth', 5, 'Normalization', 'probability');
    histogram(ampNeg, 'FaceColor', '#E69F00', 'BinWidth', 5, 'Normalization', 'probability');
    xlim([-600, 600]); 
    ylabel('Frequency');
    title(sprintf('3 Amplitudes Phase - %s',title_str));
    xlabel('Amplitudes [A.U]');
end

function r = netChange(y, framesPerHour)
% Net change beyond photobleaching, %/h: the robust (Theil-Sen) slope of
% ln(granule / background), i.e. the median of the slopes between all pairs of frames.
% NaN if more than 10% of the frames are <= 0.
tmpCheck = isfinite(y) & y > 0; % Granule must be brighter than background
if mean(tmpCheck) < 0.9
    r = NaN; 
    return; 
end
t = find(tmpCheck); 
v = log(y(tmpCheck));
slopes = [];
for i = 1:numel(t) - 1
    slopes = [slopes; (v(i+1:end) - v(i)) ./ (t(i+1:end) - t(i))]; 
end
r = median(slopes) * framesPerHour * 100;
end
