function V = summarizeGranules(M, picks)
%   Collection of granules statistics
%   picks  Indices of granules to use. ; when the input is [], defaults to 1:M.nGranules

if isempty(picks)
    picks = 1:M.nGranules;
end

% --- Gather bursts ---
fracReal = vertcat(M.fracReal{picks});      % Granule bursts, as a fraction of the granule
startReal = vertcat(M.startReal{picks});    % First frame of burst
fracCopy = vertcat(M.fracCopy{picks});      % Synthetic copy granule bursts, as a fraction of the copy granule
startCopy = vertcat(M.startCopy{picks});    % First frame of burst
amps = vertcat(M.ampAU{picks});             % Granule bursts in intensity units (a.u.)

% --- brightness: median spot mean of each granule ---
levels = M.level(picks);                                     % Median brightness per granule
V.brightnessPooled = median(levels);                         % Median over all granules
movies = M.movie(picks); 
movieList = unique(movies);
for m = 1:numel(movieList)
    movieMedians(m) = median(levels(movies == movieList(m)));
end
V.brightness = median(movieMedians);                         % Median over movies

% --- frame-to-frame noise, % of the granule ---
V.frameNoise = 100 * median(M.noise(picks));

% --- burst size in intensity units (a.u.); the mean correlates with lambda*k0 ---
V.meanBurst = mean(amps);
V.medianBurst = median(amps);

% --- burst size as % of the granule: granules, copies, and the part beyond noise ---
V.sizeGranules = 100 * median(fracReal);
V.sizeCopies = 100 * median(fracCopy);
V.sizeBeyond = V.sizeGranules - V.sizeCopies;

% --- turnover beyond noise: % of the granule lost in bursts per hour, granules minus copies ---
% (the copies' total is divided by the number of copies per granule)
hours = sum(M.hours(picks)); 
hours1 = sum(M.hours1(picks)); 
hours2 = sum(M.hours2(picks));
V.turnover = 100 * (sum(fracReal) - sum(fracCopy) / M.nCopies) / hours;
V.turnover1 = 100 * (sum(fracReal(startReal <= 180)) - sum(fracCopy(startCopy <= 180)) / M.nCopies) / hours1; % 180 = 30 minutes (IN OUR SETUP)
V.turnover2 = 100 * (sum(fracReal(startReal > 180)) - sum(fracCopy(startCopy > 180)) / M.nCopies) / hours2;   % 180 = 30 minutes (IN OUR SETUP)
V.noiseShare = (sum(fracCopy) / M.nCopies) / sum(fracReal);

% --- net change beyond photobleaching (%/h), averaged over granules ---
V.net = mean(M.granuleNet(picks), 'omitnan');
V.net1 = mean(M.granuleNet1(picks), 'omitnan');
V.net2 = mean(M.granuleNet2(picks), 'omitnan');
end
