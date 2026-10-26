function V = summarizeGranules(M, picks)
% SUMMARIZEGRANULES  All the measures of one condition, computed from a list of
% granules. This is the one place where every measure is defined.
%   M      output of burstMeasures.m (per-granule data)
%   picks  which granules to use (indices into M; repeats allowed, as in a
%          bootstrap sample). Use picks = 1:M.nGranules for the condition itself.
% Every burst here is a negative burst (material leaving the granule).

% --- collect the bursts of the picked granules (each granule's bursts are in its own list) ---
fracReal = vertcat(M.fracReal{picks});      % granule bursts, as a fraction of the granule
startReal = vertcat(M.startReal{picks});    % their first frame
fracCopy = vertcat(M.fracCopy{picks});      % bursts of the burst-free copies, same units
startCopy = vertcat(M.startCopy{picks});
amps = vertcat(M.ampAU{picks});             % granule bursts in intensity units (a.u.)

% --- brightness: median spot mean of each granule ---
levels = M.level(picks);
V.brightnessPooled = median(levels);                        % median over granules
movies = M.movie(picks); movieList = unique(movies);
movieMedians = zeros(numel(movieList), 1);
for m = 1:numel(movieList)
    movieMedians(m) = median(levels(movies == movieList(m)));
end
V.brightness = median(movieMedians);                         % median over movies

% --- frame-to-frame noise, % of the granule ---
V.frameNoise = 100 * median(M.noise(picks));

% --- burst size in intensity units (a.u.); the mean corresponds to lambda*k0 ---
V.meanBurst = mean(amps);
V.medianBurst = median(amps);

% --- burst size as % of the granule: granules, copies, and the part beyond noise ---
V.sizeGranules = 100 * median(fracReal);
V.sizeCopies = 100 * median(fracCopy);
V.sizeBeyond = V.sizeGranules - V.sizeCopies;

% --- turnover beyond noise: % of the granule lost in bursts per hour, granules minus copies ---
% (the copies' total is divided by the number of copies per granule)
cut = 180;                                                   % frame at 30 min
hours = sum(M.hours(picks)); hours1 = sum(M.hours1(picks)); hours2 = sum(M.hours2(picks));
V.turnover = 100 * (sum(fracReal) - sum(fracCopy) / M.nCopies) / hours;
V.turnover1 = 100 * (sum(fracReal(startReal <= cut)) - sum(fracCopy(startCopy <= cut)) / M.nCopies) / hours1;
V.turnover2 = 100 * (sum(fracReal(startReal > cut)) - sum(fracCopy(startCopy > cut)) / M.nCopies) / hours2;
V.noiseShare = (sum(fracCopy) / M.nCopies) / sum(fracReal);

% --- net change beyond photobleaching (%/h), averaged over granules ---
V.net = mean(M.granuleNet(picks), 'omitnan');
V.net1 = mean(M.granuleNet1(picks), 'omitnan');
V.net2 = mean(M.granuleNet2(picks), 'omitnan');
end
