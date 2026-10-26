function R = GelMetricsReport(D, opts)
% GELMETRICSREPORT  Noise-controlled burst measures for one dataset (one condition).
% This is the analysis used in the revised paper (Table I). It replaces the
% lambda/k0 fit of modifiedPoissonFit.m (now in legacy/): lambda and k0 cannot be
% estimated separately, only their product (the mean burst size) can.
%
%   R = GelMetricsReport(D)          D from gm_load_extraction.m or gm_paper_dataset.m
%   R = GelMetricsReport(D, opts)
%
% D fields
%   name        text label
%   spotMean    cell array, one trace per granule: mean intensity of the spot pixels
%               (camera offset already subtracted), as saved by Collect_intensity_data.m
%   background  cell array, same size: mean intensity of the surrounding cell /
%               dilute-phase pixels
%   spotSum     (optional) cell array: integrated granule intensity above the
%               background. If given, burst fractions, turnover, frame noise and net
%               change use it (as in vivo in the paper); otherwise spotMean - background
%               is used (as in vitro in the paper).
%   movie       (optional) movie (field of view) of every trace. Movies are the
%               independent replicates of the confidence intervals; without it every
%               granule is its own unit.
%
% opts fields (defaults = the settings of the paper)
%   window 13, SP 12, frameInterval 10 (s), minFrames 180, nBoot 2000, seed 1,
%   nCopies 1 (burst-free copies per granule; the copy-based measures are averaged
%   over them; more copies remove the random variation of a single copy), verbose true
%
% Measures (negative bursts, i.e. material leaving the granule)
%   granules, movies
%   brightness          median spot mean per granule, then median per movie, then
%                       median over movies (brightnessPooled: median over granules)
%   frameNoise          robust frame-to-frame noise, % of the granule level
%   meanBurst, medianBurst   burst amplitude in intensity units (spot mean); the mean
%                       corresponds to lambda*k0 of the Poisson model
%   sizeGranules, sizeCopies  median burst size as % of the granule, in the granules
%                       and in their burst-free copies (same trend and noise, no bursts)
%   sizeBeyond          sizeGranules - sizeCopies (% points)
%   turnover            % of the granule lost in bursts per hour, granules minus copies
%   noiseShare          share of the detected turnover reproduced by the copies
%   net                 net change beyond photobleaching, %/h: robust slope of
%                       ln(granule / background); < 0 = losing material
%   turnover1/2, net1/2 the same in the first and second 30 min; pHalvesTurnover,
%                       pHalvesNet test the change between halves (bootstrap)
% Every measure has a 95% CI (fields *_lo, *_hi): movies resampled, then granules
% within movies (nested bootstrap). Compare two conditions with GelMetricsCompare.m.

if nargin < 2, opts = struct(); end
def = struct('window', 13, 'SP', 12, 'frameInterval', 10, 'minFrames', 180, ...
             'nBoot', 2000, 'seed', 1, 'nCopies', 1, 'verbose', true);
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}), opts.(fn{k}) = def.(fn{k}); end
end
if ~isfield(D, 'name'), D.name = 'dataset'; end
S = D.spotMean(:); B = D.background(:);
hasSum = isfield(D, 'spotSum') && ~isempty(D.spotSum);
if isfield(D, 'movie') && ~isempty(D.movie), mov = D.movie(:); else, mov = (1:numel(S))'; end
keep = cellfun(@numel, S) >= opts.minFrames;
S = cellfun(@(x) double(x(:)), S(keep), 'UniformOutput', false);
B = cellfun(@(x) double(x(:)), B(keep), 'UniformOutput', false);
mov = mov(keep); n = numel(S);
if hasSum
    G = cellfun(@(x) double(x(:)), D.spotSum(keep), 'UniformOutput', false);
else
    G = cellfun(@(s, b) s - b, S, B, 'UniformOutput', false);
end
X = cellfun(@(g, b) g + b, G, B, 'UniformOutput', false);   % granule above background, plus background
L = cellfun(@numel, S);
fpm = 60 / opts.frameInterval;                  % frames per minute
tVec = (0:max(L)) / fpm;                        % minutes
perH = 100 * 60 * fpm;                          % per-frame rate -> %/h
CUT = 30 * fpm;                                 % 30 min

% burst-free copies: one bleaching shape per dataset (median of the per-trace
% cubic trends), scaled to each trace, plus white noise at the trace's own level
rng(opts.seed);
Lm = min(L); shp = zeros(n, Lm); sgm = zeros(n, 1);
for i = 1:n
    x = (1:L(i))'; trd = polyval(polyfit(x, X{i}, 3), x);
    shp(i, :) = trd(1:Lm)' / mean(trd(1:Lm));
    d = diff(X{i}); sgm(i) = 1.4826 * median(abs(d - median(d))) / sqrt(2);
end
shape = median(shp, 1); K = opts.nCopies; C = cell(n * K, 1); cOwner = zeros(n * K, 1); q = 0;
for k = 1:K
    for i = 1:n
        sh = [shape, repmat(shape(end), 1, max(L(i) - Lm, 0))]'; sh = sh(1:L(i));
        q = q + 1; C{q} = mean(X{i}) / mean(sh) * sh + sgm(i) * randn(L(i), 1); cOwner(q) = i;
    end
end

% burst detection (analysis.m, the published detector)
ws = warning('off', 'all');
[~, ~, ~, ~, ~, ~, segA] = analysis(S, B, tVec, opts.window, opts.SP, 0, 0);   % intensity units
[~, ~, ~, ~, ~, ~, segR] = analysis(X, B, tVec, opts.window, opts.SP, 0, 0);   % fraction of granule
[~, ~, ~, ~, ~, ~, segC] = analysis(C, B(cOwner), tVec, opts.window, opts.SP, 0, 0);   % copies
warning(ws);
segA = segA(segA(:, 2) < 0, :);
segR = segR(segR(:, 2) < 0 & segR(:, 6) > 0, :);
segC = segC(segC(:, 2) < 0 & segC(:, 6) > 0, :); segC(:, 1) = cOwner(segC(:, 1));

% per-granule items
it = struct('movie', num2cell(mov), 'L', num2cell(L), 'K', K);
for i = 1:n
    it(i).level = median(S{i});
    sig = G{i}; if ~hasSum, sig = S{i}; end
    d = diff(sig); it(i).noise = 1.4826 * median(abs(d - median(d))) / sqrt(2) / median(sig);
    it(i).amps = abs(segA(segA(:, 1) == i, 5));
    r = segR(segR(:, 1) == i, :); c = segC(segC(:, 1) == i, :);
    for h = 0:2
        if h == 0, sr = true(size(r, 1), 1); sc = true(size(c, 1), 1);
        elseif h == 1, sr = r(:, 3) <= CUT; sc = c(:, 3) <= CUT;
        else, sr = r(:, 3) > CUT; sc = c(:, 3) > CUT; end
        it(i).(sprintf('r%d', h)) = abs(r(sr, 5)) ./ r(sr, 6);
        it(i).(sprintf('c%d', h)) = abs(c(sc, 5)) ./ c(sc, 6);
    end
    it(i).h0 = L(i) / fpm / 60; it(i).h1 = min(L(i), CUT) / fpm / 60; it(i).h2 = max(L(i) - CUT, 0) / fpm / 60;
    y = G{i} ./ B{i};
    it(i).n0 = netRate(y, perH);
    it(i).n1 = netRate(y(1:min(CUT, end)), perH);
    if L(i) - CUT >= 60, it(i).n2 = netRate(y(CUT+1:end), perH); else, it(i).n2 = NaN; end
end

% observed values and nested bootstrap
obs = stats(it);
[umov, ~, gi] = unique(mov); nm = numel(umov);
groups = arrayfun(@(k) find(gi == k), 1:nm, 'UniformOutput', false);
keys = fieldnames(obs); bs = nan(opts.nBoot, numel(keys));
for b = 1:opts.nBoot
    pick = groups(randi(nm, nm, 1)); idx = [];
    for k = 1:nm, g = pick{k}; idx = [idx; g(randi(numel(g), numel(g), 1))]; end %#ok<AGROW>
    sb = stats(it(idx)); bs(b, :) = cellfun(@(k) sb.(k), keys);
end
R = struct('name', D.name, 'granules', n, 'movies', nm, 'bursts', size(segA, 1));
for k = 1:numel(keys)
    R.(keys{k}) = obs.(keys{k});
    q = prctile(bs(:, k), [2.5 97.5]); R.([keys{k} '_lo']) = q(1); R.([keys{k} '_hi']) = q(2);
end
dT = bs(:, strcmp(keys, 'turnover2')) - bs(:, strcmp(keys, 'turnover1'));
dN = bs(:, strcmp(keys, 'net2')) - bs(:, strcmp(keys, 'net1'));
R.pHalvesTurnover = min(1, 2 * min(mean(dT <= 0), mean(dT >= 0)));
R.pHalvesNet = min(1, 2 * min(mean(dN <= 0), mean(dN >= 0)));
R.items = it; R.opts = opts; R.boot = bs; R.bootKeys = keys;   % bootstrap draws (for CIs of differences)
if opts.verbose, printReport(R); end
end


function o = stats(it)
cat0 = @(f) vertcat(it.(f));
o.brightnessPooled = median([it.level]);
[um, ~, g] = unique([it.movie]);
o.brightness = median(arrayfun(@(k) median([it(g == k).level]), 1:numel(um)));
o.frameNoise = 100 * median([it.noise]);
a = cat0('amps'); o.medianBurst = med(a); o.meanBurst = mean(a);
r = cat0('r0'); c = cat0('c0');
o.sizeGranules = 100 * med(r); o.sizeCopies = 100 * med(c);
o.sizeBeyond = o.sizeGranules - o.sizeCopies;
hrs = [sum([it.h0]), sum([it.h1]), sum([it.h2])];
for h = 0:2
    rr = cat0(sprintf('r%d', h)); cc = cat0(sprintf('c%d', h));
    t = 100 * (sum(rr) - sum(cc) / it(1).K) / hrs(h + 1);
    nv = [it.(sprintf('n%d', h))]; nv = mean(nv(~isnan(nv)));
    if h == 0, o.turnover = t; o.net = nv; else, o.(sprintf('turnover%d', h)) = t; o.(sprintf('net%d', h)) = nv; end
end
o.noiseShare = sum(c) / it(1).K / sum(r);
end


function m = med(x)
if isempty(x), m = NaN; else, m = median(x); end
end


function r = netRate(y, perH)
% robust (Theil-Sen) slope of ln(y) per frame, in %/h; NaN if > 10% of frames <= 0
t = (1:numel(y))'; ok = isfinite(y) & y > 0;
if mean(ok) < 0.9, r = NaN; return; end
t = t(ok); v = log(y(ok)); n = numel(t);
[i, j] = find(triu(true(n), 1));
r = median((v(j) - v(i)) ./ (t(j) - t(i))) * perH;
end


function printReport(R)
f = @(k, d) sprintf(['%.' num2str(d) 'f [%.' num2str(d) 'f, %.' num2str(d) 'f]'], R.(k), R.([k '_lo']), R.([k '_hi']));
fprintf('\n%s: %d granules, %d movies, %d negative bursts\n', R.name, R.granules, R.movies, R.bursts);
fprintf('  brightness (median over movies)      %.0f   (over granules %.0f)\n', R.brightness, R.brightnessPooled);
fprintf('  frame noise (%% of granule)           %.1f\n', R.frameNoise);
fprintf('  mean burst size, lambda*k0 (a.u.)    %s\n', f('meanBurst', 0));
fprintf('  median burst size (a.u.)             %s\n', f('medianBurst', 0));
fprintf('  burst size, granules / copies (%%)    %.1f / %.1f\n', R.sizeGranules, R.sizeCopies);
fprintf('  burst size beyond noise (%% points)   %s\n', f('sizeBeyond', 1));
fprintf('  turnover beyond noise (%%/h)          %s   (noise share %.2f)\n', f('turnover', 1), R.noiseShare);
fprintf('  net change beyond bleaching (%%/h)    %s\n', f('net', 1));
fprintf('  first vs second 30 min: turnover %.1f -> %.1f (p = %.3f), net %.1f -> %.1f (p = %.3f)\n', ...
    R.turnover1, R.turnover2, R.pHalvesTurnover, R.net1, R.net2, R.pHalvesNet);
end
