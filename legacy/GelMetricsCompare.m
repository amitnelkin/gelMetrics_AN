function T = GelMetricsCompare(RA, RB, nPerm, seed)
% GELMETRICSCOMPARE  Is condition A different from condition B?
%   T = GelMetricsCompare(RA, RB)      RA, RB from GelMetricsReport.m
%   T = GelMetricsCompare(RA, RB, nPerm, seed)   (defaults 10000, 1)
% For each measure: the difference A - B and a two-sided p value from shuffling
% whole movies between the two conditions (movies are the replicates). Datasets
% without movie ids are shuffled by granule.
% Returns a table (measure, A, B, difference, p, p_AgtB) and prints it; p_AgtB is the
% one-sided p for A > B (use it only when that direction was stated in advance).
if nargin < 3, nPerm = 10000; end
if nargin < 4, seed = 1; end
rng(seed);
assert(RA.opts.nCopies == RB.opts.nCopies, 'Use the same nCopies for both conditions.');
meas = {'sizeGranules', 'sizeCopies', 'sizeBeyond', 'turnover', 'net', 'medianBurst', 'meanBurst', 'brightness', 'brightnessPooled'};
itA = RA.items; itB = RB.items;
[~, ~, ga] = unique([itA.movie]); [~, ~, gb] = unique([itB.movie]);
units = [arrayfun(@(k) itA(ga == k), 1:max(ga), 'UniformOutput', false), ...
         arrayfun(@(k) itB(gb == k), 1:max(gb), 'UniformOutput', false)];
na = max(ga); nu = numel(units);
f = @(it) cellfun(@(m) measure(it, m), meas);
obs = f(itA) - f(itB); cnt = zeros(1, numel(meas)); cntG = cnt;
for p = 1:nPerm
    o = randperm(nu);
    d = f(vertcat(units{o(1:na)})) - f(vertcat(units{o(na+1:end)}));
    cnt = cnt + (abs(d) >= abs(obs) - 1e-12); cntG = cntG + (d >= obs - 1e-12);
end
pv = (cnt + 1) / (nPerm + 1); pg = (cntG + 1) / (nPerm + 1);
T = table(meas', cellfun(@(m) RA.(m), meas)', cellfun(@(m) RB.(m), meas)', obs', pv', pg', ...
    'VariableNames', {'measure', RA.name, RB.name, 'difference', 'p', 'p_AgtB'});
fprintf('\n%s vs %s (%d vs %d movies; p from %d movie shuffles)\n', RA.name, RB.name, na, nu - na, nPerm);
disp(T);
end


function v = measure(it, m)
cat0 = @(fld) vertcat(it.(fld));
switch m
    case 'sizeGranules', v = 100 * med(cat0('r0'));
    case 'sizeCopies', v = 100 * med(cat0('c0'));
    case 'brightnessPooled', v = median([it.level]);
    case 'sizeBeyond', v = 100 * (med(cat0('r0')) - med(cat0('c0')));
    case 'turnover', v = 100 * (sum(cat0('r0')) - sum(cat0('c0')) / it(1).K) / sum([it.h0]);
    case 'net', x = [it.n0]; v = mean(x(~isnan(x)));
    case 'medianBurst', v = med(cat0('amps'));
    case 'meanBurst', v = mean(cat0('amps'));
    case 'brightness'
        [um, ~, g] = unique([it.movie]);
        v = median(arrayfun(@(k) median([it(g == k).level]), 1:numel(um)));
end
end


function m = med(x)
if isempty(x), m = NaN; else, m = median(x); end
end
