addpath('C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN');
V = load('C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN\figures_revision_local\stage\loc_ea_vitro_vars.mat');
M = V.M{2}; n = numel(M.fracReal);
mx = zeros(n, 1); nb = zeros(n, 1);
for k = 1:n, f = M.fracReal{k}; nb(k) = numel(f); if ~isempty(f), mx(k) = max(f); end; end
large = find(mx > 0.30); small = find(nb > 0 & mx <= 0.30); none = find(nb == 0);
sets = {1:n, large, small, [small; none], find(M.replicate == 2), find(M.replicate ~= 2)};
names = {'all granules', 'large-burst granules (any burst > 30%)', 'small-burst granules (bursts, all <= 30%)', 'small-burst + no-burst granules', 'movie 2 only', 'movies 1, 3, 4'};
rng(1);
for s = 1:numel(sets)
    idx = sets{s}(:); W = summarizeGranules(M, idx);
    b = zeros(2000, 1);
    for i = 1:2000, p = idx(randi(numel(idx), numel(idx), 1)); Wb = summarizeGranules(M, p); b(i) = Wb.turnover; end
    ci = prctile(b, [2.5 97.5]);
    fr = vertcat(M.fracReal{idx}); fc = vertcat(M.fracCopy{idx});
    fprintf('%-45s n=%2d granules, movies %s | bursts %d, copy bursts/copy %.1f | median size %.1f vs %.1f | turnover beyond noise %.1f [%.1f, %.1f] %%/h | net %.1f\n', ...
        names{s}, numel(idx), mat2str(unique(M.replicate(idx))'), numel(fr), numel(fc) / M.nCopies, 100 * median(fr), 100 * median(fc), W.turnover, ci(1), ci(2), W.net);
end
fprintf('stored condition value: %.1f [%.1f, %.1f]\n', M.turnover, M.turnover_lo, M.turnover_hi);
