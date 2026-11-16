V = load('C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN\figures_revision_local\stage\loc_ea_vitro_vars.mat');
disp(V.conditions);
for c = 1:3
M = V.M{c};
f = []; g = []; mv = [];
for k = 1:numel(M.fracReal)
    x = 100 * M.fracReal{k}(:);
    f = [f; x]; g = [g; k * ones(numel(x), 1)]; mv = [mv; M.replicate(k) * ones(numel(x), 1)];
end
fprintf('\n=== %s: %d bursts, %d granules, movies:', V.conditions{c}, numel(f), numel(M.fracReal)); fprintf(' %d', unique(M.replicate)); fprintf('\n');
edges = 0:2.5:60; f = min(f, 59.99);
for b = 1:numel(edges) - 1
    in = f >= edges(b) & f < edges(b + 1);
    if ~any(in), continue; end
    s = sprintf('%5.1f-%4.1f%% | %3d bursts |', edges(b), edges(b + 1), sum(in));
    for m = unique(mv(in))'
        sel = in & mv == m;
        s = [s sprintf(' movie %d: %d bursts / %d granules;', m, sum(sel), numel(unique(g(sel))))];
    end
    disp(s);
end
end
