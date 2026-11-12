function [C, BC, owner] = burstFreeCopies(X, B, nCopies)
% Burst-free copies of every trace, made as in burstMeasures.m: the median bleaching shape of
% the condition (cubic fit of each trace, divided by its mean) scaled to the trace, plus white
% noise at the trace's own frame-to-frame level. BC: the real background of each copy;
% owner: the granule each copy belongs to. Set the random seed (rng) before calling.
n = numel(X);
L = zeros(n, 1);
for k = 1:n
    L(k) = numel(X{k});
end
minLength = min(L);
shapes = zeros(n, minLength);
sigma = zeros(n, 1);
for k = 1:n
    x = (1:L(k))';
    trend = polyval(polyfit(x, X{k}, 3), x);
    shapes(k, :) = trend(1:minLength)' / mean(trend(1:minLength));
    d = diff(X{k});
    sigma(k) = 1.4826 * median(abs(d - median(d))) / sqrt(2);
end
shape = median(shapes, 1);
C = cell(n * nCopies, 1);
BC = cell(n * nCopies, 1);
owner = zeros(n * nCopies, 1);
q = 0;
for c = 1:nCopies
    for k = 1:n
        sh = [shape, repmat(shape(end), 1, L(k) - minLength)]';
        q = q + 1;
        C{q} = mean(X{k}) / mean(sh) * sh + sigma(k) * randn(L(k), 1);
        BC{q} = B{k};
        owner(q) = k;
    end
end
end


