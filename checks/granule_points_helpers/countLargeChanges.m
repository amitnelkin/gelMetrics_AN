function [nRise, nDrop] = countLargeChanges(x, threshold, maxFrames, unit)
% Rises and drops of at least `threshold` within `maxFrames` frames, counted without overlap:
% scanning forward, the first frame that reaches the threshold ends the event and the scan
% goes on from there. unit = []: the change relative to the starting value; otherwise the
% change in units of `unit` (the same amount of light for every granule).
nRise = 0;
nDrop = 0;
t = 1;
n = numel(x);
while t < n
    segment = x(t + 1:min(n, t + maxFrames));
    if isempty(unit)
        change = segment / x(t) - 1;
    else
        change = (segment - x(t)) / unit;
    end
    iRise = find(change >= threshold, 1);
    iDrop = find(change <= -threshold, 1);
    if isempty(iRise) && isempty(iDrop)
        t = t + 1;
    elseif isempty(iDrop) || (~isempty(iRise) && iRise < iDrop)
        nRise = nRise + 1;
        t = t + iRise;
    else
        nDrop = nDrop + 1;
        t = t + iDrop;
    end
end
end
