function T = gm_amplitudes_by_trace(signals, backgrounds, tVec, win, SP)
% GM_AMPLITUDES_BY_TRACE  Runs the repo's analysis.m one trace at a time, so
% every burst amplitude keeps the trace (granule) it came from. analysis.m
% treats traces independently, so pooling these gives exactly the amplitudes
% of one call on all traces. Needed for the granule bootstrap (R2 5) and used
% identically for simulated and experimental data.
% signals/backgrounds: cell arrays (experimental) or matrices (simulated,
% one trace per row). Returns struct array T(i).ampPos/ampNeg/ampVar
% (placeholders removed) and T(i).len.
if ~iscell(signals)
    signals = num2cell(signals, 2); backgrounds = num2cell(backgrounds, 2);
end
n = numel(signals);
T = repmat(struct('ampPos', [], 'ampNeg', [], 'ampVar', [], 'len', 0), n, 1);
ws = warning('off', 'all'); c = onCleanup(@() warning(ws));
for i = 1:n
    s = signals{i}; b = backgrounds{i};
    T(i).len = numel(s);
    if numel(s) < 180, continue, end              % analysis.m drops these too
    [~, ~, ap, an, av] = analysis({s(:)'}, {b(:)'}, tVec, win, SP, 0, 0);
    T(i).ampPos = ap(2:end); T(i).ampNeg = an(2:end); T(i).ampVar = av(2:end);
end
end
