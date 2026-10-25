% T2  (a) Noise floor: bursts detected in flat, burst-free signals
%         (original generateSim.m, type 1) for several stringency values.
%     (b) Bleaching asymmetry: symmetric bursts (50/50), original vs
%         corrected simulation, at low and realistic bleaching.
%
% Note: analysis.m's eventsPerSignal skips traces with zero events, so
% events per trace here is sum(eventsPerSignal) / nTracks.

clear; close all;
repoDir = '..';   % repo root when this script lives in checks/
addpath(repoDir); addpath(fullfile(repoDir, 'legacy'));
ws = warning('off', 'all');

nTracks = 300; win = 13; tVec = linspace(0, 60, 360);

% (a) Noise floor
rng(1);
[S, B] = generateSim(nTracks, 360, 1, 0.0005, 40, 2);   % original, type 1
fprintf('\n(a) Flat signals, no bursts (noise sd = 40):\n');
for SP = [10 12 14 15]
    [~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, SP, 0, 0);
    ap = ap(2:end); an = an(2:end);
    fprintf('  SP=%2d  events/trace=%.2f  pos=%d  neg=%d  median|amp|=%.1f\n', ...
        SP, sum(ev)/nTracks, numel(ap), numel(an), median(abs([ap; an])));
end

% Reference: real bursts, lambda = 2
rng(1);
[S, B] = generateSim(nTracks, 360, 3, 0.0005, 40, 2);
[~, ~, ap, an, ~, ev] = analysis(S, B, tVec, win, 15, 0, 0);
fprintf('  Reference, bursty lambda=2, SP=15: events/trace=%.2f  median pos amp=%.1f\n', ...
    sum(ev)/nTracks, median(ap(2:end)));

% (b) Bleaching asymmetry with symmetric bursts
fprintf('\n(b) Symmetric bursts (50/50), SP = 15:\n');
for fixBleach = [false true]
    for expConst = [0.0005 0.005]
        rng(2);
        [S, B] = gm_generateSim(nTracks, 3, expConst, 40, 2, fixBleach, 0.5);
        [~, ~, ap, an] = analysis(S, B, tVec, win, 15, 0, 0);
        ap = ap(2:end); an = an(2:end);
        fprintf('  fixedSim=%d expConst=%-6g (%.0f%% bleach/60 min): nPos=%d nNeg=%d meanPos=%.0f meanNeg=%.0f\n', ...
            fixBleach, expConst, 100*(1-exp(-expConst*60)), numel(ap), numel(an), mean(ap), mean(abs(an)));
    end
end
warning(ws);
