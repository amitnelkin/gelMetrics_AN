function window_duration_test(SP)
% WINDOW_DURATION_TEST  Which events does the moving-average window suppress?
% Reviewer point R2 2b. Synthetic traces (seeded), same detection as the runners
% (analysis.m) at stringency SP, moving-average windows 7..21 frames:
%  A. one negative event per trace: the level drops by `amp` linearly over
%     `dur` frames (dur 1..35), mid-trace; level 5000 + noise SD 40, background
%     1200 + noise SD 10; 150 traces per condition, 361 frames. Records the
%     detection probability (a negative burst overlapping the event) and the
%     measured / true amplitude.
%  B. two equal instantaneous drops separated by `gap` frames (2..30): number of
%     negative bursts detected around the pair (1 = merged, 2 = resolved).
% Writes analysis_outputs/robustness/window_duration.csv, window_spacing.csv.
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'analysis'));
outDir = fullfile(repoDir, 'analysis_outputs', 'robustness');
if ~exist(outDir, 'dir'), mkdir(outDir); end
tVec = (0:360) * 10/60; n = 150; L = 361; mid = 180;
wins = 7:2:21; rowsA = {}; rowsB = {};
for amp = [150 300]
    for dur = [1 3 5 9 13 17 25 35]
        rng(100 + amp + dur);
        lvl = 5000 * ones(n, L);
        ramp = min(max(((1:L) - mid) / dur, 0), 1);          % 0 -> 1 over dur frames
        S = lvl - amp * ramp + 40 * randn(n, L);
        B = 1200 + 10 * randn(n, L);
        for w = wins
            ws = warning('off', 'all');
            [~, ~, ~, ~, ~, ~, seg] = analysis(S, B, tVec, w, SP, 0, 0);
            warning(ws);
            hit = false(n,1); meas = nan(n,1);
            neg = seg(seg(:,2) < 0, :);
            for i = 1:n
                q = neg(neg(:,1) == i & neg(:,4) >= mid - w & neg(:,3) <= mid + dur + w, :);
                if ~isempty(q), hit(i) = true; meas(i) = max(abs(q(:,5))); end
            end
            % data units: (S/B - 1)*B(1) ~ S - B, so the true amplitude is amp
            rowsA(end+1,:) = {amp, dur, w, mean(hit), median(meas, 'omitnan') / amp}; %#ok<AGROW>
        end
    end
end
for gap = [2 4 6 8 10 15 20 30]
    rng(500 + gap);
    S = 5000 - 150 * ((1:L) >= mid) - 150 * ((1:L) >= mid + gap) + 40 * randn(n, L);
    B = 1200 + 10 * randn(n, L);
    for w = wins
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(S, B, tVec, w, SP, 0, 0);
        warning(ws);
        neg = seg(seg(:,2) < 0, :); nb = zeros(n,1);
        for i = 1:n
            nb(i) = sum(neg(:,1) == i & neg(:,4) >= mid - w & neg(:,3) <= mid + gap + w);
        end
        rowsB(end+1,:) = {gap, w, mean(nb == 0), mean(nb == 1), mean(nb >= 2)}; %#ok<AGROW>
    end
end
writetable(cell2table(rowsA, 'VariableNames', {'amp','durationFrames','window','pDetected','measuredOverTrue'}), ...
    fullfile(outDir, 'window_duration.csv'));
writetable(cell2table(rowsB, 'VariableNames', {'gapFrames','window','pNone','pMerged','pResolved'}), ...
    fullfile(outDir, 'window_spacing.csv'));
disp('WINDOW_TEST_DONE');
end
