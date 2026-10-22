function sim_high_sp(wins, SPs, outName)
% SIM_HIGH_SP  What raising SP costs, on simulations where the truth is known
% (RESULTS.md section 34). Two settings, each with bursts (generateSim type 3,
% lambda 2) and a burst-free twin (type 1, same noise):
%   MS-like:  noise 40 a.u. (~0.7% of the level, as MS section III), k0 = 60
%   PP7-like: noise 1100 a.u. (~18% of the 2000-10000 level, as the real PP7 spot
%             sum frame noise of 16-23%), k0 = 300 (lambda 2 -> mean burst ~10% of
%             the level, comparable to the detected in vivo burst sizes)
% 500 traces x 361 frames each, bleaching 0.0026 /min, seeded. Detector as the real
% data: analysis.m, window 13, SP 10-25, negative bursts.
% Recall: a true negative jump of the noise-free twin (same seed, noise 0) counts
% as found if it lies within a detected negative burst +- 6 frames.
% Writes analysis_outputs/detection_grid_highSP/sim_summary.csv.
% Optional (section 35): wins (default 13), SPs (default 10-25), outName, e.g.
% sim_high_sp([5 9 13 21 31], [12 15 18 21], 'sim_win_sp.csv').
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'detection_grid_highSP'); if ~exist(outDir, 'dir'), mkdir(outDir); end
sets = {'MS-like (noise 0.7%)', 40, 60; 'PP7-like (noise 18%)', 1100, 300};
if nargin < 1, wins = 13; end
if nargin < 2, SPs = [10 12 14 16 18 20 22 25]; end
if nargin < 3, outName = 'sim_summary.csv'; end
nT = 500; nF = 361; k = 0.0026; tVec = (0:360) * 10/60;
rows = cell(0, 12);
for j = 1:size(sets, 1)
    nz = sets{j,2}; k0 = sets{j,3}; seed = 3400 + j;
    rng(seed); [S, B] = generateSim(nT, nF, 3, k, nz, 2, k0);
    rng(seed); S0 = generateSim(nT, nF, 3, k, 0, 2, k0);
    rng(seed + 50); [F, BF] = generateSim(nT, nF, 1, k, nz, 2, k0);
    u = S0 ./ exp(-k * linspace(0, 60, nF)); du = diff(u, 1, 2);
    [ti, tj] = find(du < -1e-9); tsz = -du(du < -1e-9); tj = tj + 1;          % trace, frame and size of each true loss
    lvl = mean(S0, 2);
    hours = nT * nF * 10 / 3600;
    for w = wins
    for SP = SPs
        ws = warning('off', 'all');
        [~, ~, ~, ~, ~, ~, seg] = analysis(num2cell(S, 2), num2cell(B, 2), tVec, w, SP, 0, 0);
        [~, ~, ~, ~, ~, ~, segF] = analysis(num2cell(F, 2), num2cell(BF, 2), tVec, w, SP, 0, 0);
        warning(ws);
        neg = seg(seg(:,2) < 0, :); negF = segF(segF(:,2) < 0, :);
        found = false(numel(ti), 1);
        for q = 1:size(neg, 1)
            hit = ti == neg(q,1) & tj >= neg(q,3) - 6 & tj <= neg(q,4) + 6; found(hit) = true;
        end
        big = tsz >= 3 * k0;
        rows(end+1,:) = {sets{j,1}, w, SP, size(neg,1) / hours, size(negF,1) / hours, mean(found), mean(found(~big)), mean(found(big)), ...
            median(abs(neg(:,5))) / median(lvl) * 100, median(tsz) / median(lvl) * 100, numel(tsz) / hours, size(negF,1) / max(size(neg,1), 1)}; %#ok<AGROW>
        fprintf('%s window %d SP %d: real %d, burst-free %d, recall %.2f\n', sets{j,1}, w, SP, size(neg,1), size(negF,1), mean(found));
    end
    end
end
writetable(cell2table(rows, 'VariableNames', {'setting','window','SP','detected_per_granule_h','false_per_granule_h','recall','recall_1_2_units','recall_3plus_units', ...
    'median_detected_size_pct','median_true_size_pct','true_losses_per_granule_h','false_over_detected'}), fullfile(outDir, outName));
disp('SIM_HISP_DONE');
end
