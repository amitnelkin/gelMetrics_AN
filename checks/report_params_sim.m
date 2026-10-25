function report_params_sim()
% REPORT_PARAMS_SIM  The proposed report parameters on simulated datasets, with
% their true values (RESULTS.md section 24). Simulations as Fig. 13 of the
% manuscript (generateSim type 3, k0 = 60): lambda 1, 2, 4 at noise 40
% (SNR k0/noise = 1.5), lambda 2 at SNR 5 (noise 12) and 0.5 (noise 120), plus a
% burst-free control (type 1, noise 40). 200 tracks x 361 frames each, seeded.
% expConst 0.0026 /min as SimulationRunWithAnalysis.m.
% Per trace, exactly as for the real data (detection_grid.m / noise_methods_real.m):
%   negative bursts at window 13, SP 12 on the trace and on its burst-free
%   surrogate (shared bleaching shape + white noise at the trace's own level):
%   count and summed fraction shed (|amp| / level before);
%   structure function D(1), D(30) of y = signal / cubicTrend - 1;
%   net change slope of ln(signal / background) (Theil-Sen in python).
% True values from the noise-free twin: generateSim draws the same random
% numbers with noiseStd 0 (noise enters only as noiseStd*randn), so the same
% seed gives the same bursts without noise. True fraction shed = sum over
% negative jumps of |jump| / level before (bleaching divided out); true D from
% the noise-free signal.
% Writes analysis_outputs/report_params/sim_traces.csv (one row per trace) and
% sim_signals.mat (signals for the net-change slopes).
repoDir = fileparts(fileparts(mfilename('fullpath')));
addpath(repoDir); addpath(fullfile(repoDir, 'legacy')); addpath(fullfile(repoDir, 'checks'));
outDir = fullfile(repoDir, 'analysis_outputs', 'report_params'); if ~exist(outDir, 'dir'), mkdir(outDir); end
sims = {'no bursts', 1, 0, 40; 'lambda 1', 3, 1, 40; 'lambda 2', 3, 2, 40; 'lambda 4', 3, 4, 40; ...
        'lambda 2, SNR 5', 3, 2, 12; 'lambda 2, SNR 0.5', 3, 2, 120};
nT = 200; nF = 361; k = 0.0026; tVec = (0:360) * 10/60; rows = cell(0, 13); SIG = struct();
for j = 1:size(sims, 1)
    rng(500 + j); [S, B] = generateSim(nT, nF, sims{j,2}, k, sims{j,4}, sims{j,3}, 60);
    rng(500 + j); [S0, ~] = generateSim(nT, nF, sims{j,2}, k, 0, sims{j,3}, 60);
    rng(600 + j);
    Sc = num2cell(S, 2); Bc = num2cell(B, 2);
    Lm = nF; shp = zeros(nT, Lm);
    for i = 1:nT
        s = S(i,:)'; x = (1:nF)'; trd = polyval(polyfit(x, s, 3), x); shp(i,:) = trd' / mean(trd);
    end
    shape = median(shp, 1)'; Ssur = Sc;
    for i = 1:nT
        s = S(i,:)'; sgm = 1.4826 * median(abs(diff(s) - median(diff(s)))) / sqrt(2);
        Ssur{i} = (mean(s) / mean(shape) * shape + sgm * randn(nF, 1))';
    end
    res = zeros(nT, 4);
    for kind = 1:2
        Z = Sc; if kind == 2, Z = Ssur; end
        ws = warning('off', 'all'); [~, ~, ~, ~, ~, ~, seg] = analysis(Z, Bc, tVec, 13, 12, 0, 0); warning(ws);
        seg = seg(seg(:,2) < 0 & seg(:,6) > 0, :);
        res(:, (kind - 1) * 2 + 1) = accumarray(seg(:,1), 1, [nT 1]);
        res(:, (kind - 1) * 2 + 2) = accumarray(seg(:,1), abs(seg(:,5)) ./ seg(:,6), [nT 1]);
    end
    tt = 1:nF; ble = exp(-k * linspace(0, 60, nF));
    for i = 1:nT
        [D1, D30] = sfun(S(i,:)); [T1, T30] = sfun(S0(i,:));
        u = S0(i,:) ./ ble; du = diff(u); neg = du < -1e-9;
        trueFrac = sum(-du(neg) ./ u([neg false]));
        rows(end+1,:) = {sims{j,1}, sims{j,3}, sims{j,4}, i, nF, res(i,1), res(i,2), res(i,3), res(i,4), D1, D30, trueFrac, T30 - T1}; %#ok<AGROW>
    end
    SIG.(sprintf('s%d', j)) = struct('name', sims{j,1}, 'S', S, 'B', B, 'S0', S0);
    fprintf('%s done\n', sims{j,1});
end
writetable(cell2table(rows, 'VariableNames', {'dataset','lambda','noise','trace','lengthFrames','nReal','fracReal','nSurr','fracSurr', ...
    'D1','D30','trueFrac','trueD30minusD1'}), fullfile(outDir, 'sim_traces.csv'));
save(fullfile(outDir, 'sim_signals.mat'), '-struct', 'SIG');
disp('SIM_PARAMS_DONE');
end

function [D1, D30] = sfun(g)
t = 1:numel(g); trd = polyval(polyfit(t, g, 3), t); y = g ./ trd - 1;
D1 = mean((y(2:end) - y(1:end-1)).^2); D30 = mean((y(31:end) - y(1:end-30)).^2);
end
