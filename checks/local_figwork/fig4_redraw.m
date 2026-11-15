% Fig. 4 redrawn in the manuscript layout with the corrected simulation
% (whole signal bleaches) and the measured bleaching constant 0.0026/min
% (RESULTS section 3, B3). Columns: bursty, constant, increasing.
% Top: one example trace (base, + exponential, + noise); bottom: amplitudes
% detected in 1000 traces of each version (window 13, SP 15), probability.
repoDir = 'C:\Users\roeeamit-lab\Desktop\AN\GelMetrics\gelMetrics_AN';
addpath(repoDir); warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale'); warning('off', 'MATLAB:polyfit:PolyNotUnique');
outDir = fullfile(repoDir, 'figures', 'ms_revision_lab');
nTr = 1000; T = 360; noiseSD = 40; expC = 0.0005; lam = 2; win = 13; SP = 15;
tVec = linspace(0, 60, T);
types = [3 1 2]; letters = {'a', 'b'; 'c', 'd'; 'e', 'f'};
col = [65 109 172; 213 162 57; 75 155 117] / 255;        % blue, orange, green as in the MS figure
ex = [0 0; expC 0; expC noiseSD];   % versions: base, + exponential, + noise
seedEx = [101 102 103];
fig = figure('Units', 'inches', 'Position', [1 1 9.5 5.7], 'Color', 'w');
W = 0.25; x0 = [0.08 0.40 0.725];
for k = 1:3
    % example trace: same draws, three versions
    s = zeros(3, T);
    for v = 1:3
        rng(seedEx(k)); s(v, :) = generateSim(1, T, types(k), ex(v, 1), ex(v, 2), lam);
    end
    a1 = axes('Position', [x0(k) 0.80 W 0.15]); hold on; box on;
    plot(tVec, s(1, :), 'Color', col(1, :), 'LineWidth', 1.5);
    plot(tVec, s(2, :), 'Color', col(2, :), 'LineWidth', 1.5);
    xlim([0 60]);
    text(0, 1.02, ['(' letters{k, 1} ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    a2 = axes('Position', [x0(k) 0.585 W 0.15]); hold on; box on;
    plot(tVec, s(3, :), 'Color', col(3, :), 'LineWidth', 1);
    xlim([0 60]); xlabel('Time (min)');
    yl = [min(s(:)) max(s(:))]; yl = yl + [-0.05 0.05] * diff(yl);
    ylim(a1, yl); ylim(a2, yl);
    % amplitudes in 1000 traces of each version
    a3 = axes('Position', [x0(k) 0.08 W 0.33]); hold on; box on;
    for v = 1:3
        rng(10 + k);
        [S, B] = generateSim(nTr, T, types(k), ex(v, 1), ex(v, 2), lam);
        [~, ~, ap, an, av] = analysis(S, B, tVec, win, SP, 0, 0);
        amp = [reshape(ap(2:end), [], 1); reshape(an(2:end), [], 1); reshape(av(2:end), [], 1)];
        histogram(amp, 'BinWidth', 5, 'Normalization', 'probability', 'FaceColor', col(v, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        fprintf('type %d version %d: %d pos, %d neg, %d unclassified\n', types(k), v, numel(ap) - 1, numel(an) - 1, numel(av) - 1);
    end
    xlim([-600 600]); xlabel('Amplitudes (a.u.)');
    if k == 1, ylabel('Frequency'); end
    text(0, 1.02, ['(' letters{k, 2} ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'FontSize', 10);
    if k == 1
        ya = axes('Position', [x0(k) 0.585 W 0.365], 'Visible', 'off'); ylabel(ya, 'Intensity (a.u.)', 'Visible', 'on', 'FontSize', 11);
        set(ya.YLabel, 'Units', 'normalized', 'Position', [-0.2 0.5 0]);
    end
    set([a1 a2 a3], 'FontSize', 10, 'LineWidth', 1, 'TickDir', 'in');
end
% legend row between traces and histograms
la = axes('Position', [0.2 0.475 0.6 0.03], 'Visible', 'off'); hold(la, 'on');
h = gobjects(3, 1); nm = {'Base signal', '+Exponential component', '+Gaussian noise'};
for v = 1:3, h(v) = plot(la, NaN, NaN, 'Color', col(v, :), 'LineWidth', 2); end
legend(h, nm, 'Orientation', 'horizontal', 'Location', 'north', 'Box', 'off', 'FontSize', 11);
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
exportgraphics(fig, fullfile(outDir, 'Fig4.png'), 'Resolution', 300);
exportgraphics(fig, fullfile(outDir, 'Fig4.pdf'), 'ContentType', 'vector');
print(fig, fullfile(outDir, 'Fig4.svg'), '-dsvg');
