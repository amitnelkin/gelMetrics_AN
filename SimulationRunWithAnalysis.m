close all; clear; clc;

% Algorithm Parameters
movingAverageWindowSpan = 13;             % Moving average window span for signal denoising
stringencyParameter = 15;                 % Parameter for burst events detection


% Simulation parameters
simType = 3;               % simType 3 for bursts 
lambdaList = [1 2 4];      % mean number of units per burst - only relevant for simType=3
snrList = [0.5 5];         % size of one unit / frame noise
k0_sim = 60;               % intensity of one unit (a.u.)
numTracks = 200;           % granules per simulated experiment ("movie")
trackLength = 361;         % time points (60 min at 10 s)
expConst = 0.0026;         % photobleaching rate (1/min), median of the measured backgrounds - THIS SETTING IS SPECIFIC TO OUR SYSTEM SETUP AND SHOULD BE TUNED
nFits = 20;                % simulated experiments per condition
nCopies = 20;              % burst-free synthetic granule copy per granule (this is done for the first simulated "movie" only)
histSNR = 0.5;             % SNR shown in panels (d-f)

nLam = numel(lambdaList); 
nSNR = numel(snrList);
lamWS = zeros(nLam, nSNR, nFits);       % lambda from the WS fit
k0WS = zeros(nLam, nSNR, nFits);        % k0 from the WS fit
meanBurst =  zeros(nLam, nSNR, nFits);  % mean burst size from burstMeasures (a.u.)
M = cell(nLam, nSNR);                   % all measures of the first experiment
truth = cell(nLam, nSNR);               % true values of the first experiment
seeds = 1:nLam*nSNR*nFits;              
seedCounter = 1;

for si = 1:nSNR
    for li = 1:nLam
        noiseStd = k0_sim / snrList(si);
        for it = 1:nFits
            seed = seeds(seedCounter);                   % a fixed seed for every experiment
            seedCounter = seedCounter + 1;
            rng(seed);
            [signals, backgrounds] = generateSim(numTracks, trackLength, simType, expConst, noiseStd, lambdaList(li), k0_sim);
            spotMean = cell(numTracks, 1);
            background = cell(numTracks, 1);
            for i = 1:numTracks
                spotMean{i} = signals(i, :)';
                background{i} = backgrounds(i, :)';
            end

            % the copies only matter for the noise measures; the fit and the mean burst
            % size use the granule bursts only, so the later experiments
            % use 1 copy to save runtime
            if it == 1
                copies = nCopies; 
            else 
                copies = 1; 
            end

            % simulations have no movies (every granule is its own replicate) and no spot sum
            Mi = burstMeasures(spotMean, background, [], [], copies, seed, movingAverageWindowSpan, stringencyParameter, trackLength);
            amps = vertcat(Mi.ampAU{:});                            % negative bursts, a.u.
            [lamWS(li, si, it), k0WS(li, si, it)] = PoissonFit2WS(amps, 1:5, 1:1000, 10);  % Lambdas and k0 fitted by WS 
            meanBurst(li, si, it) = Mi.meanBurst;                   % Mean burst size (lambda * k0)
            if it == 1
                M{li, si} = Mi;
                truth{li, si} = trueValues(seed, numTracks, trackLength, expConst, lambdaList(li), k0_sim);
            end
        end
        fprintf('\nlambda %g, k0 %g, SNR %g (noise %.0f a.u.)\n', lambdaList(li), k0_sim, snrList(si), noiseStd);
        fprintf('  WS fit, lambda over %d experiments:   median %.1f, range %.1f - %.1f\n', nFits, ...
            median(lamWS(li, si, :)), min(lamWS(li, si, :)), max(lamWS(li, si, :)));
        fprintf('  WS fit, k0:                          median %.0f, range %.0f - %.0f\n', ...
            median(k0WS(li, si, :)), min(k0WS(li, si, :)), max(k0WS(li, si, :)));
        fprintf('  WS fit, lambda*k0:                   median %.0f, range %.0f - %.0f  (true %g)\n', ...
            median(lamWS(li, si, :) .* k0WS(li, si, :)), min(lamWS(li, si, :) .* k0WS(li, si, :)), ...
            max(lamWS(li, si, :) .* k0WS(li, si, :)), lambdaList(li) * k0_sim);
        fprintf('  mean burst size (burstMeasures):     median %.0f, range %.0f - %.0f\n', ...
            median(meanBurst(li, si, :)), min(meanBurst(li, si, :)), max(meanBurst(li, si, :)));

        Mc = M{li, si}; 
        T = truth{li, si};
        fprintf('  first experiment, %d copies per granule (true values in brackets):\n', nCopies);
        fprintf('    mean burst size (a.u.)               %.0f [%.0f, %.0f]   (true %.0f)\n', Mc.meanBurst, Mc.meanBurst_lo, Mc.meanBurst_hi, T.meanBurst);
        fprintf('    burst size, granules / copies (%%)    %.1f / %.1f   (true %.1f)\n', Mc.sizeGranules, Mc.sizeCopies, T.size);
        fprintf('    burst size beyond noise (%% points)   %.1f [%.1f, %.1f]\n', Mc.sizeBeyond, Mc.sizeBeyond_lo, Mc.sizeBeyond_hi);
        fprintf('    turnover beyond noise (%%/h)          %.1f [%.1f, %.1f]   (true %.1f)\n', Mc.turnover, Mc.turnover_lo, Mc.turnover_hi, T.turnover);
        fprintf('    net change beyond bleaching (%%/h)    %.1f [%.1f, %.1f]   (true %.1f)\n', Mc.net, Mc.net_lo, Mc.net_hi, T.net);
    end
end

%% Plotting
% Okabe-Ito pallette
% ORANGE = '#E69F00', GREEN = '#009E73', BLUE = '#0072B2', YELLOW = '#F0E442', RED = '#D55E00', PURPLE = '#CC79A7'
OKI = ["#E69F00","#009E73","#0072B2","#F0E442","#D55E00","#CC79A7"];
markers = 'osd^h';
snrColors = OKI(1:nSNR)';
f = figure();
t = tiledlayout(f, 3, 3, 'TileSpacing', 'compact', 'Padding', 'loose');

% (a) lambda from the WS fit, one dot per experiment
ax = nexttile(t); hold(ax, 'on'); box on;
p = gobjects(1, nSNR);

for si = 1:nSNR
    for li = 1:nLam
        x = lambdaList(li) + 0.15 * (si - (nSNR + 1) / 2) + linspace(-0.05, 0.05, nFits)';   % Fitted lambda + offset per noise level + offset per fitting experiment (so the dots do not overlap)
        p(si) = plot(ax, x, squeeze(lamWS(li, si, :)), 'o', 'Color', snrColors(si, :), 'MarkerSize', 2.5, 'DisplayName', sprintf('SNR = %g', snrList(si)));
    end
end

top = max(max(lambdaList) + 1, max(lamWS(:)) + 0.5);
plot(ax, [0 top], [0 top], 'k--', 'LineWidth', 0.8);
xlim(ax, [0 top]); 
ylim(ax, [0 top]);
xlabel(ax, '\lambda_{sim}'); 
ylabel(ax, '\lambda_{est}');
legend(ax, p, 'Location', 'best', 'Box', 'off', 'FontSize', 12);
text(-0.01, 1.15, '(a)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'bold');

% (b) k0 vs lambda from the same fits; dashed lines: lambda*k0 = the true value
ax = nexttile(t); hold(ax, 'on'); box on;
for li = 1:nLam
    lamLine = linspace(0.2, top, 100);
    plot(ax, lamLine, lambdaList(li) * k0_sim ./ lamLine, 'k--', 'LineWidth', 0.8);
    for si = 1:nSNR
        plot(ax, squeeze(lamWS(li, si, :)), squeeze(k0WS(li, si, :)), markers(li), 'Color', snrColors(si, :), 'MarkerSize', 3);
    end
end

xlim(ax, [0 top]); 
ylim(ax, [0 1.1 * max(k0WS(:))]);
xlabel(ax, '\lambda_{est}'); 
ylabel(ax, 'k_0_{est} (a.u.)');

for li = 1:nLam
    p(li) = plot(ax, NaN, NaN, markers(li), 'Color', 'k', 'MarkerSize', 3, 'DisplayName', sprintf('\\lambda_{sim} = %g', lambdaList(li)));
end
legend(ax, p, 'Location', 'best', 'Box', 'off', 'FontSize', 12);
text(-0.01, 1.15, '(b)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'bold');

% (c) lambda*k0: WS product and mean burst size (median and range over the experiments)
ax = nexttile(t); hold(ax, 'on'); box on;
trueProduct = lambdaList * k0_sim;
for si = 1:nSNR
    for li = 1:nLam
        ll = squeeze(lamWS(li, si, :) .* k0WS(li, si, :)); 
        mb = squeeze(meanBurst(li, si, :));
        x = trueProduct(li) + 8 * (si - (nSNR + 1) / 2);
        errorbar(ax, x - 3, median(ll), median(ll) - min(ll), max(ll) - median(ll), 'o', 'Color', snrColors(si, :), 'MarkerSize', 3, 'CapSize', 2);
        errorbar(ax, x + 3, median(mb), median(mb) - min(mb), max(mb) - median(mb), 's', 'Color', snrColors(si, :), 'MarkerFaceColor', snrColors(si, :), 'MarkerSize', 3, 'CapSize', 2);
    end
end
topTP = 1.6 * max(trueProduct);
plot(ax, [0 topTP], [0 topTP], 'k--', 'LineWidth', 0.8);
xlim(ax, [0 topTP]); 
ylim(ax, [0 topTP]);
e1 = plot(ax, NaN, NaN, 'ko', 'MarkerSize', 3, 'DisplayName','\lambda \cdot k_0');                       % Invisible plot for legend
e2 = plot(ax, NaN, NaN, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName','mean burst size'); % Invisible plot for legend
legend(ax, [e1 e2], 'Location', 'best', 'Box', 'off', 'FontSize', 12);
xlabel(ax, 'true \lambda \cdot k_0 (a.u.)'); 
ylabel(ax, 'estimated \lambda \cdot k_0 (a.u.)');
text(-0.01, 1.15, '(c)', 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'bold');

% (d-f) burst sizes of the granules and of their copies, at SNR=histSNR
histSNRi = find(snrList == histSNR, 1);
for li = 1:nLam
    ax = nexttile(t); box on;
    plotBurstHist(ax, M{li, histSNRi});
    title(ax, sprintf('\lambda_{sim}=%g, SNR=%g, N_{bursts} = %d', lambdaList(li), histSNR, numel(vertcat(M{li, histSNRi}.fracReal{:}))));
    if li == 1
        ylabel(ax, 'Frequency');
        legend(ax, {'granules', 'burst-free copies'}, 'Location', 'best', 'Box', 'off', 'FontSize', 12);
    end
    text(-0.01, 1.15, ['(' char('c' + li) ')'], 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'bold');
end

for li = nLam + 1:3
    nexttile(t); box on;
    axis off; 
end      % fewer than three lambdas: empty tiles

% (g-i) all measures of every condition, as ExperimentalAnalysis.m; black dots: true values
measures = {'meanBurst', 'turnover', 'net'};
labels = {'Mean burst size [\lambda \cdot k_0 (a.u.)]', 'Turnover beyond noise (%/h)', sprintf('Net change beyond bleaching (%%/h)')};
trueField = {'meanBurst', 'turnover', 'net'};
for j = 1:3
    ax = nexttile(t); hold(ax, 'on'); box on;
    xTicks = zeros(1, nLam); 
    
    for li = 1:nLam
        for si = 1:nSNR
            x = (li - 1) * (nSNR + 1) + si;                      % conditions grouped by lambda
            y = M{li, si}.(measures{j}); 
            lo = M{li, si}.([measures{j} '_lo']); 
            hi = M{li, si}.([measures{j} '_hi']);
            b(si) = bar(ax, x, y, 0.8, 'FaceColor', snrColors(si, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85, 'DisplayName', sprintf('SNR = %g', snrList(si)));
            errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
            p = plot(ax, x, truth{li, si}.(trueField{j}), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName', 'true value');
        end
        xTicks(li) = (li - 1) * (nSNR + 1) + (nSNR + 1) / 2;
    end

    yline(ax, 0, 'k', 'LineWidth', 0.5);
    tickText = cell(1, nLam);
    for li = 1:nLam
        tickText{li} = sprintf('\\lambda_{sim} = %g', lambdaList(li)); 
    end
    set(ax, 'XTick', xTicks, 'XTickLabel', tickText, 'XTickLabelRotation', 0);
    xlim(ax, [0.2, nLam * (nSNR + 1) - 0.2]);
    ylabel(ax, labels{j});
    text(-0.01, 1.15, ['(' char('g' + j) ')'], 'Units', 'normalized', 'FontSize', 13, 'FontWeight', 'bold');
    if j==1
       legend(ax, [b p], 'Location', 'best', 'Box', 'off', 'FontSize', 12); 
    end
end


%% ------------------------------------------------------------------ functions
function T = trueValues(seed, numTracks, trackLength, expConst, lambda, k0)
% The true values of one simulated experiment: the same granules (same seed) simulated
% without noise, before photobleaching. Burst size as burstMeasures.m measures it: the
% amplitude relative to the granule above its background, at the start of the burst.
rng(seed);
[signals, backgrounds] = generateSim(numTracks, trackLength, 3, expConst, 0, lambda, k0);
t = linspace(0, 60, trackLength);
amps = []; 
fractions = []; 
nets = zeros(numTracks, 1);
for i = 1:numTracks
    % net change as burstMeasures.m measures it: robust slope of ln(granule / background)
    ratio = (signals(i, :) - backgrounds(i, :)) ./ backgrounds(i, :);
    frames = find(ratio > 0); % Granule must be brighter than background
    v = log(ratio(frames));
    slopes = [];
    for j = 1:numel(frames) - 1
        slopes = [slopes, (v(j+1:end) - v(j)) ./ (frames(j+1:end) - frames(j))];
    end
    if numel(frames) < 0.9 * trackLength
        nets(i) = NaN;
    else
        nets(i) = median(slopes) * 360 * 100;            % %/h
    end
    % bursts
    u = signals(i, :) ./ exp(-expConst * t);             % granule without bleaching
    b = backgrounds(i, :) ./ exp(-expConst * t);         % its background without bleaching
    for j = 2:trackLength
        step = u(j) - u(j - 1);
        if step < -1e-9                                  % material leaves the granule (can't be zero due to rounding errors)
            amps(end+1, 1) = abs(step);
            fractions(end+1, 1) = abs(step) / (u(j - 1) - b(j - 1));
        end
    end
end
hours = numTracks * trackLength / 360;
T.meanBurst = mean(amps);                                % a.u.
T.size = 100 * median(fractions);                        % % of the granule
T.turnover = 100 * sum(fractions) / hours;               % %/h
T.net = mean(nets, 'omitnan');                                    % %/h 
end

function plotBurstHist(ax, M)
% Burst sizes (% of the granule) of the granules (red bars) and of their copies (black outline).
fracReal = 100 * vertcat(M.fracReal{:}); 
fracCopy = 100 * vertcat(M.fracCopy{:});
edges = 0:2.5:60;
histogram(ax, min(fracReal, 59.99), edges, 'Normalization', 'probability', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
hold(ax, 'on');
histogram(ax, min(fracCopy, 59.99), edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
xlim(ax, [0 60]); 
xlabel(ax, 'Burst size (% of granule)');
text(ax, 0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
    M.sizeGranules, M.sizeCopies, M.turnover, M.turnover_lo, M.turnover_hi), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 6.5);
end