%% PaperFigures_revised.m - Figures 6-9 of the revised manuscript (Fig. 5: SimulationRunWithAnalysis.m)
% Run one section at a time (Ctrl+Enter) or the whole script, from the repository
% folder. Each section starts from the repository data, runs the same code a user
% runs on their own data (burstMeasures.m, compareConditions.m) and draws the figure
% in the manuscript style. Figures are saved in figures/ (PNG, 300 dpi, and PDF).
%
% The noise copies use the seeds of the paper (one copy per granule), so the numbers
% are those of the manuscript. Set nCopies = 20 for values that do not depend on the
% random draw of the copies (RESULTS.md section 45).
% In vitro uses the spot sum, as in vivo (the user's decision, RESULTS.md section 53).
%
% Run time: about 1 min per figure. Figure 5 is made by SimulationRunWithAnalysis.m.

%% Settings
clearvars; close all; clc;
addpath('legacy');                % the original lambda/k0 code (only for reference)
if ~exist('figures', 'dir'), mkdir('figures'); end
makePlots = true;                 % false: compute and print the numbers only
nCopies = 1;                      % 1 = as in the paper; 20 = stable noise reference
% the random seeds of the copies used for the paper
seeds = containers.Map({'4PP7', '4PP7-FUS', '30PP7', '30PP7-FUS', '5Qb', '10Qb', 'vitro_1_1', 'vitro_10_1', 'vitro_100_1'}, ...
                       {18, 19, 20, 21, 22, 23, 18, 19, 20});
RED = [1 0.4 0.4]; BLUE = [0 0.447 0.741]; ORANGE = [0.835 0.369 0];
CPP = [0 0.447 0.698]; CFUS = [0 0.620 0.451];

%% Figure 5 - simulated granules
% Made by SimulationRunWithAnalysis.m (its default settings; set saveFigures = true to
% write figures/Fig5.png and .pdf). It runs the same analysis as the experimental data.

%% Figure 6 - in vitro granules
names = {'vitro_1_1', 'vitro_10_1', 'vitro_100_1'}; ratio = {'1:1', '10:1', '100:1'};
M6 = cell(1, 3);
for k = 1:3
    d = loadPaperDataset(names{k});
    M6{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 6 %s: median %.1f vs %.1f %%, mean burst %.0f [%.0f, %.0f], turnover %.1f [%.1f, %.1f], net %.1f [%.1f, %.1f]\n', ratio{k}, ...
        M6{k}.sizeGranules, M6{k}.sizeCopies, M6{k}.meanBurst, M6{k}.meanBurst_lo, M6{k}.meanBurst_hi, ...
        M6{k}.turnover, M6{k}.turnover_lo, M6{k}.turnover_hi, M6{k}.net, M6{k}.net_lo, M6{k}.net_hi);
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 3.7], 'Color', 'w');
    t = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:3
        ax = nexttile(t); plotBurstHist(ax, M6{k});
        title(ax, sprintf('protein:RNA %s  (n = %d)', ratio{k}, numel(vertcat(M6{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1
            lg = legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6);
            lg.ItemTokenSize = [10 6];   % short legend symbols, so the legend stays clear of the bars
        end
    end
    measures = {'meanBurst', 'turnover', 'net'};
    labels = {'Mean burst size, \lambda\cdotk_0 (a.u.)', 'Turnover beyond noise (%/h)', sprintf('Net change beyond\nbleaching (%%/h)')};
    for j = 1:3
        ax = nexttile(t); hold(ax, 'on');
        y = zeros(1, 3); lo = y; hi = y;
        for k = 1:3
            y(k) = M6{k}.(measures{j}); lo(k) = M6{k}.([measures{j} '_lo']); hi(k) = M6{k}.([measures{j} '_hi']);
            bar(ax, k, y(k), 0.6, 'FaceColor', RED, 'EdgeColor', 'none', 'FaceAlpha', 0.85);
        end
        errorbar(ax, 1:3, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
        yline(ax, 0, 'k', 'LineWidth', 0.5);
        set(ax, 'XTick', 1:3, 'XTickLabel', ratio); xlim(ax, [0.4 3.6]);
        xlabel(ax, 'protein:RNA'); ylabel(ax, labels{j});
        styleAxes(ax, char('d' + j - 1));
    end
    saveFig(f, 'Fig6');
end

%% Figure 7 - PP7 granules in vivo
names = {'4PP7', '30PP7'}; labels = {'PP7 + QCP-5x/PCP-4x', 'PP7 + PCP-30x'};
M7 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M7{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 7 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', labels{k}, ...
        M7{k}.sizeGranules, M7{k}.sizeCopies, M7{k}.turnover, M7{k}.turnover_lo, M7{k}.turnover_hi);
end
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');   % loose: the panel label (a) is not cut at the left
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M7{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M7{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    saveFig(f, 'Fig7');
end

%% Figure 8 - Qbeta granules in vivo (uses M7 from Figure 7 for the tests)
names = {'5Qb', '10Qb'}; labels = {['Q' char(946) ' + QCP-5x/PCP-4x'], ['Q' char(946) ' + QCP-10x']};
M8 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M8{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
    fprintf('Fig 8 %s: median %.1f vs %.1f %%, turnover beyond noise %.1f [%.1f, %.1f] %%/h\n', labels{k}, ...
        M8{k}.sizeGranules, M8{k}.sizeCopies, M8{k}.turnover, M8{k}.turnover_lo, M8{k}.turnover_hi);
end
compareConditions(M8{1}, M7{1}, 'Qbeta 5x', 'PP7 4x');          % Qbeta vs PP7 (text of Sec. V.B)
compareConditions(M8{2}, M7{2}, 'Qbeta 10x', 'PP7 30x');
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 2.3], 'Color', 'w');
    t = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');   % loose: the panel label (a) is not cut at the left
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M8{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M8{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    saveFig(f, 'Fig8');
end

%% Figure 9 - PP7-FUS granules in vivo (uses M7 from Figure 7 as the PP7 reference)
names = {'4PP7-FUS', '30PP7-FUS'}; labels = {'PP7-FUS + QCP-5x/PCP-4x', 'PP7-FUS + PCP-30x'};
M9 = cell(1, 2);
for k = 1:2
    d = loadPaperDataset(names{k});
    M9{k} = burstMeasures(d.spotMean, d.background, d.spotSum, d.movie, nCopies, seeds(names{k}));
end
T4 = compareConditions(M9{1}, M7{1}, 'PP7-FUS 4x', 'PP7 4x');
T30 = compareConditions(M9{2}, M7{2}, 'PP7-FUS 30x', 'PP7 30x');
pBright4 = T4.p(strcmp(T4.measure, 'brightnessPooled')); pBright30 = T30.p(strcmp(T30.measure, 'brightnessPooled'));
pBeyond4 = T4.p(strcmp(T4.measure, 'sizeBeyond')); pBeyond4one = T4.p_AgtB(strcmp(T4.measure, 'sizeBeyond'));
pBeyond30 = T30.p(strcmp(T30.measure, 'sizeBeyond'));
% (d) PP7-FUS minus PP7: total gap in median burst size, the part reproduced by the copies,
% and the part beyond noise; 95% CI from the two conditions' bootstrap samples
parts = {'sizeGranules', 'sizeCopies', 'sizeBeyond'};
gap = zeros(2, 3); gapLo = gap; gapHi = gap;
for i = 1:2
    for k = 1:3
        gap(i, k) = M9{i}.(parts{k}) - M7{i}.(parts{k});
        bootDiff = M9{i}.boot.(parts{k}) - M7{i}.boot.(parts{k});
        gapLo(i, k) = prctile(bootDiff, 2.5); gapHi(i, k) = prctile(bootDiff, 97.5);
    end
    fprintf('Fig 9(d) %s: total %.1f [%.1f, %.1f], noise %.1f [%.1f, %.1f], beyond %.1f [%.1f, %.1f]\n', labels{i}, [gap(i, :); gapLo(i, :); gapHi(i, :)]);
end
fprintf('Fig 9(c) brightness p = %.4f (4x), %.2f (30x); 9(d) beyond noise 4x: p = %.3f one-sided, %.3f two-sided; 30x: p = %.2f\n', ...
    pBright4, pBright30, pBeyond4one, pBeyond4, pBeyond30);
if makePlots
    f = figure('Units', 'inches', 'Position', [1 1 6.5 4.6], 'Color', 'w');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'loose');
    for k = 1:2
        ax = nexttile(t); plotBurstHist(ax, M9{k});
        title(ax, sprintf('%s  (n = %d)', labels{k}, numel(vertcat(M9{k}.fracReal{:}))));
        if k == 1, ylabel(ax, 'Frequency'); end
        styleAxes(ax, char('a' + k - 1));
        if k == 1, legend(ax, {'granules', 'burst-free copies'}, 'Location', 'east', 'Box', 'off', 'FontSize', 6.5); end
    end
    % (c) granule brightness (median over granules, 95% CI)
    ax = nexttile(t); hold(ax, 'on');
    Ms = {M7{1}, M9{1}, M7{2}, M9{2}}; x = [1 2 3.6 4.6]; colors = [CPP; CFUS; CPP; CFUS];
    y = zeros(1, 4); lo = y; hi = y;
    for k = 1:4
        y(k) = Ms{k}.brightnessPooled / 1000; lo(k) = Ms{k}.brightnessPooled_lo / 1000; hi(k) = Ms{k}.brightnessPooled_hi / 1000;
        bar(ax, x(k), y(k), 0.7, 'FaceColor', colors(k, :), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    errorbar(ax, x, y, y - lo, hi - y, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 3);
    set(ax, 'XTick', x, 'XTickLabel', {'PP7', 'PP7-FUS', 'PP7', 'PP7-FUS'}); ylabel(ax, 'Granule brightness (a.u./1000)');
    yt = max(hi) * 1.06; ylim(ax, [0 max(hi) * 1.3]);         % group names and p values above the bars
    text(ax, 1.5, yt, {'QCP-5x/PCP-4x', sprintf('p = %.3f', pBright4)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'FontSize', 7);
    text(ax, 4.1, yt, {'PCP-30x', sprintf('p = %.2f', pBright30)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    styleAxes(ax, 'c');
    % (d) burst size gap
    ax = nexttile(t); hold(ax, 'on');
    colors = [0.55 0.55 0.55; 0.85 0.85 0.85; CFUS]; hb = gobjects(1, 3);
    for i = 1:2
        for k = 1:3
            xx = (i - 1) * 1.1 + (k - 2) * 0.28;
            hb(k) = bar(ax, xx, gap(i, k), 0.26, 'FaceColor', colors(k, :), 'EdgeColor', 'k', 'LineWidth', 0.4);
            errorbar(ax, xx, gap(i, k), gap(i, k) - gapLo(i, k), gapHi(i, k) - gap(i, k), 'k', 'LineWidth', 0.8, 'CapSize', 3);
        end
    end
    text(ax, 0, 11.8, sprintf('p = %.3f (one-sided)\np = %.3f (two-sided)', pBeyond4one, pBeyond4), 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 6.5);
    text(ax, 1.1, 11.8, sprintf('p = %.2f', pBeyond30), 'HorizontalAlignment', 'center', 'FontSize', 6.5);
    yline(ax, 0, 'k', 'LineWidth', 0.5); ylim(ax, [-4 17]);
    set(ax, 'XTick', [0 1.1], 'XTickLabel', {'QCP-5x/PCP-4x', 'PCP-30x'}); xlim(ax, [-0.55 1.65]);
    ylabel(ax, sprintf('PP7-FUS - PP7, median burst\nsize (%% points)'));
    lg = legend(ax, hb, {'total', 'explained by noise', 'beyond noise'}, 'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 6);
    lg.ItemTokenSize = [10 6];
    styleAxes(ax, 'd');
    saveFig(f, 'Fig9');
end


%% ------------------------------------------------------------------ functions
function d = loadPaperDataset(name)
% One dataset of the paper: spot mean, background, spot sum and the movie of every trace.
if startsWith(name, 'vitro')
    folder = 'in_vitro_reextracted/same_granules_R2/';
    file = strrep(name, 'vitro', 'same');                  % e.g. vitro_1_1 -> same_1_1
    info = readtable([folder 'extraction_info_same.csv'], 'Delimiter', ',', 'TextType', 'string');
    info = info(info.dataset == file, :);
    d.spotMean = loadTraces([folder file '_data.mat']);
    d.background = loadTraces([folder file '_background.mat']);
    d.spotSum = loadTraces([folder file '_spotsum_data.mat']);
else
    folder = 'in_vivo_reextracted/after_R2_offset/';
    info = [readtable([folder 'extraction_info.csv'], 'TextType', 'string'); ...
            readtable([folder 'extraction_info_5Qb_10Qb.csv'], 'TextType', 'string')];
    info = info(info.dataset == name, :);
    d.spotMean = loadTraces([folder name '_data.mat']);
    d.background = loadTraces([folder name '_background.mat']);
    d.spotSum = loadTraces(['in_vivo_reextracted/after_R2_sums/' name '_spotsum_data.mat']);
end
d.movie = repelem((1:height(info))', info.nKept);         % the traces are stored movie by movie
end


function traces = loadTraces(file)
% The extraction stores a cell array of traces whose first cell is empty.
s = struct2cell(load(file)); traces = s{1}(:);
if isempty(traces{1}), traces = traces(2:end); end
end


function plotBurstHist(ax, M)
% Burst sizes (% of the granule) of the granules (red bars) and of their copies (black outline).
fracReal = 100 * vertcat(M.fracReal{:}); fracCopy = 100 * vertcat(M.fracCopy{:});
edges = 0:2.5:60;
histogram(ax, min(fracReal, 59.99), edges, 'Normalization', 'pdf', 'FaceColor', [1 0.4 0.4], 'EdgeColor', 'none', 'FaceAlpha', 1);
hold(ax, 'on');
histogram(ax, min(fracCopy, 59.99), edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1);
xlim(ax, [0 60]); xlabel(ax, 'Burst size (% of granule)');
text(ax, 0.97, 0.95, sprintf('median %.1f%% vs %.1f%%\nturnover beyond noise\n%.1f [%.1f, %.1f] %%/h', ...
    M.sizeGranules, M.sizeCopies, M.turnover, M.turnover_lo, M.turnover_hi), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 6.5);
end


function styleAxes(ax, panel)
% Manuscript style: boxed axes, ticks inward, Arial 8, bold panel label "(a)".
set(ax, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.6, 'Layer', 'top', ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'TitleFontWeight', 'normal');
ax.Title.FontSize = 7.5;
if isempty(ax.Title.String), title(ax, ' '); end            % leaves room for the panel label
text(ax, -0.2, 1.04, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
    'VerticalAlignment', 'bottom', 'FontName', 'Arial');
end


function saveFig(f, name)
% Saves a figure in figures/ as PNG (300 dpi) and PDF.
exportgraphics(f, fullfile('figures', [name '.png']), 'Resolution', 300);
exportgraphics(f, fullfile('figures', [name '.pdf']), 'ContentType', 'vector');
end
