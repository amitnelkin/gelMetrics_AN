function gm_style(ax, panel)
% GM_STYLE  Manuscript axes style: boxed, ticks inward, Arial 8; optional bold
% panel label, e.g. gm_style(ax, 'a') writes "(a)" above the top-left corner.
% Titles are kept at 7.5 pt, normal weight; an axes without a title gets a blank one
% so that the layout leaves room for the panel label above it.
set(ax, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.6, 'Layer', 'top', ...
    'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'TitleFontWeight', 'normal');
ax.Title.FontSize = 7.5;
if nargin > 1 && ~isempty(panel)
    if isempty(ax.Title.String), title(ax, ' ', 'FontSize', 7.5); end
    text(ax, -0.2, 1.04, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontName', 'Arial');
end
end
