function gm_style(ax, panel)
% GM_STYLE  Manuscript axes style: boxed, ticks inward, Arial 8; optional bold
% panel label, e.g. gm_style(ax, 'a') writes "(a)" above the top-left corner.
set(ax, 'Box', 'on', 'TickDir', 'in', 'FontName', 'Arial', 'FontSize', 8, 'LineWidth', 0.6, 'Layer', 'top');
if nargin > 1 && ~isempty(panel)
    text(ax, -0.14, 1.03, ['(' panel ')'], 'Units', 'normalized', 'FontWeight', 'bold', 'FontSize', 8.5, ...
        'VerticalAlignment', 'bottom', 'FontName', 'Arial');
end
end
