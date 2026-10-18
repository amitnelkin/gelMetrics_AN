function src = gm_strip_plot_lines(src)
% GM_STRIP_PLOT_LINES  Removes single-line plotting statements from a
% runner script's text (figure rendering hangs in headless MATLAB on this
% machine). Only lines whose first token is a graphics call are removed;
% computation lines are untouched.
g = 'figure|histogram|ylabel|xlabel|title|xlim|ylim|hold|yyaxis|plot|boxplot|xticklabels|legend|heatmap|colorbar|colormap|caxis|subplot|scatter|text|box';
src = regexprep(src, ['^[ \t]*(' g ')\>[^\n]*$'], '', 'lineanchors');
end
