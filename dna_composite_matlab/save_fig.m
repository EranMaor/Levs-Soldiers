function save_fig(f, cfg, name)
% SAVE_FIG  Save figure as PNG and MATLAB .fig in the results folder.
if ~exist(cfg.cacheDir, 'dir'), mkdir(cfg.cacheDir); end
try, print(f, '-dpng', '-r150', fullfile(cfg.cacheDir, [name '.png'])); catch e, warning('could not save png: %s', e.message); end
try, savefig(f, fullfile(cfg.cacheDir, [name '.fig'])); catch, end
end
