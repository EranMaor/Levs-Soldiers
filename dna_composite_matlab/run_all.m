% RUN_ALL  Reproduce Figs. 8-15 of Kaneko & Shibata (IEEE Access 2026).
%   mode = 'quick' : minutes, rough/biased curves (smoke test of the whole chain)
%   mode = 'full'  : paper-scale sample sizes (hours/days; use a parallel pool)
% Results are cached in ./results so the run can be interrupted and resumed.
if ~exist('mode', 'var'), mode = 'quick'; end
cfg = get_cfg(mode);
set_seed(0);
fprintf('== Fig. 8  (mutual information, CH-1..5)\n');          fig8_mutual_information(cfg);
fprintf('== Fig. 9  (error types)\n');                           fig9_error_types(cfg);
fprintf('== Fig. 10 (weight vectors)\n');                        fig10_weight_vectors(cfg);
fprintf('== Figs. 11-12 (inner decoder time / memory)\n');       fig11_12_complexity(cfg);
fprintf('== Figs. 13-14 (BLER of the concatenated scheme)\n');   fig13_14_bler(cfg);
fprintf('== Fig. 15 (overall time / memory)\n');                 fig15_overall_complexity(cfg);
fprintf('Done. Figures and data are in %s\n', cfg.cacheDir);
