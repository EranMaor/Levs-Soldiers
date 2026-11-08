function cfg = get_cfg(mode)
% GET_CFG  Simulation settings.  mode = 'quick' (minutes, rough curves) or 'full'.
%   The paper repeats BLER runs until 100 block errors (Table 4) and does not state the
%   number of words used for the MI curves; 'full' is much closer to that but takes
%   hours-to-days on one core (use the Parallel Computing Toolbox, parfor is used).
if nargin < 1, mode = 'quick'; end
cfg.mode = mode;
cfg.cacheDir = 'results';
cfg.n = 128;                    % inner word length
cfg.N = 2048;  cfg.dv = 2;      % outer LDPC code over GF(2^8): length, column weight
cfg.ldpcSeed = 1;
cfg.maxRounds = 20;             % inner-outer iterations (Table 4)
cfg.maxIterOuter = 20;          % FFT-SPA iterations per outer decoding (Table 4)
switch mode
    case 'quick'
        cfg.mi.nWords  = 60;                 % NB: MI is biased upward by up to ~0.3 bit here
        cfg.mi.alphas  = 0:0.2:1.4;
        cfg.mi.rlist   = [30 60 100];
        cfg.mi.rsel    = [50 100];            % depths used in Figs. 9 and 10
        cfg.bler.Ms     = [640 768 896];
        cfg.bler.rlist  = [50 100];
        cfg.bler.alphas = 0.2:0.2:1.2;
        cfg.bler.targetErrors = 5;   cfg.bler.maxBlocks = 20;   cfg.bler.batch = 4;
        cfg.bler.r14list = [50 75 100];
        cfg.cx.rlist = [50 100 150];   cfg.cx.nWords = 3;   cfg.cx.nBlocks = 2;
    case 'smoke'                              % seconds-scale wiring test, numbers meaningless
        cfg.mi.nWords = 4;  cfg.mi.alphas = [0 1];  cfg.mi.rlist = 30;  cfg.mi.rsel = 50;
        cfg.bler.Ms = 640;  cfg.bler.rlist = 50;  cfg.bler.alphas = 0.3;
        cfg.bler.targetErrors = 1;  cfg.bler.maxBlocks = 1;  cfg.bler.batch = 1;  cfg.bler.r14list = 50;
        cfg.cx.rlist = 50;  cfg.cx.nWords = 1;  cfg.cx.nBlocks = 1;
    case 'full'
        cfg.mi.nWords  = 400;
        cfg.mi.alphas  = 0:0.1:1.4;
        cfg.mi.rlist   = 30:10:100;
        cfg.mi.rsel    = [50 100];
        cfg.bler.Ms     = [512 640 768 896];
        cfg.bler.rlist  = [50 100];
        cfg.bler.alphas = 0.2:0.05:1.4;
        cfg.bler.targetErrors = 100; cfg.bler.maxBlocks = 2000; cfg.bler.batch = 12;
        cfg.bler.r14list = 50:10:100;
        cfg.cx.rlist = 50:10:150;      cfg.cx.nWords = 10;  cfg.cx.nBlocks = 10;
    otherwise
        error('mode must be ''smoke'', ''quick'' or ''full''');
end
cfg.bler.stopBLER = 0.5;        % stop sweeping alpha once BLER exceeds this
if ~exist(cfg.cacheDir, 'dir'), mkdir(cfg.cacheDir); end
end
