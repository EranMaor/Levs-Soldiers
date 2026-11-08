function set_seed(s)
% SET_SEED  Seed the random number generator (works in MATLAB and Octave).
try
    rng(s);
catch
    rand('state', s); randn('state', s); %#ok<RAND>
end
end
