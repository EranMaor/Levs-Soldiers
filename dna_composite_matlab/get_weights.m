function w = get_weights(name)
% GET_WEIGHTS  Weight vectors of Table 3 (exact fractions).
switch name
    case 'w41', w = [1/15 2/15 4/15 8/15];   % L=4, |Theta|=256, rate 8.00
    case 'w42', w = [1/8 1/8 3/8 3/8];       % L=4, |Theta|=100, rate 6.64
    case 'w21', w = [1/3 2/3];               % L=2, |Theta|=16,  rate 4.00
    otherwise, error('unknown weight vector %s', name);
end
end
