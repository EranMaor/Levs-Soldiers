function varargout = gf_helpers(op, varargin)
% GF_HELPERS  Small vectorised GF(256) utilities (dispatch on first argument).
%   Z = gf_helpers('mul', A, B)      elementwise product with broadcasting
%   V = gf_helpers('xorcols', M)     XOR-reduce the columns of M (R x C -> R x 1)
G = gf256();
switch op
    case 'mul'
        A = varargin{1}; B = varargin{2};
        % NB: indexing a vector table with a vector index keeps the table's
        % orientation, so look up through (:) and restore the index shape.
        LA = reshape(G.LOG(A(:)+1), size(A));
        LB = reshape(G.LOG(B(:)+1), size(B));
        S  = LA + LB + 1;
        Z  = reshape(G.EXP(S(:)), size(S));
        Z = Z .* ((A ~= 0) & (B ~= 0));
        varargout{1} = Z;
    case 'xorcols'
        Mx = varargin{1};
        while size(Mx, 2) > 1
            c = size(Mx, 2);  h = floor(c/2);
            Mx = [bitxor(Mx(:, 1:h), Mx(:, h+1:2*h)), Mx(:, 2*h+1:end)];
        end
        varargout{1} = Mx;
    otherwise
        error('gf_helpers: unknown op');
end
end
