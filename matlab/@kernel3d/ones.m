function obj = ones(A)
%KERNEL3D.ONES   Construct a constant kernel.
%
%   K = KERNEL3D.ONES() constructs a [1 1] kernel with value 1.
%   K = KERNEL3D.ONES(A) constructs an m x n constant block kernel, A, an
%   m x n matrix.
%
%   The returned kernel has:
%     K.eval(srcinfo, targinfo)  - returns repmat(A, nt, ns)
%     K.fmm                      - sums the density and broadcasts A*sum
%     K.getquad                  - returns the correct quadratures for 
%                                  near-field points
%
%   See also KERNEL3D, KERNEL3D.ZEROS.

if nargin < 1 || isempty(A)
    A = 1;
end
assert(isnumeric(A) && ismatrix(A), 'KERNEL3D.ONES: A must be a numeric 2D matrix.');

m = size(A, 1);
n = size(A, 2);
opdims = [m n];

obj           = kernel3d();
obj.name      = 'ones';
obj.type      = 'one';
obj.opdims    = opdims;
obj.zk        = 0;
obj.ifcomplex = 0;
obj.kernel_order = -1;
obj.iszero    = false;
obj.src_fields  = {};
obj.targ_fields = {};

obj.eval = @(s,t) repmat(A, size(t.r,2), size(s.r,2));

obj.fmm = @fmm_;
    function varargout = fmm_(eps, s, t, sigma) %#ok<INUSL>
        if isstruct(t), nt = size(t.r,2); else, nt = size(t,2); end
        blk = A * sum(reshape(sigma, n, []), 2);
        pot = repmat(blk, nt, 1);
        if nargout > 0, varargout{1} = pot; end
        if nargout > 1
            error('KERNEL3D:ones:fmm', 'Too many output arguments.');
        end
    end

obj.getquad = @(S, eps, varargin) build_near_quad(A, m, n, S, varargin);

obj.get_overs_orders = @(S,t,eps) S.norders(:);

end

function spmat = build_near_quad(A, m, n, S, args)
% Near-field quadrature over the getnear(S,targinfo) patch neighborhood.
if ~isempty(args) && isstruct(args{1}) && isfield(args{1}, 'r')
    targinfo = args{1};
else
    targinfo = S;
end

rsc = getnear(S, targinfo);
row_ptr = rsc.row_ptr;
col_ind = rsc.col_ind;

ixyzs = S.ixyzs(:);
npols = ixyzs(2:end) - ixyzs(1:end-1);
nquad = sum(npols(col_ind));

wts = S.wts(:);
wnear = zeros(1, nquad);
istart = 1;
for p = 1:numel(col_ind)
    src_inds = ixyzs(col_ind(p)):(ixyzs(col_ind(p)+1)-1);
    nelem = numel(src_inds);
    wnear(istart:istart+nelem-1) = wts(src_inds);
    istart = istart + nelem;
end

ri = kernel3d.rsc_interleave_full(m, n);
spmat = conv_rsc_to_spmat(S, row_ptr, col_ind, kron(A(:), wnear), ri);
end
