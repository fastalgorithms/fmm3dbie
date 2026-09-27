function obj = ones(m, n)
%KERNEL3D.ONES   Construct a constant kernel.
%
%   K = KERNEL3D.ONES() constructs a 1 x 1 kernel with value 1, i.e.
%   K(x,y) = 1 for all targets x and sources y. Applied to a density this
%   gives the integral of the density: (K sigma)(x) = int sigma dS.
%
%   K = KERNEL3D.ONES(M) constructs an M x M constant block kernel with
%   K(x,y) = ones(M).
%
%   K = KERNEL3D.ONES(M, N) constructs an M x N constant block kernel with
%   K(x,y) = ones(M, N).
%
%   See also KERNEL3D, KERNEL3D.ZEROS.

if nargin < 1 || isempty(m)
    m = 1;
end
if nargin < 2 || isempty(n)
    n = m;
end
assert(isnumeric(m) && isscalar(m) && m == round(m) && m > 0, ...
    'KERNEL3D:ones', 'M must be a positive integer.');
assert(isnumeric(n) && isscalar(n) && n == round(n) && n > 0, ...
    'KERNEL3D:ones', 'N must be a positive integer.');

A = builtin('ones', m, n);
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

obj.eval = @(s,t) repmat(A, size(t.r(:,:),2), size(s.r(:,:),2));

obj.fmm = @fmm_;
    function varargout = fmm_(eps, s, t, sigma) %#ok<INUSL>
        if isstruct(t), nt = size(t.r(:,:),2); else, nt = size(t(:,:),2); end
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
