function obj = zeros(opdims)
%KERNEL3D.ZEROS   Construct a zero kernel.
%
%   K = KERNEL3D.ZEROS() constructs a [1 1] zero kernel.
%   K = KERNEL3D.ZEROS(opdims) constructs a zero kernel with operator
%   dimensions opdims = [m, n].
%
%   The returned kernel has:
%     K.eval(srcinfo, targinfo)  - returns zeros(m*nt, n*ns)
%     K.fmm                      - returns zeros(m*nt, 1)
%     K.getquad                  - returns empty quadrature struct
%     K.layer_eval               - returns zeros(m*nt, 1)
%     K.iszero                   - true
%
%   This is useful as a placeholder in block-kernel matrices.

if nargin < 1
    opdims = [1 1];
end
opdims = opdims(:).';
assert(numel(opdims) == 2, 'KERNEL3D.ZEROS: opdims must have exactly 2 entries.');

obj           = kernel3d();
obj.name      = 'zero';
obj.type      = 'zero';
obj.opdims    = opdims;
obj.zk        = 0;
obj.ifcomplex = 0;
obj.iszero    = true;
obj.kernel_order = -1;
obj.src_fields = {};
obj.targ_fields = {};

m = opdims(1);
n = opdims(2);

obj.eval = @(s,t) zeros(m * size(t.r,2), n * size(s.r,2));

obj.fmm = @(eps,src,targ,sigma) zeros(m * size(targ.r,2), 1);

obj.layer_eval = @(S,sigma,targ,eps,varargin) zeros(m * size(targ.r,2), 1);

obj.getquad = @(S,eps,varargin) build_empty_quad(S, varargin);

obj.get_overs_orders = @(S,t,eps) S.norders(:);

% Zero kernels produce scalar (nker=1) zero wnear regardless of opdims.
obj.rsc_to_interleave = kernel3d.rsc_interleave_scalar();

end

function Q = build_empty_quad(S, args)
%BUILD_EMPTY_QUAD  Return an empty (zero) quadrature correction struct.
if length(args) >= 1 && ~isempty(args{1})
    targ = args{1};
    if isa(targ, 'surfer')
        ntarg = targ.npts;
    elseif isstruct(targ) && isfield(targ, 'r')
        ntarg = size(targ.r, 2);
    else
        ntarg = S.npts;
    end
else
    ntarg = S.npts;
end

Q.row_ptr  = ones(ntarg+1, 1);
Q.col_ind  = zeros(0, 1);
Q.iquad    = ones(1, 1);
Q.wnear    = zeros(1, 0);
Q.nquad    = 0;
Q.nnz      = 0;
Q.format   = 'rsc';
Q.ifcomplex = 0;
Q.kernel_order = -1;
end
