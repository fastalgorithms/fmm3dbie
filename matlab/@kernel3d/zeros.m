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
%     K.getquad                  - returns empty sparse matrix
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

obj.getquad = @(S,eps,varargin) build_empty_quad(m, n, S, varargin);

obj.get_overs_orders = @(S,t,eps) S.norders(:);

end

function spmat = build_empty_quad(m, n, S, args)
%BUILD_EMPTY_QUAD  Return an empty (zero) sparse quadrature correction.
if ~isempty(args) && ~isempty(args{1})
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
spmat = sparse(m*ntarg, n*S.npts);
end
