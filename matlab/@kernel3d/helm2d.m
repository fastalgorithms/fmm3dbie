function obj = helm2d(type, zk, coefs)
%KERNEL3D.HELM2D   Construct a Helmholtz volume potential kernel in 2D.
%
%   KERNEL3D.HELM2D(type, zk) or KERNEL3D.HELM2D(type, zk, coefs), for
%   densities on a flat surfer in the z = 0 plane, where zk is the
%   (complex) wavenumber and type is:
%      's'  - single layer,   S(x,y) = i/4 H_0^(1)(zk|x-y|)
%      'sp' - sprime,         S'(x,y) = d/dn_x S(x,y)
%      'sg' - gradient,       [d/dx S; d/dy S], [2 1], rows interleaved
%      's2trans' - [2x1]      [coefs(1)*S; coefs(2)*S'], default [1;1]
%
% See also HELM2D.KERN, HELM2D.GET_QUADRATURE_CORRECTION

if ( nargin < 2 )
    error('KERNEL3D.HELM2D: requires type and zk arguments.');
end

obj           = kernel3d();
obj.name      = 'helm2d';
obj.opdims    = [1 1];
obj.zk        = zk;
obj.ifcomplex = 1;

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.kernel_order = -1;
        obj.eval         = @(s,t) helm2d.kern(zk, s, t, 's');
        obj.fmm          = @(eps,s,t,sigma) helm2d.fmm(eps, zk, s, t, 's', sigma);

    case {'sp', 'sprime'}
        obj.type         = 'sp';
        obj.kernel_order = -1;
        obj.eval         = @(s,t) helm2d.kern(zk, s, t, 'sprime');
        obj.fmm          = @(eps,s,t,sigma) helm2d.fmm(eps, zk, s, t, 'sprime', sigma);
        obj.targ_fields  = {'n'};

    case {'s2trans'}
        if nargin < 3 || isempty(coefs), coefs = [1; 1]; end
        coefs = coefs(:);
        obj      = kernel3d.interleave([coefs(1)*kernel3d.helm2d('s', zk); ...
                                        coefs(2)*kernel3d.helm2d('sp', zk)]);
        obj.name         = 'helm2d';
        obj.type         = 's2trans';
        obj.zk           = zk;
        obj.params.coefs = coefs;
        return

    case {'sg', 'sgrad'}
        obj.type         = 'sg';
        obj.opdims       = [2 1];
        obj.kernel_order = -1;
        obj.eval         = @(s,t) helm2d.kern(zk, s, t, 'sgrad');
        obj.fmm          = @(eps,s,t,sigma) sgrad_fmm(eps, zk, s, t, sigma);

    otherwise
        error('KERNEL3D.HELM2D: unknown Helmholtz kernel type ''%s''.', type);

end

obj.getquad = @(S,eps,varargin) rsc_to_sparse( ...
    helm2d.get_quadrature_correction(S, obj.type, zk, eps, varargin{:}), S, obj.opdims(1));
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

icheck = exist(['fmm2d.' mexext], 'file');
if ( icheck ~= 3 )
    obj.fmm = [];
end

end

function spmat = rsc_to_sparse(Q, S, nker)
%RSC_TO_SPARSE  Convert helm2d getquad RSC output to a sparse matrix.
% The gradient kernel stores two rows per target, Q.wnear is (2,nquad).
if ( nker == 1 )
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
else
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
        kernel3d.rsc_interleave_full(nker, 1));
end
end

function pot = sgrad_fmm(eps, zk, s, t, sigma)
%SGRAD_FMM  Gradient of the single layer, d/dx and d/dy interleaved.
[~, grad] = helm2d.fmm(eps, zk, s, t, 's', sigma);
pot = grad(:);
end
