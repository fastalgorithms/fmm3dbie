function obj = lap2d(type)
%KERNEL3D.LAP2D   Construct a Laplace volume potential kernel in 2D.
%
%   KERNEL3D.LAP2D(type), for densities on a flat surfer in the z = 0
%   plane, where type is:
%      's'  - single layer,   S(x,y) = -log|x-y|/(2*pi)
%      'sp' - sprime,         S'(x,y) = d/dn_x S(x,y)
%      'sg' - gradient,       [d/dx S; d/dy S], [2 1], rows interleaved
%
% See also LAP2D.KERN, LAP2D.GET_QUADRATURE_CORRECTION

if ( nargin < 1 )
    error('KERNEL3D.LAP2D: missing Laplace kernel type.');
end

obj           = kernel3d();
obj.name      = 'lap2d';
obj.opdims    = [1 1];
obj.zk        = 0;
obj.ifcomplex = 0;

switch lower(type)

    case {'s', 'single'}
        obj.type         = 's';
        obj.kernel_order = -1;
        obj.eval         = @(s,t) lap2d.kern(s, t, 's');
        obj.fmm          = @(eps,s,t,sigma) lap2d.fmm(eps, s, t, 's', sigma);

    case {'sp', 'sprime'}
        obj.type         = 'sp';
        obj.kernel_order = -1;
        obj.eval         = @(s,t) lap2d.kern(s, t, 'sprime');
        obj.fmm          = @(eps,s,t,sigma) lap2d.fmm(eps, s, t, 'sprime', sigma);
        obj.targ_fields  = {'n'};

    case {'sg', 'sgrad'}
        obj.type         = 'sg';
        obj.opdims       = [2 1];
        obj.kernel_order = -1;
        obj.eval         = @(s,t) lap2d.kern(s, t, 'sgrad');
        obj.fmm          = @(eps,s,t,sigma) sgrad_fmm(eps, s, t, sigma);

    otherwise
        error('KERNEL3D.LAP2D: unknown Laplace kernel type ''%s''.', type);

end

obj.getquad = @(S,eps,varargin) rsc_to_sparse( ...
    lap2d.get_quadrature_correction(S, obj.type, eps, varargin{:}), S, obj.opdims(1));
obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, obj.zk, obj.kernel_order);

icheck = exist(['fmm2d.' mexext], 'file');
if ( icheck ~= 3 )
    obj.fmm = [];
end

end

function spmat = rsc_to_sparse(Q, S, nker)
%RSC_TO_SPARSE  Convert lap2d getquad RSC output to a sparse matrix.
% The gradient kernel stores two rows per target, Q.wnear is (2,nquad).
if ( nker == 1 )
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
else
    spmat = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ...
        kernel3d.rsc_interleave_full(nker, 1));
end
end

function pot = sgrad_fmm(eps, s, t, sigma)
%SGRAD_FMM  Gradient of the single layer, d/dx and d/dy interleaved.
[~, grad] = lap2d.fmm(eps, s, t, 's', sigma);
pot = grad(:);
end
