function obj = iceflex(type, nu, rts, ejs)
%KERNEL3D.ICEFLEX   Construct a flexural surface wave kernel in 3D.
%
%   KERNEL3D.ICEFLEX(type, nu, rts, ejs)
%
%   NU is the Poisson ratio; RTS and EJS are the dispersion roots and their
%   partial fraction residues from SURFWAVE.FLEX.FIND_ROOTS_FLEX.
%
%   Scalar kernels, [1 1]:
%
%      'gs'          - single layer of G_S
%      'gphi'        - single layer of G_phi
%      'gphi_bilap'  - bilaplacian applied to G_phi
%      's3d_gphi'    - 3D single layer applied to G_phi
%      's3d'         - Laplace 3D single layer
%      's3d_sum'     - s3d + s3d_gphi
%
%   Volume-to-boundary kernels, [2 1]:
%
%      'gs_v2b'      - G_S traces on the boundary
%      'gphi_v2b'    - G_phi traces on the boundary
%      'v2b'         - -gamma/4*gs_v2b + 1/2*gphi_v2b, the combined field of
%                      the polynya solve (one Fortran near-quadrature call)
%
%   The Fortran near quadrature (getnearquad_flex_bcs) integrates
%   Q(g) = -g/4*G_S + 1/2*G_phi for the gamma g stored in zpars(12). Each
%   v2b type gets its own weights from this by linearity: c_s*G_S + c_p*G_phi
%   is 2*c_p*Q(-2*c_s/c_p) for c_p ~= 0 (one call), and G_S alone is
%   2*(Q(-2) - Q(0)) (two calls). So scaling and adding the types gives
%   correct near weights.
%
%
% See also SURFWAVE.FLEX.KERN, KERNEL3D.GRAVITY, KERNEL3D.CAPILLARY

if ( nargin < 4 )
    error('KERNEL3D.ICEFLEX: need type, nu, rts and ejs.');
end

rts = rts(:);  ejs = ejs(:);

% zpars = [rts(5); ejs(5); alpha; gamma; nu], with alpha and gamma
% recovered from the roots of alpha z^5 + gamma z - 1.
alpha = 1/prod(rts);
gamma = alpha*sum(prod(nchoosek(rts,4), 2));
if ( abs(imag(alpha)) > 1e-10*abs(alpha) || abs(imag(gamma)) > 1e-10*max(abs(gamma),1) )
    error('KERNEL3D.ICEFLEX: rts do not come from a real (alpha, gamma) pair.');
end
alpha = real(alpha);  gamma = real(gamma);

zpars = complex([rts; ejs; alpha; gamma; nu]);

obj           = kernel3d();
obj.name      = 'flexural';
obj.ifcomplex = 1;
obj.zk        = max(abs(rts));
obj.params.nu    = nu;
obj.params.rts   = rts;
obj.params.ejs   = ejs;
K = @(ktype) @(s,t) surfwave.flex.kern(s, t, ktype, nu, rts, ejs);

switch lower(type)

    case {'gs'}
        obj = scalar_kern(obj, 'gs',       K('gs_s'),      -1, zpars, 1);
    case {'gphi'}
        obj = scalar_kern(obj, 'gphi',     K('gphi_s'),    -1, zpars, 2);
    case {'gphi_bilap', 'bilap'}
        obj = scalar_kern(obj, 'gphi_bilap', K('gphi_bilap'), 1, zpars, 4);
    case {'s3d_gphi'}
        obj = scalar_kern(obj, 's3d_gphi', K('s3d_gphi'),  -1, zpars, 5);
    case {'s3d'}
        obj = scalar_kern(obj, 's3d', ...
                  @(s,t) surfwave.flex.lap3dkern(s.r, t.r), -1, zpars, 6);
    case {'s3d_sum', 's3d_plus_s3d_gphi'}
        obj = scalar_kern(obj, 's3d_sum', ...
                  @(s,t) surfwave.flex.lap3dkern(s.r, t.r) + ...
                         surfwave.flex.kern(s, t, 's3d_gphi', nu, rts, ejs), ...
                  -1, zpars, 7);

    case {'gs_v2b'}
        obj = eval_kern(obj, 'gs_v2b',   K('gs_v2b'),   [2 1], 1);
        obj.getquad = v2b_getquad_handle(zpars, kernel3d.rsc_interleave_full(2,1), 1, 0);
        obj.src_fields = [];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n', 'd', 'd2', 'kappa'};
    case {'gphi_v2b'}
        obj = eval_kern(obj, 'gphi_v2b', K('gphi_v2b'), [2 1], 1);
        obj.getquad = v2b_getquad_handle(zpars, kernel3d.rsc_interleave_full(2,1), 0, 1);
        obj.src_fields = [];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n', 'd', 'd2', 'kappa'};
    case {'v2b'}
        Kgs_v2b = K('gs_v2b');  Kgp_v2b = K('gphi_v2b');
        obj = eval_kern(obj, 'v2b', ...
                  @(s,t) (-gamma/4)*Kgs_v2b(s,t) + (1/2)*Kgp_v2b(s,t), [2 1], 1);
        obj.getquad = v2b_getquad_handle(zpars, kernel3d.rsc_interleave_full(2,1), -gamma/4, 1/2);
        obj.src_fields = [];
        obj.kernel_order = -1;
        obj.targ_fields  = {'n', 'd', 'd2', 'kappa'};
    otherwise
        error('KERNEL3D.ICEFLEX: unknown flexural kernel type ''%s''.', type);

end

obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, ...
                                     obj.zk, obj.kernel_order);

end


function obj = scalar_kern(obj, tname, evalfn, korder, zpars, iker)
obj.type         = tname;
obj.opdims       = [1 1];
obj.kernel_order = korder;
obj.eval         = evalfn;
ri = kernel3d.rsc_interleave_full(1, 1);
obj.getquad = @(S,eps,varargin) flex_getquad(S, eps, zpars, iker, ri, varargin{:});
end


function obj = eval_kern(obj, tname, evalfn, opdims, korder)
obj.type         = tname;
obj.opdims       = opdims;
obj.kernel_order = korder;
obj.eval         = evalfn;
obj.getquad      = [];
obj.src_fields   = {'n', 'd'};
obj.targ_fields  = {'n', 'd', 'd2'};
end


function h = v2b_getquad_handle(zpars, ri, cs, cp)
% near quadrature for cs*gs_v2b + cp*gphi_v2b
h = @(S,eps,varargin) flex_v2b_getquad(S, eps, zpars, ri, cs, cp, varargin{:});
end


function spmat = flex_getquad(S, eps, zpars, iker, ri, targinfo, opts)
%FLEX_GETQUAD  Near-quadrature for a scalar flexural volume kernel.
if ( nargin < 6 || isempty(targinfo) ), targinfo = S;        end
if ( nargin < 7 || isempty(opts) ),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);
if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end
nquad = rsc.iquad(end) - 1;
[patch_id, uvs_targ] = targ_ids(targinfo);

wnear = surfwave.flex.getnearquad_flex(S.npatches, norders, ixyzs, iptype, ...
    S.npts, srccoefs, srcvals, targinfo, patch_id, uvs_targ, eps, 1, ...
    rsc.nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, zpars, ...
    nquad, iker);

spmat = pack(S, rsc, wnear, nquad, ri);
end


function spmat = flex_v2b_getquad(S, eps, zpars, ri, cs, cp, targinfo, opts)
%FLEX_V2B_GETQUAD  Near-quadrature for cs*gs_v2b + cp*gphi_v2b.
%  The Fortran routine integrates Q(g) = -g/4*G_S + 1/2*G_phi, g = zpars(12);
%  the requested combination is built from Q by linearity (see header).
if ( nargin < 7 || isempty(targinfo) )
    error('KERNEL3D.ICEFLEX: v2b kernels require a boundary targinfo.');
end
if ( nargin < 8 || isempty(opts) ), opts = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);
if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end
nquad = rsc.iquad(end) - 1;
[patch_id, uvs_targ] = targ_ids(targinfo);

Q = @(g) surfwave.flex.getnearquad_v2b_flex(S.npatches, norders, ixyzs, ...
    iptype, S.npts, srccoefs, srcvals, targinfo, patch_id, uvs_targ, eps, ...
    1, rsc.nnz, rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, ...
    [zpars(1:11); g; zpars(13:end)], nquad);
if ( cp ~= 0 )
    wnear = (2*cp) * Q(-2*cs/cp);
else
    wnear = (2*cs) * (Q(-2) - Q(0));
end

spmat = pack(S, rsc, wnear, nquad, ri);
end


function [patch_id, uvs_targ] = targ_ids(targinfo)
ntarg = size(extract_targ_array(targinfo), 2);
patch_id = [];
if ( isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id') )
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
end
if ( isempty(patch_id) )
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end
end


function spmat = pack(S, rsc, wnear, nquad, ri)
if ( size(wnear, 1) == nquad && size(wnear, 2) ~= nquad )
    wnear = wnear.';
end
spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);
end

