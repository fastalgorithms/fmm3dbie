function obj = radpp(green, rmax, type, varargin)
%KERNEL3D.RADPP   Kernel built from a piecewise polynomial fit of a radial
%   Green's function and its derivatives.
%
%   K = KERNEL3D.RADPP(green, rmax, type)
%   K = KERNEL3D.RADPP(green, rmax, type, 'nu', nu, 'ipv', ipv, ...)
%   K = KERNEL3D.RADPP(green, rmax, wfun, 'nderiv', nd, 'opdims', [m 1], ...)
%   K = KERNEL3D.RADPP(green, rmax)            % type 's'
%   K = KERNEL3D.RADPP(green, rmax, opts)      % type 's', options struct
%
%   In the last form OPTS holds the fit options of RADPP_FIT and may also
%   set rmin, zk and kernel_order (the old KERNEL3D.RADCHEB interface).
%
%   GREEN is a function handle with the usual calling sequence
%
%      [val, grad, hess, third] = green(src, targ)
%
%   for src (dim,ns) and targ (dim,nt). G must be radial; it is only
%   sampled on the x axis. Sums can be formed with RADPP_SUM.
%
%   GREEN may also be a scalar kernel3d object. Its zk and kernel_order are
%   then used unless given.
%
%   RMAX is the largest distance needed, or [rmin rmax]. If rmin is not
%   given it is 1e-12. Distances outside [rmin rmax] give NaN in K.eval.
%
%   The kernel is
%
%      K_j(x,y) = sum_c w_jc(x) D_c G(x-y),   j = 1..opdims(1)
%
%   with weights w that depend only on the target and the derivative
%   tensor components D_c G (see RADPP_TENSORS), for 'dim' = 3
%
%      [G | G_x G_y G_z | H_xx H_xy H_xz H_yy H_yz H_zz |
%       T_xxx T_xxy T_xxz T_xyy T_xyz T_xzz T_yyy T_yyz T_yzz T_zzz]
%
%   and for 'dim' = 2
%
%      [G | G_x G_y | H_xx H_xy H_yy | T_xxx T_xxy T_xyy T_yyy]
%
%   Each row is computed by its own call to the GGQ near quadrature, with
%   the contraction inside the integrand.
%
%   Built-in types (targinfo fields: n, d, d2 as in chunkie):
%      's'            ('single', 'slp')            G
%      'sp'           ('sprime')                   G_{n_x}
%      'clamped_bc'   ('clamped', 'clamped_plate_bcs')
%                                                  [G; G_{n_x}]
%   dim = 2 only:
%      'supported_bc' ('supported', 'supported_plate_bcs')
%                                                  [G; G_{n_x n_x} + nu G_{t_x t_x}]
%      'free_bc'      ('free', 'free_plate_bcs')   [moment; shear], as in
%                                                  chnk.flex2d.kern 'free_plate_bcs'
%
%   Custom kernels: pass a function handle WFUN in place of TYPE,
%
%      W = wfun(targinfo)     % (nt, ncomp, opdims(1))
%
%   with ncomp = 1, 4, 10, 20 (dim 3) or 1, 3, 6, 10 (dim 2) for
%   nderiv = 0, 1, 2, 3, and set 'nderiv',
%   'opdims' and 'targ_fields'.
%
%   Name-value options:
%      'dim'          - 3 (default) or 2, see above
%      'nu'           - Poisson ratio (supported_bc, free_bc)
%      'ipv'          - ipv for the GGQ quadrature, scalar or one per row.
%                       Default 0 (weakly singular), except 1 for the
%                       derivative rows of 'sp' and 'clamped_bc' in 3D.
%                       Use 1 for principal value or hypersingular kernels
%      'opts'         - fit options, see RADPP_FIT (eps default 1e-13)
%      'nderiv'       - custom kernels: highest derivative used (0..3)
%      'opdims'       - custom kernels: [m 1] (default [1 1])
%      'targ_fields'  - custom kernels: target fields used by wfun
%      'name'         - custom kernels: name (default 'custom')
%      'zk'           - wavenumber for oversampling (default 0)
%      'kernel_order' - for oversampling (default from nderiv)
%      'ifcomplex'    - custom kernels: whether the weights are complex
%                       (default true)
%
%   Example:
%      green = @(s,t) chnk.flex2d.hkdiffgreen(zk, s, t);
%      K = kernel3d.radpp(green, 10, 'free_bc', 'nu', 0.3, 'dim', 2);
%
% See also RADPP_FIT, RADPP_TENSORS

if nargin < 2 || isempty(rmax)
    error('KERNEL3D.RADPP: need green and rmax.');
end
if nargin < 3 || isempty(type), type = 's'; end

% radcheb-style call: radpp(green, rmax, opts) with an options struct
% (fit options plus rmin, zk, kernel_order) and type 's'
if isstruct(type)
    fo = type;
    type = 's';
    if isfield(fo, 'rmin') && numel(rmax) == 1, rmax = [fo.rmin rmax]; end
    for f = {'zk', 'kernel_order'}
        if isfield(fo, f{1}), varargin = [varargin, f, {fo.(f{1})}]; end
    end
    fo = rmfield(fo, intersect(fieldnames(fo), {'rmin','zk','kernel_order'}));
    varargin = [varargin, {'opts', fo}];
end

o = struct('nu', [], 'ipv', [], 'opts', struct(), 'nderiv', [], ...
           'opdims', [1 1], 'targ_fields', {{}}, 'name', 'custom', ...
           'zk', [], 'kernel_order', [], 'ifcomplex', true, 'dim', 3);
o = parse_options(o, varargin);

valueonly = isa(green, 'kernel3d');
if valueonly
    [green, zk0, korder0] = kernel3d_green(green);
    if isempty(o.zk), o.zk = zk0; end
    if isempty(o.kernel_order), o.kernel_order = korder0; end
end
if isempty(o.zk), o.zk = 0; end

rmin = 1e-12;
if numel(rmax) == 2
    rmin = rmax(1);
    rmax = rmax(2);
end

% the kernel: weights, derivative order, shape
if isa(type, 'function_handle')
    if isempty(o.nderiv)
        error('KERNEL3D.RADPP: custom kernels need ''nderiv''.');
    end
    spec = struct('name', o.name, 'nderiv', o.nderiv, ...
                  'nrow', o.opdims(1), 'targ_fields', {o.targ_fields}, ...
                  'wfun', type);
else
    spec = builtin_kernel(type, o.nu, o.dim);
end
nderiv = spec.nderiv;
nrow   = spec.nrow;
if valueonly && nderiv > 0
    error('KERNEL3D.RADPP: a kernel3d green only gives values, nderiv must be 0.');
end
if ~any(o.dim == [2 3])
    error('KERNEL3D.RADPP: dim must be 2 or 3.');
end
if o.dim == 3
    ncomp = [1 4 10 20];
else
    ncomp = [1 3 6 10];
end
ncomp  = ncomp(nderiv+1);

ipv = o.ipv(:).';
if isempty(ipv)
    % 3D derivative rows use the principal value rules, as helm3d 'sp'
    ipv = 0;
    if isfield(spec, 'ipv3') && o.dim == 3, ipv = spec.ipv3; end
end
if isscalar(ipv), ipv = repmat(ipv, 1, nrow); end
if numel(ipv) ~= nrow
    error('KERNEL3D.RADPP: need one ipv per row (%d).', nrow);
end

kernel_order = o.kernel_order;
if isempty(kernel_order), kernel_order = min(nderiv, 2) - 1; end

% fit the on-axis profiles
fitopts = o.opts;
if ~isfield(fitopts, 'eps'), fitopts.eps = 1e-13; end
[ipars, dpars, info] = radpp_fit(@(r) profiles(green, r, nderiv), ...
                                 rmin, rmax, fitopts);
ipars = [ipars; nderiv; ncomp; o.dim];

obj              = kernel3d();
obj.name         = 'radpp';
obj.type         = spec.name;
obj.opdims       = [nrow 1];
obj.zk           = o.zk;
obj.kernel_order = kernel_order;
obj.ifcomplex    = ~isreal(info.coefs) || ~isreal(o.nu) || ...
                   (isa(type, 'function_handle') && o.ifcomplex);
obj.targ_fields  = spec.targ_fields;
obj.src_fields   = {};

obj.params = struct('green', green, 'rmin', rmin, 'rmax', rmax, ...
    'nderiv', nderiv, 'ncomp', ncomp, 'dim', o.dim, 'nu', o.nu, 'ipv', ipv, ...
    'wfun', spec.wfun, 'ipars', ipars, 'dpars', dpars, ...
    'breaks', info.breaks, 'coefs', info.coefs, 'err', info.err, ...
    'fitopts', fitopts);

obj.params.targs = @(t) radpp_targs(t, spec.wfun);   % for radpp_kern

p = obj.params;
obj.eval = @(s,t) radpp_eval_kern(s, t, p);
obj.fmm  = [];

ri = kernel3d.rsc_interleave_full(nrow, 1);
ifcomplex = obj.ifcomplex;
obj.getquad = @(S, eps, varargin) radpp_getquad(S, eps, p, ifcomplex, ...
                                                ri, varargin{:});
obj.get_overs_orders = @(S, t, eps) kernel3d.kernel3d_getnear_overs(S, ...
                           t, eps, obj.zk, obj.kernel_order);

end


%%% built-in kernels

function spec = builtin_kernel(type, nu, dim)

needs_nu = {'supported_bc', 'supported', 'supported_plate_bcs', ...
            'free_bc', 'free', 'free_plate_bcs'};
if any(strcmpi(type, needs_nu)) && isempty(nu)
    error('KERNEL3D.RADPP: type ''%s'' needs ''nu''.', type);
end
if any(strcmpi(type, needs_nu)) && dim ~= 2
    error('KERNEL3D.RADPP: type ''%s'' needs ''dim'' = 2.', type);
end

switch lower(type)
    case {'s', 'single', 'slp'}
        spec = kspec('s', 0, 1, {}, @(t) ones(size(t.r,2), 1));
    case {'sp', 'sprime'}
        spec = kspec('sp', 1, 1, {'n'}, @(t) w_sp(t, dim));
        spec.ipv3 = 1;
    case {'clamped_bc', 'clamped', 'clamped_plate_bcs'}
        spec = kspec('clamped_bc', 1, 2, {'n'}, @(t) w_clamped(t, dim));
        spec.ipv3 = [0 1];
    case {'supported_bc', 'supported', 'supported_plate_bcs'}
        spec = kspec('supported_bc', 2, 2, {'n', 'd'}, ...
                     @(t) w_supported(t, nu));
    case {'free_bc', 'free', 'free_plate_bcs'}
        spec = kspec('free_bc', 3, 2, {'n', 'd', 'd2'}, @(t) w_free(t, nu));
    otherwise
        error('KERNEL3D.RADPP: unknown type ''%s''.', type);
end

end


function spec = kspec(name, nderiv, nrow, targ_fields, wfun)
spec = struct('name', name, 'nderiv', nderiv, 'nrow', nrow, ...
              'targ_fields', {targ_fields}, 'wfun', wfun);
end


function W = w_sp(t, dim)
% G_{n_x}
n = normals(t, dim);
W = [zeros(size(n,1),1), n];
end


function W = w_clamped(t, dim)
% [G; G_{n_x}]
n = normals(t, dim);
z = zeros(size(n,1), 1);
W = cat(3, [1+z, 0*n], [z, n]);
end


function n = normals(t, dim)
% unit normals (nt,dim), padded with zeros if t.n has fewer rows
n = zeros(size(t.n,2), dim);
m = min(dim, size(t.n,1));
n(:,1:m) = t.n(1:m,:).';
end


function W = w_supported(t, nu)
% [G; G_{n_x n_x} + nu G_{t_x t_x}]
[n, tau, ~] = frame(t);
z = zeros(size(n,1), 1);
W = cat(3, [1+z, z, z, z, z, z], ...
           [z, z, z, sym2(n,n) + nu*sym2(tau,tau)]);
end


function W = w_free(t, nu)
% chnk.flex2d.kern 'free_plate_bcs'
%   moment = G_{nn} + nu G_{tt}
%   shear  = G_{nnn} + (2-nu) G_{ntt} + (1-nu) kappa (G_{tt} - G_{nn})
[n, tau, kappa] = frame(t);
z = zeros(size(n,1), 1);
moment = [z, z, z, sym2(n,n) + nu*sym2(tau,tau), z, z, z, z];
shear  = [z, z, z, (1-nu)*kappa.*(sym2(tau,tau) - sym2(n,n)), ...
          sym3(n,n,n) + (2-nu)*sym3(n,tau,tau)];
W = cat(3, moment, shear);
end


function [n, tau, kappa] = frame(t)
% unit normal, unit tangent (nt,2) and signed curvature (nt,1) at the
% targets, from chunkie fields n, d, d2
n = t.n(1:2,:).';
tau = [];  kappa = [];
if isfield(t, 'd') || isprop(t, 'd')
    d = t.d(1:2,:).';
    ds = sqrt(sum(d.^2, 2));
    tau = d./ds;
    if isfield(t, 'd2') || isprop(t, 'd2')
        d2 = t.d2(1:2,:).';
        kappa = (d(:,1).*d2(:,2) - d2(:,1).*d(:,2))./ds.^3;
    end
end
end


function w = sym2(u, v)
% weights of [H_xx H_xy H_yy] for the second derivative along u then v
w = [u(:,1).*v(:,1), u(:,1).*v(:,2) + u(:,2).*v(:,1), u(:,2).*v(:,2)];
end


function w = sym3(u, v, s)
% weights of [T_xxx T_xxy T_xyy T_yyy] for the third derivative along u,
% v, s
ux = u(:,1);  uy = u(:,2);
vx = v(:,1);  vy = v(:,2);
sx = s(:,1);  sy = s(:,2);
w = [ux.*vx.*sx, ...
     ux.*vx.*sy + ux.*vy.*sx + uy.*vx.*sx, ...
     ux.*vy.*sy + uy.*vx.*sy + uy.*vy.*sx, ...
     uy.*vy.*sy];
end


%%% fit, eval and quadrature

function [green, zk, kernel_order] = kernel3d_green(K)
% value-only Green's function from a scalar, position-only kernel3d
if ~isequal(K.opdims, [1 1])
    error('KERNEL3D.RADPP: only scalar kernel3d objects can be used as green.');
end
if ~isempty(K.src_fields) || ~isempty(K.targ_fields)
    error(['KERNEL3D.RADPP: kernel depends on source or target fields, ' ...
           'so it is not a function of r alone.']);
end
pt3 = @(x) struct('r', [x; zeros(3-size(x,1), size(x,2))]);
green = @(s,t) K.eval(pt3(s), pt3(t));
zk = K.zk;
kernel_order = K.kernel_order;
end


function f = profiles(green, r, nderiv)
% on-axis profiles [G | G_x | G_xx G_yy | G_xxx G_xyy], (numel(r), *)
r = r(:);
src  = [0; 0];
targ = [r.'; zeros(1, numel(r))];
out = cell(1, nderiv+1);
[out{:}] = green(src, targ);
f = out{1}(:);
if nderiv >= 1, f = [f, out{2}(:,1,1)]; end
if nderiv >= 2, f = [f, out{3}(:,1,1), out{3}(:,1,3)]; end
if nderiv >= 3, f = [f, out{4}(:,1,1), out{4}(:,1,3)]; end
end


function submat = radpp_eval_kern(srcinfo, targinfo, p)
% (nrow*nt, ns), rows interleaved
src  = srcinfo.r;
targ = targinfo.r;
ns = size(src, 2);  nt = size(targ, 2);

dx = targ(1,:).' - src(1,:);
dy = targ(2,:).' - src(2,:);
dz = 0;
if size(src,1) > 2 && size(targ,1) > 2
    dz = targ(3,:).' - src(3,:);
end
dz = dz + zeros(nt, ns);

comps = radpp_tensors(dx, dy, dz, p.breaks, p.coefs, p.nderiv, p.dim);
comps = reshape(comps, nt, ns, p.ncomp);

W = p.wfun(targinfo);
nrow = size(W, 3);
submat = zeros(nrow*nt, ns);
for j = 1:nrow
    Kj = 0;
    for c = 1:p.ncomp
        Kj = Kj + W(:,c,j).*comps(:,:,c);
    end
    submat(j:nrow:end, :) = Kj;
end
end


function spmat = radpp_getquad(S, eps, p, ifcomplex, ri, targinfo, opts)
% near quadrature: target locations followed by the weights of each row

if nargin < 6 || isempty(targinfo), targinfo = S;        end
if nargin < 7 || isempty(opts),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);

if isfield(opts, 'rsc') && ~isempty(opts.rsc)
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end
nquad = rsc.iquad(end) - 1;

targs = radpp_targs(targinfo, p.wfun);
ntarg = size(targs, 2);

patch_id = [];
if isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id')
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
end
if isempty(patch_id)
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end

wnear = getnearquad_radpp(S.npatches, norders, ixyzs, iptype, S.npts, ...
    srccoefs, srcvals, targs, patch_id, uvs_targ, eps, p.ipv, rsc.nnz, ...
    rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, nquad, p.dpars, ...
    p.ipars);
if ~ifcomplex
    wnear = real(wnear);
end

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end


function targs = radpp_targs(targinfo, wfun)
% target array for radpp_kern: locations, then the complex weights of
% each row as interlaced real/imaginary parts, (3 + 2*ncomp*nrow, nt)
r = targinfo.r;
ntarg = size(r, 2);
xyz = zeros(3, ntarg);
xyz(1:size(r,1), :) = r;

W = wfun(targinfo);                            % (nt, ncomp, nrow)
W = permute(W, [2 3 1]);                       % (ncomp, nrow, nt)
wri = zeros([2, size(W)]);
wri(1,:,:,:) = real(W);
wri(2,:,:,:) = imag(W);
targs = [xyz; reshape(wri, [], ntarg)];
end


function o = parse_options(o, args)
if mod(numel(args), 2) ~= 0
    error('KERNEL3D.RADPP: options must be name-value pairs.');
end
for i = 1:2:numel(args)
    name = lower(args{i});
    if ~isfield(o, name)
        error('KERNEL3D.RADPP: unknown option ''%s''.', args{i});
    end
    o.(name) = args{i+1};
end
end
