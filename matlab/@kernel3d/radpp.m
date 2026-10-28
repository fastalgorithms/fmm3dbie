function obj = radpp(green, rmax, type, varargin)
%KERNEL3D.RADPP   Kernel built from a piecewise polynomial fit of a radial
%   Green's function and its derivatives.
%
%   K = KERNEL3D.RADPP(green, rmax, type)
%   K = KERNEL3D.RADPP(green, rmax, type, coefs)  % or nu, plate kernels
%   K = KERNEL3D.RADPP(green, rmax, type, 'nu', nu, 'ipv', ipv, ...)
%   K = KERNEL3D.RADPP(green, rmax, wfun, 'nderiv', nd, 'opdims', [m 1], ...)
%   K = KERNEL3D.RADPP(green, rmax)            % type 's'
%   K = KERNEL3D.RADPP(green, rmax, opts)      % type 's', options struct
%
%   In the last form OPTS holds the fit options of RADPP_FIT and may also
%   set rmin, zk and kernel_order (the old KERNEL3D.RADCHEB interface).
%
%   GREEN is a radial function handle in the helm3d.green format, or the
%   equivalent 2D format. Sums can be formed with RADPP_SUM.
%
%   GREEN may also be a scalar kernel3d object. Its zk and kernel_order are
%   then used unless given.
%
%   RMAX is the largest distance needed, or [rmin rmax]. If rmin is not
%   given it is 1e-12. Distances outside [rmin rmax] give NaN in K.eval.
%
%   The kernel is
%
%      K_j(x,y) = sum_c w_jc(x) D_c G(x-y) + sum_c v_jc(x) d_{n_y} D_c G(x-y)
%
%   for j = 1...opdims(1), with weights w, v that depend only on the target
%   (v only for kernels with source derivatives) and the derivative tensor
%   components D_c G (see RADPP_ROTATE), for 'dim' = 3
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
%   Built-in types (names as in helm3d; targinfo fields n and, for the
%   plate kernels, d, d2; srcinfo field n), with c = 'coefs' = [alpha beta]:
%      's'            ('single', 'slp')            G
%      'd'            ('double', 'dlp')            G_{n_y}
%      'c'            ('comb', 'combined')         alpha G + beta G_{n_y}
%      'sp'           ('sprime')                   G_{n_x}
%      'dp'           ('dprime')                   G_{n_x n_y}
%      'cp'           ('combprime', 'cprime')      alpha G_{n_x} + beta G_{n_x n_y}
%      'grad'         ('sgrad', 'sg')              grad_x G, opdims [dim 1]
%      'dgrad'        ('dg')                       grad_x G_{n_y}
%      'cgrad'        ('cg')                       alpha grad_x G + beta grad_x G_{n_y}
%      'hess'         ('shess')                    [G_xx G_xy G_xz G_yy G_yz G_zz]
%                                                  (dim 2: [G_xx G_xy G_yy])
%      'trans_sys'    ('transmission_sys', 'trans', 'transmission', 'tsys')
%                                                  [c11 D, c12 S; c21 D', c22 S'],
%                                                  coefs 2x2 (default ones)
%      'trans_rep'    ('transmission_rep', 'trep') [c1 D, c2 S]
%      'trep_prime'   ('trepp', 'trans_rep_prime') [c1 D', c2 S']
%      'trep_grad'    ('trepg', 'trans_rep_grad')  [c1 grad_x D, c2 grad_x S]
%   dim = 3 only, on Cartesian vector densities, opdims [3 3]:
%      'mfie'         ('nxcurls', 'n_cross_curl_s') n_x x curl S[J]
%      'curls'                                     curl S[J] = grad G x J
%      'clamped_bc'   ('clamped', 'clamped_plate_bcs')
%                                                  [G; G_{n_x}]
%      's2trans'      ('s2t')                      [c1 G; c2 G_{n_x}]
%      'd2trans'      ('d2t')                      [G_{n_y}; G_{n_x n_y}]
%      'c2trans'      ('c2t', 'c2tr')              ['c'; 'cp']
%   dim = 2 only:
%      'supported_bc' ('supported', 'supported_plate_bcs')
%                                                  [G; G_{n_x n_x} + nu G_{t_x t_x}]
%      'free_bc'      ('free', 'free_plate_bcs')   [moment; shear]
%
%   Custom kernels: pass a function handle WFUN in place of TYPE,
%
%      W = wfun(targinfo)     % (nt, ncomp, opdims(1)*opdims(2))
%
%   with the block entries in column-major order%
%   with ncomp = 1, 4, 10, 20 (dim 3) or 1, 3, 6, 10 (dim 2) for
%   nderiv = 0, 1, 2, 3, and set 'nderiv', 'opdims' and 'targ_fields'.
%   'srcwfun' gives the weights v in the same form; nderiv + 1 is then at
%   most 3.
%
%   Name-value options:
%      'dim'          - 3 (default) or 2, see above
%      'nu'           - Poisson ratio (supported_bc, free_bc)
%      'coefs'        - [alpha beta] for 'c', 'cp', 'c2trans' (default [1 1])
%      'ipv'          - ipv for the GGQ quadrature, scalar or one per entry:
%                       0 weakly singular, 1 principal value, 2
%                       hypersingular. By default this is
%                       detected automatically from the behaviour on a
%                       test surface
%      'srcwfun'      - custom kernels: source derivative weights
%      'opts'         - fit options, see RADPP_FIT (eps default 1e-13)
%      'nderiv'       - custom kernels: highest derivative used (0..3)
%      'opdims'       - custom kernels: [m n] (default [1 1])
%      'targ_fields'  - custom kernels: target fields used by wfun
%      'name'         - custom kernels: name (default 'custom')
%      'zk'           - wavenumber for oversampling (default 0), the
%                       largest if several
%      'kernel_order' - for oversampling (default max(ipv) - 1)
%      'ifcomplex'    - custom kernels: whether the weights are complex
%                       (default true)
%      'fmm'          - FMM for green, U = fmm(eps, srcuse, pg, targ, pgt)
%                       e.g. @(eps,s,pg,t,pgt) hfmm3d(eps, zk, s, pg, t, pgt).
%
%   Sums of radpp kernels of one built-in type are radpp kernels refitted
%   on the sum of their greens, on the distances common to them; multiples
%   scale the fit.
%
%   Example:
%      green = @(s,t) helm3d.green(zk, s, t);
%      fmm = @(eps,s,pg,t,pgt) hfmm3d(eps, zk, s, pg, t, pgt);
%      K = kernel3d.radpp(green, 10, 'cp', [1 1i], 'fmm', fmm);
%
% See also RADPP_FIT, RADPP_ROTATE, RADPP_FMM

% kernel arithmetic: kernel3d.radpp('summable', KS) and
% kernel3d.radpp('sum', KS, CS) = sum_k CS(k)*KS{k}; and
% kernel3d.radpp('fit', K, S, targinfo, rsc) for radgreen kernels
if ischar(green)
    switch green
        case 'summable'
            obj = summable(rmax);
        case 'fit'
            obj = radgreen_fit(rmax.params, type, varargin{:});
        otherwise
            obj = ppsum(rmax, type);
    end
    return
end

% green mode (from KERNEL3D.RADGREEN): green is a struct with the exact
% profiles, eval is exact and the fit is built in getquad
gmode = isstruct(green) && isfield(green, 'profiles');

if nargin < 2 || (isempty(rmax) && ~gmode)
    error('KERNEL3D.RADPP: need green and rmax.');
end
if nargin < 3 || isempty(type), type = 's'; end

% a prebuilt kernel spec (from kernel arithmetic on green kernels)
isspec = isstruct(type) && isfield(type, 'wfun');

% radcheb-style call: radpp(green, rmax, opts) with an options struct
% (fit options plus rmin, zk, kernel_order) and type 's'
if isstruct(type) && ~isspec
    fo = type;
    type = 's';
    if isfield(fo, 'rmin') && numel(rmax) == 1, rmax = [fo.rmin rmax]; end
    for f = {'zk', 'kernel_order'}
        if isfield(fo, f{1}), varargin = [varargin, f, {fo.(f{1})}]; end
    end
    fo = rmfield(fo, intersect(fieldnames(fo), {'rmin','zk','kernel_order'}));
    varargin = [varargin, {'opts', fo}];
end

% positional parameter after the type, as in helm3d: coefs, or nu for
% the plate kernels
npos = find(cellfun(@ischar, varargin), 1) - 1;
if isempty(npos), npos = numel(varargin); end
if npos > 1
    error('KERNEL3D.RADPP: one positional parameter (coefs or nu) after the type.');
end
if npos == 1
    if ~ischar(type)
        error('KERNEL3D.RADPP: positional coefs or nu need a built-in type.');
    end
    pname = 'coefs';
    if any(strcmpi(type, {'supported_bc', 'supported', 'supported_plate_bcs', ...
                          'free_bc', 'free', 'free_plate_bcs'}))
        pname = 'nu';
    end
    if any(strcmpi(varargin(2:2:end), pname))
        error('KERNEL3D.RADPP: %s given twice.', pname);
    end
    varargin = [{pname}, varargin];
end
if ischar(type) && any(strcmpi(varargin(1:2:end), 'coefs')) && ...
        ~any(strcmpi(type, {'c', 'comb', 'combined', 'cp', 'combprime', ...
            'cprime', 'c2trans', 'c2t', 'c2tr', 'trans_sys', 'transmission_sys', ...
            'trans', 'transmission', 'tsys', 'trans_rep', 'transmission_rep', ...
            'trep', 's2trans', 's2t', 'cgrad', 'cg', 'trep_prime', 'trepp', 'trep_p', 'trans_rep_prime', ...
            'trans_rep_p', 'trep_grad', 'trepg', 'trep_g', 'trans_rep_grad', ...
            'trans_rep_g'}))
    error('KERNEL3D.RADPP: type ''%s'' takes no coefs.', type);
end

o = struct('nu', [], 'ipv', [], 'opts', struct(), 'nderiv', [], ...
           'opdims', [1 1], 'targ_fields', {{}}, 'name', 'custom', ...
           'zk', [], 'kernel_order', [], 'ifcomplex', true, 'dim', 3, ...
           'coefs', [1 1], 'srcwfun', [], 'fmm', [], 'fit', []);
o = parse_options(o, varargin);

valueonly = isa(green, 'kernel3d');
if valueonly
    [green, zk0, korder0] = kernel3d_green(green);
    if isempty(o.zk), o.zk = zk0; end
    if isempty(o.kernel_order), o.kernel_order = korder0; end
end
if isempty(o.zk), o.zk = 0; end

rmin = 1e-12;
if gmode, rmin = []; end
if numel(rmax) == 2
    rmin = rmax(1);
    rmax = rmax(2);
end

% the kernel: weights, derivative order, shape
if isspec
    spec = type;
elseif isa(type, 'function_handle')
    if isempty(o.nderiv)
        error('KERNEL3D.RADPP: custom kernels need ''nderiv''.');
    end
    spec = kspec(o.name, o.nderiv, o.opdims, o.targ_fields, type);
    spec.vfun = o.srcwfun;
else
    spec = builtin_kernel(type, o.nu, o.dim, o.coefs);
end
nderiv = spec.nderiv;
nrow   = size(spec.ent, 1);          % one GGQ call per block entry
isrc   = ~isempty(spec.vfun);
nfit   = nderiv + isrc;
if nfit > 3
    error('KERNEL3D.RADPP: at most 3 derivatives, target and source together.');
end
if valueonly && nfit > 0
    error('KERNEL3D.RADPP: a kernel3d green only gives values, no derivatives.');
end
if ~any(o.dim == [2 3 13])          % 13: tensor green, from radgreen
    error('KERNEL3D.RADPP: dim must be 2 or 3.');
end
if o.dim == 13 && isrc
    error('KERNEL3D.RADPP: no source derivatives with a tensor green.');
end
if gmode
    green.check(nfit, spec.name);
end
if o.dim == 3
    ncomp = [1 4 10 20];
elseif o.dim == 13
    ncomp = [6 24];
else
    ncomp = [1 3 6 10];
end
ncomp  = ncomp(nderiv+1);

ipv = o.ipv(:).';
ipvauto = isempty(ipv);
if ipvauto
    % from the singularity of each entry of the kernel itself (its weights
    % applied to the green's derivatives) on a test surface, flat in 2D
    if gmode
        prof0 = @(r) green.profiles(r, nfit);
    else
        prof0 = @(r) profiles(green, r, nfit);
    end
    p0 = struct('nderiv', nderiv, 'isrc', isrc, 'dim', o.dim, ...
                'gmode', true, 'profun', prof0, ...
                'wfun', spec.wfun, 'vfun', spec.vfun, 'opdims', spec.opdims, ...
                'ent', spec.ent, 'ncomp', ncomp);
    ipv = surface_ipv(p0, spec.targ_fields, o.dim == 2);
end
if isscalar(ipv), ipv = repmat(ipv, 1, nrow); end
if numel(ipv) ~= nrow
    error('KERNEL3D.RADPP: need one ipv per entry (%d).', nrow);
end

kernel_order = o.kernel_order;
if gmode
    if ~isempty(kernel_order)
        error('KERNEL3D.RADGREEN: kernel_order follows ipv; set ''ipv'' instead.');
    end
    kernel_order = max(ipv) - 1;
elseif isempty(kernel_order)
    kernel_order = max(ipv) - 1;
end

% fit the on-axis profiles (in green mode, at getquad time)
fitopts = o.opts;
if ~isfield(fitopts, 'eps'), fitopts.eps = 1e-13; end
ipost = [nderiv; ncomp; o.dim; isrc];
if gmode
    profun = @(r) green.profiles(r, nfit);
    ipars = [];  dpars = [];
    info = struct('breaks', [], 'coefs', profun([0.5; 1; 2]), 'err', []);
else
    profun = @(r) profiles(green, r, nfit);
    if isempty(o.fit)
        [ipars, dpars, info] = radpp_fit(profun, rmin, rmax, fitopts);
    else
        % a fit given by the kernel arithmetic, see ppsum
        ipars = o.fit.ipars;  dpars = o.fit.dpars;  info = o.fit;
    end
    ipars = [ipars; ipost];
end

obj              = kernel3d();
obj.name         = 'radpp';
if gmode, obj.name = 'radgreen'; end
obj.type         = spec.name;
obj.opdims       = spec.opdims;
obj.zk           = max(abs(o.zk));      % a plate pair gives its largest
obj.kernel_order = kernel_order;
obj.ifcomplex    = double(~isreal(info.coefs) || ~isreal(o.nu) || ...
                   ~isreal(o.coefs) || ...
                   ((isa(type, 'function_handle') || isspec) && o.ifcomplex));
obj.targ_fields  = spec.targ_fields;
obj.src_fields   = {};
if isrc, obj.src_fields = {'n'}; end

obj.params = struct('green', green, 'rmin', rmin, 'rmax', rmax, ...
    'nderiv', nderiv, 'ncomp', ncomp, 'dim', o.dim, 'isrc', isrc, ...
    'nu', o.nu, 'kcoefs', o.coefs, 'ipv', ipv, ...
    'wfun', spec.wfun, 'vfun', spec.vfun, 'opdims', spec.opdims, ...
    'ent', spec.ent, 'ipars', ipars, 'dpars', dpars, ...
    'breaks', info.breaks, 'coefs', info.coefs, 'err', info.err, ...
    'fitopts', fitopts, 'spec', spec, 'gmode', gmode, ...
    'profun', [], 'ipost', ipost, 'kernel_order', kernel_order, ...
    'fmm', o.fmm);
if gmode, obj.params.profun = profun; end

obj.params.targs = @(t) radpp_targs(t, spec.wfun, spec.vfun);  % for radpp_kern
% mex availability
obj.params.usemex = exist('fmm3dbie_routs', 'file') == 3;

p = obj.params;
obj.eval = @(s,t) radpp_eval_kern(s, t, p);
obj.fmm  = [];
if ~isempty(o.fmm)
    if o.dim == 13 || nderiv > 2
        error('KERNEL3D.RADPP: ''fmm'' needs a scalar green and at most 2 derivatives.');
    end
    fmmfun = o.fmm;
    obj.fmm = @(eps, s, t, sigma) radpp_fmm(fmmfun, p, eps, s, t, sigma);
end

ri = kernel3d.rsc_interleave_full(spec.opdims(1), spec.opdims(2), nrow, ...
        struct('ker_id', num2cell(1:nrow), ...
               'row_id', num2cell(spec.ent(:,1).'), ...
               'col_id', num2cell(spec.ent(:,2).'), ...
               'coef',   num2cell(ones(1,nrow))));
ifcomplex = obj.ifcomplex;
if gmode
    obj.getquad = @(S, eps, varargin) radgreen_getquad(S, eps, p, ...
                                          ifcomplex, ri, varargin{:});
else
    obj.getquad = @(S, eps, varargin) radpp_getquad(S, eps, p, ifcomplex, ...
                                                    ri, varargin{:});
end
obj.get_overs_orders = @(S, t, eps) kernel3d.kernel3d_getnear_overs(S, ...
                           t, eps, obj.zk, obj.kernel_order);

end


%%% built-in kernels

function spec = builtin_kernel(type, nu, dim, coefs)

if dim == 13
    spec = dyad_kernel(type);
    return
end

needs_nu = {'supported_bc', 'supported', 'supported_plate_bcs', ...
            'free_bc', 'free', 'free_plate_bcs'};
if any(strcmpi(type, needs_nu)) && isempty(nu)
    error('KERNEL3D.RADPP: type ''%s'' needs ''nu''.', type);
end
if any(strcmpi(type, needs_nu)) && dim ~= 2
    error('KERNEL3D.RADPP: type ''%s'' needs ''dim'' = 2.', type);
end

al = coefs(1);  be = coefs(2);
one  = @(t) ones(size(t.r(:,:), 2), 1);
zero = @(t) zeros(size(t.r(:,:), 2), 1);
dn   = @(t) [zero(t), normals(t, dim)];        % weights of G_{n_x}
z1   = @(t) zeros(size(t.r(:,:), 2), 1+dim);

switch lower(type)
    case {'s', 'single', 'slp'}
        spec = kspec('s', 0, 1, {}, one, []);
    case {'d', 'double', 'dlp'}
        spec = kspec('d', 0, 1, {}, zero, one);
    case {'c', 'comb', 'combined'}
        spec = kspec('c', 0, 1, {}, @(t) al*one(t), @(t) be*one(t));
    case {'sp', 'sprime'}
        spec = kspec('sp', 1, 1, {'n'}, dn, []);
    case {'dp', 'dprime'}
        spec = kspec('dp', 1, 1, {'n'}, z1, dn);
    case {'cp', 'combprime', 'cprime'}
        spec = kspec('cp', 1, 1, {'n'}, @(t) al*dn(t), @(t) be*dn(t));
    case {'grad', 'sgrad', 'sg'}
        spec = kspec('grad', 1, dim, {}, @(t) w_grad(t, dim), []);
    case {'dgrad', 'dg'}
        spec = kspec('dgrad', 1, dim, {}, @(t) 0*w_grad(t, dim), ...
                     @(t) w_grad(t, dim));
    case {'cgrad', 'cg'}
        spec = kspec('cgrad', 1, dim, {}, @(t) al*w_grad(t, dim), ...
                     @(t) be*w_grad(t, dim));
    case {'hess', 'shess'}
        nh = dim*(dim+1)/2;
        spec = kspec('hess', 2, nh, {}, @(t) w_hess(t, dim), []);
    case {'trans_sys', 'transmission_sys', 'trans', 'transmission', 'tsys'}
        % entries (1,1) D, (2,1) D', (1,2) S, (2,2) S', column-major
        c = coefs;
        if isequal(size(c), [1 2]) && all(c == [1 1]), c = ones(2); end
        if ~isequal(size(c), [2 2])
            error('KERNEL3D.RADPP: trans_sys coefs must be 2x2.');
        end
        spec = kspec('trans_sys', 1, [2 2], {'n'}, ...
                     @(t) cat(3, z1(t), z1(t), [c(1,2)*one(t), zeros(size(t.r(:,:), 2), dim)], ...
                              c(2,2)*dn(t)), ...
                     @(t) cat(3, [c(1,1)*one(t), zeros(size(t.r(:,:), 2), dim)], c(2,1)*dn(t), ...
                              z1(t), z1(t)));
    case {'trans_rep', 'transmission_rep', 'trep'}
        spec = kspec('trans_rep', 0, [1 2], {}, ...
                     @(t) cat(3, zero(t), be*one(t)), ...
                     @(t) cat(3, al*one(t), zero(t)));
    case {'trep_prime', 'trepp', 'trep_p', 'trans_rep_prime', 'trans_rep_p'}
        spec = kspec('trep_prime', 1, [1 2], {'n'}, ...
                     @(t) cat(3, z1(t), be*dn(t)), ...
                     @(t) cat(3, al*dn(t), z1(t)));
    case {'trep_grad', 'trepg', 'trep_g', 'trans_rep_grad', 'trans_rep_g'}
        % entries (j,1) c1 D_{x_j}, (j,2) c2 S_{x_j}, column-major
        spec = kspec('trep_grad', 1, [dim 2], {}, ...
                     @(t) cat(3, 0*w_grad(t, dim), be*w_grad(t, dim)), ...
                     @(t) cat(3, al*w_grad(t, dim), 0*w_grad(t, dim)));
    case {'clamped_bc', 'clamped', 'clamped_plate_bcs'}
        spec = kspec('clamped_bc', 1, 2, {'n'}, ...
                     @(t) cat(3, [one(t), zeros(size(t.r(:,:), 2), dim)], dn(t)), []);
    case {'s2trans', 's2t'}
        spec = kspec('s2trans', 1, 2, {'n'}, ...
                     @(t) cat(3, [al*one(t), zeros(size(t.r(:,:), 2), dim)], be*dn(t)), []);
    case {'d2trans', 'd2t'}
        spec = kspec('d2trans', 1, 2, {'n'}, @(t) cat(3, z1(t), z1(t)), ...
                     @(t) cat(3, [one(t), zeros(size(t.r(:,:), 2), dim)], dn(t)));
    case {'c2trans', 'c2t', 'c2tr'}
        spec = kspec('c2trans', 1, 2, {'n'}, ...
                     @(t) al*cat(3, [one(t), zeros(size(t.r(:,:), 2), dim)], dn(t)), ...
                     @(t) be*cat(3, [one(t), zeros(size(t.r(:,:), 2), dim)], dn(t)));
    case {'mfie', 'nxcurls', 'n_cross_curl_s'}
        if dim ~= 3, error('KERNEL3D.RADPP: ''mfie'' needs dim = 3.'); end
        spec = kspec('mfie', 1, [3 3], {'n'}, @w_mfie, []);
    case {'curls'}
        if dim ~= 3, error('KERNEL3D.RADPP: ''curls'' needs dim = 3.'); end
        spec = kspec('curls', 1, [3 3], {}, @w_curls, []);
    case {'supported_bc', 'supported', 'supported_plate_bcs'}
        spec = kspec('supported_bc', 2, 2, {'n', 'd'}, ...
                     @(t) w_supported(t, nu), []);
    case {'free_bc', 'free', 'free_plate_bcs'}
        spec = kspec('free_bc', 3, 2, {'n', 'd', 'd2'}, ...
                     @(t) w_free(t, nu), []);
    otherwise
        error('KERNEL3D.RADPP: unknown type ''%s''.', type);
end

end


function spec = kspec(name, nderiv, opdims, targ_fields, wfun, vfun)
% opdims: m or [m n]; the weights hold the block entries in column-major
% order
if nargin < 6, vfun = []; end
if isscalar(opdims), opdims = [opdims 1]; end
[i, j] = ndgrid(1:opdims(1), 1:opdims(2));
spec = struct('name', name, 'nderiv', nderiv, 'opdims', opdims, ...
              'ent', [i(:) j(:)], 'targ_fields', {targ_fields}, ...
              'wfun', wfun, 'vfun', vfun);
end


function W = w_hess(t, dim)
% Hessian: entry j has weight 1 on the j-th second derivative component
nt = size(t.r(:,:), 2);
nh = dim*(dim+1)/2;
W = zeros(nt, 1+dim+nh, nh);
for j = 1:nh
    W(:, 1+dim+j, j) = 1;
end
end


function ipv = surface_ipv(p, targ_fields, flat)
% ipv per block entry from the kernel on the paraboloid z = (x^2+3y^2)/2,
% or the plane z = 0 if flat (in-plane kernels, normals in the plane),
% target at the origin: |K| ~ rho^-q in the geodesic distance rho (log
% or bounded: q = 0) gives ipv 0 for q = 0, for q = 1 0 if isotropic or
% flat and 1 if not, 1 for q = 2 and 2 for q >= 3. Target fields other
% than n, d, d2 and mean_curv are filled with random values.
if flat
    a = 0;  b = 0;  Q = eye(3);
    u = [cos(0.7); sin(0.7); 0];
    targ = struct('r', [0;0;0], 'n', u, 'd', [-u(2); u(1); 0], 'd2', -u);
else
    a = 1;  b = 3;
    % a generic orientation, so fixed-axis weights (grad, hess) are not
    % aligned with the normal
    c1 = cos(0.7);  s1 = sin(0.7);  c2 = cos(1.1);  s2 = sin(1.1);
    Q = [1 0 0; 0 c1 -s1; 0 s1 c1]*[c2 0 s2; 0 1 0; -s2 0 c2];
    targ = struct('r', [0;0;0], 'n', Q*[0;0;1], 'mean_curv', -(a+b)/2);
end
fs = setdiff(targ_fields, fieldnames(targ));
for i = 1:numel(fs)
    targ.(fs{i}) = randn(3, 1);
end
rho = 10.^-(2:0.5:5);  th = 0.1 + (0:11)*pi/6;
[T, R] = ndgrid(th, rho);
x = R(:).'.*cos(T(:).');  y = R(:).'.*sin(T(:).');
src.r = Q*[x; y; (a*x.^2 + b*y.^2)/2];
src.n = [-a*x; -b*y; ones(size(x))];  src.n = Q*(src.n./vecnorm(src.n));
if flat, src.n = repmat([cos(1.1); sin(1.1); 0], 1, numel(x)); end

try
    A = radpp_eval_kern(src, targ, p);          % (m, n*ns)
catch ME
    error(['KERNEL3D.RADGREEN: could not evaluate the kernel on the test ' ...
           'surface to choose ipv (%s); pass ''ipv''.'], ME.message);
end
m = p.opdims(1);  n = p.opdims(2);
ne = size(p.ent, 1);
ipv = zeros(1, ne);
% largest entry at each distance: an entry below 1e-10 of it everywhere
% is zero up to rounding (e.g. the curl of a gradient)
big = max(reshape(abs(A), m*n*numel(th), numel(rho)), [], 1);
for e = 1:ne
    v = abs(A(p.ent(e,1), p.ent(e,2):n:end));
    v = reshape(v, numel(th), numel(rho));
    mx = max(v, [], 1);
    if all(mx <= 1e-10*big), continue; end
    sl = diff(log(mx + realmin))./diff(log(rho));
    q = max(0, ceil(-median(sl(end-2:end)) - 0.25));
    if q == 1 && flat
        ipv(e) = 0;
    elseif q == 1
        w = v(:,end)*rho(end);
        ipv(e) = (max(w) - min(w)) > 0.05*max(w);
    else
        ipv(e) = min(max(q - 1, 0), 2);
    end
end
end


function spec = dyad_kernel(type)
% kernels of an isotropic tensor green G_ij on Cartesian densities,
% opdims [3 3]; components G_ij then d_k G_ij (see RADPP_ROTATE)
switch lower(type)
    case {'s', 'single'}
        spec = kspec('s', 0, [3 3], {}, @(t) w_dyad(t, 's'), []);
    case {'efie', 'nxs'}
        spec = kspec('efie', 0, [3 3], {'n'}, @(t) w_dyad(t, 'efie'), []);
    case {'curl', 'curls'}
        spec = kspec('curl', 1, [3 3], {}, @(t) w_dyad(t, 'curl'), []);
    case {'mfie', 'nxcurls', 'n_cross_curl_s'}
        spec = kspec('mfie', 1, [3 3], {'n'}, @(t) w_dyad(t, 'mfie'), []);
    otherwise
        error('KERNEL3D.RADPP: unknown type ''%s'' for a tensor green.', type);
end
end


function W = w_dyad(t, type)
% (G J)_i = G_ij J_j, (n x G J)_i = eps_iab n_a G_bj,
% (curl G J)_i = eps_iab d_a G_bj, (n x curl G J)_i = eps_icd n_c eps_dab d_a G_bj
sidx = [1 2 3; 2 4 5; 3 5 6];
gidx = @(b, j, a) 6 + 3*(sidx(b,j) - 1) + a;
nt = size(t.r(:,:), 2);
nd = 0;  if any(strcmp(type, {'curl', 'mfie'})), nd = 1; end
W = zeros(nt, 6 + 18*nd, 9);
if isfield(t, 'n') || isprop(t, 'n'), n = t.n(1:3,:).'; end
for j = 1:3
    for i = 1:3
        e = i + 3*(j-1);
        switch type
            case 's'
                W(:, sidx(i,j), e) = 1;
            case 'efie'
                for a = 1:3, for b = 1:3
                    if levi(i,a,b) ~= 0
                        W(:, sidx(b,j), e) = W(:, sidx(b,j), e) + levi(i,a,b)*n(:,a);
                    end
                end, end
            case 'curl'
                for a = 1:3, for b = 1:3
                    if levi(i,a,b) ~= 0
                        W(:, gidx(b,j,a), e) = W(:, gidx(b,j,a), e) + levi(i,a,b);
                    end
                end, end
            case 'mfie'
                for c = 1:3, for d = 1:3
                    if levi(i,c,d) == 0, continue; end
                    for a = 1:3, for b = 1:3
                        if levi(d,a,b) ~= 0
                            W(:, gidx(b,j,a), e) = W(:, gidx(b,j,a), e) ...
                                + levi(i,c,d)*levi(d,a,b)*n(:,c);
                        end
                    end, end
                end, end
        end
    end
end
end


function W = w_mfie(t)
% (n x (grad G x J))_i = sum_c (n_c G_i - d_ic n.grad G) J_c, entries
% (i,c) in column-major order, weights on [G | G_x G_y G_z]
n = t.n(1:3,:).';
nt = size(n, 1);
W = zeros(nt, 4, 9);
for c = 1:3
    for i = 1:3
        e = i + 3*(c-1);
        W(:, 1+i, e) = n(:, c);
        if i == c, W(:, 2:4, e) = W(:, 2:4, e) - n; end
    end
end
end


function W = w_curls(t)
% (grad G x J)_i = sum_c eps_ibc G_b J_c, entries (i,c) column-major
nt = size(t.r(:,:), 2);
W = zeros(nt, 4, 9);
for c = 1:3
    for i = 1:3
        b = 6 - i - c;
        if i ~= c && b >= 1 && b <= 3
            W(:, 1+b, i + 3*(c-1)) = levi(i, b, c);
        end
    end
end
end


function s = levi(i, j, k)
s = (j - i)*(k - i)*(k - j)/2;
end


function W = w_grad(t, dim)
% grad_x: row j has weight 1 on G_{x_j}, (nt, 1+dim, dim)
nt = size(t.r(:,:), 2);
W = zeros(nt, 1+dim, dim);
for j = 1:dim
    W(:, 1+j, j) = 1;
end
end


function n = normals(t, dim)
% unit normals (nt,dim), padded with zeros if t.n has fewer rows
n = zeros(size(t.n(:,:), 2), dim);
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
% free plate boundary conditions
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
% targets, from the target fields n, d, d2
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
% on-axis profiles [G | G_x | G_xx G_yy | G_xxx G_xyy], (numel(r), *).
% green is called with points (3,n), or (2,n) if that fails; hess and
% third may be full, (nt,ns,dim,dim) and (nt,ns,dim,dim,dim) for dim 2
% or 3, or packed, (nt,ns,3) = xx xy yy and (nt,ns,4) = xxx xxy xyy yyy
r = r(:);
n = numel(r);
out = cell(1, nderiv+1);
try
    [out{:}] = green(zeros(3,1), [r.'; zeros(2, n)]);
catch
    [out{:}] = green(zeros(2,1), [r.'; zeros(1, n)]);
end
f = out{1}(:);
if nderiv >= 1, f = [f, out{2}(:,1,1)]; end
if nderiv >= 2
    if size(out{3}, 4) > 1                        % full (nt,ns,dim,dim)
        f = [f, out{3}(:,1,1,1), out{3}(:,1,2,2)];
    else                                          % packed (nt,ns,3)
        f = [f, out{3}(:,1,1), out{3}(:,1,3)];
    end
end
if nderiv >= 3
    if size(out{4}, 5) > 1                        % full (nt,ns,dim,dim,dim)
        f = [f, out{4}(:,1,1,1,1), out{4}(:,1,1,2,2)];
    else                                          % packed (nt,ns,4)
        f = [f, out{4}(:,1,1), out{4}(:,1,3)];
    end
end
end


function submat = radpp_eval_kern(srcinfo, targinfo, p)
% (m*nt, n*ns), block entries interleaved
src  = srcinfo.r(:,:);
targ = targinfo.r(:,:);
ns = size(src, 2);  nt = size(targ, 2);

dx = targ(1,:).' - src(1,:);
dy = targ(2,:).' - src(2,:);
dz = 0;
if size(src,1) > 2 && size(targ,1) > 2
    dz = targ(3,:).' - src(3,:);
end
dz = dz + zeros(nt, ns);

nfit = p.nderiv + p.isrc;
if p.gmode
    % exact profiles from the green
    r = sqrt(dx(:).^2 + dy(:).^2 + dz(:).^2);
    comps = radpp_rotate(dx, dy, dz, p.profun(r), nfit, p.dim);
else
    % the fit
    if p.usemex
        comps = radpp_eval_derivs(dx, dy, dz, p.ipars, p.dpars);
    else
        r = sqrt(dx(:).^2 + dy(:).^2 + dz(:).^2);
        prof = radpp_eval(r, p.breaks, p.coefs);
        comps = radpp_rotate(dx, dy, dz, prof, nfit, p.dim);
    end
end
comps = reshape(comps, nt, ns, []);

W = p.wfun(targinfo);
nrow = size(W, 3);
m = p.opdims(1);  n = p.opdims(2);

if p.isrc
    % d_{n_y} D_c G = -sum_i n_y,i D_{c+i} G
    up = radpp_up(p.dim);
    ny = zeros(p.dim, ns);
    mn = min(p.dim, size(srcinfo.n, 1));
    ny(1:mn, :) = srcinfo.n(1:mn, :);
    dcomps = zeros(nt, ns, p.ncomp);
    for c = 1:p.ncomp
        for i = 1:p.dim
            dcomps(:,:,c) = dcomps(:,:,c) - ny(i,:).*comps(:,:,up(c,i));
        end
    end
    V = p.vfun(targinfo);
end

submat = zeros(m*nt, n*ns);
for j = 1:nrow
    Kj = 0;
    for c = 1:p.ncomp
        Kj = Kj + W(:,c,j).*comps(:,:,c);
        if p.isrc
            Kj = Kj + V(:,c,j).*dcomps(:,:,c);
        end
    end
    ie = p.ent(j,1);  je = p.ent(j,2);
    submat(ie:m:end, je:n:end) = submat(ie:m:end, je:n:end) + Kj;
end
end


function up = radpp_up(dim)
% up(c,i): index of D_c G differentiated once more along x_i, for the
% components of order 0...2 in the order of RADPP_ROTATE
idx = {zeros(1,0)};
for ord = 1:3
    for c = 1:numel(idx)
        if numel(idx{c}) ~= ord-1, continue; end
        last = 1;
        if ord > 1, last = idx{c}(end); end
        for i = last:dim
            idx{end+1} = [idx{c} i]; 
        end
    end
end
key = cellfun(@(v) sprintf('%d', v), idx, 'UniformOutput', false);
nlow = 1 + dim + dim*(dim+1)/2;
up = zeros(nlow, dim);
for c = 1:nlow
    for i = 1:dim
        up(c,i) = find(strcmp(key, sprintf('%d', sort([idx{c} i]))));
    end
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

targs = radpp_targs(targinfo, p.wfun, p.vfun);
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

if ifcomplex
    nq = @getnearquad_radpp;
else
    nq = @getnearquad_radpp_d;
end
wnear = nq(S.npatches, norders, ixyzs, iptype, S.npts, ...
    srccoefs, srcvals, targs, patch_id, uvs_targ, eps, p.ipv, rsc.nnz, ...
    rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, nquad, p.dpars, ...
    p.ipars);

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end


function spmat = radgreen_getquad(S, eps, p, ifcomplex, ri, targinfo, opts)
% fit the profiles over the distances this near quadrature needs, then
% use the radpp quadrature

if nargin < 6 || isempty(targinfo), targinfo = S;        end
if nargin < 7 || isempty(opts),     opts     = struct(); end
if ~isfield(opts, 'rsc') || isempty(opts.rsc)
    opts.rsc = getnear(S, targinfo);
end

p = radgreen_fit(p, S, targinfo, opts.rsc);

spmat = radpp_getquad(S, eps, p, ifcomplex, ri, targinfo, opts);
end


function p = radgreen_fit(p, S, targinfo, rsc)
% params p of a radgreen kernel with ipars, dpars fitted over the near
% field distances of S at targinfo
if nargin < 4 || isempty(rsc), rsc = getnear(S, targinfo); end
[rmin, rmax] = radgreen_range(S, targinfo, rsc);
G = p.green;
if numel(G.parts) == 1
    [ipars, dpars] = radpp_fit(p.profun, rmin, rmax, p.fitopts);
else
    % a sum of greens: fit each (their singular parts may cancel only up
    % to rounding) and add the expansions
    nfit = p.nderiv + p.isrc;
    fits = cell(1, numel(G.parts));
    for k = 1:numel(G.parts)
        [~, ~, fits{k}] = radpp_fit(@(r) G.c(k)*G.parts{k}.profiles(r, nfit), ...
                                    rmin, rmax, p.fitopts);
    end
    [ipars, dpars] = radpp_addfits(fits);
end
p.ipars = [ipars; p.ipost];
p.dpars = dpars;
end


function [ipars, dpars, info] = radpp_addfits(fits)
% sum of radpp_fit expansions on the union of their panels, on the
% distances common to them, in the vpp format
breaks = unique(cell2mat(cellfun(@(f) f.breaks(:), fits(:), ...
                                 'UniformOutput', false)));
lo = max(cellfun(@(f) f.breaks(1), fits));
hi = min(cellfun(@(f) f.breaks(end), fits));
breaks = unique([lo; breaks(breaks > lo & breaks < hi); hi]);
[~, n, m] = size(fits{1}.coefs);
a = breaks(1:end-1);  nbin = numel(a);
mid = (a + breaks(2:end))/2;
cf = zeros(nbin, n, m);
for k = 1:numel(fits)
    ib = discretize(mid, fits{k}.breaks);
    d = fits{k}.coefs(ib, end:-1:1, :);       % d(:,p+1,:): (r-b)^p
    del = a - fits{k}.breaks(ib);             % re-centre at a: r-b = (r-a)+del
    e = zeros(size(d));
    for q = 0:n-1
        for pp = q:n-1
            e(:,q+1,:) = e(:,q+1,:) + nchoosek(pp, q)*del.^(pp-q).*d(:,pp+1,:);
        end
    end
    cf = cf + e(:, end:-1:1, :);
end
cp = permute(cf, [3 2 1]);                    % (m,n,nbin)
cfv = zeros(2*m, n, nbin);
cfv(1:2:end,:,:) = real(cp);
cfv(2:2:end,:,:) = imag(cp);
ipars = [nbin; n; 1; nbin+2; 2*m];
dpars = [breaks; cfv(:)];
info = struct('breaks', breaks, 'coefs', cf);
end


function tf = summable(Ks)
% Compatible fitted kernels, with identity terms ignored.
tf = ~isempty(Ks) && all(cellfun(@(K) isequal(K.opdims,Ks{1}.opdims),Ks));
if ~tf, return; end
Ks = radpp_sum_parts(Ks);
if isempty(Ks), return; end
tf = all(cellfun(@(K) strcmp(K.name, 'radpp'), Ks));
if ~tf, return; end
p = Ks{1}.params;
if numel(Ks) == 1
    tf = ~p.gmode;
    return
end
same = @(K) ~K.params.gmode && ~strcmp(K.params.spec.name, 'custom') && ...
    strcmp(K.params.spec.name, p.spec.name) && ...
    isequal(K.params.kcoefs, p.kcoefs) && isequal(K.params.nu, p.nu) && ...
    K.params.dim == p.dim && isequal(K.params.ipv, p.ipv) && ...
    size(K.params.coefs, 2) == size(p.coefs, 2) && ...
    size(K.params.coefs, 3) == size(p.coefs, 3);
tf = all(cellfun(same, Ks)) && ...
    max(cellfun(@(K) K.params.breaks(1), Ks)) < ...
    min(cellfun(@(K) K.params.breaks(end), Ks));
end


function K = ppsum(Ks, cs)
% sum_k cs(k)*Ks{k}: a sum is refitted on the green sum_k cs(k)*green_k
% over the distances common to the fits; a multiple scales the fit
opdims = Ks{1}.opdims;
[Ks,cs,diagonal] = radpp_sum_parts(Ks,cs);
if isempty(Ks)
    K = kernel3d.zeros(opdims);
    if ~isempty(diagonal), K = kernel3d.eye(zeros(opdims)); end
    K.diag = diagonal;
    return
elseif numel(Ks) == 1 && cs == 1
    K = Ks{1};
    K.diag = diagonal;
    return
end
p = Ks{1}.params;
greens = cellfun(@(K) K.params.green, Ks, 'UniformOutput', false);
% the fmm of the sum of the greens
fmm = {};
if all(cellfun(@(Kk) ~isempty(Kk.params.fmm), Ks))
    fmms = cellfun(@(Kk) Kk.params.fmm, Ks, 'UniformOutput', false);
    fmm = {'fmm', radpp_sum(cs, fmms{:})};
end
args = {p.spec, 'dim', p.dim, 'coefs', p.kcoefs, 'nu', p.nu, ...
        'opts', p.fitopts, 'zk', max(cellfun(@(K) K.zk, Ks)), fmm{:}, ...
        'ifcomplex', any(cellfun(@(K) logical(K.ifcomplex), Ks)) || ~isreal(cs)};
if numel(Ks) == 1
    info = struct('breaks', p.breaks, 'coefs', cs*p.coefs, 'err', abs(cs)*p.err);
    [info.ipars, info.dpars] = radpp_addfits({info});
    K = kernel3d.radpp(radpp_sum(cs, greens{:}), info.breaks([1 end]), args{:}, ...
            'ipv', p.ipv, 'kernel_order', p.kernel_order, 'fit', info);
else
    rmin = max(cellfun(@(K) K.params.breaks(1), Ks));
    rmax = min(cellfun(@(K) K.params.breaks(end), Ks));
    K = kernel3d.radpp(radpp_sum(cs, greens{:}), [rmin rmax], args{:});
end
% otherwise the sum of the fmms of the kernels, if they have them
if isempty(fmm) && all(cellfun(@(Kk) isa(Kk.fmm, 'function_handle'), Ks))
    fmms = cellfun(@(Kk) Kk.fmm, Ks, 'UniformOutput', false);
    K.fmm = radpp_sum(cs, fmms{:});
end
K.diag = diagonal;
end


function [rmin, rmax] = radgreen_range(S, targinfo, rsc)
% distances from each target to its near patches: rmax bounds them with
% the patch radii, rmin is a small fraction of the smallest patch
c = S.cms;
R = reshape(S.rads, 1, []);

tr = zeros(3, size(targinfo.r(:,:), 2));
tr(1:size(targinfo.r,1), :) = targinfo.r(:,:);
rp  = double(rsc.row_ptr(:));
ci  = double(rsc.col_ind(:));
tid = reshape(repelem((1:numel(rp)-1).', diff(rp)), [], 1);

rmin = 1e-10*min(R);
rmax = 1;
if ~isempty(ci)
    d = vecnorm(tr(:,tid) - c(:,ci)) + reshape(R(ci), 1, []);
    rmax = 1.1*max(d);
end
rmax = max(rmax, 10*rmin);
end


function targs = radpp_targs(targinfo, wfun, vfun)
% target array for radpp_kern: locations, then the complex weights w (and
% v) of each row as interlaced real/imaginary parts,
% (3 + 2*ncomp*(1+isrc)*nrow, nt)
r = targinfo.r(:,:);
ntarg = size(r, 2);
xyz = zeros(3, ntarg);
xyz(1:size(r,1), :) = r;

W = wfun(targinfo);                            % (nt, ncomp, nrow)
if ~isempty(vfun)
    W = cat(2, W, vfun(targinfo));             % (nt, 2*ncomp, nrow)
end
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
