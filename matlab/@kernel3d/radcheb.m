function obj = radcheb(fkern,rmax,opts)
%KERNEL3D.RADCHEB   Adaptive Chebyshev interpolant of a radial kernel in 3D.
%
%   KERNEL3D.RADCHEB(fkern, rmax) or KERNEL3D.RADCHEB(fkern, rmax, opts)
%   builds an adaptive piecewise Chebyshev interpolant of a kernel which
%   depends only on r = |x-y|.
%
%   FKERN should be a function handle with the usual kernel calling
%   sequence fkern(srcinfo,targinfo), or a kernel3d object. It must be
%   scalar valued, and may not depend on any source or target fields other
%   than the positions.
%
%   RMAX is the largest radius of interest, or [rmin, rmax]. If rmin is
%   not given it is taken from opts.rmin, which defaults to 1e-12.
%
%   OPTS is an options struct:
%      opts.rmin    - smallest radius of interest (default 1e-12)
%      opts.norder  - terms per panel (default 16)
%      opts.eps     - build tolerance (default 1e-12)
%      opts.nlevmax - max subdivision levels (default 50)
%      opts.maxsub  - max number of panels (default 10000)
%      opts.kernel_order - kernel order for oversampling. Inherited from
%                     fkern if it is a kernel3d, otherwise -1
%      opts.zk      - wavenumber for oversampling purposes. Inherited from
%                     fkern if it is a kernel3d, otherwise 0
%
% See also RADCHEB_FIT, RADCHEB_EVAL

if nargin < 1
    error('KERNEL3D.RADCHEB: missing kernel.');
end
if nargin < 2 || isempty(rmax)
    error('KERNEL3D.RADCHEB: missing rmax.');
end
if nargin < 3 || isempty(opts), opts = struct(); end

rmin = 1e-12;
if isfield(opts,'rmin'), rmin = opts.rmin; end
if numel(rmax) == 2
    rmin = rmax(1);
    rmax = rmax(2);
end

kernel_order = -1;
zk = 0;

if isa(fkern,'kernel3d')
    if ~isequal(fkern.opdims,[1 1])
        error('KERNEL3D.RADCHEB: only scalar kernels are supported.');
    end
    if ~isempty(fkern.src_fields) || ~isempty(fkern.targ_fields)
        error(['KERNEL3D.RADCHEB: kernel depends on source or target ' ...
               'fields, so it is not a function of r alone.']);
    end
    feval = fkern.eval;
    if ~isempty(fkern.kernel_order), kernel_order = fkern.kernel_order; end
    if ~isempty(fkern.zk), zk = fkern.zk; end
elseif isa(fkern,'kernel')
    if ~isequal(fkern.opdims,[1 1])
        error('KERNEL3D.RADCHEB: only scalar kernels are supported.');
    end
    feval = fkern.eval;
elseif isa(fkern,'function_handle')
    feval = fkern;
else
    error('KERNEL3D.RADCHEB: kernel not of a supported type.');
end

if isfield(opts,'kernel_order'), kernel_order = opts.kernel_order; end
if isfield(opts,'zk'), zk = opts.zk; end

f = @(r) radcheb_probe(feval,r);

fv = f([rmin; rmax]);
if numel(fv) ~= 2
    error('KERNEL3D.RADCHEB: only scalar kernels are supported.');
end

[ipars,dpars,info] = radcheb_fit(f,rmin,rmax,opts);

obj           = kernel3d();
obj.name      = 'radcheb';
obj.type      = 'radcheb';
obj.opdims    = [1 1];
obj.zk        = zk;
obj.ifcomplex = ~isreal(info.coefs);
obj.kernel_order = kernel_order;

obj.params.fkern  = fkern;
obj.params.rmin   = rmin;
obj.params.rmax   = rmax;
obj.params.opts   = opts;
obj.params.ipars  = ipars;
obj.params.dpars  = dpars;
obj.params.breaks = info.breaks;
obj.params.coefs  = info.coefs;
obj.params.nbin   = info.nbin;
obj.params.err    = info.err;

breaks = info.breaks;
coefs  = info.coefs;
obj.eval = @(s,t) radcheb_kern(s,t,breaks,coefs);
obj.fmm  = [];

ifcomplex = obj.ifcomplex;
obj.getquad = @(S,eps,varargin) radcheb_getquad(S, eps, ipars, dpars, ...
                  ifcomplex, kernel3d.rsc_interleave_full(1,1), varargin{:});

obj.get_overs_orders = @(S,t,eps) kernel3d.kernel3d_getnear_overs(S, t, eps, ...
                                     obj.zk, obj.kernel_order);

end


function val = radcheb_probe(feval,r)
%  evaluate the underlying kernel at radii r, with the source at the
%  origin and the targets along the x axis

r = r(:);
s = []; s.r = zeros(3,1);
t = []; t.r = [r.'; zeros(2,numel(r))];
val = feval(s,t);
val = val(:);

end


function val = radcheb_kern(srcinfo,targinfo,breaks,coefs)

dim = size(srcinfo.r,1);
src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

r = 0;
for i = 1:dim
    r = r + (targ(i,:).' - src(i,:)).^2;
end
r = sqrt(r);

val = radcheb_eval(r,breaks,coefs);

end


function spmat = radcheb_getquad(S,eps,ipars,dpars,ifcomplex,ri,targinfo,opts)
%  near quadrature correction for a radcheb kernel, obtained by handing
%  the precomputed vpp expansion to the fortran quadrature routine

if nargin < 7 || isempty(targinfo), targinfo = S;        end
if nargin < 8 || isempty(opts),     opts     = struct(); end

[srcvals, srccoefs, norders, ixyzs, iptype, ~] = extract_arrays(S);
npatches = S.npatches;
npts     = S.npts;

if ( isfield(opts, 'rsc') && ~isempty(opts.rsc) )
    rsc = opts.rsc;
else
    rsc = getnear(S, targinfo);
end

ntarg = size(extract_targ_array(targinfo), 2);
if isfield(targinfo, 'patch_id') || isprop(targinfo, 'patch_id')
    patch_id = targinfo.patch_id;
    uvs_targ = targinfo.uvs_targ;
else
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end
if isempty(patch_id)
    patch_id = -ones(ntarg, 1);
    uvs_targ = zeros(2, ntarg);
end

nquad = rsc.iquad(end) - 1;

wnear = getnearquad_radcheb(npatches, norders, ixyzs, iptype, npts, ...
    srccoefs, srcvals, targinfo, patch_id, uvs_targ, eps, 1, rsc.nnz, ...
    rsc.row_ptr, rsc.col_ind, rsc.iquad, rsc.rfac0, nquad, dpars, ipars, ...
    ifcomplex);

spmat = conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, wnear, ri);

end
