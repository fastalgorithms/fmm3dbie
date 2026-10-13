function obj = kernel_eval(type, nu, rts, ejs, opts)
%SURFWAVE.FLEX.KERNEL_EVAL  Chunkie kernel object (with splitinfo) for the
%   surfwave.flex boundary-to-volume free-plate evaluation kernels.
%
%   obj = surfwave.flex.kernel_eval('free_plate_gs_eval', nu, rts, ejs)
%   obj = surfwave.flex.kernel_eval('free_plate_gphi_eval', nu, rts, ejs)
%
%   [1 3] kernel, columns [G_ny, (1+nu)/2*G_tauy, G]. By the moment
%   conditions on (rts,ejs) (sum ejs.*rts.^q = 0 for q < 4, = 1/alpha at
%   q = 4). Only 'free_plate_gs_eval' gets a splitinfo.
%
%   obj = surfwave.flex.kernel_eval(type, nu, rts, ejs, opts)
%      opts.rmax    - build a radial Chebyshev representation of G and its
%                     gradient (surfwave.radcheb_green) and evaluate through
%                     it; rmax is the largest source-target distance (or
%                     [rmin rmax]). No splitinfo is attached in that case.
%      opts.radcheb - options passed to kernel3d.radcheb (eps default 1e-13)
%      opts.R       - alternatively a prebuilt surfwave.radcheb_green struct

if ~any(strcmpi(type, {'free_plate_gs_eval', 'free_plate_gphi_eval'}))
    error('surfwave:flex:kernel_eval', 'unknown type ''%s''', type);
end

obj = kernel();
obj.name       = 'flexural';
obj.type       = type;
obj.opdims     = [1 3];
if nargin < 5, opts = []; end
R = [];
if isfield(opts, 'R') && ~isempty(opts.R)
    R = opts.R;
elseif isfield(opts, 'rmax') && ~isempty(opts.rmax)
    ro = struct('eps', 1e-13);
    if isfield(opts, 'radcheb'), ro = opts.radcheb; end
    if strcmpi(type, 'free_plate_gs_eval')
        vg = @(s,t) surfwave.flex.gsflex(rts, ejs, s, t);
    else
        vg = @(s,t) surfwave.flex.gphiflex(rts, ejs, s, t);
    end
    R = surfwave.radcheb_green(vg, opts.rmax, ro);
end
obj.params.R   = R;
obj.eval       = @(s,t) surfwave.flex.kern(s, t, type, nu, rts, ejs, R);
% obj.src_fields = {'n', 'd'};

if strcmpi(type, 'free_plate_gs_eval') && isempty(R)
    alpha = 1/sum(ejs(:).*rts(:).^4);
    obj.sing = 'log';
    obj.splitinfo = [];
    obj.splitinfo.type     = {[0 0 0 0], [1 0 0 0]};
    obj.splitinfo.action   = {'r', 'r'};
    obj.splitinfo.functions = @(s, t) gs_eval_split(s, t, nu, rts, ejs, alpha);
elseif strcmpi(type, 'free_plate_gs_eval')
    obj.sing = 'log';          % radcheb evaluation: no splitinfo
else
    obj.sing = 'smooth';
end

end


function f = gs_eval_split(s, t, nu, rts, ejs, alpha)
%GS_EVAL_SPLIT  Log split for free_plate_gs_eval, f{1} = K + log(r)*B/(2*pi),
%   B fit directly to the r -> 0 asymptotics of surfwave.flex.gsflex.

K = surfwave.flex.kern(s, t, 'free_plate_gs_eval', nu, rts, ejs);

ns = size(s.r, 2);  nt = size(t.r, 2);

srcnorm = s.n;  srctang = s.d;
nx = repmat(srcnorm(1,:), nt, 1);  ny = repmat(srcnorm(2,:), nt, 1);
dx = repmat(srctang(1,:), nt, 1);  dy = repmat(srctang(2,:), nt, 1);
ds = sqrt(dx.*dx + dy.*dy);
taux = dx./ds;  tauy = dy./ds;

xt = repmat(t.r(1,:).', 1, ns);  yt = repmat(t.r(2,:).', 1, ns);
xs = repmat(s.r(1,:), nt, 1);    ys = repmat(s.r(2,:), nt, 1);
ux = xt - xs;  uy = yt - ys;
rr2 = ux.*ux + uy.*uy;

Bn   = 2*(ux.*nx + uy.*ny)/alpha;
Btau = (1+nu)/2*(2*(ux.*taux + uy.*tauy)/alpha);
Bval = -rr2/alpha;

B = zeros(nt, 3*ns);
B(:,1:3:end) = Bn;
B(:,2:3:end) = Btau;
B(:,3:3:end) = Bval;

rr = repelem(sqrt(rr2), 1, 3);

f = cell(2,1);
f{1} = K + log(rr).*B/(2*pi);
f{2} = B;

end
