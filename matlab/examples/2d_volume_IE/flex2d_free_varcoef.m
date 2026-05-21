%
%  flex2d_free_varcoef
%
%  Solves the variable-coefficient flexural wave equation
%
%    a(r) \Delta^2 u + b(r) \Delta u + c(r) u = f    in \Omega
%
%  with free plate boundary conditions on \partial\Omega,
%  where the coefficients a, b, c vary spatially (variable thickness).
%
%  All four operator blocks are assembled as dense matrices using
%  near-correction + smooth oversampled quadrature. The coupled
%  volume/boundary system is solved with gmres.
%

%% Geometry

S = geometries.disk([],[],[4 4 4],6);

[~,pcoefs] = bump4(S);

a = pcoefs{1};
b = 0; % can only handle one wavenumber (so far)
c = pcoefs{3};
nu = pcoefs{end};

zk1 = sqrt((-b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((-b - sqrt(b^2 + 4*a*c)) / (2*a));
zpars_f = [zk1, zk2];

cparams = []; cparams.ta = pi/16; cparams.tb = cparams.ta + 2*pi;
chnkr = chunkerfuncuni(@(t) ellipse(t),16,cparams);
chnkr = sort(chnkr);

bcoefs = bump4(chnkr,1);

chnkr = makedatarows(chnkr,3);
chnkr.data(1,:) = bcoefs{1};
chnkr.data(2,:) = bcoefs{2};
chnkr.data(3,:) = bcoefs{3};

kappa = signed_curvature(chnkr);
kappa = kappa(:);

%% chnkrdim: 3-D lift of the 2-D chunker

chnkrdim = [];
loc_field = fieldnames(chnkr)';
for i = 1:length(loc_field)
    if mod(numel(chnkr.(loc_field{i})), chnkr.npt) == 0
        chnkrdim.(loc_field{i}) = chnkr.(loc_field{i})(:,:);
    end
end
chnkrdim.r   = chnkrdim.r;
chnkrdim.d   = chnkrdim.d; 
chnkrdim.d2  = chnkrdim.d2;
chnkrdim.n   = chnkrdim.n;
chnkrdim.wts = chnkrdim.wts(:);

[~, patch_inds, uvsloc] = get_closest_pts(S, chnkrdim);
chnkrdim.patch_id = patch_inds;
chnkrdim.uvs_targ = uvsloc;
chnkrdim.kappa = kappa;

%% Target info for v2v (geom nodes -> geom nodes)

[icoefs,pcoefs,H,alpha] = bump4(S);

targ_geom = [];
targ_geom.r        = S.r;
targ_geom.patch_id = S.patch_id;
targ_geom.uvs_targ = S.uvs_targ;
targ_geom.icoefs   = icoefs;

%% Kernels

dpars   = nu;
skern   = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 's',              nu);
v2vkern = @(s,t) v2v_kern(zk1,s,t);
repkern = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 'free_plate_eval',nu);
rhskern = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_bcs_var',  nu);
fkern1  = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_var',      nu);
double  = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

%% Quadrature

opts  = []; opts.sing  = 'log';
opts2 = []; opts2.sing = 'pv';

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double,  opts);
H = chunkermat(chnkr, hilbert, opts2);

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H  + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H  ...
    + 2*((1+nu)/2)^2*(bcoefs{2}(:)./bcoefs{1}(:)).*D*D - diag(bcoefs{2}(:)./bcoefs{1}(:)).*(1/8)*(1+nu).^2 + 1/2*diag(bcoefs{2}(:)./bcoefs{1}(:));
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

sysmat(1:2:end,1:2:end) = sysmat(1:2:end,1:2:end) - (1-nu)*(sysmat(1:2:end,2:2:end).*(bcoefs{3}./bcoefs{1}))*H;
sysmat(2:2:end,1:2:end) = sysmat(2:2:end,1:2:end) - (1-nu)*(sysmat(2:2:end,2:2:end).*(bcoefs{3}./bcoefs{1}))*H;

Djump = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];
Djump = kron(eye(chnkr.npt), Djump);

sysb2b = Djump + sysmat;
fprintf('%5.2e s : time to assemble b2b matrix\n', toc(start))

densmap = zeros(3*chnkr.npt, 2*chnkr.npt);
densmap(1:3:end, 1:2:end) = eye(chnkr.npt);
densmap(3:3:end, 2:2:end) = eye(chnkr.npt);
densmap(2:3:end, 1:2:end) = -H;

eps = 1e-12;
Q = getnear(S); Q.targinfo = S; Q.wavenumber = zk1; Q.kernel_order = 0;
novers   = get_oversampling_parameters(S, Q, eps);
[geom_over, xinterp] = S.oversample(novers);
nearquad_opts = []; nearquad_opts.nover = novers;

getnearquad_v2v = @(varargin) flex2d.getnearquad_var(varargin{1:17}, dpars, zpars_f, varargin{18}, 'v2v');
getnearquad_s   = @(varargin) flex2d.getnearquad(    varargin{1:17}, dpars, zpars_f, varargin{18}, 'v2v');
getnearquad_v2b = @(varargin) flex2d.getnearquad(    varargin{1:17}, dpars, zpars_f, varargin{18}, 'free_var');

start = tic;
v2v_near = getnearquad_kern(S, v2vkern, eps, getnearquad_v2v, targ_geom, nearquad_opts);
v2vkern_eval = kernel(v2vkern);
v2v_smth = (v2vkern_eval.eval(geom_over, targ_geom) .* geom_over.wts(:).') * xinterp;
flex_v2v = v2v_near + v2v_smth;
fprintf('%5.2e s : time to assemble v2v matrix\n', toc(start))

start = tic;
v2b_near = getnearquad_kern(S, rhskern, eps, getnearquad_v2b, chnkrdim, nearquad_opts);
rhskern_eval = kernel(rhskern);
v2b_smth = (rhskern_eval.eval(geom_over, chnkrdim) .* geom_over.wts(:).') * xinterp;
flex_v2b = v2b_near + v2b_smth;
fprintf('%5.2e s : time to assemble v2b matrix\n', toc(start))

start = tic;
opts_cor = []; opts_cor.corrections = 0;

repkerns = cell(7,1);
repkerns{1} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapgx', nu);
repkerns{2} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapgy', nu);
repkerns{3} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapg',  nu);
repkerns{4} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gyy',   nu);
repkerns{5} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gxx',   nu);
repkerns{6} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gxy',   nu);
repkerns{7} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval',       nu);

flex_b2v = 0;
for ii = 1:7
    flex_b2v = flex_b2v + icoefs(ii,:).' .* chunkerkernevalmat(chnkr, repkerns{ii}, S.r(1:2,:), opts_cor) * densmap;
end
fprintf('%5.2e s : time to assemble b2v matrix\n', toc(start))

%% Assemble full system
%   [  a*I + flex_v2v ,   flex_b2v  ]   [ eta  ]   [ rhs_vol ]
%   [    flex_v2b     ,   sysb2b    ] * [ zeta ] = [ rhs_bc  ]

nvol  = S.npts;
nbdry = 2 * chnkr.npt;

sysmat = [a*eye(nvol) + flex_v2v,  flex_b2v;
          flex_v2b,                 sysb2b  ];

%% Right-hand side (manufactured solution u = sin(x)*sin(y))

rhs_vol = (4*a + 2*b - c) .* sin(S.r(1,:)) .* sin(S.r(2,:)) ...
    - icoefs(1,:) .* 2 .* cos(S.r(1,:)) .* sin(S.r(2,:))   ...
    - icoefs(2,:) .* 2 .* sin(S.r(1,:)) .* cos(S.r(2,:))   ...
    - icoefs(3,:) .* 2 .* sin(S.r(1,:)) .* sin(S.r(2,:))   ...
    - icoefs(4,:) .*      sin(S.r(1,:)) .* sin(S.r(2,:))   ...
    - icoefs(5,:) .*      sin(S.r(1,:)) .* sin(S.r(2,:))   ...
    + icoefs(6,:) .*      cos(S.r(1,:)) .* cos(S.r(2,:))   ...
    + icoefs(7,:) .*      sin(S.r(1,:)) .* sin(S.r(2,:));

hess = zeros(chnkr.npt, 3);
hess(:,1) = -sin(chnkr.r(1,:)) .* sin(chnkr.r(2,:));
hess(:,2) =  cos(chnkr.r(1,:)) .* cos(chnkr.r(2,:));
hess(:,3) = -sin(chnkr.r(1,:)) .* sin(chnkr.r(2,:));

third = zeros(chnkr.npt, 4);
third(:,1) = -cos(chnkr.r(1,:)) .* sin(chnkr.r(2,:));
third(:,2) = -sin(chnkr.r(1,:)) .* cos(chnkr.r(2,:));
third(:,3) = -cos(chnkr.r(1,:)) .* sin(chnkr.r(2,:));
third(:,4) = -sin(chnkr.r(1,:)) .* cos(chnkr.r(2,:));

nx   = chnkr.n(1,:).';
ny   = chnkr.n(2,:).';
dx   = chnkr.d(1,:).';
dy   = chnkr.d(2,:).';
ds   = sqrt(dx.*dx + dy.*dy);
taux = dx ./ ds;
tauy = dy ./ ds;

firstbc = (hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
    nu .* (hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy));

secondbc = (third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) + ...
    third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny)) + ...
    (2-nu).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
    third(:,3).*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + third(:,4).*(tauy.*tauy.*ny)) + ...
    (1-nu).*kappa.*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy - ...
    hess(:,1).*nx.*nx - hess(:,2).*(2*nx.*ny) - hess(:,3).*ny.*ny) + ...
    (bcoefs{2}.'./bcoefs{1}.').*firstbc + ...
    2*(1-nu)*(bcoefs{3}.'./bcoefs{1}.').*(nx.*taux.*hess(:,1) + (nx.*tauy+ny.*taux).*hess(:,2) + ny.*tauy.*hess(:,3));

rhs_bc = zeros(2*chnkr.npt, 1);
rhs_bc(1:2:end) = firstbc;
rhs_bc(2:2:end) = secondbc;

rhs = [rhs_vol(:); rhs_bc];

start = tic;
sol = gmres(sysmat, rhs, [], 1e-10, size(sysmat,1));
fprintf('%5.2e s : time for gmres\n', toc(start))

%% Evaluate solution

[v2v_s_near, norderup] = getnearquad_kern(S, skern, eps, getnearquad_s, targ_geom, novers);

skern_eval = kernel(skern);
v2v_s_smth = (skern_eval.eval(geom_over, targ_geom) .* geom_over.wts(:).') * xinterp;

flex_s_v2v = v2v_s_near + v2v_s_smth;

opts_cor = []; opts_cor.corrections = 0;
flex_s_b2v = chunkerkernevalmat(chnkr, repkern, S.r(1:2,:), opts_cor) * densmap;

u     = [flex_s_v2v, flex_s_b2v] * sol;
ref_u = sin(S.r(1,:)) .* sin(S.r(2,:));
err   = abs(u - ref_u(:)) / max(abs(ref_u(:)));

fprintf('max relative error in u: %.4e\n', max(err));

figure(1); clf
plot(S, log10(S.patch_max(err)));
colorbar; title('log_{10} pointwise error in u');

%% Local functions

function out = v2v_kern(zk,srcinfo,targinfo)

src  = srcinfo.r;
targ = targinfo.r;

[val,~,hess,third] = chnk.flex2d.hkdiffgreen(zk,src,targ);
val   = 1/(2*zk^2)*val;
hess  = 1/(2*zk^2)*hess;
third = 1/(2*zk^2)*third;

hessxx = hess(:,:,1);
hessxy = hess(:,:,2);
hessyy = hess(:,:,3);

lap      = hessxx + hessyy;
gradlapx = third(:,:,1) + third(:,:,3);
gradlapy = third(:,:,2) + third(:,:,4);

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);
xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt - xs;
dy = yt - ys;
r2 = dx.*dx + dy.*dy;

eulgam = 0.5772156649015328606065120900824;

val(r2<1e-16)    = 1/(2*zk^2)/(2*pi)*(log(2/zk)-log(2/(1i*zk)));
hessxx(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
hessxy(r2<1e-16) = 0;
hessyy(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
lap(r2<1e-16)    = 2/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
gradlapx(r2<1e-16) = 0;
gradlapy(r2<1e-16) = 0;

out = targinfo.icoefs(1,:).'.*gradlapx ...
    + targinfo.icoefs(2,:).'.*gradlapy ...
    + targinfo.icoefs(3,:).'.*lap      ...
    + targinfo.icoefs(4,:).'.*hessyy   ...
    + targinfo.icoefs(5,:).'.*hessxx   ...
    + targinfo.icoefs(6,:).'.*hessxy   ...
    + targinfo.icoefs(7,:).'.*val;

end
