% flex2d_free_var.m
%
% Variable thickness version of flex2d_free_surferfun_dense.
% 
%
% Builds all four operator blocks (v2v, v2b, b2v, b2b) as explicit dense
% matrices using the same near-correction + smooth-oversampled-part
% strategy as test_surfer_vs_surferfun, then solves the coupled
% volume/boundary system with backslash.
%
% Blocks:
%   v2v: geom nodes -> geom nodes   (skern,      wrap_getnearquad 'v2v')
%   v2b: geom nodes -> chnkr points (rhskern,    wrap_getnearquad 'free')
%   b2v: chnkr pts  -> geom nodes   (repkern,    chunkerkernevalmat)
%   b2b: chnkr pts  -> chnkr points (biekern,    chunkermat + diagval)

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

% %% Build surferfun geometry
% [bdist, normals, align_type, bdry_ids, ibdry_type] = get_boundary_patches(S, chnkr);
% 
% iptype_sfun = S.iptype;
% % iptype_sfun(bdry_ids) = ibdry_type;
% 
% S_rot = rotate_patches(S, bdist, normals, align_type);
% geom  = surferfun(S_rot, iptype_sfun);
geom = S;
%% Build chnkrdim
chnkrdim = [];
loc_field = fieldnames(chnkr)';
for i = 1:length(loc_field)
    if mod(numel(chnkr.(loc_field{i})), chnkr.npt) == 0
        chnkrdim.(loc_field{i}) = chnkr.(loc_field{i})(:,:);
    end
end
chnkrdim.r   = [chnkrdim.r;  zeros(3-chnkr.dim, chnkr.npt)];
chnkrdim.d   = [chnkrdim.d;  zeros(3-chnkr.dim, chnkr.npt)];
chnkrdim.d2  = [chnkrdim.d2; zeros(3-chnkr.dim, chnkr.npt)];
chnkrdim.du  = chnkrdim.d;
chnkrdim.n   = [chnkrdim.n;  zeros(3-chnkr.dim, chnkr.npt)];
chnkrdim.wts = chnkrdim.wts(:);

[~, patch_inds, uvsloc] = get_closest_pts(geom, chnkrdim);
chnkrdim.patch_id = patch_inds;
chnkrdim.uvs_targ = uvsloc;
chnkrdim.uvs_targ(2,:) = 0;
chnkrdim.kappa = kappa;

%% Target info for v2v (geom nodes -> geom nodes)

[icoefs,pcoefs,H,alpha] = bump4(geom);

% V = eval_gauss(geom.r);
% icoefs = [zeros(2,geom.npts); icoefs(3:7,:)];
% icoefs = [zeros(6,geom.npts); ones(1,geom.npts)];

targ_geom = [];
targ_geom.r        = geom.r;
targ_geom.patch_id = geom.patch_id;
targ_geom.uvs_targ = geom.uvs_targ;
targ_geom.icoefs = icoefs;

%% Kernels
dpars   = nu;
skern   = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 's',              nu);
v2vkern   = @(s,t) v2v_kern(zk1,s,t);
repkern = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 'free_plate_eval',nu);
rhskern = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_bcs_var', nu);
fkern1 =  @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_var', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

%% b2b boundary integral system (chunkermat, same as flex2d_free_surferfun)
fprintf('Building dense b2b ...\n'); tic;

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H  + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H  ... 
    + 2*((1+nu)/2)^2*(bcoefs{2}(:)./bcoefs{1}(:)).*D*D - diag(bcoefs{2}(:)./bcoefs{1}(:)).*(1/8)*(1+nu).^2 + 1/2*diag(bcoefs{2}(:)./bcoefs{1}(:));
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

sysmat(1:2:end,1:2:end) = sysmat(1:2:end,1:2:end) - (1-nu)*(sysmat(1:2:end,2:2:end).*(bcoefs{3}./bcoefs{1}))*H;
sysmat(2:2:end,1:2:end) = sysmat(2:2:end,1:2:end) - (1-nu)*(sysmat(2:2:end,2:2:end).*(bcoefs{3}./bcoefs{1}))*H ;

D = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

sys =  D + sysmat;

densmap = zeros(3*chnkr.npt, 2*chnkr.npt);
densmap(1:3:end, 1:2:end) = eye(chnkr.npt);
densmap(3:3:end, 2:2:end) = eye(chnkr.npt);
densmap(2:3:end, 1:2:end) = -H;

sysb2b =  D + sysmat;
fprintf('  done (%.2f s)\n', toc);

%% Oversampling parameters
eps = 1e-12;
Q = getnear(S_rot); Q.targinfo = S_rot; Q.wavenumber = zk1; Q.kernel_order = 0;
novers   = get_oversampling_parameters(S_rot, Q, eps);
norderup = max(novers) - S_rot.norders(1);
norderup = min(norderup, 11 - S_rot.norders(1));
norderup = max(norderup, 1);
nover    = geom.norders(1) + norderup;

[geom_over, xinterp] = geom.oversample(nover);
[S_over,    xinterp_surfer] = oversample(S_rot, nover);

%% Near-quad function handles
getnearquad_v2v = @(varargin) flex2d.getnearquad_var(varargin{1:17}, dpars, zpars_f, varargin{18}, 'v2v');
getnearquad_s = @(varargin) flex2d.getnearquad(varargin{1:17}, dpars, zpars_f, varargin{18}, 'v2v');
getnearquad_v2b = @(varargin) flex2d.getnearquad(varargin{1:17}, dpars, zpars_f, varargin{18}, 'free_var');

%% v2v block: geom -> geom
fprintf('Building v2v near correction ...\n'); tic;
[v2v_near, norderup] = wrap_getnearquad(geom, v2vkern, eps, getnearquad_v2v, targ_geom, norderup);
fprintf('  done (%.2f s)\n', toc);

v2vkern_eval = kernel(v2vkern);
fprintf('Building v2v smooth part ...\n'); tic;
v2v_smth = (v2vkern_eval.eval(geom_over, targ_geom) .* geom_over.wts(:).') * xinterp;
fprintf('  done (%.2f s)\n', toc);

flex_v2v = v2v_near + v2v_smth;

%% v2b block: geom -> chnkr
fprintf('Building v2b near correction ...\n'); tic;
[v2b_near, ~] = wrap_getnearquad(geom, rhskern, eps, getnearquad_v2b, chnkrdim, norderup);
fprintf('  done (%.2f s)\n', toc);

rhskern_eval = kernel(rhskern);
fprintf('Building v2b smooth part ...\n'); tic;
v2b_smth = (rhskern_eval.eval(geom_over, chnkrdim) .* geom_over.wts(:).') * xinterp;
fprintf('  done (%.2f s)\n', toc);

flex_v2b = v2b_near + v2b_smth;

%% b2v block: chnkr -> geom  (chunkerkernevalmat, no PCFFT needed)
fprintf('Building b2v ...\n'); tic;
opts_cor = []; opts_cor.corrections = 0; 

repkerns = cell(7,1);
repkerns{1} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapgx',nu);
repkerns{2} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapgy',nu);
repkerns{3} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_lapg',nu);
repkerns{4} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gyy',nu);
repkerns{5} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gxx',nu);
repkerns{6} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval_gxy',nu);
repkerns{7} = @(s,t) flex2d.kern(zpars_f, s, t, 'free_plate_eval',nu);

flex_b2v = 0;
for ii = 1:7
flex_b2v = flex_b2v + icoefs(ii,:).' .* chunkerkernevalmat(chnkr, repkerns{ii}, geom.r(1:2,:), opts_cor) * densmap;
end
fprintf('  done (%.2f s)\n', toc);

%% Assemble full system matrix
%   [  a*I + V.*flex_v2v ,   V.*flex_b2v  ]   [ eta  ]   [ rhs_vol ]
%   [      flex_v2b      ,     sysb2b      ] * [ zeta ] = [ rhs_bc  ]

nvol  = geom.npts;
nbdry = 2 * chnkr.npt;

A_vv = spdiags(a * ones(nvol,1), 0, nvol, nvol) + flex_v2v;
% A_vv = spdiags(a * ones(nvol,1), 0, nvol, nvol) + spdiags(V(:), 0, nvol, nvol) * flex_v2v;
A_vb = flex_b2v;
A_bv = full(flex_v2b);
A_bb = sysb2b;

sysmat = [A_vv, full(A_vb); full(A_bv), full(A_bb)];

%% Right-hand side (manufactured solution u = sin(x)*sin(y))
rhs_vol = (4*a + 2*b - c) .* sin(geom.r(1,:)) .* sin(geom.r(2,:)) ...
    - icoefs(1,:) .* 2 .* cos(geom.r(1,:)) .* sin(geom.r(2,:)) ...     
    - icoefs(2,:) .* 2 .* sin(geom.r(1,:)) .* cos(geom.r(2,:)) ...     
    - icoefs(3,:) .* 2 .* sin(geom.r(1,:)) .* sin(geom.r(2,:)) ...     
    - icoefs(4,:) .* sin(geom.r(1,:)) .* sin(geom.r(2,:)) ... 
    - icoefs(5,:) .* sin(geom.r(1,:)) .* sin(geom.r(2,:)) ... 
    + icoefs(6,:) .* cos(geom.r(1,:)) .* cos(geom.r(2,:)) ...
    + icoefs(7,:) .* sin(geom.r(1,:)) .* sin(geom.r(2,:)) ;

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

rhs_bc = zeros(2*chnkr.npt, 1);
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

rhs_bc(1:2:end) = firstbc;
rhs_bc(2:2:end) = secondbc;

rhs = [rhs_vol(:); rhs_bc];

%% Solve
fprintf('Solving %d x %d system with backslash ...\n', size(sysmat,1), size(sysmat,2)); tic;
sol = sysmat \ rhs;
fprintf('  done (%.2f s)\n', toc);

%% eval block: geom -> geom
fprintf('Building v2v eval near correction ...\n'); tic;
[v2v_s_near, norderup] = wrap_getnearquad(geom, skern, eps, getnearquad_s, targ_geom, norderup);
fprintf('  done (%.2f s)\n', toc);

skern_eval = kernel(skern);
fprintf('Building v2v eval smooth part ...\n'); tic;
v2v_s_smth = (skern_eval.eval(geom_over, targ_geom) .* geom_over.wts(:).') * xinterp;
fprintf('  done (%.2f s)\n', toc);

flex_s_v2v = v2v_s_near + v2v_s_smth;

%% eval block: chnkr -> geom  (chunkerkernevalmat, no PCFFT needed)
fprintf('Building b2v eval ...\n'); tic;
opts_cor = []; opts_cor.corrections = 0;
flex_s_b2v = chunkerkernevalmat(chnkr, repkern, geom.r(1:2,:), opts_cor) * densmap;
fprintf('  done (%.2f s)\n', toc);

%% Error
eta = sol(1:nvol);
u = [flex_s_v2v, flex_s_b2v]*sol;
ref_u = sin(geom.r(1,:)) .* sin(geom.r(2,:));
err = abs(u - ref_u(:)) / max(abs(ref_u(:)));

fprintf('max relative error in u: %.4e\n', max(err));

figure(1); clf
plot(geom, log10(geom.patch_max(err)));
colorbar; title('log_{10} pointwise error in u');

%% Local functions

function out = v2v_kern(zk,srcinfo,targinfo)

src = srcinfo.r;
targ = targinfo.r;

[val,~,hess,third] = chnk.flex2d.hkdiffgreen(zk,src,targ);
val = 1/(2*zk^2)*val;
% grad = 1/(2*zk^2)*grad;
hess = 1/(2*zk^2)*hess;
third = 1/(2*zk^2)*third;

hessxx = hess(:,:,1);
hessxy = hess(:,:,2);
hessyy = hess(:,:,3);

lap = hessxx+hessyy;
gradlapx = third(:,:,1) + third(:,:,3);
gradlapy = third(:,:,2) + third(:,:,4);

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

dx = xt-xs;
dy = yt-ys;

dx2 = dx.*dx;
dy2 = dy.*dy;
r2 = dx2 + dy2;

eulgam = 0.5772156649015328606065120900824;

val(r2<1e-16) = 1/(2*zk^2)/(2*pi)*(log(2/zk)-log(2/(1i*zk)));
hessxx(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
hessxy(r2<1e-16) = 0;
hessyy(r2<1e-16) = 1/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
lap(r2<1e-16) = 2/(2*zk^2)*(zk^2/(4*pi)*(log(zk)-1i*pi/2+eulgam-0.5-log(2))-(1i*zk)^2/(4*pi)*(log(1i*zk)-1i*pi/2+eulgam-0.5-log(2)));
gradlapx(r2<1e-16) = 0;
gradlapy(r2<1e-16) = 0;

out = targinfo.icoefs(1,:).'.*gradlapx;
out = out + targinfo.icoefs(2,:).'.*gradlapy;
out = out + targinfo.icoefs(3,:).'.*lap;
out = out + targinfo.icoefs(4,:).'.*hessyy;
out = out + targinfo.icoefs(5,:).'.*hessxx;
out = out + targinfo.icoefs(6,:).'.*hessxy;
out = out + targinfo.icoefs(7,:).'.*val;

end
