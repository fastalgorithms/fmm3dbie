% ice_floe.m
%
% Solve the ice floe problem (gravity outside and flexural gravity inside)
%   gs_v2v  : S -> S   (gravity single-layer)
%   flex_v2v: S -> S   (flex single-layer, 'v2v')
%   flex_v2b: S -> chnkr  (flex, 'free_plate_bcs')
%   flex_b2v: chnkr -> S  (flex, 'free_plate_eval')
%   flex_b2b: chnkr -> chnkr (BIE system)
%
% The coupled system is:
%   h = 2*mu + Gs*mu - flex_v2v*eta - flex_b2v*zeta   = rhs_grav
%   f = -g/2*Gs*mu  + a*eta + b*(flex_v2v*eta + flex_b2v*zeta) = rhs_vol
%   g =  flex_v2b*eta + sysb2b*zeta                    = rhs_bc
%
% Solved as a dense linear system with gmres.

%% Parameters

% PDE coefficients
alpha = 1;
gamma = 1;   % gravity wavenumber
beta  = 0.1;

% flex coefficients
a  = 1;
b  = 0.;
c  = -1/pi;
nu = 0.3;

% Manufactured solution (plane wave)
theta    = pi;  % propagation angle

% Tolerances
eps_quad = 1e-8;
eps_gmres = 1e-8;

gamma_pw = gamma;

zk1 = sqrt((-b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((-b - sqrt(b^2 + 4*a*c)) / (2*a));
zpars_f = [zk1, zk2];
gravpars = [alpha, -beta+c/a, gamma];
zpars = gamma;

fprintf('alpha = %f, beta = %f, gamma = %f\n', alpha, beta, gamma)

%% Geometry
cparams = []; cparams.maxchunklen = 0.2;
chnkr = chunkerfunc(@(t) starfish(t,5,0.1),cparams);
[S, chnkr] = triangulate_chunker_interior(chnkr, 6, 0.2);

chnkr = sort(chnkr);

kappa = signed_curvature(chnkr);
kappa = kappa(:);

fprintf('npatches = %d, nchunks = %d\n', S.npatches, chnkr.nch);

%% Build chnkrdim

chnkr_info = [];
loc_field = fieldnames(chnkr)';
for i = 1:length(loc_field)
    if mod(numel(chnkr.(loc_field{i})), chnkr.npt) == 0
        chnkr_info.(loc_field{i}) = chnkr.(loc_field{i})(:,:);
    end
end
chnkr_info.r   = chnkr_info.r;
% chnkr_info.d   = chnkr_info.d./vecnorm(chnkr_info.d);
chnkr_info.d2  = chnkr_info.d2;
chnkr_info.n   = chnkr_info.n;
chnkr_info.wts = chnkr_info.wts(:);

[~, patch_inds, uvsloc] = get_closest_pts(S, chnkr_info);
chnkr_info.patch_id = patch_inds;
chnkr_info.uvs_targ = uvsloc;
% triangulate_chunker_interior aligns the edge patches
chnkr_info.uvs_targ(2,:) = 0;
chnkr_info.kappa = kappa;

targ_S = [];
targ_S.r        = S.r;
targ_S.patch_id = S.patch_id;
targ_S.uvs_targ = S.uvs_targ;

%% Kernels

dpars    = nu;
skern    = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 's',              nu);
repkern  = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 'free_plate_eval',nu);
rhskern  = @(s,t) chnk.flex2d.kern(zpars_f, s, t, 'free_plate_bcs', nu);
biekern  = @(s,t) free_plate_kern(zpars_f, s, t, nu);

ivpp     = 1;
iker     = 0;
S3d_scal = zpars*4;
gs_kern  = @(s,t) surfwave.gravity.kern(zpars, 1, s, t, 'gs_s');

%% b2b boundary integral system

double_kern  = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert_kern = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts_log = []; opts_log.sing = 'log';
opts_pv  = []; opts_pv.sing  = 'pv';

D = chunkermat(chnkr, double_kern,  opts_log);
H = chunkermat(chnkr, hilbert_kern, opts_pv);
flexsys_D = 2*((1+nu)/2)^2 * D * D;

diagval = -[-1/2 + (1/8)*(1+nu)^2, 0; 0, 1/2];
diagval = kron(eye(chnkr.npt), diagval);
diagval(1:2:end,1:2:end) = diagval(1:2:end,1:2:end) + flexsys_D;

flex_b2b_val = chunkermat(chnkr, biekern, opts_log);
npt = chnkr.npt;
densmap = zeros(3*npt, 2*npt);
densmap(1:3:end, 1:2:end) = eye(npt);
densmap(3:3:end, 2:2:end) = eye(npt);
densmap(2:3:end, 1:2:end) = -H;

sysb2b = diagval + flex_b2b_val * densmap;

%% Oversampling parameters

Q = getnear(S); Q.targinfo = S; Q.wavenumber = zk1; Q.kernel_order = 0;
novers = get_oversampling_parameters(S, Q, 1e2*eps_quad);
novers = max(novers,S.norders+1);
[S_over, xinterp] = S.oversample(novers);
nearquad_opts = []; nearquad_opts.nover = novers;
%% Near-quad function handles

getnearquad_grav = @(varargin) surfwave.gravity.getnearquad_gravity(varargin{1:17}, zpars, varargin{18}, iker, ivpp, S3d_scal);
getnearquad_v2v  = @(varargin) flex2d.getnearquad(varargin{1:17}, dpars, zpars_f, varargin{18}, 'v2v');
getnearquad_v2b  = @(varargin) flex2d.getnearquad(varargin{1:17}, dpars, zpars_f, varargin{18}, 'free');

%% Gravity block

fprintf('Building gravity blocks ...\n'); tic;
gs_v2v_near = getnearquad_kern(S, gs_kern, eps_quad, getnearquad_grav, targ_S, nearquad_opts);
gs_v2v_smth = (gs_kern(S_over, targ_S) .* S_over.wts(:).') * xinterp;

gs_v2v = gs_v2v_near + gs_v2v_smth;

fprintf('  done (%.2f s)\n', toc);

%% Flex v2v block

fprintf('Building flex ...\n'); 
flex_v2v_near = getnearquad_kern(S, skern, eps_quad, getnearquad_v2v, targ_S, nearquad_opts);
flex_v2v_smth = (skern(S_over, targ_S) .* S_over.wts(:).') * xinterp;

flex_v2v = flex_v2v_near + flex_v2v_smth;

flex_v2b_near = getnearquad_kern(S, rhskern, eps_quad, getnearquad_v2b, chnkr_info, nearquad_opts);
flex_v2b_smth = (rhskern(S_over, chnkr_info) .* S_over.wts(:).') * xinterp;

flex_v2b = flex_v2b_near + flex_v2b_smth;

flex_b2v = chunkerkernevalmat(chnkr, repkern, S.r(1:2,:)) * densmap;
fprintf('  done (%.2f s)\n', toc);

%% Assemble dense system matrix
%
% Unknowns: [mu (nvol); eta (nvol); zeta (2*npt)]
%   mu   = gravity surface density
%   eta  = flex volume density
%   zeta = flex boundary density pair
%
% Equations (matching apply_grav):
%   2*mu + Gs*mu - (flex_v2v*eta + flex_b2v*zeta)      = rhs_grav
%   -g/2*Gs*mu + a*eta + b*(flex_v2v*eta + flex_b2v*zeta) = rhs_vol
%   flex_v2b*eta + sysb2b*zeta                          = rhs_bc

nvol  = S.npts;
nbdry = 2*chnkr.npt;

g_val = gravpars(3);   % gamma
b_val = gravpars(2);   % -beta + c/a
a_val = gravpars(1);   % alpha

% Row block 1: gravity equation (mu unknowns)
A_hmu  = 2*eye(nvol) + gs_v2v;
A_heta = -flex_v2v;
A_hzeta = -flex_b2v;

% Row block 2: volume flex equation (eta unknowns)
A_fmu  = -g_val/2 * gs_v2v;
A_feta = a_val*eye(nvol) + b_val*flex_v2v;
A_fzeta = b_val*flex_b2v;

% Row block 3: boundary equation (zeta unknowns)
A_gmu   = zeros(nbdry, nvol);
A_geta  = flex_v2b;
A_gzeta = sysb2b;

sysmat = [A_hmu,  A_heta,  A_hzeta; ...
          A_fmu,  A_feta,  A_fzeta; ...
          A_gmu,  A_geta,  A_gzeta];

fprintf('System size: %d x %d\n', size(sysmat,1), size(sysmat,2));

%% RHS and solve (accuracy test: gamma=0 so manufactured solution decouples)

[rhs, ref_u] = ice_floe_rhs(S, chnkr, kappa, nu, alpha, beta, gamma_pw, theta);

gravpars_acc = gravpars;
gravpars_acc(3) = 0;
g_acc = gravpars_acc(3);

A_fmu_acc  = -g_acc/2 * gs_v2v;
sysmat_acc = [A_hmu,        A_heta,  A_hzeta; ...
              A_fmu_acc,    A_feta,  A_fzeta; ...
              A_gmu,        A_geta,  A_gzeta];

fprintf('Solving with gmres (accuracy test, gamma=0) ...\n'); tic;
sol_acc = gmres(sysmat_acc, rhs, [], eps_gmres, size(sysmat_acc,1));
fprintf('  done (%.2f s)\n', toc);

etazeta_acc = sol_acc(nvol+1:end);
eta_acc     = sol_acc(nvol+1:2*nvol);
zeta_acc    = sol_acc(2*nvol+1:end);
u = flex_v2v*eta_acc + flex_b2v*zeta_acc;

err = abs(u + ref_u) / norm(ref_u, inf);
fprintf('Max relative error in u: %.4e\n', max(err));

figure(10); clf
subplot(1,2,1)
plot(S, real(ref_u)); colorbar; title('u_{true}'); 
subplot(1,2,2)
plot(S, S.patch_max(log10(err))); colorbar; title('log_{10} error');

figure(11); clf
plot(S, log10(surf_fun_error(S, u))); colorbar
title('log_{10} resolution'); 

%% Full solve (gamma != 0)

fprintf('Solving with gmres (full system) ...\n'); tic;
sol = gmres(sysmat, rhs, [], eps_gmres, size(sysmat,1));
fprintf('  done (%.2f s)\n', toc);

eta  = sol(nvol+1:2*nvol);
zeta = sol(2*nvol+1:end);
u_s  = flex_v2v*eta + flex_b2v*zeta;
u_in = ref_u;

figure(20); clf
plot(S, real(u_s + u_in)); colorbar;
title('Re(u_{scattered} + u_{in})')

figure(21); clf
plot(S, log10(surf_fun_error(S, u_s))); colorbar; 
title('log_{10} resolution')

%% Local functions

function submat = free_plate_kern(zk, s, t, nu)
submat_0 = chnk.flex2d.kern(zk, s, t, 'free_plate', nu);
submat = zeros(2, size(t.r(:,:),2), 3, size(s.r(:,:),2));
submat_0 = reshape(submat_0, 4, size(t.r(:,:),2), 2, size(s.r(:,:),2));
submat(:,:,1,:) = submat_0(1:2,:,1,:);
submat(:,:,2,:) = submat_0(3:4,:,1,:);
submat(:,:,3,:) = submat_0(1:2,:,2,:);
submat = reshape(submat, 2*size(t.r(:,:),2), 3*size(s.r(:,:),2));
end
