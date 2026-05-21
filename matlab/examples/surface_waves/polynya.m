% polynya.m
%
% Solve the polynya problem (flexural gravity waves, starfish geometry).
% All volume and mixed operators are assembled as dense matrices using
% oversampled smooth quadrature + near-field correction, following the
% same pattern as ice_floe.m.
%
% Unknowns: [mu (vol); rho (2*npt bdry)]
%
% System:
%   v2v(mu) + b2v(rho)   = rhs_vol
%   v2b(mu) + b2b(rho)   = rhs_bdry

cparams = []; cparams.maxchunklen = 0.3;
chnkr = chunkerfunc(@(t) starfish(t,5,0.),cparams);
[S, chnkr] = triangulate_chunker_interior(chnkr, 6, 0.3);

chnkr = sort(chnkr);

kappa = signed_curvature(chnkr);
kappa = kappa(:);

novers = 10;
ivpp = 0;

nu = 0.33;

maxrad = max(S.r(:));

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
wireframe(S,struct('wfill',0,'nsign',-1))
plot(chnkr,'x-')
view(0,90)
drawnow

omega = 0.124;

alpha = 23;
gamma = 0.5;
nondimfac = 0.00016; % 1/m

[rts,ejs] = surfwave.flex.find_roots(alpha,gamma);
zpars = [rts;ejs;alpha;gamma;nu];

zk_ext = rts(1)

chnkr = chnkr.move([0;0],[0;0],0,nondimfac);
chnkr = sort(chnkr);

chnkr2 = chnkr;
chnkr2.rstor = [chnkr2.r; 0*chnkr2.r(1,:,:)];
chnkr2.nstor = [chnkr2.n; 0*chnkr2.n(1,:,:)];
chnkr2.dstor = [chnkr2.d; 0*chnkr2.d(1,:,:)];
chnkr2.d2stor = [chnkr2.d2; 0*chnkr2.d2(1,:,:)];

targinfo = []; targinfo.r = S.r;

irhs1 = sin(1.1*S.r(1,:)) - 1i*cos(1.3*(S.r(2,:))) + 1/pi*S.r(1,:);
figure(2); clf
plot(S,real(irhs1))
title('rhs1')
view(0,90)
drawnow

%% Setup

gs_kern            = @(s,t) surfwave.flex.kern(s,t,'gs_s',nu,rts,ejs);
gp_kern            = @(s,t) surfwave.flex.kern(s,t,'gphi_s',nu,rts,ejs);
s3d_gp_kern        = @(s,t) surfwave.flex.kern(s,t,'s3d_gphi',nu,rts,ejs);
s3d_kern           = @(s,t) surfwave.flex.lap3dkern(s.r,t.r);
s3d_plus_s3d_gp_kern = @(s,t) s3d_kern(s,t) + s3d_gp_kern(s,t);
bilap_gp_kern      = @(s,t) 2/alpha*s3d_kern(s,t) ...
                           + 2/alpha*s3d_gp_kern(s,t) - gamma/alpha*gp_kern(s,t);

ikerngs   = @(s,t) surfwave.flex.kern(s, t, 'free_plate_gs_eval',   nu, rts, ejs);
ikerngphi = @(s,t) surfwave.flex.kern(s, t, 'free_plate_gphi_eval', nu, rts, ejs);

gs_v2b_kern   = @(s,t) surfwave.flex.kern(s, t, 'gs_v2b',   nu, rts, ejs);
gphi_v2b_kern = @(s,t) surfwave.flex.kern(s, t, 'gphi_v2b', nu, rts, ejs);
ext_rep_kern  = @(s,t) -gamma/4*surfwave.flex.kern(s,t,'gs_s',  nu,rts,ejs) ...
                       +1/2   *surfwave.flex.kern(s,t,'gphi_s',nu,rts,ejs);
v2b_bc_kern   = @(s,t) -gamma/4*gs_v2b_kern(s,t) + 1/2*gphi_v2b_kern(s,t);

eps = 1e-9;

[S_over, xinterp] = S.oversample(novers);

targ_S = []; targ_S.r = S.r;

getnearquad_s_v2v   = @(varargin) surfwave.flex.getnearquad_flex( ...
    varargin{1:17}, zpars, varargin{18}, 1, ivpp);
getnearquad_p_v2v   = @(varargin) surfwave.flex.getnearquad_flex( ...
    varargin{1:17}, zpars, varargin{18}, 2, ivpp);
getnearquad_s3dp_v2v = @(varargin) surfwave.flex.getnearquad_flex( ...
    varargin{1:17}, zpars, varargin{18}, 7, ivpp);

nearquad_opts = []; nearquad_opts.nover = novers;

%% Quadrature

t2 = tic;

fprintf('Gs vol-to-vol quadrature ...           '); tic;
gs_v2v_near = getnearquad_kern(S, gs_kern, eps, getnearquad_s_v2v, targ_S, nearquad_opts);
gs_v2v_smth = (gs_kern(S_over, targ_S) .* S_over.wts(:).') * xinterp;
gs_v2v      = gs_v2v_near + gs_v2v_smth;
tgs = toc; fprintf('%.2f s\n', tgs);

fprintf('Gphi vol-to-vol quadrature ...         '); tic;
gp_v2v_near = getnearquad_kern(S, gp_kern, eps, getnearquad_p_v2v, targ_S, nearquad_opts);
gp_v2v_smth = (gp_kern(S_over, targ_S) .* S_over.wts(:).') * xinterp;
gp_v2v      = gp_v2v_near + gp_v2v_smth;
tgp = toc; fprintf('%.2f s\n', tgp);

fprintf('S3d+S3d*Gphi vol-to-vol quadrature ... '); tic;
s3dp_v2v_near = getnearquad_kern(S, s3d_plus_s3d_gp_kern, eps, getnearquad_s3dp_v2v, targ_S, nearquad_opts);
s3dp_v2v_smth = (s3d_plus_s3d_gp_kern(S_over, targ_S) .* S_over.wts(:).') * xinterp;
s3dp_v2v      = s3dp_v2v_near + s3dp_v2v_smth;
ts3d = toc; fprintf('%.2f s\n', ts3d);

bilap_gp_mat = 2/alpha*s3dp_v2v - gamma/alpha*gp_v2v;

hilbert = @(s,t) surfwave.flex.kern(s,t,'hilb',nu,rts,ejs);
opts_pv = []; opts_pv.sing = 'pv';
H = chunkermat(chnkr, hilbert, opts_pv);

npt = chnkr.npt;
densmap = zeros(3*npt, 2*npt);
densmap(1:3:end,1:2:end) = eye(npt);
densmap(3:3:end,2:2:end) = eye(npt);
densmap(2:3:end,1:2:end) = H;

opts_corr = []; opts_corr.eps = eps; 

fprintf('Bdry-to-vol ... '); tic;
gs_b2v = chunkerkernevalmat(chnkr, ikerngs,   S.r(1:2,:), opts_corr);
gp_b2v = chunkerkernevalmat(chnkr, ikerngphi, S.r(1:2,:), opts_corr);
tb2v = toc; fprintf('%.2f s\n', tb2v);

gs_b2v = gs_b2v * densmap;
gp_b2v = gp_b2v * densmap;

b2v_mat = 1/2*gs_b2v - gp_b2v;

targinfo2 = [];
targinfo2.r     = [chnkr.r(:,:); zeros(1,chnkr.npt)];
targinfo2.n     = [chnkr.n(:,:); zeros(1,chnkr.npt)];
targinfo2.d     = [chnkr.d(:,:); zeros(1,chnkr.npt)];
targinfo2.d2    = [chnkr.d2(:,:); zeros(1,chnkr.npt)];
targinfo2.kappa = kappa(:);
[~, patch_inds, uvsloc, dists] = get_closest_pts(S, targinfo2);
targinfo2.patch_id  = patch_inds;
targinfo2.uvs_targ  = uvsloc;
targinfo2.uvs_targ(2,:) = 0;

chnkr_info2 = []; chnkr_info2.r = chnkr.r(1:2,:);
chnkr_info2.n = chnkr.n(1:2,:);

fprintf('Vol-to-bdry quadrature ...             '); tic;
getnearquad_v2b = @(varargin) surfwave.flex.getnearquad_v2b_flex( ...
    varargin{1:17}, zpars, varargin{18});
v2b_near = getnearquad_kern(S, kernel(v2b_bc_kern), eps, getnearquad_v2b, targinfo2, nearquad_opts);
v2b_smth = (v2b_bc_kern(S_over, targinfo2) .* S_over.wts(:).') * xinterp;
v2b_mat  = v2b_near + v2b_smth;
tv2b = toc; fprintf('%.2f s\n', tv2b);

opts_log  = []; opts_log.sing  = 'log';
opts_smth = []; opts_smth.quad = 'native'; opts_smth.sing = 'smooth';

fkern1   = @(s,t) surfwave.flex.kern(s, t, 'free plate first part',    nu, rts, ejs);
fkern1bh = @(s,t) surfwave.flex.kern(s, t, 'free plate first part bh', nu, rts, ejs);
fkern2bh = @(s,t) surfwave.flex.kern(s, t, 'free plate hilbert bh',    nu, rts, ejs);
fkern2   = @(s,t) surfwave.flex.kern(s, t, 'free plate hilbert',       nu, rts, ejs);
double   = @(s,t) chnk.lap2d.kern(s,t,'d',nu,rts,ejs);

fprintf('Bdry-to-bdry quadrature ...            '); tic;
sysmat1   = chunkermat(chnkr, fkern1,   opts_log);
sysmat1bh = chunkermat(chnkr, fkern1bh, opts_smth);
sysmat2bh = chunkermat(chnkr, fkern2bh, opts_smth);
sysmat2   = chunkermat(chnkr, fkern2,   opts_log);
D         = chunkermat(chnkr, double,   opts_log);

sysmat1bh(isnan(sysmat1bh)) = 0;
sysmat1bh(2:2:end,1:2:end)  = sysmat1bh(2:2:end,1:2:end) + (2/alpha)*diag((-3+3*nu)/(8*pi)*kappa.^2.*chnkr.wts(:));
sysmat1 = sysmat1 + sysmat1bh;

sysmat2bh(isnan(sysmat2bh)) = 0;
sysmat2 = sysmat2 + sysmat2bh;

sysmat2(1:2:end,1:2:end) = sysmat2(1:2:end,1:2:end)*H - 2*((1+nu)/2)^2*D*D*(2/alpha);
sysmat2(2:2:end,1:2:end) = sysmat2(2:2:end,1:2:end)*H;

sysmat_b2b = sysmat1 + sysmat2;

jumpmat = kron(eye(chnkr.npt), (1/alpha)*[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2]);
idsysmat = jumpmat + sysmat_b2b/2;

tb2b = toc; fprintf('  done (%.2f s)\n', tb2b);

tpre = toc(t2);
fprintf('\nPreprocessing summary:\n');
fprintf('  tgs        = %.2f s\n', tgs);
fprintf('  tgp        = %.2f s\n', tgp);
fprintf('  ts3d       = %.2f s\n', ts3d);
fprintf('  tb2v_corr  = %.2f s\n', tb2v_corr);
fprintf('  tb2v       = %.2f s\n', tb2v);
fprintf('  tv2b       = %.2f s\n', tv2b);
fprintf('  tb2b       = %.2f s\n', tb2b);
fprintf('  tpre total = %.2f s\n', tpre);

%% Forming RHS

rhsbdry = zeros(2*chnkr.npt,1);
rhs1 = [irhs1.'; rhsbdry];

%% Assembling matrix operators and solving

v2v_mat = 1/2*eye(S.npts) - gamma/4*gs_v2v + 1/2*gp_v2v - alpha/2*bilap_gp_mat;

% Full dense system matrix
A = [v2v_mat,    b2v_mat; ...
     v2b_mat,    idsysmat];

fprintf('System size: %d x %d\n', size(A,1), size(A,2));

fprintf('Solving with gmres ...\n'); tic;
[sol1,flag1,relres1,iter1,resvec1] = gmres(A, rhs1, [], 1e3*eps, 4000, [], []);
tsolve = toc; fprintf('  done (%.2f s)\n', tsolve);

if flag1
    fprintf('GMRES failed with relative residual %f.\n', relres1);
else
    fprintf('GMRES succeeded in %d iterations.\n', iter1(2));
end

mu  = sol1(1:S.npts);
rho = sol1(S.npts+1:end);

gs_b2v_rho = gs_b2v * rho;
gp_b2v_rho = gp_b2v * rho;

phizint1 = 1/2*mu - gamma/4*gs_v2v*mu + 1/2*gp_v2v*mu + 1/2*gs_b2v_rho;
phiint1  = alpha/2*bilap_gp_mat*mu + gp_b2v_rho;

figure(3); clf
plot(S,real(phiint1)); colorbar

figure(2); clf
plot(S,log10(abs(surf_fun_error(S,phiint1)))); colorbar
