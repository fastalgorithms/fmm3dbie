%
%  capillary_ext_dirichlet_disk
%
%  Solves the exterior Dirichlet problem for the capillary-gravity equation
%  in the exterior of a disk. The incident wave is a plane wave with
%  wavenumber zk_ext.
%
%  The unknown densities are:
%    mu  - volume density on S (the disk)
%    rho - boundary density on chnkr (the disk boundary)
%
%  After solving, the script plots phi in the surfer (interior) and on a
%  Cartesian grid outside the surfer boundary.
%
%  Requires chunkie

%% Parameters

nplot = 120;
nplot = 60;

alpha = 0.1;
beta  = -0.1;

[rts, ejs] = surfwave.capillary.find_roots(alpha, beta);
zpars  = [rts; ejs];
zk_ext = rts(1);

%% Geometry

norder = 6;
rad    = 5;

S = rad*geometries.disk([], [], [3 2 4], norder);

nch          = 4*4;
cparams.ta   = pi/nch;
cparams.tb   = 2*pi + cparams.ta;
chnkr        = rad*chunkerfuncuni(@(t) ellipse(t), nch, cparams);
chnkr        = sort(chnkr);

figure(1); clf
plot(S, rand(S.npatches, 1))
hold on
wireframe(S,struct('nsign',-1))
plot(chnkr, 'k-', 'LineWidth', 2)
title('Geometry: disk and boundary')

%% Right-hand side (incident plane wave)

rhs_vol = -(zk_ext - 1) * exp(1i*zk_ext*S.r(1,:));
rhs_bdy = -zk_ext .* exp(1i*zk_ext*chnkr.r(1,:));

%% Matrix assembly

eps = 1e-8;

% Volume-to-volume
tic
Gsv2v      = surfwave.capillary.matgen(S, 'gs',      zpars, eps).';
Gphiv2v    = surfwave.capillary.matgen(S, 'gphi',    zpars, eps).';
lapGphiv2v = surfwave.capillary.matgen(S, 'lapgphi', zpars, eps).';
fprintf('%5.2e s : v2v matrices\n', toc)

% Volume-to-boundary
tic
v2bmats  = surfwave.capillary.gsevalmat(S, chnkr, zpars, eps);
Gsv2b    = reshape(v2bmats(1,:), [S.npts, chnkr.npt]).';
Gphiv2b  = reshape(v2bmats(2,:), [S.npts, chnkr.npt]).';
fprintf('%5.2e s : v2b matrices\n', toc)

% Boundary-to-boundary
fkern_gs_d  = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, 'gs_d');
tic
Dgsb2b = chunkermat(chnkr, fkern_gs_d);
fprintf('%5.2e s : b2b matrix\n', toc)

% Boundary-to-volume
fkern_gphi_d = @(s,t) surfwave.capillary.kern(rts, ejs, s, t, 'gphi_d');
tic
Dgsb2v   = chunkerkernevalmat(chnkr, fkern_gs_d,   S.r(1:2,:));
Dgphib2v = chunkerkernevalmat(chnkr, fkern_gphi_d, S.r(1:2,:));
fprintf('%5.2e s : b2v matrices\n', toc)

%% Exterior evaluation targets

L     = 2 * max(vecnorm(chnkr.r(:,:)));
x1    = linspace(-L, L, nplot);
[xx, yy] = meshgrid(x1, x1);
in    = chunkerinterior(chnkr, {x1, x1});
out   = ~in;

targs_ext      = [];
targs_ext.r    = [xx(out).'; yy(out).'; 0*xx(out).'];
ntarg          = sum(out(:));

tic
evalmats     = surfwave.capillary.gsevalmat(S, targs_ext, zpars, eps);
lapGphiemats = surfwave.capillary.lapgphievalmat(S, targs_ext, zpars, eps);
Gsv2ext      = reshape(evalmats(1,:),    [S.npts, ntarg]).';
Gphiv2ext    = reshape(evalmats(2,:),    [S.npts, ntarg]).';
lapGphiv2ext = reshape(lapGphiemats(1,:),[S.npts, ntarg]).';
fprintf('%5.2e s : exterior eval matrices\n', toc)

%% Linear system

T = [1/2*eye(S.npts),              zeros(S.npts, chnkr.npt); ...
     zeros(chnkr.npt, S.npts),     1/(2*alpha)*eye(chnkr.npt)];

K = [1/4*(beta+abs(beta))*Gsv2v + (1-abs(beta))/2*Gphiv2v + alpha/2*lapGphiv2v, ...
     1/2*Dgsb2v - Dgphib2v; ...
     1/4*(beta+abs(beta))*Gsv2b + 1/2*Gphiv2b, ...
     1/2*Dgsb2b];

lhs = T + K;
rhs = [rhs_vol, rhs_bdy].';

tic
sol = gmres(lhs, rhs, [], 1e-10, 400);
fprintf('%5.2e s : GMRES solve\n', toc)

mu  = sol(1:S.npts);
rho = sol(S.npts+1:end);

%% Evaluate phi in the interior (surfer)

phiint = -alpha/2*lapGphiv2v*mu + abs(beta)/2*Gphiv2v*mu ...
         + Dgphib2v*rho + exp(1i*zk_ext*S.r(1,:)).';

%% Evaluate phi in the exterior (grid)

Dgsb2ext   = chunkerkerneval(chnkr, fkern_gs_d,   rho, targs_ext.r(1:2,:));
Dgphib2ext = chunkerkerneval(chnkr, fkern_gphi_d, rho, targs_ext.r(1:2,:));

phiext_vals = -alpha/2*lapGphiv2ext*mu + abs(beta)/2*Gphiv2ext*mu ...
              + Dgphib2ext + exp(1i*zk_ext*targs_ext.r(1,:)).';

phiext = NaN(size(xx));
phiext(out) = phiext_vals;

%% Plot

figure(2); clf
plot(S, real(phiint));
hold on
surf(xx, yy, xx*0, real(phiext), 'EdgeColor', 'none')
view(0, 90)
axis equal
plot(chnkr, 'k-', 'LineWidth', 2)
colorbar
title('\phi  (real part)')
