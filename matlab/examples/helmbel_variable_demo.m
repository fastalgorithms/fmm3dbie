% Solve the Helmholtz-Beltrami equation with a variable wavenumber on an
% ellipsoid, against an analytic solution, using dense GMRES

%% Setup geometry and right hand side
S = geometries.ellipsoid([1,1,1.5], 4*[1,1,1], [0;0;0], 8);

zk2fun = @(x) 0 + 4*(1.-x(1,:)).^2;     % wavenumber squared as a function of x

% Refine where the local wavenumber is large
dirs = eye(3);
fpw = @(S) exp(1i*sqrt(zk2fun(S.r)).*(dirs.'*S.r));
S = S.resolve_fun(1e-5, fpw,S.norders(1)-1,3);

x_source = [0.15; 0.25; 0.35];   % interior source point
charge = 1;

[uex, rhs] = belpde.belpde_source(S, x_source, charge, @(x) sqrt(zk2fun(x)));

kparam  = kernel3d('bel','klb');
krem    = kernel3d('bel','rlb');

zk2f = @(t) reshape(zk2fun(t.r), 1, 1, []);
kvar = krem + zk2f.*kparam;

eps = 1e-9;

%% Build the system matrix and solve
tic;
Rmat = surfermat(S, kvar, eps);
Amat = Rmat + eye(size(Rmat));
Kmat = surfermat(S, kparam, eps);
tbuild = toc

tic;
sigma = gmres(Amat, rhs, [], 1e-12, 200);
tsolve = toc

u = Kmat*sigma;

err = norm((u-uex).*sqrt(S.wts))/norm(uex.*sqrt(S.wts))

%% Plot the solution and the error
figure(1); clf
subplot(1,2,1)
plot(S, real(u));
colorbar
title('variable Helmholtz-Beltrami solution')

subplot(1,2,2)
plot(S, log10(patch_max(S, abs(u-uex))));
colorbar
title('log_{10} patchwise error')

%% Solve again with a Gaussian right hand side

x0 = [0; 0; 1.5];    % centre of the Gaussian, at the top of the ellipsoid
sig0 = 0.3;
f = exp(-vecnorm(S.r - x0).^2/(2*sig0^2)).';

ug = Kmat*gmres(Amat, f, [], 1e-12, 200);

figure(2); clf
subplot(1,2,1)
plot(S, f);
colorbar
title('Gaussian right hand side')

subplot(1,2,2)
plot(S, real(ug));
colorbar
title('variable Helmholtz-Beltrami solution')
