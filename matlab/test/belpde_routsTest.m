%
% This file tests the Laplace-Beltrami and Helmholtz-Beltrami
% parametrix/remainder quadrature on an ellipsoid, against an
% analytic solution built from belpde_source
%
%
run ../startup.m

S = geometries.ellipsoid([1,1,1.5], 2*[1,1,1], [0;0;0], 8);
wts = S.wts;

x_source = [0.15; 0.25; 0.35];
charge = 1;
eps = 1e-9;

%% Laplace-Beltrami

[uex, rhs] = belpde.belpde_source(S, x_source, charge, 0);

Amat = surfermat(S, kernel3d('bel','rlb'), eps);
Amat = Amat + eye(size(Amat));
Kmat = surfermat(S, kernel3d('bel','klb'), eps);

sigma = gmres(Amat, rhs, [], 1e-12, 200);
u = Kmat*sigma;
u = u - sum(u.*wts)/sum(wts);

err1 = norm((u-uex).*sqrt(wts))/norm(uex.*sqrt(wts));
fprintf('Error in Laplace-Beltrami solution=%d\n',err1);
assert(err1 < 1e-3, 'Laplace-Beltrami error too large');

%% Helmholtz-Beltrami

zk = 2.1;
[uex, rhs] = belpde.belpde_source(S, x_source, charge, zk);

Amat = surfermat(S, kernel3d('bel','rhb',zk), eps);
Amat = Amat + eye(size(Amat));
Kmat = surfermat(S, kernel3d('bel','khb',zk), eps);

sigma = gmres(Amat, rhs, [], 1e-12, 200);
u = Kmat*sigma;

err2 = norm((u-uex).*sqrt(wts))/norm(uex.*sqrt(wts));
fprintf('Error in Helmholtz-Beltrami solution=%d\n',err2);
assert(err2 < 1e-3, 'Helmholtz-Beltrami error too large');

%% Variable wavenumber, using the Laplace-Beltrami parametrix

zk2fun = @(x) 4 + x(3,:).^2;

[uex, rhs] = belpde.belpde_source(S, x_source, charge, @(x) sqrt(zk2fun(x)));

kparam = kernel3d('bel','klb');
zk2f = @(t) reshape(zk2fun(t.r), 1, 1, []);

Amat = surfermat(S, kernel3d('bel','rlb') + zk2f.*kparam, eps);
Amat = Amat + eye(size(Amat));
Kmat = surfermat(S, kparam, eps);

sigma = gmres(Amat, rhs, [], 1e-12, 200);
u = Kmat*sigma;

err3 = norm((u-uex).*sqrt(wts))/norm(uex.*sqrt(wts));
fprintf('Error in variable wavenumber solution=%d\n',err3);
assert(err3 < 1e-3, 'Variable wavenumber Beltrami error too large');
