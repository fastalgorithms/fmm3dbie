%
% This file demonstrates the Maxwell pec solver with the
% NRCCIE representations
%
%

% Generate Geometry
radii = [1,0.3,.3];
S = geometries.wobbletorus(radii,5,[1,1,1],[8,10],7);
S = geometries.ellipsoid([1,1,1.5],[1,1,1], [0;0;0], 8);
[srcvals,~,~,~,~,wts] = extract_arrays(S);
[~, npts] = size(srcvals);

% Set wavenumber and representation parameter
zk = 7;
alpha = 1.0;

% generate RHS
eps = 1e-7;
dir = [pi/3,0]; pol = [1;0.2];
[Einc, Hinc] = em3d.planewave(zk, dir, pol, S);

% Solve
zpars = complex([zk, alpha]);
opts = [];
opts.eps_gmres = 1e-10;
tic;
[densities] = em3d.pec.solver(S, Einc, Hinc, eps, zk, alpha, opts);

% evaluate scattered field
targ_info = [];
npts = 200;
xx = linspace(-3,3,npts);
[X,Y,Z] = meshgrid(xx,xx,0);
targ_info.r = [X(:).';Y(:).';Z(:).'];
[E, H] = em3d.pec.eval(S, densities, targ_info, eps, zk, alpha);

[E_ex, H_ex] = em3d.planewave(zk, dir, pol, targ_info);

u_in = [E; H];
u_ex = [E_ex; H_ex];
utot = u_in + u_ex;

% plot densities
figure(1);clf; t = tiledlayout('flow'); t.TileSpacing = 'tight';
labs = {"$J_x$","$J_y$","$J_z$","$\rho$"};
for i = 1:size(densities,1)
    nexttile()
    plot(S,real(densities(i,:)))
    title(labs{i},'Interpreter','latex')
    set(gca,'fontsize',14)
end


% plot electric field
figure(3);clf
idx = 1;
h = pcolor(X,Y,reshape(real(utot(idx,:)),npts,npts));hold on
plot(S,real(S.n(idx,:).*densities(4,:))) % E on the surface is n*rho
colorbar
title('$\Re E_x$','Interpreter','latex')
set(gca,'fontsize',14)

figure(4);clf
h = pcolor(X,Y,reshape(vecnorm(utot(1:3,:),2,1),npts,npts));hold on
plot(S,abs(densities(4,:))) % E on the surface is n*rho
colorbar
title('$|E|$','Interpreter','latex')

set(gca,'fontsize',14)