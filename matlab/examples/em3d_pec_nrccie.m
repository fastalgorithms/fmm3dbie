% EM3D_PEC_NRCCIE  PEC scattering via NRCCIE using surfermat.
%
% Solves  n x (E + E_inc) = 0  on the surface of a torus,
% then evaluates the scattered field on a 2D slice.

zk = 7;
alpha = 0.5;
eps   = 1e-6;

% S = geometries.torus(1, 0.4, [3,3,3], [0;0;0], 8, 1);
% S = geometries.ellipsoid([1,1,1.5],[3,3,3], [0;0;0], 8);
S = geometries.ellipsoid([1,1,1.5],[1,1,1], [0;0;0], 8);

% % --- build system matrix (nrccie-bc kernel + 1/2 identity) ---
K   = kernel3d.em3d('nrccie-bc', zk, alpha);
tic;
Q = K.getquad(S,eps);
toc;
tic;
Q = K.getquad(S,eps,S);
toc;
tic;
A   = surfermat(S, K, eps);
A   = A + 0.5*eye(size(A));
toc;
% 
% % --- incident plane wave and nrccie RHS ---
% [einc, hinc] = em3d.planewave(zk, [pi/3; pi/4], [1; 0], S);
% rhs = nrccie_rhs(S, einc, hinc, alpha);
% 
% % --- solve ---
% sol = A \ rhs;
% 
% % --- recover Cartesian J and rho from orthonormal density ---
% [Jvec, rho] = nrccie_to_cartesian(S, sol);
% dens_eval = [Jvec; rho];   % (4, npts): [Jx; Jy; Jz; rho]

[einc, hinc] = em3d.planewave(zk, [pi/3; pi/4], [1; 0], S);

opts = [];
opts.eps_gmres = 1e-10;
tic;
[dens_eval] = em3d.pec.solver(S, einc, hinc, eps, zk, alpha, opts);
toc;

%%
% --- evaluate scattered field on a 2D slice ---
nplot = 200;
xx = linspace(-3, 3, nplot);
[XX, YY] = meshgrid(xx, xx);
targs.r = [XX(:).'; YY(:).'; zeros(1, nplot^2)];

Keval = kernel3d.em3d('nrccie-eval', zk);
opts.usematlab = 1;
tic;
Escat = surferkerneval(S, Keval, dens_eval(:), targs, eps,opts);
toc;
Escat = reshape(Escat, 6, []);

[Einc, Hinc] = em3d.planewave(zk, [pi/3; pi/4], [1; 0], targs);
Etot = Einc + Escat(1:3, :);
Htot = Hinc + Escat(4:6, :);

figure(1); clf
h = pcolor(XX, YY, reshape(vecnorm([Etot;Htot]), nplot, nplot));
h.EdgeColor = 'none';
colorbar; title('Re(E_x) total');
hold on; plot(S, zeros(S.npts,1)); hold off


% =========================================================================
% Helpers
% =========================================================================

function rhs = nrccie_rhs(S, einc, hinc, alpha)
% NRCCIE_RHS  Build the NRCCIE right-hand side in the orthonormal frame.
%   RHS rows 1-2: (n x H_inc - alpha*nxnxE_inc) . {ru, rv}
%   RHS row  3:    n . E_inc
n = S.n(:,:);   % (3, npts)

nxH = cross(n, hinc, 1);
ndotE = sum(n .* einc, 1);
nxnxE = ndotE .* n - einc;

zvec = nxH - alpha * nxnxE;   % (3, npts) Cartesian RHS for rows 1-2

du = S.du(:,:);
ru = du ./ vecnorm(du);
rv = cross(n ./ vecnorm(n), ru, 1);

rhs1 = sum(ru .* zvec, 1);   % (1, npts)
rhs2 = sum(rv .* zvec, 1);
rhs3 = ndotE;

rhs = [rhs1; rhs2; rhs3];
rhs = rhs(:);   % interleaved (3*npts, 1)
end


function [Jvec, rho] = nrccie_to_cartesian(S, sol)
% NRCCIE_TO_CARTESIAN  Convert orthonormal density [j_ru;j_rv;rho] to
%   Cartesian current Jvec (3,npts) and scalar rho (1,npts).
npts = S.npts;
sol3 = reshape(sol, 3, npts);
j_ru = sol3(1,:);
j_rv = sol3(2,:);
rho  = sol3(3,:);

n  = S.n(:,:);
du = S.du(:,:);
ru = du ./ vecnorm(du);
rv = cross(n ./ vecnorm(n), ru, 1);

Jvec = bsxfun(@times, ru, j_ru) + bsxfun(@times, rv, j_rv);   % (3, npts)
end
