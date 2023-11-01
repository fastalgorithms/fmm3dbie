%
% This file illustrates the application of 
% S_{k}[n \times \nabla_{\Gamma} \Delta_{\Gamma}^{-1}[Ynm]] 
% on the sphere.
%
% Note as things are setup, k cannot be set to 0 or
% very small quantities, the laplace version of routines
% will be made available soon.
%
% The computation of these quantities naturally arises in the
% computation of the surface currents for the generalized debye
% formulation.
%
% \nabla_{\Gamma} is the surface gradient
% Delta_{\Gamma}^{-1} is the surface Laplacian inverse computed using
% dan's code, \bn is the normal to the surface, and Ynm is the spherical
% harmonic of order n and degree l

addpath(genpath('~/git/fmm3dbie/matlab'))
run ~/git/surface-hps/setup.m
%% initialize surface of the sphere

iref = 1;
domtmp = surfacemesh.sphere(9,iref);
[S,dom] = surfer.surfacemesh_to_surfer(domtmp);

plot(dom);

%% Construct the inverse of the Laplace Beltrami operator
pdo = [];
pdo.lap = 1;
n = normal(dom);
L = surfaceop(dom, pdo);
L.rankdef = true;
tic, L.build(); t1 = toc;
fprintf('Time taken to build surface Laplacian inverse=%d\n',t1);

%% Set boundary data and verify Delta^{\Gamma}^{-1} Y_{nm} = -Y_{nm}/(n(n+1))
ndeg = 2;

f = spherefun.sphharm(ndeg,1);

frhs = surfacefun(@(x,y,z) f(x,y,z),dom);
L.rhs = frhs;
u = solve(L);
err1 = norm(u + frhs./(ndeg*(ndeg+1)));
fprintf('Error in Laplace beltami solve=%d\n',err1);

% Compute n \times surface grad
gradu = grad(u);
gradu = cross(n, gradu);

% Convert it to flat array
gradu_flat = complex(surfacefun_to_array(gradu,dom,S));

%% Setup computation of S_{k}[J] using 

% Set wavenumber and zpars for single layer evaluation
zk = 1.1;
zpars = complex([zk, 1.0, 0.0]);
rfac = 1j*zk*jn(ndeg)*hn(ndeg);

% multipliers required for analytic solution
jn = @(n) sqrt(pi/2/zk)*besselj(n+0.5,zk);
hn = @(n) sqrt(pi/2/zk)*besselh(n+0.5,1,zk);


eps = 1e-8;

% Precompute quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
Q = helm3d.dirichlet.get_quadrature_correction(S,zpars, ...
   eps,S,opts_quad);

Skj = zeros(size(gradu_flat));
for j=1:3
    Skj(:,j) = helm3d.dirichlet.eval(S,zpars,gradu_flat(:,j),eps,S,Q);
end


Skj_ex = gradu_flat.*rfac;

wts3 = repmat(S.wts,[1,3]);
err1 = sqrt(sum(abs(Skj_ex-Skj).^2.*wts3)/(sum(abs(Skj_ex).^2.*wts3)));
fprintf('error in computed gradient in surfer land=%d\n',err1);
