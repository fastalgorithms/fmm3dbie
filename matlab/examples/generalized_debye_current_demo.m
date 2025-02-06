%
% This file illustrates the application of 
% \nabla_{\Gamma} \Delta_{\Gamma}^{-1} S_{k}[Ynm] 
% and \bn \times \nabla_{\Gamma} \Delta_{\Gamma}^{-1} S_{k}[Ynm]
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
%% initialize surface of the sphere

iref = 1;
domtmp = surfacemesh.sphere(9,iref);
[S,dom] = surfer.surfacemesh_to_surfer(domtmp);

plot(dom);

%% Construct the inverse of the Laplace Beltrami operator
pdo = [];
pdo.lap = 1;
n = normal(dom);
L = surfaceop(dom,pdo);
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

%% Evaluate rhs on surfer object

% Note that the rhs for surfer objects is evaluated at a 
% different set of nodes
rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = f(S.r(1,:)./rr,S.r(2,:)./rr,S.r(3,:)./rr);
rhs = rhs(:);

%% Now test conversion from flat array to surfacefun

frhs2 = array_to_surfacefun(rhs,dom,S);
err1 = norm(frhs-frhs2)/norm(frhs);
fprintf('Error in interpolant =%d\n',err1);

%% Setup computation of S_{k}[Ynm] using 

% Set wavenumber and zpars for single layer evaluation
zk = 1.1;
zpars = complex([zk, 1.0, 0.0]);


% multipliers required for analytic solution
jn = sqrt(pi/2/zk)*besselj(ndeg+0.5,zk);
hn = sqrt(pi/2/zk)*besselh(ndeg+0.5,1,zk);

eps = 1e-12;

% Precompute quadrature corrections
opts_quad = [];
opts_quad.format='rsc';
rep_pars = [1.0; 0.0];
Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, ...
   S,opts_quad);

opts_eval = [];
opts_eval.precomp_quadrature = Q;
% Now evaluate potential
p = helm3d.dirichlet.eval(S,rhs,S,eps,zk,rep_pars,opts_eval); 

% The exact potential
p_ex = rhs*jn*hn*1j*zk;

[srcvals,~,~,~,~,wts] = extract_arrays(S); 

err1 = norm((p-p_ex).*sqrt(wts))/norm(rhs.*sqrt(wts));


fprintf('Error in single layer potential with precomp corr=%d\n',err1);

%% Now evaluate \Delta_{\Gamma}^{-1} S_{k} [Ynm]

prhs = array_to_surfacefun(p,dom,S);
L.rhs = prhs;
u2 = solve(L);
u2ex = -frhs*(jn*hn*1j*zk)/(ndeg*(ndeg+1));
err1 = norm(u2-u2ex)/norm(frhs);
fprintf("Error in evaluating R^{-1} S_{k}[Ynm]=%d\n",err1);


%% Compute surface grad \nabla_{\Gamma} \Delta_{\Gamma}^{-1} S_{k}[Ynm]

% first in surfacefun land
gradu2 = grad(u2);

% convert surfacefunv to flat array
gradu = surfacefun_to_array(gradu2,dom,S);

% Next compute the gradient in surfer land

rhs_ex = -p_ex/(ndeg*(ndeg+1));
gradu_ex = get_surface_grad(S,rhs_ex);
gradu_ex = gradu_ex.';

wts3 = repmat(wts,[1,3]);
err1 = sqrt(sum(abs(gradu_ex-gradu).^2.*wts3)/(sum(abs(gradu_ex).^2.*wts3)));
fprintf('error in computed gradient in surfer land=%d\n',err1);

% Convert it back to a surfacefunv and compare
gradu2_ex = array_to_surfacefun(gradu_ex,dom,S);
err1 = norm(norm(gradu2-gradu2_ex))./norm(norm(gradu2_ex));
fprintf('error in computed gradient in surfacefun land=%d\n',err1);





