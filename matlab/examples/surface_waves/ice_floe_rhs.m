function [rhs, ref_u] = ice_floe_rhs(geom, chnkr, kappa, nu, alpha, beta, gamma_pw, theta)
% ICE_FLOE_RHS  RHS and reference solution for a plane-wave manufactured solution.
%
% [rhs, ref_u] = ice_floe_rhs(geom, chnkr, kappa, nu, alpha, beta, a, b, gamma_pw, theta)
%
% Assembles the full right-hand side vector [rhs_grav; rhs_vol; rhs_bc]
% for the coupled gravity/flexural ice-floe system with manufactured
% solution u = exp(i*gamma_pw*(cos(theta)*x + sin(theta)*y)), and returns
% the reference solution evaluated at the volume nodes.
%
% Inputs:
%   geom      - surferfun object (volume discretization)
%   chnkr     - chunker object (boundary discretization)
%   kappa     - signed curvature at boundary nodes (chnkr.npt x 1)
%   nu        - Poisson ratio
%   alpha     - volume equation coefficient (coefficient of eta)
%   beta      - volume equation coefficient (coefficient of flex solution)
%   gamma_pw  - plane wave frequency (wavenumber magnitude)
%   theta     - plane wave propagation angle (radians)
%
% Outputs:
%   rhs    - right-hand side vector of length geom.npts + geom.npts + 2*chnkr.npt
%   ref_u  - reference solution at volume nodes (geom.npts x 1)

%% Reference solution and its derivatives at volume and boundary nodes
[ref_u, ~, ~, ~] = planewave2(gamma_pw, geom.r(1:2,:), theta);

[~, ~, hess, third] = planewave2(gamma_pw, chnkr.r(1:2,:), theta);

% Volume part
% Applying alpha*Delta^2 - beta to u = exp(i*k.r):
%   Delta u = -gamma_pw^2 * u  =>  a*Delta^2 u + b*Delta u = (a*gamma_pw^4 - b*gamma_pw^2) * u
rhs_grav = zeros(geom.npts, 1);
rhs_vol  = (alpha*gamma_pw^4 - beta) .* ref_u;

% Boundary part: free-plate BCs
% Normal and tangent vectors
nx   = chnkr.n(1,:).';
ny   = chnkr.n(2,:).';
dx   = chnkr.d(1,:).';
dy   = chnkr.d(2,:).';
ds   = sqrt(dx.*dx + dy.*dy);
taux = dx ./ ds;
tauy = dy ./ ds;

% Bending moment BC (M_nn + nu*M_tt = 0)
rhs_bc = zeros(2*chnkr.npt, 1);
rhs_bc(1:2:end) = ...
    (hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
    nu.*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy));

% Kirchhoff shear BC (V_n = 0)
rhs_bc(2:2:end) = ...
    (third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) + ...
     third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny)) + ...
    (2-nu).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) + ...
             third(:,3).*(2*taux.*tauy.*ny + tauy.*tauy.*nx) + third(:,4).*(tauy.*tauy.*ny)) + ...
    (1-nu).*kappa(:).*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy ...
                     - hess(:,1).*nx.*nx    - hess(:,2).*(2*nx.*ny)      - hess(:,3).*ny.*ny);

rhs = -[rhs_grav; rhs_vol; rhs_bc];

end


function [r1,grad,hess,third] = planewave2(k,r,theta)

k1 = k*cos(theta);
k2 = k*sin(theta);

r1 = exp(1i*(k1*r(1,:)+k2*r(2,:))).';

gx = r1.*(1i*k1);
gy = r1.*(1i*k2);
grad = [gx, gy];

hxx = r1.*(1i*k1).^2; 
hxy = r1.*(1i*k1)*(1i*k2); 
hyy = r1.*(1i*k2).^2;
hess = [hxx, hxy, hyy];

txxx = r1.*(1i*k1).^3;
txxy = r1.*(1i*k1).^2*(1i*k2); 
txyy = r1.*(1i*k1)*(1i*k2).^2; 
tyyy = r1.*(1i*k2).^3; 
third = [txxx, txxy, txyy, tyyy];

end