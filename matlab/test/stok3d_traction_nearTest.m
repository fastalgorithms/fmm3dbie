%
% This file tests the near quadrature for the Stokes traction kernels
%
%     S' = kernel3d('stok3d', 'strac')
%     D' = kernel3d('stok3d', 'ctrac', [0; 1])
%
% at targets close to the unit sphere. For a Stokes flow u inside, with
% surface traction f = sigma(u) n (n the OUTWARD normal, mu = 1), the
% representation u = S[f] - D[u] gives for the traction at x inside, with
% normal m,
%
%     t(x) = S'[f](x) - D'[u](x).
%
% The flow is that of a Stokeslet outside the sphere. The targets lie at
% distance d inside the surface, with random normals.
%

run ../startup.m

eps = 1e-8;
d = 1e-2;
S = geometries.sphere(1, 3, [0;0;0], 12, 1);

% Stokeslet outside the unit sphere
src = [];  src.r = [1.6; -0.9; 1.1];  f0 = [0.7; -1.1; 0.5];
stoklet_u = @(t) reshape(sum(stok3d.kern(src, t, 's') .* ...
    reshape(f0, [1,1,3]), 3), 3, []);
stoklet_t = @(t) reshape(sum(stok3d.kern(src, t, 'sprime') .* ...
    reshape(f0, [1,1,3]), 3), 3, []);
u = stoklet_u(S);
f = stoklet_t(S);

% targets at distance d inside the sphere, with random normals
rng(1);
nt = 20;
x = randn(3, nt);  x = (1 - d) * x ./ vecnorm(x);
m = randn(3, nt);  m = m ./ vecnorm(m);
targ = struct('r', x, 'n', m);

tS = surferkerneval(S, kernel3d('stok3d', 'strac'), f(:), targ, eps);
tD = surferkerneval(S, kernel3d('stok3d', 'ctrac', [0; 1]), u(:), targ, eps);
t = reshape(tS - tD, 3, []);

tex = stoklet_t(targ);
err = max(vecnorm(t - tex)) / max(vecnorm(tex));
fprintf('traction Green''s identity at distance %.0e: rel err %.2e\n', d, err);
assert(err < 1e-4, 'traction Green''s identity failed (%.3e)', err);
