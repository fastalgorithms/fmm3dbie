%
% This file tests the Stokes single layer
% potenial on the sphere
%
%
run ../startup.m

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
tic, [srcvals,~,~,~,~,wts] = extract_arrays(S); toc;


%% Now test the solver + kernel evaluation routines

xyz_in = [0.3;0.5;0.1];
xyz_out = [1.3;-5.2;0.1];
src_info = [];
src_info.r = xyz_in;
src_info.v = [1;-2;3];
rhs = sum(stok3d.kern(src_info,S,'s').*reshape(src_info.v,[1,1,3,1]),3);  

eps = 1e-7;
dpars = [1,1];
sig = stok3d.velocity.solver(S,rhs,eps,dpars);

targ_info = [];
targ_info.r = xyz_out;
targ_info.patch_id = -1;
dat = stok3d.kern(S,targ_info,'c',dpars(1),dpars(2));
pot = sum(dat.*(reshape(sig.*wts.',[1,1,3,S.npts])),[3,4]);
pot_eval = stok3d.velocity.eval(S,sig,targ_info,eps,dpars);
pot_ex = sum(stok3d.kern(src_info,targ_info,'s').*reshape(src_info.v,[1,1,3]),3);
fprintf('Error in iterative solver=%d\n',vecnorm(pot-pot_ex)/vecnorm(pot_ex));
fprintf('Error in stok3d_eval=%d\n',vecnorm(pot_eval-pot_ex)/vecnorm(pot_ex));

%% Near quadrature for single- and double-layer traction

%
% This file tests the near quadrature for the Stokes traction kernels
%
%     S' = kernel3d('stok3d', 'strac')
%     D' = kernel3d('stok3d', 'dtrac')
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
tD = surferkerneval(S, kernel3d('stok3d', 'dtrac'), u(:), targ, eps);
t = reshape(tS - tD, 3, []);

tex = stoklet_t(targ);
err = max(vecnorm(t - tex)) / max(vecnorm(tex));
fprintf('traction Green''s identity at distance %.0e: rel err %.2e\n', d, err);
assert(err < 1e-4, 'traction Green''s identity failed (%.3e)', err);
