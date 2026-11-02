% Test the near quadrature corrections of the package routines,
% get_quadrature_correction with opts.corrections:
%   - through the getquad of the kernel3d kernels, against the near
%     quadrature minus smooth_sparse_quad
%   - as the precomputed quadrature of the evaluators and solvers,
%     against their add and subtract path

%% Now run the tests

test_kernel_corrections();
test_eval_solver_corrections();


function test_kernel_corrections()
% getquad with opts.corrections equals getquad minus smooth_sparse_quad.

S    = geometries.sphere(1, 1, [0;0;0], 4, 1);
eps  = 1e-7;
tol  = 1e-11;
zk   = 1.1;
nvec = 5 + mod(1:S.npatches, 3).';
toff = struct('r', 1.05*S.r(:,1:40), 'n', S.n(:,1:40), ...
    'du', S.du(:,1:40), 'dv', S.dv(:,1:40));

kerns = {kernel3d.lap3d('s'), kernel3d.lap3d('d'), kernel3d.lap3d('sp'), ...
         kernel3d.lap3d('dp'), kernel3d.lap3d('c', [1; 0.5]), ...
         kernel3d.lap3d('cp', [1; 0.5]), ...
         kernel3d.helm3d('s', zk), kernel3d.helm3d('d', zk), ...
         kernel3d.helm3d('sp', zk), kernel3d.helm3d('dp', zk), ...
         kernel3d.helm3d('c', zk, [1; -1i]), ...
         kernel3d.helm3d('cp', zk, [1; -1i]), ...
         kernel3d.stok3d('s'), kernel3d.stok3d('d'), kernel3d.stok3d('sp'), ...
         kernel3d.stok3d('c', [1; 0.5]), ...
         kernel3d.em3d('nrccie-bc', zk, 1.0), kernel3d.em3d('nrccie-eval', zk)};
for k = 1:numel(kerns)
    check_kernel(kerns{k}, S, {S, toff}, {7, nvec}, eps, tol);
end

% the near quadrature of stok3d('cp') is not that of its eval (the
% traction of the double layer differs), so its getquad does not
% compute the corrections
K = kernel3d.stok3d('cp', [1; 0.5]);
assert(~K.getquad_sub, 'kernel corrections: stokes ctrac flag');

% Beltrami kernels, on surface only
kerns = {kernel3d.belpde('klb'), kernel3d.belpde('rlb'), ...
         kernel3d.belpde('khb', zk), kernel3d.belpde('rhb', zk)};
for k = 1:numel(kerns)
    check_kernel(kerns{k}, S, {S}, {7, nvec}, eps, tol);
end

end


function check_kernel(K, S, targs, novers, eps, tol)

label = [K.name ' ' K.type];
assert(K.getquad_sub, 'kernel corrections: %s flag', label);
ri = kernel3d.rsc_interleave_full(K.opdims(1), K.opdims(2));
for t = targs
    rsc = getnear(S, t{1});
    A   = K.getquad(S, eps, t{1});
    for nv = novers
        ws   = smooth_sparse_quad(K, t{1}, S, rsc.row_ptr, rsc.col_ind, nv{1});
        Aref = A - conv_rsc_to_spmat(S, rsc.row_ptr, rsc.col_ind, ws, ri);
        Asub = K.getquad(S, eps, t{1}, ...
            struct('corrections', true, 'nover', nv{1}));
        err  = full(max(abs(Asub(:) - Aref(:)))) / full(max(abs(Aref(:))));
        assert(err < tol, 'kernel corrections: %s, error %g', label, err);
    end
end

end


function test_eval_solver_corrections()
% The evaluators and solvers given the corrections, computed with the
% default oversampling, reproduce their add and subtract path. The
% solvers iterate to a small residual, so that the solutions differ by
% the rounding of the two paths, which these coarse discretizations
% amplify by up to 1e5 (a 1e-14 change of the data changes the
% solutions by up to 3e-9), while an inconsistent correction changes
% them by the quadrature error at least.

S     = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps   = 1e-7;
zk    = 1.1;
tole  = 1e-10;
tols  = 1e-7;
og    = struct('eps_gmres', 1e-12);
qc    = @(Q) struct('eps_gmres', 1e-12, 'quadrature_correction', Q);
src   = struct('r', [0.3; 0.5; 0.1]);
idx   = 1:7:S.npts;
toff  = struct('r', S.r(:,idx) + 0.02*S.n(:,idx), 'n', S.n(:,idx));
oc    = struct('corrections', true);

% Laplace, Dirichlet
dpars = [1, 1];
rhs = lap3d.kern(src, S, 's');
Q = lap3d.dirichlet.get_quadrature_correction(S, eps, dpars, S, oc);
s0 = lap3d.dirichlet.solver(S, rhs, eps, dpars, og);
s1 = lap3d.dirichlet.solver(S, rhs, eps, dpars, ...
    qc(Q));
check(s0, s1, tols, 'lap3d.dirichlet.solver');
Q = lap3d.dirichlet.get_quadrature_correction(S, eps, dpars, toff, oc);
p0 = lap3d.dirichlet.eval(S, s0, toff, eps, dpars);
p1 = lap3d.dirichlet.eval(S, s0, toff, eps, dpars, ...
    struct('precomp_quadrature', Q));
check(p0, p1, tole, 'lap3d.dirichlet.eval');

% Laplace, Neumann
rhs = lap3d.kern(src, S, 'sprime');
Q = lap3d.neumann.get_quadrature_correction(S, eps, S, ...
    struct('corrections', true, 'rep', 'bc'));
s0 = lap3d.neumann.solver(S, rhs, eps, og);
s1 = lap3d.neumann.solver(S, rhs, eps, qc(Q));
check(s0, s1, tols, 'lap3d.neumann.solver');
Q = lap3d.neumann.get_quadrature_correction(S, eps, toff, ...
    struct('corrections', true, 'rep', 'eval'));
p0 = lap3d.neumann.eval(S, s0, toff, eps);
p1 = lap3d.neumann.eval(S, s0, toff, eps, struct('precomp_quadrature', Q));
check(p0, p1, tole, 'lap3d.neumann.eval');

% Helmholtz, Dirichlet
rep_pars = [-1j*zk, 1];
rhs = helm3d.kern(zk, src, S, 's');
Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, S, oc);
s0 = helm3d.dirichlet.solver(S, rhs, eps, zk, rep_pars, og);
s1 = helm3d.dirichlet.solver(S, rhs, eps, zk, rep_pars, ...
    qc(Q));
check(s0, s1, tols, 'helm3d.dirichlet.solver');
Q = helm3d.dirichlet.get_quadrature_correction(S, eps, zk, rep_pars, toff, oc);
p0 = helm3d.dirichlet.eval(S, s0, toff, eps, zk, rep_pars);
p1 = helm3d.dirichlet.eval(S, s0, toff, eps, zk, rep_pars, ...
    struct('precomp_quadrature', Q));
check(p0, p1, tole, 'helm3d.dirichlet.eval');

% Helmholtz, Neumann
alpha = 1.0;
rhs = helm3d.kern(zk, src, S, 'sprime');
Q = helm3d.neumann.get_quadrature_correction(S, eps, zk, alpha, S, ...
    struct('corrections', true, 'rep', 'rpcomb-bc'));
d0 = helm3d.neumann.solver(S, rhs, eps, zk, alpha, og);
d1 = helm3d.neumann.solver(S, rhs, eps, zk, alpha, ...
    qc(Q));
check(d0, d1, tols, 'helm3d.neumann.solver');
Q = helm3d.neumann.get_quadrature_correction(S, eps, zk, alpha, toff, ...
    struct('corrections', true, 'rep', 'rpcomb-eval'));
p0 = helm3d.neumann.eval(S, d0, toff, eps, zk, alpha);
p1 = helm3d.neumann.eval(S, d0, toff, eps, zk, alpha, ...
    struct('precomp_quadrature', Q));
check(p0, p1, tole, 'helm3d.neumann.eval');

% Helmholtz, impedance
zlams = complex(1 + 0.3*cos(S.r(1,:).'));
rhs = helm3d.kern(zk, src, S, 'sprime') + ...
    1j*zk.*zlams.*helm3d.kern(zk, src, S, 's');
Q = helm3d.impedance.get_quadrature_correction(S, eps, zk, alpha, S, ...
    struct('corrections', true, 'rep', 'rpcomb-bc'));
d0 = helm3d.impedance.solver(S, zlams, rhs, eps, zk, alpha, og);
d1 = helm3d.impedance.solver(S, zlams, rhs, eps, zk, alpha, ...
    qc(Q));
check(d0, d1, tols, 'helm3d.impedance.solver');
Q = helm3d.impedance.get_quadrature_correction(S, eps, zk, alpha, toff, ...
    struct('corrections', true, 'rep', 'rpcomb-eval'));
p0 = helm3d.impedance.eval(S, d0, toff, eps, zk, alpha);
p1 = helm3d.impedance.eval(S, d0, toff, eps, zk, alpha, ...
    struct('precomp_quadrature', Q));
check(p0, p1, tole, 'helm3d.impedance.eval');

% Helmholtz, transmission
zks = [1.1; 2.1];
rep_params = [1; 1; 1; 1];
rhs = complex(zeros(2, S.npts));
rhs(1,:) = helm3d.kern(zks(1), src, S, 's');
rhs(2,:) = helm3d.kern(zks(1), src, S, 'sprime');
Q = helm3d.transmission.get_quadrature_correction(S, eps, zks, rep_params, ...
    S, struct('corrections', true, 'rep', 'trans-bc'));
d0 = helm3d.transmission.solver(S, rhs, eps, zks, rep_params, og);
d1 = helm3d.transmission.solver(S, rhs, eps, zks, rep_params, ...
    qc(Q));
check(d0, d1, tols, 'helm3d.transmission.solver');

% Stokes, velocity and traction
dpars = [1, 1];
sv  = struct('r', src.r, 'v', [1; -2; 3]);
rhs = sum(stok3d.kern(sv, S, 's').*reshape(sv.v, [1,1,3,1]), 3);
Q = stok3d.velocity.get_quadrature_correction(S, eps, dpars, S, oc);
s0 = stok3d.velocity.solver(S, rhs, eps, dpars, og);
s1 = stok3d.velocity.solver(S, rhs, eps, dpars, ...
    qc(Q));
check(s0, s1, tols, 'stok3d.velocity.solver');
Q = stok3d.velocity.get_quadrature_correction(S, eps, dpars, toff, oc);
p0 = stok3d.velocity.eval(S, s0, toff, eps, dpars);
p1 = stok3d.velocity.eval(S, s0, toff, eps, dpars, ...
    struct('precomp_quadrature', Q));
check(p0, p1, tole, 'stok3d.velocity.eval');
Q = stok3d.traction.get_quadrature_correction(S, eps, [], toff, oc);
p0 = stok3d.traction.eval(S, s0, toff, eps);
p1 = stok3d.traction.eval(S, s0, toff, eps, struct('precomp_quadrature', Q));
check(p0, p1, tole, 'stok3d.traction.eval');

% Maxwell, perfect conductor and dielectric
alpha = 1.0;
se = struct('r', [3.1; 1.0; 1.0], 'edips', [1; 2; -1] + 1i*[0.5; -1; 0.2], ...
    'hdips', [-1; 0.5; 2] + 1i*[0.3; 0.1; -0.4]);
[Einc, Hinc] = em3d.incoming_sources(zk, se, S, 'ehd');
Q = em3d.pec.get_quadrature_correction(S, eps, zk, alpha, S, ...
    struct('corrections', true, 'rep', 'nrccie-bc'));
d0 = em3d.pec.solver(S, Einc, Hinc, eps, zk, alpha, og);
d1 = em3d.pec.solver(S, Einc, Hinc, eps, zk, alpha, ...
    qc(Q));
check(d0, d1, tols, 'em3d.pec.solver');
Q = em3d.pec.get_quadrature_correction(S, eps, zk, alpha, toff, ...
    struct('corrections', true, 'rep', 'nrccie-eval'));
[E0, H0] = em3d.pec.eval(S, d0, toff, eps, zk, alpha);
[E1, H1] = em3d.pec.eval(S, d0, toff, eps, zk, alpha, ...
    struct('precomp_quadrature', Q));
check([E0; H0], [E1; H1], tole, 'em3d.pec.eval');

om = 1.1;
rep_params = complex([1.4, 1, 1, 1]);
[Einc, Hinc] = em3d.incoming_sources(om, se, S, 'ehd');
Q = em3d.dielectric.get_quadrature_correction(S, eps, om, rep_params, oc);
d0 = em3d.dielectric.solver(S, Einc, Hinc, eps, om, rep_params, og);
d1 = em3d.dielectric.solver(S, Einc, Hinc, eps, om, rep_params, ...
    qc(Q));
check(d0, d1, tols, 'em3d.dielectric.solver');

end


function check(a, b, tol, label)

err = max(abs(a(:) - b(:))) / max(abs(a(:)));
assert(err < tol, '%s with corrections: error %g', label, err);

end
