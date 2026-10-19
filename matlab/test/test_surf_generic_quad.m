% Test the name-dispatched near quadrature: getnearquad_mex, the kernel
% registry, the kernel3d getquad wiring and the Stokes 'sparse' format.

%% Now run the tests

test_registry();
test_matches_wrappers();
test_formats();
test_errors();
test_stokes_sparse();


function test_registry()
% Every kernel3d getquad name is registered, with the right real/complex flag.

names = list_quad_kernels();
reals = {'lap3d_s','lap3d_d','lap3d_sp','lap3d_c','lap3d_cp','lap3d_cp_pv', ...
         'belpde_klb','belpde_rlb'};
cplxs = {'helm3d_s','helm3d_d','helm3d_sp','helm3d_c','helm3d_cp', ...
         'helm3d_cp_pv','belpde_khb','belpde_rhb'};
assert(all(ismember([reals, cplxs], names)), 'registry: missing kernel');

S = geometries.sphere(1, 1, [0;0;0], 4, 1);
t = struct('r', 1.5*S.r(:,1:5));
for k = 1:numel(reals)
    if startsWith(reals{k}, 'belpde'), continue; end
    Q = getnearquad_mex(S, reals{k}, 1e-6, struct('dpars',[1;1]), t);
    assert(Q.ifcomplex == 0 && isreal(Q.wnear), 'registry: %s not real', reals{k});
end
for k = 1:numel(cplxs)
    if startsWith(cplxs{k}, 'belpde'), continue; end
    Q = getnearquad_mex(S, cplxs{k}, 1e-6, struct('zpars',[1.1;1;1]), t);
    assert(Q.ifcomplex == 1, 'registry: %s not complex', cplxs{k});
end

end


function test_matches_wrappers()
% kernel3d getquad (now through getnearquad_mex) reproduces the
% kernel-specific get_quadrature_correction wrappers exactly.

S   = geometries.sphere(1, 2, [0;0;0], 4, 1);
eps = 1e-6;
zk  = 1.3;
o1  = struct('iprime', 1);
old = @(Q) conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear);
err = @(A,B) full(max(abs(A(:)-B(:)))/max(abs(B(:))));
tol = 1e-14;

c = [1.2; -0.4];
cases = { ...
  kernel3d.lap3d('s'),    old(lap3d.dirichlet.get_quadrature_correction(S,eps,[1;0])), 'lap s'; ...
  kernel3d.lap3d('d'),    old(lap3d.dirichlet.get_quadrature_correction(S,eps,[0;1])), 'lap d'; ...
  kernel3d.lap3d('c',c),  old(lap3d.dirichlet.get_quadrature_correction(S,eps,c)), 'lap c'; ...
  kernel3d.lap3d('sp'),   old(lap3d.dirichlet.get_quadrature_correction(S,eps,[1;0],S,o1)), 'lap sp'; ...
  kernel3d.lap3d('dp'),   old(lap3d.dirichlet.get_quadrature_correction(S,eps,[0;1],S,o1)), 'lap dp'; ...
  kernel3d.lap3d('cp',c), old(lap3d.dirichlet.get_quadrature_correction(S,eps,c,S,o1)), 'lap cp'; ...
  kernel3d.lap3d('cp',[1;0]), old(lap3d.dirichlet.get_quadrature_correction(S,eps,[1;0],S,o1)), 'lap cp pv'; ...
  kernel3d.helm3d('s',zk),   old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1;0])), 'helm s'; ...
  kernel3d.helm3d('d',zk),   old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0;1])), 'helm d'; ...
  kernel3d.helm3d('c',zk,[1;-0.5i]), old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1;-0.5i])), 'helm c'; ...
  kernel3d.helm3d('sp',zk),  old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1;0],S,o1)), 'helm sp'; ...
  kernel3d.helm3d('dp',zk),  old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[0;1],S,o1)), 'helm dp'; ...
  kernel3d.helm3d('cp',zk,c), old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,c,S,o1)), 'helm cp'; ...
  kernel3d.helm3d('cp',zk,[1;0]), old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,[1;0],S,o1)), 'helm cp pv'};

for k = 1:size(cases,1)
    A = cases{k,1}.getquad(S, eps);
    assert(err(A, cases{k,2}) < tol, 'getquad mismatch: %s', cases{k,3});
end

% off-surface targets
t = struct('r', 1.05*S.r(:,1:40));
A = kernel3d.helm3d('c',zk,c).getquad(S, eps, t);
B = old(helm3d.dirichlet.get_quadrature_correction(S,eps,zk,c,t));
assert(err(A, B) < tol, 'getquad mismatch: helm c off-surface');

end


function test_formats()
% rsc, csc and sparse outputs describe the same quadrature.

S = geometries.sphere(1, 2, [0;0;0], 4, 1);
p = struct('zpars', [1.2; 1; 0.5]);
Qr = getnearquad_mex(S, 'helm3d_c', 1e-6, p);
Qc = getnearquad_mex(S, 'helm3d_c', 1e-6, p, [], struct('format','csc'));
Qs = getnearquad_mex(S, 'helm3d_c', 1e-6, p, [], struct('format','sparse'));

assert(isequal(Qr.wnear, Qc.wnear), 'formats: csc wnear differs');
A = conv_rsc_to_spmat(S, Qr.row_ptr, Qr.col_ind, Qr.wnear);
assert(isequal(size(Qs.spmat), [S.npts S.npts]), 'formats: sparse size');
assert(full(max(abs(A(:) - Qs.spmat(:)))) == 0, 'formats: sparse differs');

% csc lists the same (target, patch) pairs as rsc
trsc = repelem((1:S.npts)', diff(Qr.row_ptr));
pcsc = repelem((1:S.npatches)', diff(Qc.col_ptr));
assert(isequal(sortrows([trsc, Qr.col_ind(:)]), sortrows([Qc.row_ind(:), pcsc])), ...
    'formats: csc structure differs');

% off-surface targets give a rectangular matrix
t  = struct('r', 1.05*S.r(:,1:40));
Qr = getnearquad_mex(S, 'helm3d_c', 1e-6, p, t);
Qs = getnearquad_mex(S, 'helm3d_c', 1e-6, p, t, struct('format','sparse'));
A  = conv_rsc_to_spmat(S, Qr.row_ptr, Qr.col_ind, Qr.wnear);
assert(isequal(size(Qs.spmat), [40 S.npts]), 'formats: off-surface sparse size');
assert(full(max(abs(A(:) - Qs.spmat(:)))) == 0, 'formats: off-surface sparse differs');

% real kernels give real sparse matrices
Ql = getnearquad_mex(S, 'lap3d_s', 1e-6, [], [], struct('format','sparse'));
assert(isreal(Ql.spmat), 'formats: real kernel gave complex sparse');

end


function test_errors()
% Unknown names and missing options raise MATLAB errors (no Fortran stop).

S = geometries.sphere(1, 2, [0;0;0], 4, 1);

errored = false;
try
    getnearquad_mex(S, 'no_such_kernel', 1e-6);
catch
    errored = true;
end
assert(errored, 'errors: unknown kernel did not error');

errored = false;
try
    getnearquad_mex(S, 'lap3d_s', 1e-6, [], [], struct('corrections', true));
catch
    errored = true;
end
assert(errored, 'errors: missing nover did not error');

errored = false;
try
    getnearquad_mex(S, 'lap3d_s', 1e-6, [], [], ...
        struct('corrections', true, 'nover', [5 6]));
catch
    errored = true;
end
assert(errored, 'errors: bad nover did not error');

end


function test_stokes_sparse()
% format 'sparse' in the Stokes wrappers uses the symmetric 3x3 layout.

S    = geometries.sphere(1, 2, [0;0;0], 3, 1);
eps  = 1e-6;
ri   = kernel3d.rsc_interleave_symmetric3();
err  = @(A,B) full(max(abs(A(:)-B(:)))/max(abs(B(:))));

dpars = [1; 0.5];
Q  = stok3d.velocity.get_quadrature_correction(S, eps, dpars);
Qs = stok3d.velocity.get_quadrature_correction(S, eps, dpars, S, struct('format','sparse'));
A  = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ri);
assert(isequal(size(Qs.spmat), [3*S.npts 3*S.npts]), 'stokes velocity: sparse size');
assert(err(Qs.spmat, A) == 0, 'stokes velocity: sparse differs');

Q  = stok3d.traction.get_quadrature_correction(S, eps, []);
Qs = stok3d.traction.get_quadrature_correction(S, eps, [], S, struct('format','sparse'));
A  = conv_rsc_to_spmat(S, Q.row_ptr, Q.col_ind, Q.wnear, ri);
assert(err(Qs.spmat, A) == 0, 'stokes traction: sparse differs');

end

