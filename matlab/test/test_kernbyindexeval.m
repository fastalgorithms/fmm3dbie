% Test byindex.kernbyindexeval against surferkernevalmat.

run ../startup.m

rng(42);

%% Now run the tests

test_single_surfer_laplace();
test_two_surfer_helmholtz();
test_extra_targ_field();
test_oversampled_decomposition();
test_two_surfer_oversampled_decomposition();


function test_single_surfer_laplace()
% Single surfer, Laplace SLP: kernbyindexeval matches surferkernevalmat.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 1, [0;0;0], 4, 1);
kern = kernel3d('l', 's');

rng(1);
ntarg = 50;
targs.r = 3 * (2*rand(3,ntarg) - 1);
i = sort(randperm(ntarg,  20)).';
j = sort(randperm(S.npts, 30)).';

check_eval('single-surfer Laplace', i, j, S, kern, targs, eps, tol);

end


function test_two_surfer_helmholtz()
% Two surfers, Helmholtz: kernbyindexeval matches surferkernevalmat.

eps = 1e-10;  tol = 1e-9;

zk = 1.3 + 0.1i;
S1 = geometries.sphere(2, 1, [0;0;0], 3, 1);
S2 = geometries.sphere(2, 1, [6;0;0], 3, 1);
kern = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

rng(2);
ntarg = 40;
targs.r = [3*(2*rand(3,ntarg/2)-1), [6;0;0] + 3*(2*rand(3,ntarg/2)-1)];
i = sort(randperm(ntarg,          15)).';
j = sort(randperm(S1.npts+S2.npts, 25)).';

check_eval('two-surfer Helmholtz', i, j, [S1,S2], kern, targs, eps, tol);

end


function test_extra_targ_field()
% Kernel declaring a target field beyond the geometry: the targets carry
% mean_curv as an ntarg x 1 column, which must be sliced by point index.

eps = 1e-10;  tol = 1e-9;

S    = geometries.ellipsoid([1,1,1.2], [2,2,2], [], 5);
kern = kernel3d('bel', 'rlb');

ntarg = 120;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
targs = [];
targs.r = (2 + rand(1,ntarg)) .* [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.n = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.mean_curv = rand(ntarg,1);

i = randi(ntarg, 60, 1);
j = randi(S.npts, 60, 1);

check_eval('extra targ field', i, j, S, kern, targs, eps, tol);

end


function test_oversampled_decomposition()
% Targets close to the surface force oversampling. The full matrix is the
% smooth rule on the oversampled surfer times the interpolation matrix,
% plus corrections=1 output from surferkernevalmat.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 1, [0;0;0], 4, 1);
kern = kernel3d('l', 's');

rng(4);
ntarg = 40;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
targs.r = (1 + 0.05*rand(1,ntarg)) .* [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];

[Qc, objover] = surferkernevalmat(S, kern, targs, eps, struct('corrections', 1));
Afull = surferkernevalmat(S, kern, targs, eps);
Sover   = objover{1}{1};
xinterp = objover{2}{1};

assert(Sover.npts > S.npts, 'expected the targets to force oversampling');

i = sort(randperm(ntarg,  15)).';
j = sort(randperm(S.npts, 20)).';

Aover = byindex.kernbyindexeval(i, 1:size(xinterp,1), Sover, kern, targs);
Atest = Aover * xinterp(:, j) + Qc(i, j);
err = norm(Atest - Afull(i,j), 'fro') / norm(Afull(i,j), 'fro');
assert(err < tol, 'oversampled decomposition: err %.2e > tol %.2e', err, tol);

end


function test_two_surfer_oversampled_decomposition()
% Same as test_oversampled_decomposition with two surfers and a kernel
% per surfer: the oversampled surfers are concatenated and the
% interpolation matrices are block diagonal.

eps = 1e-10;  tol = 1e-9;

zk = 1.1;
S1 = geometries.sphere(1, 1, [0;0;0], 4, 1);
S2 = geometries.sphere(1, 1, [3;0;0], 4, 1);
kern = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

rng(5);
ntarg = 40;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
dirs = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.r = [(1 + 0.05*rand(1,ntarg/2)) .* dirs(:,1:ntarg/2), ...
           [3;0;0] + (1 + 0.05*rand(1,ntarg/2)) .* dirs(:,ntarg/2+1:end)];

[Qc, objover] = surferkernevalmat([S1,S2], kern, targs, eps, struct('corrections', 1));
Afull = surferkernevalmat([S1,S2], kern, targs, eps);
Sover   = [objover{1}{:}];
xinterp = blkdiag(objover{2}{:});

i = sort(randperm(ntarg, 15)).';
j = sort(randperm(S1.npts+S2.npts, 25)).';

Aover = byindex.kernbyindexeval(i, 1:size(xinterp,1), Sover, kern, targs);
Atest = Aover * xinterp(:, j) + Qc(i, j);
err = norm(Atest - Afull(i,j), 'fro') / norm(Afull(i,j), 'fro');
assert(err < tol, 'two-surfer oversampled decomposition: err %.2e > tol %.2e', err, tol);

end


function check_eval(label, i, j, surfers, kern, targs, eps, tol)
% Helper: compare kernbyindexeval against full surferkernevalmat without
% oversampling.
opts = []; opts.ifoversamp = 0;
opts_ns = opts; opts_ns.nonsmoothonly = true;
Qns   = surferkernevalmat(surfers, kern, targs, eps, opts_ns);
Afull = surferkernevalmat(surfers, kern, targs, eps, opts);
Atest = byindex.kernbyindexeval(i, j, surfers, kern, targs, Qns);
err   = norm(Atest - Afull(i,j), 'fro') / norm(Afull(i,j), 'fro');
assert(err < tol, '%s: err %.2e > tol %.2e', label, err, tol);
end
