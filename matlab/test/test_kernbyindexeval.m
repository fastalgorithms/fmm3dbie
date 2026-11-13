% Test byindex.kernbyindexeval against surferkernevalmat.

run ../startup.m

rng(42);

%% Now run the tests

test_two_surfer_helmholtz();
test_mixed_opdims();
test_extra_targ_field();
test_two_surfer_oversampled_decomposition();
test_target_groups();
test_flam_ifmm_oversampled();


function test_two_surfer_helmholtz()
% Two surfers, Helmholtz: kernbyindexeval matches surferkernevalmat.

eps = 1e-8;  tol = 1e-9;

zk = 1.3 + 0.1i;
S1 = geometries.sphere(2, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(2, 2, [6;0;0], 3, 1);
kern = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

ntarg = 40;
targs.r = [3*(2*rand(3,ntarg/2)-1), [6;0;0] + 3*(2*rand(3,ntarg/2)-1)];
i = sort(randperm(ntarg,           15)).';
j = sort(randperm(S1.npts+S2.npts, 25)).';

check_eval('two-surfer Helmholtz', i, j, [S1,S2], kern, targs, eps, tol);

end


function test_mixed_opdims()
% Source surfers with opdims [3 1] and [3 3]: kernbyindexeval matches
% the smooth rule.

tol = 1e-12;

S1 = geometries.sphere(1, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(1, 2, [3;0;0], 3, 1);
kern = [grad(kernel3d('l', 's')), kernel3d('stok', 'd')];

ntarg = 30;
targs.r = [1.5; 0; 0] + 3*(rand(3, ntarg) - 0.5);

w1 = S1.wts(:).';
w2 = repelem(S2.wts(:).', 1, 3);
Aref = [kern(1).eval(S1, targs).*w1, kern(2).eval(S2, targs).*w2];

i = randi(3*ntarg, 40, 1);
j = randi(S1.npts + 3*S2.npts, 40, 1);
err = norm(byindex.kernbyindexeval(i, j, [S1, S2], kern, targs) - Aref(i,j), 'fro') ...
    / norm(Aref(i,j), 'fro');
assert(err < tol, 'mixed opdims: %.2e', err);

end


function test_extra_targ_field()
% Kernel declaring a target field beyond the geometry: the targets carry
% mean_curv as an ntarg x 1 column, which must be sliced by point index.

eps = 1e-8;  tol = 1e-9;

S    = geometries.ellipsoid([1,1,1.2], [2,2,2], [], 4);
kern = kernel3d('bel', 'rlb');

ntarg = 60;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
targs = [];
targs.r = (2 + rand(1,ntarg)) .* [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.n = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.mean_curv = rand(ntarg,1);

i = randi(ntarg, 40, 1);
j = randi(S.npts, 40, 1);

check_eval('extra targ field', i, j, S, kern, targs, eps, tol);

end


function test_two_surfer_oversampled_decomposition()
% Targets close to the surfaces force oversampling: the smooth rule on the
% oversampled surfers times the interpolation, plus corrections.

eps = 1e-8;  tol = 1e-7;

zk = 1.1;
S1 = geometries.sphere(1, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(1, 2, [3;0;0], 3, 1);
kern = [kernel3d('h', 's', zk), kernel3d('h', 'd', zk)];

ntarg = 40;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
dirs = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs.r = [(1 + 0.05*rand(1,ntarg/2)) .* dirs(:,1:ntarg/2), ...
           [3;0;0] + (1 + 0.05*rand(1,ntarg/2)) .* dirs(:,ntarg/2+1:end)];

[Qc, objover] = surferkernevalmat([S1,S2], kern, targs, eps, struct('corrections', 1));
Afull = surferkernevalmat([S1,S2], kern, targs, eps);
Sover   = [objover{1}{:}];
xinterp = blkdiag(objover{2}{:});
assert(sum([Sover.npts]) > S1.npts + S2.npts, ...
    'expected the targets to force oversampling');

i = sort(randperm(ntarg, 15)).';
j = sort(randperm(S1.npts+S2.npts, 25)).';

Aover = byindex.kernbyindexeval(i, 1:size(xinterp,1), Sover, kern, targs);
Atest = Aover * xinterp(:, j) + Qc(i, j);
err = norm(Atest - Afull(i,j), 'fro') / norm(Afull(i,j), 'fro');
assert(err < tol, 'two-surfer oversampled decomposition: err %.2e > tol %.2e', err, tol);

end


function test_target_groups()
% Surfers as target groups with a 2x2 kernel array with mixed opdims:
% matches the on-surface kernbyindex, which uses the same smooth rule.

S1 = geometries.sphere(1, 2, [0;0;0], 3, 1);
S2 = geometries.sphere(1, 2, [3;0;0], 3, 1);
srfrs = [S1, S2];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('stok', 's');
kerns(1,2) = grad(kernel3d('l', 's'));
kerns(2,1) = kernel3d.zeros([1 3]);
kerns(2,2) = kernel3d('l', 'sp');

n = 3*S1.npts + S2.npts;
i = randi(n, 80, 1);
j = randi(n, 80, 1);

Aref = byindex.kernbyindex(i, j, srfrs, kerns);
err = norm(byindex.kernbyindexeval(i, j, srfrs, kerns, srfrs) - Aref, 'fro') ...
    / norm(Aref, 'fro');
assert(err < 1e-14, 'target groups: %.2e', err);

end


function test_flam_ifmm_oversampled()
% The off-surface ifmm usage in byindex.flamfuns. Only run if FLAM is on
% the path.

if exist('ifmm', 'file') ~= 2
    fprintf('FLAM not found, skipping test_flam_ifmm_oversampled\n');
    return
end

eps = 1e-4;  tol = 1e-3;

S    = geometries.sphere(1, 2, [0;0;0], 3, 1);
kern = kernel3d('h', 's', 2.0);

ntarg = 300;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
dirs = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
targs = [];
targs.r = (1.02 + 1.5*rand(1,ntarg)) .* dirs;

funs = byindex.flamfuns(S, kern, eps, targs);
opts_flam = struct('proxybylevel', contains(fileread(which('ifmm')), 'proxybylevel'));
F = ifmm(funs.Afun, funs.rx, funs.cx, 32, eps, funs.pxyfun, opts_flam);

Afull = surferkernevalmat(S, kern, targs, eps);
sig = randn(S.npts, 1);
u   = Afull*sig;
err = norm(ifmm_mv(F, funs.xinterp*sig, funs.Afun) + funs.Q*sig - u)/norm(u);
assert(err < tol, 'ifmm oversampled: err %.2e > tol %.2e', err, tol);

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
