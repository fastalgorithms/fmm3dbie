% Test byindex.proxy_points and byindex.proxyfun.
%
% The proxy functions are called with the arguments FLAM would pass, and
% the skeletons they give are checked against the interactions they
% stand in for. FLAM is not needed.

run ../startup.m

rng(42);

%% Now run the tests

test_proxy_points_calling_sequences();
test_proxy_points_rejects();
test_rskelf_single_surfer();
test_rskelf_addtrans();
test_ifmm_kernel_array();
test_ifmm_stokes();
test_ifmm_offsurface();
test_proxybylevel();


function test_proxy_points_calling_sequences()
% All calling sequences, and the shape of the output.

eps  = 1e-8;
kern = kernel3d('l', 's');

p1 = byindex.proxy_points(kern, 1, eps);
p2 = byindex.proxy_points(kern, 1, eps, 1);
p3 = byindex.proxy_points(kern, 1, eps, struct('nshell_max', 2, 'ifaddtrans', true));
p4 = byindex.proxy_points(kern, 1, eps, [], struct('nproxy_max', 4000));

assert(p1.nshell == 1 && p2.nshell == 1 && p2.nshell_max == 1, ...
    'Laplace should need one shell');
assert(p3.nshell_max == 2 && p3.ifaddtrans, 'opts not used as 4th argument');
assert(p4.nshell_max == 20, 'default nshell_max should be 20');
npxy = p1.nproxy*p1.nshell;
assert(isequal(size(p1.r), [3 npxy]) && isequal(size(p1.n), [3 npxy]) && ...
    numel(p1.w) == npxy, 'wrong proxy sizes');
assert(abs(p1.radius - 1.5) < 1e-12, 'default radius should be 1.5 box sizes');
assert(abs(sum(p1.w) - 4*pi*1.5^2) < 1e-10, 'weights should sum to the sphere area');

% Helmholtz also needs one shell
ph = byindex.proxy_points(kernel3d('h', 's', 1.2), 1, eps);
assert(ph.nshell == 1, 'Helmholtz should need one shell');

% Tighter tolerances need more proxy points
p5 = byindex.proxy_points(kern, 1, 1e-4);
assert(p5.nproxy < p1.nproxy, 'nproxy should grow as eps decreases');

end


function test_proxy_points_rejects()
% Kernels needing fields other than n at the proxies, and proxy surfaces
% that would contain non-neighbor boxes, are rejected.

caught = false;
try
    byindex.proxy_points(kernel3d('bel', 'rlb'), 1, 1e-6);
catch
    caught = true;
end
assert(caught, 'mean_curv at proxy targets should raise an error');

caught = false;
try
    byindex.proxy_points(kernel3d('l', 's'), 1, 1e-6, struct('crad', 2));
catch
    caught = true;
end
assert(caught, 'crad > sqrt(3) should raise an error');

end


function test_rskelf_single_surfer()
% rskelf calling sequence: proxy-based skeleton gives a low-rank
% approximation of the src->targ block, and nbr keeps only the points
% inside the proxy surface.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('l', 's');
ctr  = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kern, l, eps);

dist    = vecnorm(S.r - ctr);
srcids  = find(dist < 0.5);
targids = find(dist > 2.0);
nbrin   = find(dist > 0.5 & dist < 1.4);

src2targ = byindex.kernbyindex(targids, srcids, S, kern);

[Kpxy, nbr] = byindex.proxyfun(S, kern, [], proxy, S.r, srcids, ...
    [nbrin, targids], l, ctr);
assert(isequal(sort(nbr), sort(nbrin)), 'nbr filtered incorrectly');
assert(size(Kpxy, 1) == size(proxy.r, 2), 'wrong number of proxy rows');
check_lr('rskelf single surfer', src2targ, Kpxy, tol);

end


function test_rskelf_addtrans()
% Non-symmetric kernel with proxy.ifaddtrans: the skeleton compresses both
% the columns and the rows of slf.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 10, [0;0;0], 6, 1);
kern = kernel3d('l', 'd');
ctr  = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kern, l, eps, struct('ifaddtrans', true));

dist    = vecnorm(S.r - ctr);
srcids  = find(dist < 0.5);
targids = find(dist > 2.0);

src2targ = byindex.kernbyindex(targids, srcids, S, kern);
targ2src = byindex.kernbyindex(srcids, targids, S, kern);

Kpxy = byindex.proxyfun(S, kern, [], proxy, S.r, srcids, targids, l, ctr);
npxy = size(proxy.r, 2);
assert(size(Kpxy, 1) == 2*npxy, 'ifaddtrans should give %d rows, got %d', ...
    2*npxy, size(Kpxy, 1));

check_lr('rskelf addtrans, columns', src2targ,    Kpxy, tol);
check_lr('rskelf addtrans, rows',    targ2src.', Kpxy, tol);

end


function test_ifmm_kernel_array()
% ifmm calling sequence on two surfers with a 2x2 Laplace kernel array:
% rc = 'c' compresses sources on surfer 1, rc = 'r' compresses targets on
% surfer 1.

eps = 1e-10;  tol = 1e-9;

S1 = geometries.sphere(2, 8, [0;0;0], 5, 1);
S2 = geometries.sphere(2, 8, [8;0;0], 5, 1);
srfrs = [S1, S2];

kerns(2,2) = kernel3d();
kerns(1,1) = kernel3d('l', 's');
kerns(2,1) = kernel3d('l', 'sp');
kerns(1,2) = kernel3d('l', 'd');
kerns(2,2) = kernel3d('l', 'dp');

ctr = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kerns, l, eps);
x = [S1.r, S2.r];

near = find(vecnorm(S1.r - ctr) < 0.5);
far  = [find(vecnorm(S1.r - ctr) > 2.0), S1.npts + (1:S2.npts)];

% columns
src2targ = byindex.kernbyindex(far, near, srfrs, kerns);
[Kpxy, nbr] = byindex.proxyfun(srfrs, kerns, [], proxy, 'c', x, x, near, far, l, ctr);
assert(isempty(nbr), 'far points should be removed from nbr');
assert(size(Kpxy, 1) == 2*size(proxy.r, 2), 'one proxy copy per block row');
check_lr('ifmm kernel array, rc = c', src2targ, Kpxy, tol);

% rows
targ2src = byindex.kernbyindex(near, far, srfrs, kerns);
Kpxy = byindex.proxyfun(srfrs, kerns, [], proxy, 'r', x, x, near, far, l, ctr);
assert(size(Kpxy, 2) == 2*size(proxy.r, 2), 'one proxy copy per block column');
check_lr('ifmm kernel array, rc = r', targ2src.', Kpxy.', tol);

end


function test_ifmm_stokes()
% opdims = [3 3]: Stokes single layer, both rc.

eps = 1e-8;  tol = 1e-7;

S    = geometries.sphere(2, 8, [0;0;0], 5, 1);
kern = kernel3d('stok', 's');
ctr  = [2; 0; 0];  l = 1;
proxy = byindex.proxy_points(kern, l, eps);
x = repelem(S.r, 1, 3);

near = find(vecnorm(S.r - ctr) < 0.5);
far  = find(vecnorm(S.r - ctr) > 2.0);
near = reshape(3*near + (-2:0).', 1, []);
far  = reshape(3*far  + (-2:0).', 1, []);

src2targ = byindex.kernbyindex(far, near, S, kern);
Kpxy = byindex.proxyfun(S, kern, [], proxy, 'c', x, x, near, far, l, ctr);
assert(size(Kpxy, 1) == 3*size(proxy.r, 2), 'wrong number of proxy rows');
check_lr('ifmm stokes, rc = c', src2targ, Kpxy, tol);

targ2src = byindex.kernbyindex(near, far, S, kern);
Kpxy = byindex.proxyfun(S, kern, [], proxy, 'r', x, x, near, far, l, ctr);
check_lr('ifmm stokes, rc = r', targ2src.', Kpxy.', tol);

end


function test_ifmm_offsurface()
% Off-surface targets: rc = 'c' compresses sources near a box, rc = 'r'
% compresses targets in a box.

eps = 1e-10;  tol = 1e-9;

S    = geometries.sphere(2, 8, [0;0;0], 6, 1);
kern = kernel3d('h', 's', 1.2);
l = 1;
proxy = byindex.proxy_points(kern, l, eps);

rng(6);
ntarg = 200;
th = acos(2*rand(1,ntarg)-1);  ph = 2*pi*rand(1,ntarg);
dirs = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];

% columns: sources on S near ctr, targets far from it
ctr   = [2; 0; 0];
targs = [];
targs.r = ctr + (2 + rand(1,ntarg)) .* dirs;
srcids  = find(vecnorm(S.r - ctr) < 0.5);
src2targ = byindex.kernbyindexeval((1:ntarg).', srcids, S, kern, targs);
[Kpxy, nbr] = byindex.proxyfun(S, kern, targs, proxy, 'c', targs.r, S.r, ...
    srcids, 1:ntarg, l, ctr);
assert(isempty(nbr), 'far targets should be removed from nbr');
check_lr('offsurface, rc = c', src2targ, Kpxy, tol);

% rows: targets in a box at the origin, all of S at distance 2
ctr   = [0; 0; 0];
targs = [];
targs.r = 0.5*rand(1,ntarg) .* dirs;
targ2src = byindex.kernbyindexeval((1:ntarg).', (1:S.npts).', S, kern, targs);
Kpxy = byindex.proxyfun(S, kern, targs, proxy, 'r', targs.r, S.r, ...
    1:ntarg, 1:S.npts, l, ctr);
assert(isequal(size(Kpxy), [ntarg, size(proxy.r, 2)]), 'wrong rc = r size');
check_lr('offsurface, rc = r', targ2src.', Kpxy.', tol);

end


function test_proxybylevel()
% With one argument, proxyfun returns the proxy function for that level,
% rebuilt for a box of size width/2^(lvl-1).

eps = 1e-6;  tol = 1e-5;

S    = geometries.sphere(2, 8, [0;0;0], 6, 1);
kern = kernel3d('h', 's', 4);
width = 4;
proxy = byindex.proxy_points(kern, width, eps);
pxyfun = @(varargin) byindex.proxyfun(S, kern, [], proxy, varargin{:});

lvl = 3;
pxyfunlvl = pxyfun(lvl);
assert(isa(pxyfunlvl, 'function_handle'), 'expected a function handle for lvl');

plvl = byindex.proxy_points(kern, width/2^(lvl-1), eps);
assert(plvl.nproxy < proxy.nproxy, 'smaller boxes should need fewer proxies');

ctr = [2; 0; 0];  l = width/2^(lvl-1);
dist    = vecnorm(S.r - ctr);
srcids  = find(dist < l/2);
targids = find(dist > 1.6*l);

src2targ = byindex.kernbyindex(targids, srcids, S, kern);
Kpxy = pxyfunlvl(S.r, srcids, targids, l, ctr);
assert(size(Kpxy, 1) == size(plvl.r, 2), 'level proxy function should use plvl');
check_lr('proxybylevel', src2targ, Kpxy, tol);

end


function check_lr(label, src2targ, Kpxy, tol)
% Helper: verify the interpolative decomposition of Kpxy (column-pivoted
% QR, as in FLAM) also compresses src2targ.
[~, R, jpvt] = qr(Kpxy, 0);
d  = abs(diag(R));
k  = max(sum(d > tol/10 * d(1)), 1);
T  = R(1:k, 1:k) \ R(1:k, k+1:end);
sk = jpvt(1:k);
rd = jpvt(k+1:end);
err = norm(src2targ(:, rd) - src2targ(:, sk)*T, 'fro') / norm(src2targ, 'fro');
assert(err < tol, '%s: err %.2e > tol %.2e', label, err, tol);
end
