% Test kernel3d.eye and the kernel3d.diag identity term in surfermat and
% surfermatapply.

%% Now run the tests

test_eye_scalar();
test_eye_block();
test_eye_algebra();


function test_eye_scalar()
% surfermat(kern + kernel3d.eye(a)) == a*I + surfermat(kern), dense and
% corrections + surfermatapply.

tol = 1e-12;

S   = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
eps = 1e-10;
a   = -0.5;

kern  = kernel3d('h','d',1.1);
keye  = kern + kernel3d.eye(a);

A1 = surfermat(S, kern, eps) + a*eye(S.npts);
A2 = surfermat(S, keye, eps);
assert(norm(A1 - A2, 'fro') < tol*norm(A1, 'fro'), 'eye: dense surfermat mismatch');

[cors, objover] = surfermat(S, keye, eps, struct('corrections', 1));
sigma = randn(S.npts, 1);
pot   = surfermatapply(S, keye, sigma, eps, objover, cors);
assert(norm(pot - A2*sigma) < 1e-6*norm(A2*sigma), 'eye: surfermatapply mismatch');

end


function test_eye_block()
% Block identity terms: a 3x3 Stokes jump and an interleaved matrix of
% kernels with an off-diagonal identity block.

tol = 1e-12;

S   = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
eps = 1e-10;
n   = S.npts;

kstok = kernel3d('stok','d');
A1 = surfermat(S, kstok, eps) + 0.5*eye(3*n);
A2 = surfermat(S, kstok + kernel3d.eye(0.5*eye(3)), eps);
assert(norm(A1 - A2, 'fro') < tol*norm(A1, 'fro'), 'eye: stokes block mismatch');

a = 0.5; c = 0.3;
kd = kernel3d('l','d');
kz = kernel3d.zeros();

kerns(2,2) = kernel3d();
kerns(1,1) = kd;  kerns(1,2) = kz;  kerns(2,1) = kz;  kerns(2,2) = kd;
Asm = surfermat(S, kernel3d(kerns), eps);

kerns(1,1) = kd + kernel3d.eye(a);
kerns(2,1) = kernel3d.eye(c);
kerns(2,2) = kd + kernel3d.eye(a);
Ad = surfermat(S, kernel3d(kerns), eps);

Aref = Asm + kron(eye(n), [a 0; c a]);
assert(norm(Ad - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: interleaved block mismatch');

end


function test_eye_algebra()
% diag is carried through scalar/matrix multiplication, negation,
% subtraction, division, and function-handle eye.

tol = 1e-12;

S   = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
eps = 1e-10;
n   = S.npts;

kern = kernel3d('l','s');
Km   = surfermat(S, kern, eps);
I    = eye(n);

A = surfermat(S, 2*(kern + kernel3d.eye(0.5)), eps);
assert(norm(A - (2*Km + I), 'fro') < tol*norm(A, 'fro'), 'eye: scalar times mismatch');

A = surfermat(S, -(kern - kernel3d.eye(2)), eps);
assert(norm(A - (-Km + 2*I), 'fro') < tol*norm(A, 'fro'), 'eye: minus/uminus mismatch');

A = surfermat(S, (kern + kernel3d.eye(3))/3, eps);
assert(norm(A - (Km/3 + I), 'fro') < tol*norm(A, 'fro'), 'eye: mrdivide mismatch');

% variable coefficient, left multiplied by a function of the target
hfun = @(t) reshape(1 + t.r(1,:).^2, 1, 1, []);
A = surfermat(S, hfun*(kern + kernel3d.eye(@(t) reshape(t.r(3,:), 1, 1, []))), eps);
h = 1 + S.r(1,:).^2;
Aref = h(:).*(Km + diag(S.r(3,:)));
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: variable coefficient mismatch');

end
