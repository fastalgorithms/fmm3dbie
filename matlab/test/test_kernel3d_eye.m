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
% Block identity terms: a 3x3 Stokes jump, a 3x3 matrix valued function
% handle (pages and stacked), and a matrix of kernels with an off-diagonal
% identity block.

tol = 1e-12;

S   = slicesurfer(geometries.ellipsoid([1,1,1.1],[3,3,3],[],4), 1);
eps = 1e-10;
n   = S.npts;

kstok = kernel3d('stok','d');
A1 = surfermat(S, kstok, eps) + 0.5*eye(3*n);
A2 = surfermat(S, kstok + kernel3d.eye(0.5*eye(3)), eps);
assert(norm(A1 - A2, 'fro') < tol*norm(A1, 'fro'), 'eye: stokes block mismatch');

% matrix valued function handle: tangential projector P(x) = I - n n^T,
% returned as [3 x 3 x nt] pages, and as (3*nt x 3) stacked by target
Pfun   = @(t) eye(3) - pagemtimes(reshape(t.n(:,:),3,1,[]), 'none', ...
                                  reshape(t.n(:,:),3,1,[]), 'transpose');
Pstack = @(t) reshape(permute(Pfun(t), [1 3 2]), [], 3);

nrm  = S.n(:,:);
Pref = zeros(3*n);
for k = 1:n
    ii = 3*(k-1) + (1:3);
    Pref(ii,ii) = eye(3) - nrm(:,k)*nrm(:,k).';
end
Aref = surfermat(S, kstok, eps) + Pref;

A3 = surfermat(S, kstok + kernel3d.eye(Pfun), eps);
A4 = surfermat(S, kstok + kernel3d.eye(Pstack), eps);
assert(norm(A3 - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: matrix handle (pages) mismatch');
assert(norm(A4 - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: matrix handle (stacked) mismatch');

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

% diag handle using other fields
A = surfermat(S, kern + kernel3d.eye(@(t) 1 + t.mean_curv(:)), eps);
Aref = Km + I + diag(S.mean_curv(:));
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: mean curvature diag mismatch');

end
