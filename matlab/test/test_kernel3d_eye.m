% Test kernel3d.eye and the kernel3d.diag identity term in surfermat and
% surfermatapply.

%% Now run the tests

test_eye_scalar();
test_eye_block();
test_eye_algebra();
test_eye_fundamental_forms();


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
%
% Also test various possible shapes for diagonal functions

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
Pside  = @(t) reshape(Pfun(t), 3, []);
A5 = surfermat(S, kstok + kernel3d.eye(Pside), eps);
assert(norm(A5 - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: matrix handle (side by side) mismatch');
assert(norm(A3 - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: matrix handle (pages) mismatch');
assert(norm(A4 - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: matrix handle (stacked) mismatch');

% non-square blocks: the normal as a 3x1 block (side by side == S.n, and
% stacked) and as a 1x3 block (stacked == S.n.', and side by side)
Nref = zeros(3*n, n);
for k = 1:n
    Nref(3*(k-1) + (1:3), k) = nrm(:,k);
end
A = surfermat(S, kernel3d.eye(@(t) t.n(:,:)), eps);
assert(norm(A - Nref, 'fro') < tol*norm(Nref, 'fro'), 'eye: 3x1 side by side mismatch');
A = surfermat(S, kernel3d.eye(@(t) reshape(t.n(:,:), [], 1)), eps);
assert(norm(A - Nref, 'fro') < tol*norm(Nref, 'fro'), 'eye: 3x1 stacked mismatch');
A = surfermat(S, kernel3d.eye(@(t) reshape(t.n(:,:), 3, 1, [])), eps);
assert(norm(A - Nref, 'fro') < tol*norm(Nref, 'fro'), 'eye: 3x1 pages mismatch');
A = surfermat(S, kernel3d.eye(@(t) t.n(:,:).'), eps);
assert(norm(A - Nref.', 'fro') < tol*norm(Nref, 'fro'), 'eye: 1x3 stacked mismatch');
A = surfermat(S, kernel3d.eye(@(t) reshape(t.n(:,:), 1, [])), eps);
assert(norm(A - Nref.', 'fro') < tol*norm(Nref, 'fro'), 'eye: 1x3 side by side mismatch');
A = surfermat(S, kernel3d.eye(@(t) reshape(t.n(:,:), 1, 3, [])), eps);
assert(norm(A - Nref.', 'fro') < tol*norm(Nref, 'fro'), 'eye: 1x3 pages mismatch');

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
hfun = @(t) 1 + t.r(1,:).^2.';
A = surfermat(S, hfun*(kern + kernel3d.eye(@(t) t.r(3,:))), eps);
h = 1 + S.r(1,:).^2;
Aref = h(:).*(Km + diag(S.r(3,:)));
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: variable coefficient mismatch');

% scalar left multiplier: row, column and pages all give the same matrix
Aref = h(:).*(Km + diag(S.r(3,:)));
kd3  = kern + kernel3d.eye(@(t) t.r(3,:));
A = surfermat(S, (@(t) 1 + t.r(1,:).^2)*kd3, eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar row multiplier mismatch');
A = surfermat(S, (@(t) (1 + t.r(1,:).^2).')*kd3, eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar column multiplier mismatch');
A = surfermat(S, (@(t) reshape(1 + t.r(1,:).^2, 1, 1, []))*kd3, eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar pages multiplier mismatch');

% scalar diag handle: row, column and pages all give the same matrix
Aref = Km + diag(S.r(3,:));
A = surfermat(S, kern + kernel3d.eye(@(t) t.r(3,:)), eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar row mismatch');
A = surfermat(S, kern + kernel3d.eye(@(t) t.r(3,:).'), eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar column mismatch');
A = surfermat(S, kern + kernel3d.eye(@(t) reshape(t.r(3,:), 1, 1, [])), eps);
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: scalar pages mismatch');

% diag handle using other fields
A = surfermat(S, kern + kernel3d.eye(@(t) 1 + t.mean_curv(:)), eps);
Aref = Km + I + diag(S.mean_curv(:));
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: mean curvature diag mismatch');

end


function test_eye_fundamental_forms()
% diag handle built from the surfer's fundamental forms: Gauss curvature

tol = 1e-12;

S   = geometries.ellipsoid([1,1,1.1],[3,3,3],[],4);
eps = 1e-10;

kern = kernel3d('l','s');
A = surfermat(S, kern + kernel3d.eye(@gauss_curv), eps);
Aref = surfermat(S, kern, eps) + diag(gauss_curv(S));
assert(norm(A - Aref, 'fro') < tol*norm(Aref, 'fro'), 'eye: gauss curvature diag mismatch');

% the same identities carried through every kernel operation. The diag
% handles read surfer-only data (fundamental forms, mean curvature).
Km = surfermat(S, kern, eps);
g  = gauss_curv(S);  H = S.mean_curv(:);  h = 1 + S.r(1,:).'.^2;
Dg = diag(g);        DH = diag(H);
kg = kern + kernel3d.eye(@gauss_curv);
kh = kernel3d.eye(@(t) t.mean_curv(:));
chk = @(K, Aref, msg) assert(norm(surfermat(S, K, eps) - Aref, 'fro') ...
    < tol*norm(Aref, 'fro'), ['eye forms: ' msg ' mismatch']);

chk(kg + kh, Km + Dg + DH, 'plus');
chk(kg - kh, Km + Dg - DH, 'minus');
chk(kh - kg, DH - Km - Dg, 'minus (reversed)');
chk(-kg,     -(Km + Dg),   'uminus');
chk(kg/3,    (Km + Dg)/3,  'mrdivide');
chk(2*kg,    2*(Km + Dg),  'scalar mtimes');
chk(kg.*2,   2*(Km + Dg),  'scalar times');

% elementwise and constant matrix products, changing opdims
chk([1;2].*kg, kron(Km + Dg, [1;2]), 'column times');
chk(kg.*[1 2], kron(Km + Dg, [1 2]), 'row times');
chk([1;2]*kg,  kron(Km + Dg, [1;2]), 'matrix left mtimes');
chk(kg*[1 2],  kron(Km + Dg, [1 2]), 'matrix right mtimes');

% function handle multipliers, left (targets) and right (sources)
hfun = @(t) 1 + t.r(1,:).^2;
chk(hfun*kg, h.*(Km + Dg),   'handle left mtimes');
% (identity only on the right: with oversampling, K*h and Km*diag(h)
% differ at the level of the discretization error)
chk(kernel3d.eye(@gauss_curv)*hfun, diag(g.*h), 'handle right mtimes');

% multiplier that itself reads a non-geometry field: it must be declared
% in targ_fields so the smooth evaluation receives it
kmc = (@(t) t.mean_curv)*kg;
kmc.targ_fields = {'mean_curv'};
chk(kmc, H.*(Km + Dg), 'mean curvature multiplier');

% interleaved blocks
kerns(2,2) = kernel3d();
kerns(1,1) = kg;  kerns(1,2) = kernel3d.zeros();
kerns(2,1) = kh;  kerns(2,2) = kg;
chk(kernel3d(kerns), kron(Km + Dg, eye(2)) + kron(DH, [0 0; 1 0]), 'interleave');

% sanity check of the forms themselves: K = 1/R^2 on a sphere
R  = 2;
Sp = geometries.ellipsoid([R,R,R],[3,3,3],[],12);
assert(norm(gauss_curv(Sp) - 1/R^2, inf) < 1e-5, 'eye: gauss curvature of sphere');

end

function K = gauss_curv(t)
I  = cat(3, t.ffform{:});
II = cat(3, t.sfform{:});
detI  = I(1,1,:).*I(2,2,:)   - I(1,2,:).^2;
detII = II(1,1,:).*II(2,2,:) - II(1,2,:).^2;
K = detII(:) ./ detI(:);
end
