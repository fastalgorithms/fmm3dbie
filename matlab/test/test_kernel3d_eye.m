% Test kernel3d.eye and the kernel3d.diag identity term in surfermat and
% surfermatapply.

%% Now run the tests

test_eye_scalar();
test_eye_block();
test_eye_algebra();
test_eye_radial();


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

end


function test_eye_radial()
% Identity shifts must preserve radial quadrature and FMM representations.

g = @lap3d.green;
K = kernel3d(g, 'sp', 'fmm', @lfmm3d);
I = kernel3d.eye;
J = K - 0.5*I;
src = struct('r', [0 0.2; 0.1 -0.3; 0.4 0.5]);
targ = struct('r', [1 1.2; 0.2 0.4; -0.1 -0.2], ...
              'n', [1 0; 0 1; 0 0]);
checks = {J, K, -0.5; I+K, K, 1; I-K, -K, 1; ...
          J+0.5*I, K, 0; J+J, 2*K, -1; J-J, K-K, 0; ...
          2*J, 2*K, -1; -J, -K, 0.5; J/2, K/2, -0.25};
for k = 1:size(checks,1)
    A = checks{k,1}; B = checks{k,2}; d = checks{k,3};
    assert(strcmp(A.name,'radgreen'), 'Identity shift lost radial kernel family');
    assert(norm(A.eval(src,targ)-B.eval(src,targ),'fro') < 1e-12);
    assert(norm(A.params.profun([0.4;0.8])-B.params.profun([0.4;0.8]),'fro') < 1e-12);
    assert(norm(A.params.wfun(targ)-B.params.wfun(targ),'fro') < 1e-12);
    assert(norm(A.diag(targ)-d*ones(2,1)) < 1e-12);
end
assert(strcmp(J.type,'sp'));
assert(isequal(J.params.ipv,K.params.ipv));
assert(isa(J.fmm,'function_handle'));

pairs = {K+I, I+K; ...
         K+(K+I), (K+K)+I; ...
         K+(I+K), (K+I)+K; ...
         K+(I+I), (K+I)+I; ...
         I+I, 2*I; ...
         K+(I-I), K};
for k = 1:size(pairs,1)
    check_same(pairs{k,1},pairs{k,2},src,targ);
end
A = K+J;  assert(strcmp(A.name,'radgreen'));
A = K+kernel3d.zeros([1 1]);  assert(strcmp(A.name,'radgreen'));
caught = false;
try, K + kernel3d.eye(eye(2)); catch, caught = true; end
assert(caught, 'dimension mismatch not caught');

% interleave: identity terms on radial blocks, and identity or zero
% blocks, keep one radial kernel
K2 = kernel3d(g, 's', 'fmm', @lfmm3d);
K3 = 2*K;
Z  = kernel3d.zeros([1 1]);
blocks = {[K, I; K2, K3], [J, K2; K2, K3+I], [K, Z; I, K3]};
for k = 1:numel(blocks)
    B = blocks{k};
    A = kernel3d.interleave(B);
    assert(strcmp(A.name,'radgreen'), 'interleave lost radial kernel family');
    check_interleave(A, B, src, targ);
end

% The same identity operation must preserve fitted radial kernels too.
R = kernel3d.radpp(g,[1e-6 4],'sp');
Q = R-0.5*I;
assert(strcmp(Q.name,'radpp'));
assert(isequal(Q.params.ipars,R.params.ipars));
assert(isequal(Q.params.dpars,R.params.dpars));
A = R+Q;  assert(strcmp(A.name,'radpp'));
check_same(R+I,I+R,src,targ);
check_same(R+(R+I),(R+R)+I,src,targ);
Rc = kernel3d.radpp(g,[1e-6 4],@(t) [zeros(size(t.r,2),1),t.n.'], ...
    'nderiv',1,'targ_fields',{'n'},'ifcomplex',false);
A = Rc+I;  assert(strcmp(A.name,'radpp'));
check_same(Rc+I,I+Rc,src,targ);
check_same(2*(Rc+I),2*Rc+2*I,src,targ);

% interleaved [K, I; K2, K3] against the blocks
S = slicesurfer(geometries.sphere(1,2,[],4,1),1);
eps = 1e-7;
A = surfermat(S,kernel3d.interleave([K, I; K2, K3]),eps);
A11 = surfermat(S,K,eps);  A21 = surfermat(S,K2,eps);  A22 = surfermat(S,K3,eps);
B = zeros(2*S.npts);
B(1:2:end,1:2:end) = A11;  B(1:2:end,2:2:end) = eye(S.npts);
B(2:2:end,1:2:end) = A21;  B(2:2:end,2:2:end) = A22;
assert(norm(A-B,'fro') < 1e-8*norm(B,'fro'));

S = slicesurfer(geometries.sphere(1,2,[],4,1),1);
eps = 1e-7;
A = surfermat(S,K,eps);
B = surfermat(S,J,eps);
assert(norm(B-(A-0.5*eye(S.npts)),'fro') < 1e-12*max(1,norm(A,'fro')));
[C,over] = surfermat(S,J,eps,struct('corrections',true));
sigma = cos(S.r(1,:)).';
u = surfermatapply(S,J,sigma,eps,over,C);
assert(norm(u-B*sigma) < 1e-6*max(1,norm(B*sigma)));
end

function check_same(A,B,src,targ)
assert(strcmp(A.name,B.name));
assert(strcmp(A.type,B.type));
assert(isequal(A.opdims,B.opdims));
assert(norm(A.eval(src,targ)-B.eval(src,targ),'fro') < 1e-12);
DA=zeros(size(targ.r,2)*A.opdims(1),A.opdims(2)); DB=DA;
if ~isempty(A.diag), DA=A.diag(targ); end
if ~isempty(B.diag), DB=B.diag(targ); end
assert(norm(DA-DB,'fro') < 1e-12);
if ~isempty(A.fmm) && ~isempty(B.fmm)
    sigma=ones(size(src.r,2)*A.opdims(2),1);
    assert(norm(A.fmm(1e-10,src,targ,sigma)-B.fmm(1e-10,src,targ,sigma)) < 1e-12);
end
end

function check_interleave(A,B,src,targ)
% A = interleave(B) against the blocks of B
[m,n] = size(B);
ns = size(src.r,2);  nt = size(targ.r,2);
EA = A.eval(src,targ);
DA = zeros(m*nt,n);
if ~isempty(A.diag), DA = A.diag(targ); end
DA = reshape(DA,m,nt,n);
for k = 1:m
    for l = 1:n
        assert(norm(EA(k:m:end,l:n:end)-B(k,l).eval(src,targ),'fro') < 1e-12);
        Dkl = zeros(nt,1);
        if ~isempty(B(k,l).diag), Dkl = B(k,l).diag(targ); end
        assert(norm(reshape(DA(k,:,l),[],1)-Dkl) < 1e-12);
    end
end
sigma = ones(n*ns,1);
pot = A.fmm(1e-10,src,targ,sigma);
ref = zeros(m*nt,1);
for k = 1:m
    for l = 1:n
        ref(k:m:end) = ref(k:m:end) + B(k,l).fmm(1e-10,src,targ,ones(ns,1));
    end
end
assert(norm(pot-ref) < 1e-10*max(1,norm(ref)));
end
