function test_kernel3d_radial_eye(do_matrix)
% Identity shifts must preserve radial quadrature and FMM representations.
% Pass false to run only the algebra checks, without compiled routines.
if nargin < 1, do_matrix = true; end

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
    check_same(pairs{k,1},pairs{k,2},src,targ,do_matrix);
end
assert(strcmp((K+J).name,'radgreen'));
assert(strcmp((K+kernel3d.zeros([1 1])).name,'radgreen'));
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
    check_interleave(A, B, src, targ, do_matrix);
end

% The same identity operation must preserve fitted radial kernels too.
R = kernel3d.radpp(g,[1e-6 4],'sp');
Q = R-0.5*I;
assert(strcmp(Q.name,'radpp'));
assert(isequal(Q.params.ipars,R.params.ipars));
assert(isequal(Q.params.dpars,R.params.dpars));
assert(strcmp((R+Q).name,'radpp'));
check_same(R+I,I+R,src,targ,do_matrix);
check_same(R+(R+I),(R+R)+I,src,targ,do_matrix);
Rc = kernel3d.radpp(g,[1e-6 4],@(t) [zeros(size(t.r,2),1),t.n.'], ...
    'nderiv',1,'targ_fields',{'n'},'ifcomplex',false);
assert(strcmp((Rc+I).name,'radpp'));
check_same(Rc+I,I+Rc,src,targ,do_matrix);
check_same(2*(Rc+I),2*Rc+2*I,src,targ,do_matrix);

if ~do_matrix, return; end
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

function check_same(A,B,src,targ,do_fmm)
assert(strcmp(A.name,B.name));
assert(strcmp(A.type,B.type));
assert(isequal(A.opdims,B.opdims));
if do_fmm || ~strcmp(A.name,'radpp')
    assert(norm(A.eval(src,targ)-B.eval(src,targ),'fro') < 1e-12);
else
    assert(norm(A.params.coefs(:)-B.params.coefs(:)) < ...
           1e-12*max(1,norm(B.params.coefs(:))));
end
DA=zeros(size(targ.r,2)*A.opdims(1),A.opdims(2)); DB=DA;
if ~isempty(A.diag), DA=A.diag(targ); end
if ~isempty(B.diag), DB=B.diag(targ); end
assert(norm(DA-DB,'fro') < 1e-12);
if do_fmm && ~isempty(A.fmm) && ~isempty(B.fmm)
    sigma=ones(size(src.r,2)*A.opdims(2),1);
    assert(norm(A.fmm(1e-10,src,targ,sigma)-B.fmm(1e-10,src,targ,sigma)) < 1e-12);
end
end

function check_interleave(A,B,src,targ,do_fmm)
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
if ~do_fmm, return; end
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
