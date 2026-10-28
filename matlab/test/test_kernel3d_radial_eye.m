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
assert(kernel3d.radgreen('summable',{K,I}));
assert(kernel3d.radgreen('summable',{I,K}));
assert(kernel3d.radgreen('summable',{K,J}));
assert(~kernel3d.radgreen('summable',{K,kernel3d.eye(eye(2))}));

% The same identity operation must preserve fitted radial kernels too.
R = kernel3d.radpp(g,[1e-6 4],'sp');
Q = R-0.5*I;
assert(strcmp(Q.name,'radpp'));
assert(isequal(Q.params.ipars,R.params.ipars));
assert(isequal(Q.params.dpars,R.params.dpars));
assert(kernel3d.radpp('summable',{R,I}));
assert(kernel3d.radpp('summable',{I,R}));
assert(kernel3d.radpp('summable',{R,Q}));
check_same(R+I,I+R,src,targ,do_matrix);
check_same(R+(R+I),(R+R)+I,src,targ,do_matrix);
Rc = kernel3d.radpp(g,[1e-6 4],@(t) [zeros(size(t.r,2),1),t.n.'], ...
    'nderiv',1,'targ_fields',{'n'},'ifcomplex',false);
assert(kernel3d.radpp('summable',{Rc,I}));
check_same(Rc+I,I+Rc,src,targ,do_matrix);
check_same(2*(Rc+I),2*Rc+2*I,src,targ,do_matrix);

if ~do_matrix, return; end
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
