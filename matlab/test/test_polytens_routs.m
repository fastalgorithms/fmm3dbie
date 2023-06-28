addpath '../';

%% Legendre tests
f = @(u,v) cos(3.2*u + 4.1*v);
norder = 14;
rnodes = polytens.lege_nodes(norder);
rwts = polytens.lege_weights(norder);
plot(rnodes(1,:),rnodes(2,:),'kx');


umat = polytens.lege_vals2coefs(norder,rnodes);
vmat = polytens.lege_coefs2vals(norder,rnodes);

Q = integral2(f,-1,1,-1,1);
fvals = f(rnodes(1,:),rnodes(2,:));
Q2 = fvals*rwts';
rerr = norm(Q2-Q)/norm(Q);
fprintf('Error in integral (legendre)=%d\n',rerr);

fvals = fvals.';
fcoefs = umat*fvals;
semilogy(abs(fcoefs),'k.');

uvtest = zeros(2,1);
uvtest(1) = 0.3;
uvtest(2) = 0.4;
pp = polytens.lege_pols(norder,uvtest);

fex = f(uvtest(1),uvtest(2));
ftest = fcoefs.'*pp;
rerr = norm(ftest-fex)/norm(fex);
fprintf('Error in interpolant (legendre)=%d\n',rerr);


%% Chebyshev tests
f = @(u,v) cos(3.2*u + 4.1*v);
norder = 14;
rnodes = polytens.cheb_nodes(norder);
rwts = polytens.cheb_weights(norder);
plot(rnodes(1,:),rnodes(2,:),'kx');


umat = polytens.cheb_vals2coefs(norder,rnodes);
vmat = polytens.cheb_coefs2vals(norder,rnodes);

Q = integral2(f,-1,1,-1,1);
fvals = f(rnodes(1,:),rnodes(2,:));
Q2 = fvals*rwts';
rerr = norm(Q2-Q)/norm(Q);
fprintf('Error in integral (Cheb)=%d\n',rerr);

fvals = fvals.';
fcoefs = umat*fvals;
semilogy(abs(fcoefs),'k.');

uvtest = zeros(2,1);
uvtest(1) = 0.3;
uvtest(2) = 0.4;
pp = polytens.cheb_pols(norder,uvtest);

fex = f(uvtest(1),uvtest(2));
ftest = fcoefs.'*pp;
rerr = norm(ftest-fex)/norm(fex);
fprintf('Error in interpolant (cheb)=%d\n',rerr);





