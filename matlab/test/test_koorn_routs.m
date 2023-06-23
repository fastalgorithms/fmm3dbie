addpath '../';
f = @(u,v) sin(3.2*u + 4.1*v);
norder = 14;
rnodes = koorn.rv_nodes(norder);
rwts = koorn.rv_weights(norder);
plot(rnodes(1,:),rnodes(2,:),'kx');
umat = koorn.vals2coefs(norder,rnodes);
vmat = koorn.coefs2vals(norder,rnodes);

vmax = @(u) 1-u;
Q = integral2(f,0,1,0,vmax);
fvals = f(rnodes(1,:),rnodes(2,:));
Q2 = fvals*rwts';
rerr = norm(Q2-Q)/norm(Q);
fprintf('Error in integral=%d\n',rerr);

fvals = fvals.';
fcoefs = umat*fvals;
semilogy(abs(fcoefs),'k.');

uvtest = zeros(2,1);
uvtest(1) = 0.3;
uvtest(2) = 0.4;
pp = koorn.pols(norder,uvtest);

fex = f(uvtest(1),uvtest(2));
ftest = fcoefs.'*pp;
rerr = norm(ftest-fex)/norm(fex);
fprintf('Error in interpolant=%d\n',rerr);







