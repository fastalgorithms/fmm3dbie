%%%%% a basic test code for constructing operators for transmission 
%%%%% problems. Not supposed to be used more generally. Compares with 
%%%%% analytic solution (exterior only...)

eps1 = 1.4;
eps2 = 1.0;
epss = [eps1,eps2];

mu1 = 1;
mu2 = 1;
mus = [mu1,mu2];

norder = 5;
occ = 200;
zk = 4.0;
rank_or_tol = 1E-4;
m = 10;
m2= 12;


a = 1.0;
b = 1;
c = 1.0;
ifc = 0;
rmax = min(0.5,2*pi/(abs(zk)*max(sqrt(epss.*mus))));

% For triangles
S = surfer.ellipsoid([1,1,1],rmax,norder);


% Pick 'm' points inside object
xyz_in = zeros(3, m);

rng(42);

xyz_in = randn(3,m);
xyz_in = (1-rmax)*xyz_in.*(rand(1,m)./vecnorm(xyz_in,2));
xyz_in = diag([a,b,c])*xyz_in;

% Pick 'm' points outside surface to the object
xyz_out = randn(3, m2);
rads = (0.5*rand(1,m2)+1.5);
xyz_out = xyz_out.*(rads./vecnorm(xyz_out,2));
xyz_out = diag([a,b,c])*xyz_out;

q = rand(m, 1)-0.5 + 1j*(rand(m,1)-0.5);


% Quadrature points
x = S.r;
x = repelem(x, 1,2);
area = S.wts';
N = S.npts;

% Complex parameterizations
zpars = complex([zk;eps1;mu1;eps2;mu2]);


eps = rank_or_tol;

% Compute the quadrature corrections
tic, qcorr = helm_trans_near_corr(S,zpars,eps); tquad =toc;
P = zeros(N,1);
w = whos('qcorr');
fprintf('quad: %10.4e (s) / %6.2f (MB)\n',tquad,w.bytes/1e6)

% Set system matrix
opts = [];
opts.l2scale = true;
Afun_use = @(i, j) helm3d.transmission.matgen(i, j, S, zpars, P, qcorr,opts);

zk0 = zpars(1)*sqrt(zpars(2)*zpars(3));
zk1 = zpars(2)*sqrt(zpars(4)*zpars(5));
zks = [zk0,zk1];
% Set proxy function
pxyfun_use = @(xx, slf, nbr, proxy, l, ctr) ...
    helm3d.transmission.proxyfun(xx, slf, nbr, proxy, l, ctr,zks,S);

% Factorize the matrix
opts = struct('verb', 1, 'symm','n', 'zk', zk);

tic, F = srskelf_asym_new(Afun_use, x, occ, rank_or_tol, pxyfun_use, opts); tfac = toc;
w = whos('F');
fprintf([repmat('-',1,80) '\n'])
fprintf('mem: %6.4f (GB)\n',w.bytes/1048576/1024)


zk0 = zpars(1)*sqrt(zpars(4)*zpars(5));
deps0 = zpars(4);
zk  = zpars(1)*sqrt(zpars(2)*zpars(3));
deps  = zpars(2);

% Generate the data due to a collection of point charges
targ = [];
targ.r = xyz_in;
u0  = helm3d.kern(zk,targ,S,'s');
du0 = helm3d.kern(zk,targ,S,'sprime');

B1 =  (u0*q).*sqrt(area).';
B2 = -(du0*q/deps).*sqrt(area).';
Btmp = zeros(size(B1').*[1,2]);
Btmp(1:2:end) = B1;
Btmp(2:2:end) = B2;
B = Btmp.';

% Solve for surface density
tic, X = srskelf_sv_nn(F, B); tsolve = toc;

% A = Afun_use(1:2*N, 1:2*N);
% X_test = A \ B;
X_test = X;

 dens1 = X_test(1:2:end)./sqrt(area.');
 dens2 = X_test(2:2:end)./sqrt(area.');
 
 targ.r = xyz_out;
 G  = helm3d.kern(zk,S,targ,'s').';
 dG = helm3d.kern(zk,S,targ,'d').'; 
 u0comp  =  sum(G.*dens2.*(area.'),1);
 u0dcomp =  sum(dG.*dens1.*(area.'),1);
 u0 = deps^2*u0comp + deps*u0dcomp;
 
 nvecs = zeros(size(xyz_in));
 u0true = helm3d.green(zk,xyz_in,xyz_out);
 u0true = (u0true)*q;
 errs = u0.' - u0true
 norm(errs)/norm(u0true)
