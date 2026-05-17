%
%  flex2d_free_bie_new_bcs
%
%  Solves the flexural plate equation with variable-coefficient free
%  boundary conditions on a starfish domain, using a boundary integral
%  equation formulation (boundary-only, no volume term).
%  The right-hand side and reference solution are generated from
%  point sources placed outside the domain.
%

%% Geometry and problem definition

iseed = 8675309;
rng(iseed);

a = 0.093517375528373;
b = 0;
c = 1;
nu = 0.3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));
zk = [zk1 zk2];

narms = 3;
amp = 0.25;
cparams = []; cparams.nover = 1; cparams.maxchunklen = 4.0/max(abs(zk));
pref = []; pref.k = 16;
start = tic;
chnkr = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref);
opts = []; opts.nover = 2;
chnkr = refine(chnkr, opts);
fprintf('%5.2e s : time to build geometry\n', toc(start))

% Variable coefficient data on boundary
a_var  = 1 + chnkr.r(1,:).^2;
dadn   = chnkr.r(1,:);
dadt   = chnkr.r(2,:);
chnkr  = makedatarows(chnkr, 3);
chnkr.data(1,:) = a_var;
chnkr.data(2,:) = dadn;
chnkr.data(3,:) = dadt;

% Point sources outside domain and interior targets
ns = 3;
ts = 2*pi*rand(ns, 1);
sources = 3.0*starfish(ts, narms, amp);
strengths = randn(ns, 1);

nt = 10;
ts = 2*pi*rand(nt, 1);
targets = starfish(ts, narms, amp) .* repmat(rand(1,nt), 2, 1);

%% Quadrature

kern1 = @(s,t) chnk.flex2d.kern(zk, s, t, 's');
kern2 = @(s,t) flex2d.kern(zk, s, t, 'free_plate_bcs_var', nu);
fkern1 = @(s,t) flex2d.kern(zk, s, t, 'free_plate_var', nu);
double  = @(s,t) chnk.lap2d.kern(s, t, 'd');
hilbert = @(s,t) chnk.lap2d.kern(s, t, 'hilb');

opts  = []; opts.sing  = 'log';
opts2 = []; opts2.sing = 'pv';

start = tic;
sysmat1 = chunkermat(chnkr, fkern1, opts);
D = chunkermat(chnkr, double,  opts);
H = chunkermat(chnkr, hilbert, opts2);

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H ...
    + 2*((1+nu)/2)^2*(dadn(:)./a_var(:)).*D*D - diag(dadn(:)./a_var(:)).*(1/8)*(1+nu).^2 + 0.5*diag(dadn(:)./a_var(:));
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);
sysmat(1:2:end,1:2:end) = sysmat(1:2:end,1:2:end) - (1-nu)*(sysmat(1:2:end,2:2:end).*(dadt./a_var))*H;
sysmat(2:2:end,1:2:end) = sysmat(2:2:end,1:2:end) - (1-nu)*(sysmat(2:2:end,2:2:end).*(dadt./a_var))*H;

Djump = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];
Djump = kron(eye(chnkr.npt), Djump);
sys = Djump + sysmat;
fprintf('%5.2e s : time to assemble matrix\n', toc(start))

%% Right hand side and solve

% RHS from exterior point sources evaluated on boundary
srcinfo = []; srcinfo.r = sources;
rhs = kern2(srcinfo, chnkr) * strengths;

% Reference solution at interior targets from point sources
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = targets;
utarg = kern1(srcinfo, targinfo) * strengths;

start = tic;
sol = gmres(sys, rhs, [], 1e-10, 200);
fprintf('%5.2e s : time for dense gmres\n', toc(start))

% Evaluate at targets and compute error
dens_comb = zeros(3*chnkr.npt, 1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = -H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

ikern = @(s,t) flex2d.kern(zk, s, t, 'free_plate_eval_var', nu);
Dsol = chunkerkerneval(chnkr, ikern, dens_comb, targets);

relerr = norm(utarg - Dsol, 'fro') / (sqrt(chnkr.nch)*norm(utarg, 'fro'));
fprintf('relative error: %5.2e\n', relerr)

assert(relerr < 1e-8);
