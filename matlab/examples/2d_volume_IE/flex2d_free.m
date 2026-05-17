% solving a \Delta^2 u - b \Delta u - c u + V u = 0 w/ free BCs

S = geometries.disk([],[],[4 4 4],6);

a = 1;
b = 0.7;
c = 1/pi;
% c = 4*a+2*b;
nu = .3;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

% zk=1;
zk = [zk1,zk2];

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t),cparams);
chnkr = sort(chnkr);

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

V = eval_gauss(S.r);
%%

% volume to volume part 
eps = 1e-10;

start = tic; 

A = flex2d.v2v_matgen(S,zk,eps);
l11 = a*eye(S.npts) + (V).*A;

t1 = toc(start);

fprintf('%5.2e s : time to assemble v2v matrix\n',t1)


% volume to boundary part 

%%

targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
targinfo.d = [chnkr.d(:,:);0*chnkr.d(1,:)];
targinfo.du = targinfo.d;
targinfo.d2 = [chnkr.d2(:,:);0*chnkr.d2(1,:)];
kappa = signed_curvature(chnkr);
targinfo.kappa = kappa(:);
[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, targinfo);
targinfo.patch_id = patch_inds;
targinfo.uvs_targ = uvsloc;
%%

start = tic; 
v2b_supp = flex2d.v2b_matgen_supp2(S,zk,nu,targinfo,eps,targinfo.patch_id,targinfo.uvs_targ);
v2b_free = flex2d.v2b_matgen_free2(S,zk,nu,targinfo,eps,targinfo.patch_id,targinfo.uvs_targ);
l21 = zeros(2*chnkr.npt,S.npts);
l21(1:2:end,:) = v2b_supp;
l21(2:2:end,:) = v2b_free;
% %%
% 
% dpars = nu;
% getnearquad = @(varargin) flex2d.getnearquad(varargin{1:17},dpars,zk,varargin{18},'free');
% % [flex_v2b_cor,norderup] = edge_quad_cor_sub(geom, rhskern, eps, getnearquad, geom.chnkrdim, norderup);
% rhskern = @(s,t) chnk.flex2d.kern(zk,s,t,'free_plate_bcs',nu);
% [flex_v2b_cor,norderup] = wrap_getnearquad(S, rhskern, eps, getnearquad, targinfo);
% [S_over,xinterp] = oversample(S,S.norders + norderup);
% l21_2 = flex_v2b_cor + (rhskern(S_over,targinfo).*S_over.wts(:).')*xinterp;
% [norm(l21 - l21_2), norm(l21(1:2:end,:) - l21_2(1:2:end,:)),...
% norm(l21(2:2:end,:) - l21_2(2:2:end,:))]
% 
% [norm(l21*V - l21_2*V), norm(l21(1:2:end,:)*V - l21_2(1:2:end,:)*V),...
% norm(l21(2:2:end,:)*V - l21_2(2:2:end,:)*V)]

t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)
%%

% boundary to boundary part 

% build system matrix

fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate', nu);        % build the desired kernel
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) - sysmat1(3:4:end,1:2:end)*H  + 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) - sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

D = -[-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

l22 =  D + sysmat;

t1 = toc(start);

fprintf('%5.2e s : time to assemble b2b matrix\n',t1)

% boundary to volume part 
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval', nu);
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);

start = tic; 
l12 = zeros(S.npts,2*chnkr.npt);
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);
l12(:,1:2:end) = V.*(b2v(:,1:3:end) - b2v(:,2:3:end)*H);
l12(:,2:2:end) = V.*b2v(:,3:3:end);
t1 = toc(start);
fprintf('%5.2e s : time to assemble b2v matrix\n',t1)


%% form system matrix and rhs 
lhs = [l11, l12; 
    l21, l22];


rhs_vol = (4*a+2*b-c+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
% rhs_vol = (4-zk^2+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_bc = zeros(chnkr.npt*2,1);
% rhs_bc(1:2:end) = -(1+nu)*sin(chnkr.r(1,:)).*sin(chnkr.r(2,:)) + ...
%     2*(1-nu)*cos(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(1,:).*chnkr.n(2,:);

hess = zeros(chnkr.npt,3);
hess(:,1) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
hess(:,2) = cos(chnkr.r(1,:)).*cos(chnkr.r(2,:));
hess(:,3) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));

third = zeros(chnkr.npt,4);
third(:,1) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:));
third(:,2) = -sin(chnkr.r(1,:)).*cos(chnkr.r(2,:));
third(:,3) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:));
third(:,4) = -sin(chnkr.r(1,:)).*cos(chnkr.r(2,:));

nx = chnkr.n(1,:).'; 
ny = chnkr.n(2,:).';

dx = chnkr.d(1,:).';
dy = chnkr.d(2,:).';

ds = sqrt(dx.*dx+dy.*dy);
taux = (dx./ds);                                                                       % normalization
tauy = (dy./ds);

rhs_bc(1:2:end) = ((hess(:,1).*(nx.*nx) + hess(:,2).*(2*nx.*ny) + hess(:,3).*(ny.*ny)) + ...
           nu.*(hess(:,1).*(taux.*taux) + hess(:,2).*(2*taux.*tauy) + hess(:,3).*(tauy.*tauy)));

rhs_bc(2:2:end) = ((third(:,1).*(nx.*nx.*nx) + third(:,2).*(3*nx.*nx.*ny) +...
       third(:,3).*(3*nx.*ny.*ny) + third(:,4).*(ny.*ny.*ny))  + ...
        (2-nu).*(third(:,1).*(taux.*taux.*nx) + third(:,2).*(taux.*taux.*ny + 2*taux.*tauy.*nx) +...
        third(:,3).*(2*taux.*tauy.*ny+ tauy.*tauy.*nx) +...
        + third(:,4).*(tauy.*tauy.*ny)) + ...
        (1-nu).*kappa(:).*(hess(:,1).*taux.*taux + hess(:,2).*(2*taux.*tauy) + hess(:,3).*tauy.*tauy+...
        -(hess(:,1).*nx.*nx + hess(:,2).*(2*nx.*ny) + hess(:,3).*ny.*ny)));

rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;

% solve 

start = tic; 
sol = gmres(lhs,rhs,[],eps,100); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)    

dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(S.npts+1:2:end);
dens_comb(2:3:end) = -H*sol(S.npts+1:2:end);
dens_comb(3:3:end) = sol(S.npts+2:2:end);

ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);
start1 = tic;
u = A*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    dens_comb, S.r(1:2,:));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
err = abs(u - ref_u(:)) / max(abs(u));


figure(2); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); 
colorbar
%%
figure(3);clf
subplot(1,3,1)
plot(S,real(u));colorbar
subplot(1,3,2)
plot(S,real(ref_u));colorbar
subplot(1,3,3)
plot(S,real(u-ref_u));colorbar

%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end

