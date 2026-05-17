% solving a \Delta^2 u - b \Delta u - c u + V u = 0 w/ clamped BCs

S = geometries.disk([],[],[4 4 4],6);

a = 1;
b = 0.7;
c = 1/pi;
% c = 4*a+2*b;

zk1 = sqrt((- b + sqrt(b^2 + 4*a*c)) / (2*a));
zk2 = sqrt((- b - sqrt(b^2 + 4*a*c)) / (2*a));

zk = [zk1 zk2];

cparams = []; cparams.maxchunklen = 4/max(abs(zk));
chnkr = chunkerfunc(@(t) ellipse(t),cparams);
chnkr = sort(chnkr);

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
wireframe(S,struct('wfill',0))
view(0,90)

V = 0*eval_gauss(S.r);

% volume to volume part 

start = tic; 

A = flex2d.v2v_matgen(S,zk,1e-8);
l11 = a*eye(S.npts) + V.*A;

t1 = toc(start);

fprintf('%5.2e s : time to assemble v2v matrix\n',t1)



% boundary to volume part 
fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
targetinfo = [];
targetinfo.r = S.r(1:2,:);
targetinfo.n = S.n(1:2,:);


start = tic; 
b2v = chunkerkernevalmat(chnkr,fkern,targetinfo);
l12 = V.*b2v;
t1 = toc(start);
fprintf('%5.2e s : time to assemble b2v matrix\n',t1)

% volume to boundary part 


targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];

start = tic; 
v2b_dir = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
v2b_neu = flex2d.v2b_matgen_neu(S,zk,targinfo,1e-8);
l21 = zeros(2*chnkr.npt,S.npts);
l21(1:2:end,:) = v2b_dir;
l21(2:2:end,:) = v2b_neu;

getnearquad = @(varargin) flex2d.getnearquad(varargin{1:17},[],zk,varargin{18},'clamped');
% [flex_v2b_cor,norderup] = edge_quad_cor_sub(geom, rhskern, eps, getnearquad, geom.chnkrdim, norderup);
rhskern = @(s,t) chnk.flex2d.kern(zk,s,t,'clamped_plate_bcs');
[flex_v2b_cor,norderup] = wrap_getnearquad(S, rhskern, eps, getnearquad, targinfo);
[S_over,xinterp] = oversample(S,S.norders + norderup);
l21 = flex_v2b_cor + (rhskern(S_over,targinfo).*S_over.wts(:).')*xinterp;

t1 = toc(start);
fprintf('%5.2e s : time to assemble v2b matrix\n',t1)


% boundary to boundary part 

fkern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate');
kappa = signed_curvature(chnkr);
kappa = kappa(:);

opts = [];
opts.sing = 'log';

start = tic;
b2b = chunkermat(chnkr,fkern, opts);
l22 = b2b + 0.5*eye(2*chnkr.npt);
l22(2:2:end,1:2:end) = l22(2:2:end,1:2:end) - kappa.*eye(chnkr.npt);

t1 = toc(start);
fprintf('%5.2e s : time to assemble b2b matrix\n',t1)

%%
% rhs = - a \Delta^2 ui + b \Delta ui + c ui - V ui 

% form system matrix and rhs 
lhs = [l11, l12; 
    l21, l22];


rhs_vol = (-4*a-2*b+c-V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_bc = zeros(chnkr.npt*2,1);
rhs_bc(1:2:end) = -sin(chnkr.r(1,:)).*sin(chnkr.r(2,:));
rhs_bc(2:2:end) = -cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          - sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; 



rhs = zeros(S.npts+2*chnkr.npt,1);
rhs(1:S.npts) = rhs_vol;
rhs(S.npts+1:end) = rhs_bc;


% solve 

start = tic; 
sol = gmres(lhs,rhs,[],1e-10,100); 
t1 = toc(start);
fprintf('%5.2e s : time for dense gmres\n',t1)   



ikern = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');
start1 = tic;
u = A*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), S.r(1:2,:));
u = real(u);
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

ref_u = (sin(S.r(1,:)).*sin(S.r(2,:))).';
err = abs(u + ref_u(:)) / max(abs(u));


figure; clf
scatter(S.r(1,:),S.r(2,:),8,log10(err));
title('error')
colorbar

return

%%

nx = 50;
x = linspace(-1,1,nx);
y = linspace(-1,1,nx);
[x,y] = ndgrid(x,y); 
x = x(:).';
y = y(:).';
idx = x.^2+y.^2<1;
xi = x(idx);
yi = y(idx);
targinfo = [];
targinfo.r = [xi; yi; zeros(size(xi))];
targinfo.n = zeros(size(targinfo.r));
targinfo.n(3,:) = 1;


A2 = flex2d.v2b_matgen_dir(S,zk,targinfo,1e-8);
u = A2*sol(1:S.npts)+chunkerkerneval(chnkr, ikern,...
    sol(S.npts+1:end), targinfo.r(1:2,:));
u = real(u);

U = NaN(size(x));
U(idx) = u;
U = reshape(U, [nx,nx]);

figure; clf
imagesc(U); 
set(gca, 'Color', 'w');
colorbar




%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end

