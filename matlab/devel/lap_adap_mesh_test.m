clear all;close all
run /home/fryklund/Projects/research/adapmesh3d/startup.m
% Two torii test to get 4 digits of accuracy in the solution
nu = 1;
nv = 3;
a = 1;
d = 0.05;

rmajor = 1;
rminor = 0.5;
rwave = 0.0;

radii2 = [rmajor,rminor,rwave];
radii = [rminor,rmajor,rwave];
nuv = [nu,nv];
norder = 4;
nosc = 1;
scales = [1,1,1];
iptype_pat = 1;
iptyperef = 1;
tol = 1e-3;
tolref = 1e-4;

lptype = 0;

%S = geometries.wobbletorus(radii, 5, [], [], 8, 1)

%WT = geometries.wobbletorus(radii,nosc,scales,nuv,norder,iptype_pat)
S = geometries_adap.wobblytorus_adap(radii,nosc,scales,nuv,norder,iptype_pat,tol,lptype);

[srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
rho = sqrt(srcvals(1,:).^2 + srcvals(2,:).^2).';
cosv = (rho - rmajor) ./ rminor;

H = (rmajor + 2*rminor.*cosv) ./ (2*rminor.*(rmajor + rminor.*cosv));
erra = norm((S.mean_curv - H).*sqrt(S.wts))/S.area;

fprintf('error in resolving curvature = %d\n', erra);
%%
nuvwt = [16,16];
WT = geometries.wobbletorus(radii2,1,scales,nuvwt,norder,iptype_pat);
[srcvalswt,srccoefswt,norderswt,ixyzswt,iptypewc,wtswc] = extract_arrays(WT);
rhowt = sqrt(srcvalswt(1,:).^2 + srcvalswt(2,:).^2).';
cosvwt = (rhowt - rmajor) ./ rminor;

Hwt = (rmajor + 2*rminor.*cosvwt) ./ (2*rminor.*(rmajor + rminor.*cosvwt));
errawt = norm((WT.mean_curv - Hwt).*sqrt(WT.wts))/WT.area;

fprintf('error in resolving curvature = %d\n', errawt);
%%
Sref = geometries_adap.wobblytorus_adap(radii,nosc,scales,2*nuv,norder,...
    iptyperef,tolref);
[srcvalsref,srccoefsref,nordersref,ixyzsref,iptyperef,wtsref] = ...
    extract_arrays(Sref);

rhoref = sqrt(srcvalsref(1,:).^2 + srcvalsref(2,:).^2).';
cosvref = (rhoref - rmajor) ./ rminor;

Href = (rmajor + 2*rminor.*cosvref) ./ (2*rminor.*(rmajor + rminor.*cosvref));
erraref = norm((Sref.mean_curv - Href).*sqrt(Sref.wts))/Sref.area;

fprintf('error in resolving curvature = %d\n', erraref);
%%
tol_st = 1e-7;
norder_st = 14;
ST = geometries_adap.stellarator_adap(nuv,norder_st,iptype_pat,tol_st,[],[],10);
STuni = geometries.stellarator();
%%
%
dpars(1) = 5.0d0;
dpars(2) = 1.9d0;
dpars(3) = 2.0d0*acos(-1.0d0)/3.0d0;
dpars(4) = 2.0d0;


coefs = dpars;
tol_ftube = 1e-4;
norder_ftube = 16;
nuvF = [1 4];
SFtube = geometries_adap.fourier_curve_tube_adap(coefs,nuvF,norder_ftube,iptype_pat,tol_ftube);

%%
tol_tube = 1e-4;
norder_tube = 6;
rmaj = 0.25;
rmin = 0.13;
polfq = 3;
torfq = 2;
rcross = 0.07;

radii = [rmaj,rmin,rcross];
freq = [polfq,torfq];
Stube = geometries_adap.torus_curve_tube_adap(radii,freq,nuv,norder_tube,iptype_pat,tol_tube)


