
test_surferinterior0()
test_surferinterior1()
function test_surferinterior0()
S = geometries.startorus([1,0.5,0.05], 5, [], [5,10], 6, 11);
% figure(1); clf
% plot(S)
% hold on
% wireframe(S,struct('wfill',0))

nplot = 80;
xx = linspace(0,2,nplot);
yy = linspace(-1,1,nplot);
[X,Y] = meshgrid(xx,yy);
targs = []; targs.r = [0*X(:).';X(:).';Y(:).'];

in = surferinterior(S,targs);

% figure(3);clf
% wireframe(S,struct('wfill',0))
% hold on
% scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],in,'.'); colorbar()
% hold off

hs = 10.^(-3:-1);
hs = [-hs,hs];
targs = [];
targs.r = S.r(:,1) + hs.*S.n(:,1);
in = surferinterior(S,targs);
% hold on
% scatter3(targs.r(1,in),targs.r(2,in),targs.r(3,in),'r.');
% scatter3(targs.r(1,~in),targs.r(2,~in),targs.r(3,~in),'b.');
% hold off
% axis equal

% Points shifted inward (negative hs) should be interior and vice versa
assert(all(in(hs<0)),  'surferinterior: inward points should be interior');
assert(all(~in(hs>0)), 'surferinterior: outward points should be exterior');


end




function test_surferinterior1()

% on a disk
S = geometries.disk([],[],[3 3 5],10,1);

rng(1)
ntarg = 200;
bdryrs = randn(3,ntarg);
bdryrs(3,:) = 0;
bdryrs = bdryrs./vecnorm(bdryrs);

intrs = randn(3,ntarg);
intrs(3,:) = 0;
intrs = rand(1,ntarg).*intrs./vecnorm(intrs);

rs = [bdryrs,intrs];

tol = 1e-12;
opts = []; opts.tol = tol; 
[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, struct('r',rs),opts);

% figure(1);clf
% wireframe(S),hold on
% scatter(rs(1,:),rs(2,:),[],log10(dists));colorbar

assert(norm(flags)==0)

assert(max(dists)<10*tol)
assert(norm(rs - sxyz)<10*tol*norm(rs))

% on a curved surface
abc = [1,1.2,0.8];
S1 = geometries.ellipsoid(abc,4*[1,1,1],[],6);
S2 = geometries.ellipsoid(abc,4*[1,1,1],[],8);
S1 = slicesurfer(S1,find(S1.cms(3,:)>0));
S2 = slicesurfer(S2,find(S2.cms(3,:)<0));
S = merge([S1,S2]);

rs = randn(3,ntarg);

rs = rs./vecnorm(rs./abc.');

tol = 1e-10;
opts = []; opts.tol = tol; 
[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, rs);

% figure(3);clf
% wireframe(S,struct('wfill',1))
% hold on
% scatter3(rs(1,:),rs(2,:),rs(3,:),[],log10(dists))
% colorbar

assert(norm(flags)==0)
assert(max(dists)<1e-5)
assert(norm(rs - sxyz)<1e-5*norm(rs))

% off surface

r_check = NaN*rs;
for p = 1:S.npatches
    idx = find(patch_inds == p);
    if isempty(idx), continue; end
    pols_p = koorn.pols(S.norders(p), uvsloc(:,idx));
    r_check(:,idx) = S.srccoefs{p}(1:3,:) * pols_p;
end
assert(norm(r_check - sxyz)<1e-14*norm(r_check))

S = geometries.ellipsoid(abc,8*[1,1,1],[],6,11);
patch_ids = randi(S.npatches,1,ntarg);
uvs = 2*rand(2,ntarg)-1;
rs = zeros(3,ntarg);
ns = zeros(3,ntarg);
for p = 1:S.npatches
    idx = find(patch_ids == p);
    if isempty(idx), continue; end
    pols_p = polytens.lege.pols(S.norders(p), uvs(:,idx));
    rs(:,idx) = S.srccoefs{p}(1:3,:) * pols_p;

    us = S.srccoefs{p}(4:6,:) * pols_p;
    vs = S.srccoefs{p}(7:9,:) * pols_p;
    ns(:,idx) = cross(us,vs);
end
ns = ns./vecnorm(ns);
rs = rs+1e-2*ns;

[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, rs);
assert(norm(patch_ids(:) - patch_inds(:))==0)
assert(norm(dists-1e-2)<1e-10)

assert(norm(flags)==0)


end
