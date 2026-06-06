
% test_surferinterior0()
test_surferinterior1()
function test_surferinterior0()
S = geometries.startorus([1,0.5,0.05], 5, [], [5,10], 6, 11);
% figure(1); clf
% plot(S)
% hold on
% wireframe(S,struct('wfill',0))

nplot = 100;
xx = linspace(0,2,nplot);
yy = linspace(-1,1,nplot);
[X,Y] = meshgrid(xx,yy);
targs = []; targs.r = [0*X(:).';X(:).';Y(:).'];

in = surferinterior(S,targs);

figure(3);clf
wireframe(S,struct('wfill',0))
hold on
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],in,'.'); colorbar()
hold off

hs = 10.^(-3:-1);
hs = [-hs,hs];
targs = [];
targs.r = S.r(:,1) + hs.*S.n(:,1);
in = surferinterior(S,targs);
hold on
scatter3(targs.r(1,in),targs.r(2,in),targs.r(3,in),'r.');
scatter3(targs.r(1,~in),targs.r(2,~in),targs.r(3,~in),'b.');
hold off
axis equal

% Points shifted inward (negative hs) should be interior and vice versa
assert(all(in(hs<0)),  'surferinterior: inward points should be interior');
assert(all(~in(hs>0)), 'surferinterior: outward points should be exterior');


end




function test_surferinterior1()

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

figure(1);clf
wireframe(S),hold on
scatter(rs(1,:),rs(2,:),[],log10(dists));colorbar

assert(norm(flags)==0)

assert(max(dists)<10*tol)
assert(norm(rs - sxyz)<10*tol*norm(rs))

abc = [1,1.2,0.8];
S = geometries.ellipsoid(abc,4*[1,1,1],[],6);


rs = randn(3,ntarg);

rs = rs./vecnorm(rs./abc.');

tol = 1e-10;
opts = []; opts.tol = tol; 
[sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, rs);

figure(3);clf
wireframe(S,struct('wfill',1))
hold on
scatter3(rs(1,:),rs(2,:),rs(3,:),[],log10(dists))
colorbar

assert(norm(flags)==0)
assert(max(dists)<1e-5)
assert(norm(rs - sxyz)<1e-5*norm(rs))

r_check = NaN*rs;
for p = 1:S.npatches
    idx = find(patch_inds == p);
    if isempty(idx), continue; end
    pols_p = koorn.pols(S.norders, uvsloc(:,idx));
    r_check(:,idx) = S.srccoefs{p}(1:3,:) * pols_p;
end
assert(norm(r_check - sxyz)<1e-14)
end
