
S = geometries.startorus([1,0.5,0.05], 5, [], [5,10], 8, 11);
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
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],in,'.');colorbar
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






