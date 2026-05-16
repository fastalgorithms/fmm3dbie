
S = geometries.startorus([1,0.5,0.05], 5, [], [5,10], 8, 11);
% S = geometries.stellarator([5,15],8,11);
figure(1); clf
plot(S)
hold on
wireframe(S,struct('wfill',0))


%%
nplot = 100;
xx = linspace(0,2,nplot);
yy = linspace(-1,1,nplot);
[X,Y] = meshgrid(xx,yy);
targs = []; targs.r = [0*X(:).';X(:).';Y(:).'];

in = surferinterior(S,targs);

figure(2);clf
h = pcolor(X,Y,reshape(in+0,size(X))); h.EdgeColor = 'none'; colorbar;
title('val = 1 indicates interior','FontSize',14)

figure(3);clf
wireframe(S,struct('wfill',0))
hold on
scatter3(targs.r(1,:),targs.r(2,:),targs.r(3,:),[],in,'.');colorbar
hold off


%%

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


% %%
% 
% figure(10);clf
% h = plot(S);
% h.FaceColor = 'b';
% 
% figure(11);clf
% h = plot(S);
% h.FaceColor = 'b';
% material('dull')
% camlight('headlight')





