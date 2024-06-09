% example of acoustic sound-hard (Dirichet BC) scattering from ellipsoid
% Barnett 6/7/24
clear; format long g

ord = 6;          % panel order
refine = 1;       % integer refinement factor (increase 1 to 2 to check conv)
S = geometries.ellipsoid([1.5,0.5,1.0], ceil(refine*[5 2 3]), [], ord, 11);
k = 10;           % wavenumber
dir = [-1,4,1]; dir = dir/norm(dir);  % inc direction (unit norm vec)
uinc = helm3d.planewave(k,dir,S.r);

dx = 0.03; L=3; g = -L:dx:L; n = numel(g); z0=0.0; [xx yy] = ndgrid(g,g);
sheet.r = [xx(:) yy(:) z0+0*xx(:)]';  % sheet of test points (3*n^2)
pt.r = [0.5;1.0;1.0];                 % single test pt (col vec)

figure; subplot(1,2,1); S.plot(real(uinc));
hold on; plot3([0 dir(1)], [0 dir(2)], [0 dir(3)], 'k-','linewidth',2) 
uincsh = helm3d.planewave(k,dir,sheet.r);
surf(xx,yy,z0+0*xx, real(reshape(uincsh,[n n]))); shading interp;
h = S.plot_nodes; h.MarkerEdgeColor='black'; h.SizeData=1;  % overlay nodes
xlabel('x'); ylabel('y'); zlabel('z');
plot3(pt.r(1,:),pt.r(2,:),pt.r(3,:),'k.','markersize',10);
colormap(jet(256)); colorbar; caxis(1.5*[-1,1]);
%scatter3(sheet.r(1,:),sheet.r(2,:),sheet.r(3,:),1,real(uincsh), '.');
title('direction, test pt, and u^{inc} on S and sheet'); drawnow

tol = 1e-4;
rhs = -uinc;
fprintf('solving...\n'); tic
dens = helm3d.solver(S, 'dir', rhs, tol, k);         % solve (uses CFIE)
toc, fprintf('eval soln...\n'); tic
uscpt = helm3d.eval(S, 'dir', dens, pt, tol, k)       % eval u scatt @ pt
uscsh = helm3d.eval(S, 'dir', dens, sheet, tol, k);   % eval u scatt @ sheet
toc
subplot(1,2,2); S.plot(zeros(S.npts,1)); hold on;
surf(xx,yy,z0+0*xx, real(reshape(uincsh+uscsh,[n n]))); shading interp;
colormap(jet(256)); colorbar; caxis(1.5*[-1,1]);
title('u^{tot} on sheet');

%exportgraphics(gcf,'helm_scatt_Dir.png')
