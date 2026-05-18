%DRIVER_FELIPE_2 Minimal CAD-skeleton multiscale mesher visualization.

driverFile = mfilename('fullpath');
examplesDir = fileparts(driverFile);
matlabDir = fileparts(examplesDir);
repoRoot = fileparts(matlabDir);

run(fullfile(matlabDir, 'startup.m'));

fname = fullfile(repoRoot, 'geometries', 'meshes', 'cylinder.gidmsh');
fcad = fullfile(repoRoot, 'geometries', 'meshes', 'cylinder_skeleton.txt');

opts = struct();
opts.nrefine = 1;
opts.nquad = 16;
opts.rlam = 5;
opts.fcad = fcad;

Sall = multiscale_mesher(fname, 4, opts);
S = Sall{end};

figure('Name', 'driver_felipe_2: CAD skeleton cylinder');
plot(S);
axis equal;
view(3);
colorbar;
title('cylinder.gidmsh smoothed with exact CAD skeleton', ...
    'Interpreter', 'none');
xlabel('x');
ylabel('y');
zlabel('z');
