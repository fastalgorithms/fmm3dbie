iertot = 0;
filename = "sphere.msh";
[ier,geo] = readmsh(filename);
iertot = iertot + ier;
filename = "lens_r00.msh";
[ier,geo] = read_gmsh_v2(filename);
iertot = iertot + ier;
filename = "prism_50.gidmsh";
[ier,geo] = readgidmsh(filename);
iertot = iertot + ier;
filename = "cuboid_a1_b2_c1p3.tri";
[ier,geo] = readtri(filename);
iertot = iertot + ier;
filename = "cow_new_gmshv4.msh";
[ier,geo] =read_gmsh_v4(filename)
iertot = iertot + ier;
scatter3(geo.points(1,:),geo.points(2,:),geo.points(3,:));
filename = "lens_r00_gmshv4.msh";
[ier,geo] =read_gmsh_v4(filename)
iertot = iertot + ier;

iertot