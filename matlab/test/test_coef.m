S1 = geometries.sphere(1, 2, [0;0;0],4,1);
S2 = geometries.sphere(1, 2, [3;0;0],5,1);
S3 = geometries.sphere(1, 2, [0;3;0],5,11);
S4 = geometries.sphere(1, 2, [0;0;3],5,12);
S = merge([S1,S2,S3,S4]);