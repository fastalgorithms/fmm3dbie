% Tester for surfer scaling
run ../startup.m

% Unit sphere in origin
R = 1;
S = geometries.sphere(R, 6, [0;0;0]);

% Error in surface area before scaling
aerr_sphere = abs(4*pi*R^2-area(S));
fprintf('Error in area before scaling=%d\n',aerr_sphere);

% Scale to prolate spheroid using vector [a,b,c]
a = 2.1; b = 2.1; c = 3.4;
S_spheroid = S.scale([a,b,c]);

% Validate spheroidal surface area
ecc = sqrt(1-a^2/c^2);
aref = 2*pi*a^2*(1+c/a/ecc*asin(ecc));
aerr_spheroid = abs(aref-area(S_spheroid));
fprintf('Error in area after vector scaling=%d\n',aerr_spheroid);

% Scale to ellipsoid using vector [a,b,c]
a = 0.1; b = 2.1; c = 3.4;
S_ellipsoid = S.scale([a,b,c]);

% Validate scaling factors of coordinates in each direction
assert(max(S_ellipsoid.r(1,:))/max(S.r(1,:))==a);
assert(max(S_ellipsoid.r(2,:))/max(S.r(2,:))==b);
assert(max(S_ellipsoid.r(3,:))/max(S.r(3,:))==c);

% Scalar scaling of sphere using scalar a
a = 1.37;
S_scalar_scaled = S.scale(a);
aerr_sphere_scaled = abs(4*pi*(R*a)^2-area(S_scalar_scaled));
fprintf('Error in area after scalar scaling=%d\n',aerr_sphere_scaled);

