% Test the flat geometries (geometries.disk and geometries.square):
%
%   - the area is correct
%   - the normals are (0,0,+-1) and agree with du x dv
%   - a polynomial is integrated correctly
%   - interpolate_data reproduces a smooth function
%

norder  = 8;
iptypes = [1, 11, 12];
iorts   = [1, -1];
names   = {'disk', 'square'};

% reference area, and integral of x^2 y^2 + x^2
arearef = [pi, 4];
polyref = [pi/24 + pi/4, 4/9 + 4/3];
tols    = [1e-9, 1e-12];

% targets for the interpolation test
uvs_targ = [0.1 0.8; 0.3 0.5; 0.7 0.1; 0.24 0.31].';

rng(319);

for igeo = 1:2
    for iptype = iptypes
        for iort = iorts

            if igeo == 1
                S = geometries.disk([], [], [3 3 3], norder, iptype, iort);
            else
                S = geometries.square([6 5], norder, iptype, iort);
            end
            tol = tols(igeo);
            str = sprintf('%-6s iptype = %2d, iort = %2d:', ...
                names{igeo}, iptype, iort);

            % flat, in the z = 0 plane
            assert(norm(S.r(3,:), inf) < 1e-14, '%s not in z = 0', str);

            % area
            err = abs(area(S) - arearef(igeo))/arearef(igeo);
            fprintf('%s area err    %5.2e\n', str, err);
            assert(err < tol, '%s area error %5.2e', str, err);

            % normals are (0,0,sign(iort)) and agree with du x dv
            nref = repmat([0; 0; sign(iort)], 1, S.npts);
            err = norm(S.n - nref, inf);
            assert(err < 1e-13, '%s normals not (0,0,%d)', str, sign(iort));
            dn = cross(S.du, S.dv);
            dn = dn./vecnorm(dn);
            err = norm(dn - nref, inf);
            fprintf('%s normal err  %5.2e\n', str, err);
            assert(err < 1e-12, '%s du x dv does not match normal', str);

            % integral of a polynomial
            x = S.r(1,:).'; y = S.r(2,:).';
            val = sum((x.^2.*y.^2 + x.^2).*S.wts(:));
            err = abs(val - polyref(igeo))/polyref(igeo);
            fprintf('%s integral err %5.2e\n', str, err);
            assert(err < tol, '%s integral error %5.2e', str, err);

            % interpolation of a smooth function
            ipatchids = randi(S.npatches, [4,1]);
            dens = exp(-x.^2 - y.^2);
            xs = S.interpolate_data(x, ipatchids, uvs_targ);
            ys = S.interpolate_data(y, ipatchids, uvs_targ);
            vals = S.interpolate_data(dens, ipatchids, uvs_targ);
            err = norm(vals(:) - exp(-xs(:).^2 - ys(:).^2), inf);
            fprintf('%s interp err  %5.2e\n', str, err);
            assert(err < 1e-8, '%s interpolation error %5.2e', str, err);

        end
    end
end
