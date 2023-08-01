%
%  Additional dependency: surface-hps
%
run ~/git/surface-hps/setup.m
addpath(genpath('~/git/fmm3dbie/matlab'))


% Test surface_mesh conversion script

norder = 5;
f = @(u,v) cos(3.2*u + 4.1*v);

dom = surfacemesh.sphere(norder+1,1);

[n,~] = size(dom.x{1});

n2 = n;
xleg = legpts(n2);
xcheb = chebpts(n,[-1,1]);
eval = barymat(xleg,xcheb);

x2 = eval * dom.x{5} * eval.';
% x2 = x2(:);



xcoeffs = surfacefun.vals2coeffs(dom.x);
xcoeffs{5};

uvs = polytens.cheb_nodes(norder,2);
umat = polytens.cheb_vals2coefs(norder,uvs);
xuse = dom.x{5}';
xuse = xuse(:);
xcoeffs2 = umat*xuse;
xcoeffs2 = reshape(xcoeffs2,[n,n]);
xcoeffs2test = xcoeffs2';



fprintf('error in vals to coefs=%d\n',norm(xcoeffs2test-xcoeffs{5}));


uvs2 = polytens.lege_nodes(n2-1);
vmat = polytens.cheb_coefs2vals(norder,uvs2);

xuse = xcoeffs2;
xuse = xuse(:);
x3 = vmat*xuse;
x3 = reshape(x3,[n2,n2]);

fprintf('error in interpolant=%d\n',norm(x3-x2'));

%% Test normals
rndom = normal(dom);
xyzu  = surfacefunv(surfacefun(dom.xu,dom),surfacefun(dom.yu,dom),surfacefun(dom.zu,dom));
xyzv  = surfacefunv(surfacefun(dom.xv,dom),surfacefun(dom.yv,dom),surfacefun(dom.zv,dom));

rndom2 = cross(xyzu,xyzv);
rndom2 = rndom2./norm(rndom2);
vv = norm(rndom-rndom2);
err1 = norm(vv.vals{1});
fprintf('Error in normals before fix (should be large)=%d\n',err1);

x = dom.x;
y = dom.y;
z = dom.z;

x{1} = x{1}';
y{1} = y{1}';
z{1} = z{1}';
dom2 = surfacemesh(x,y,z);

rndom = normal(dom2);
xyzu  = surfacefunv(surfacefun(dom2.xu,dom2),surfacefun(dom2.yu,dom),surfacefun(dom2.zu,dom));
xyzv  = surfacefunv(surfacefun(dom2.xv,dom2),surfacefun(dom2.yv,dom),surfacefun(dom2.zv,dom));
rndom2 = cross(xyzu,xyzv);
rndom2 = rndom2./norm(rndom2);
vv2 = norm(rndom-rndom2);
err2 = norm(vv2.vals{1});
fprintf('Error in normals after fix=%d\n',err2);
