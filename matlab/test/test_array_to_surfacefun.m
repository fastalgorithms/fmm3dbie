% choose domain type
run ../startup.m

if ~exist('surfacemesh'); return; end;
geomtype = 'sphere';
funtype = 'spharm';

iref = 0;
if(strcmpi(geomtype,'sphere'))
    domtmp = surfacemesh.sphere(7,3+iref);
elseif(strcmpi(geomtype,'stellarator'))
    nref = 2^iref;
    nu = 8*nref;
    nv = 24*nref;
    domtmp = surfacemesh.stellarator(12,nu,nv);    
end
[S,dom] = surfer.surfacemesh_to_surfer(domtmp);

% Choose function type

if(strcmpi(funtype,'spharm') && strcmpi(geomtype,'sphere'))
    ndeg = 2;
    f = spherefun.sphharm(ndeg,1);
else
    f = @(x,y,z) 1.0./sqrt((x-2.2).^2 + (y-1.1).^2 + z.^2);
end

frhs = surfacefun(@(x,y,z) f(x,y,z),dom);

rhs = f(S.r(1,:),S.r(2,:),S.r(3,:));
rhs = rhs(:);

frhs2 = array_to_surfacefun(rhs,dom,S);
err1 = norm(frhs-frhs2)/norm(frhs);
fprintf('Error in interpolant = %d\n',err1);


