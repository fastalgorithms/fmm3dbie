function [obj] = ellipsoid(abc,rmax,norder,nref)
% Generate the triangulation of an ellipsoid
% by considering it an axissymmetric problem
%
% All input arguments are optional, default values are indicated in
% brackets
% 
% Input
%  - abc: double(3) semi-axes along coordinate directions of the ellipsoid
%           ([1;1;1] == unit sphere)
%  - rmax: max patch diameter (0.25)
%  - norder: order of discretization on each patch (3)
%  - nref: additional refinement of patches (0)
    if(nargin < 1)
        abc = [1;1;1];
    end
    
    if(nargin < 2)
        rmax = 0.25;
    end
    
    if(nargin < 3)
        norder = 3;
    end
    
    if(nargin < 4)
        nref = 0;
    end
   
    
    a = abc(1);
    b = abc(2);
    c = abc(3);
    nthet = ceil(2*c/rmax)*2^nref;
    hthet = pi/nthet;
    
    ts = zeros(2,nthet);
    ts(1,:) = (0:(nthet-1))*hthet;
    ts(2,:) = (1:nthet)*hthet;
    tuse = ts(1,:);
    for i=1:nthet
        if(abs(ts(1,i)-pi/2)>abs(ts(2,i)-pi/2)), tuse(i) = ts(2,i); end
    end
    tt = (ts(1,:)-pi/2).*(ts(2,:)-pi/2);
    tuse(tt<=0) = pi/2;
    
    
    alpha = a*sin(tuse);
    beta = b*sin(tuse);
    hh = ((alpha-beta)./(alpha+beta)).^2;
    ellip_p = pi*(alpha+beta).*(1 + 3*hh./(10+sqrt(4-3*hh)));
    nphis = ceil(ellip_p/rmax)*2^nref;
    npatches = sum(nphis)*2;
    
    tri = [];
    verts = [];
    istart = 0;
    for i=1:nthet
        [uu,vv] = meshgrid(0:1,0:nphis(i));
        uu = uu(:)*hthet + (i-1)*hthet;
        vv = vv(:)/nphis(i)*2*pi;
        tri0 = delaunay(uu,vv)+istart;
        n = length(uu);
        istart = istart+n;
        
        tri = [tri; tri0];
        verts = [verts, [uu, vv, zeros(size(uu))]'];
    end

    % collect all list of vertices
    v1s = verts(1:3,tri(:,1));
    v2s = verts(1:3,tri(:,2));
    v3s = verts(1:3,tri(:,3));

    % Generate rv nodes on a single patch
    rnodes = koorn.rv_nodes(norder);
 
    % compute the projection along with all derivatives
    ru  = rnodes(1,:)';
    rv = rnodes(2,:)';
    np = length(ru);

    x = repmat(v1s(1,:),[np,1]) + ru*(v2s(1,:)-v1s(1,:)) + ...
         rv*(v3s(1,:)-v1s(1,:));

    y = repmat(v1s(2,:),[np,1]) + ru*(v2s(2,:)-v1s(2,:)) + ...
         rv*(v3s(2,:)-v1s(2,:));

    z = repmat(v1s(3,:),[np,1]) + ru*(v2s(3,:)-v1s(3,:)) + ...
         rv*(v3s(3,:)-v1s(3,:));


    dxdu0 = v2s(1,:)-v1s(1,:);
    dydu0 = v2s(2,:)-v1s(2,:);
    dzdu0 = v2s(3,:)-v1s(3,:);

    dxdv0 = v3s(1,:)-v1s(1,:);
    dydv0 = v3s(2,:)-v1s(2,:);
    dzdv0 = v3s(3,:)-v1s(3,:);


    dxdu = repmat(dxdu0,[np,1]);
    dydu = repmat(dydu0,[np,1]);
    dzdu = repmat(dzdu0,[np,1]);

    dxdv = repmat(dxdv0,[np,1]);
    dydv = repmat(dydv0,[np,1]);
    dzdv = repmat(dzdv0,[np,1]);

    
    npts = npatches*length(ru);
    srcvals = zeros(12,npts);

    xuse = a*sin(x).*cos(y);
    yuse = b*sin(x).*sin(y);
    zuse = c*cos(x);

    srcvals(1,:) = xuse(:)';
    srcvals(2,:) = yuse(:)';
    srcvals(3,:) = zuse(:)';

    dxduuse = a*(cos(x).*cos(y).*dxdu - sin(x).*sin(y).*dydu);
    dyduuse = b*(cos(x).*sin(y).*dxdu + sin(x).*cos(y).*dydu);
    dzduuse = -c*sin(x).*dxdu;

    srcvals(4,:) = dxduuse(:)';
    srcvals(5,:) = dyduuse(:)';
    srcvals(6,:) = dzduuse(:)';

    dxdvuse = a*(cos(x).*cos(y).*dxdv - sin(x).*sin(y).*dydv);
    dydvuse = b*(cos(x).*sin(y).*dxdv + sin(x).*cos(y).*dydv);
    dzdvuse = -c*sin(x).*dxdv;

    srcvals(7,:) = dxdvuse(:)';
    srcvals(8,:) = dydvuse(:)';
    srcvals(9,:) = dzdvuse(:)';

    dn = cross(srcvals(4:6,:),srcvals(7:9,:));
%   srcvals(10:12,:) = dn./repmat(vecnorm(dn,1),[3,1]);
    srcvals(10:12,:) = dn./repmat(vecnorm(dn),[3,1]);

    obj = surfer(npatches,norder,srcvals); 

end
