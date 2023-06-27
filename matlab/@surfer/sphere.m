function [obj] = sphere(norder,nu,nref)
% Generate the triangulation of a sphere
% through projecting out a cube discretized with
% nu \times nu patches in each direction, refined
% by nref levels of refinement.
%
%
% Input
%  - norder: order of discretization on each patch
%  - nu: number of patches in u,v direction on each face of cube
%  - nref: additional refinement of patches in each direction
    if(nargin < 1)
        norder = 3;
    end
    
    if(nargin < 2)
        nu = 2;
    end
    
    if(nargin < 3)
        nref = 0;
    end
    
    nuuse = nu*2^nref;
    nverts = 6*(nuuse+1)*(nuuse+1);
    verts = zeros(3,nverts);

    [uu,vv] = meshgrid(0:nuuse);
    uu = uu/(nuuse)*2-1;
    vv = vv/(nuuse)*2-1;
    uu = uu(:);
    vv = vv(:);
    nv0 = length(uu);

    tri0 = delaunay(uu,vv);

    [nt,~] = size(tri0);

    npatches = 6*nt;
    tri = zeros(npatches,3);

    % +z face
    verts(1,1:nv0) = uu;
    verts(2,1:nv0) = vv;
    verts(3,1:nv0) = 1;
    tri(1:nt,:) = tri0;

    % -z face
    istart = nv0+1;
    iind = (istart):(istart+nv0-1);
    itind = (nt+1):(2*nt);
    verts(1,iind) = uu;
    verts(2,iind) = vv;
    verts(3,iind) = -1;
    tri(itind,:) = fliplr(tri0) + istart-1;

    % + x face
    istart = 2*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (2*nt+1):(3*nt);
    verts(1,iind) = 1;
    verts(2,iind) = uu;
    verts(3,iind) = vv;
    tri(itind,:) = tri0 + istart-1;


    % -x face
    istart = 3*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (3*nt+1):(4*nt);
    verts(1,iind) = -1;
    verts(2,iind) = uu;
    verts(3,iind) = vv;
    tri(itind,:) = fliplr(tri0) + istart-1;


    % + y face
    istart = 4*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (4*nt+1):(5*nt);
    verts(1,iind) = vv;
    verts(2,iind) = 1;
    verts(3,iind) = uu;
    tri(itind,:) = tri0 + istart-1;


    % -y face
    istart = 5*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (5*nt+1):(6*nt);
    verts(1,iind) = vv;
    verts(2,iind) = -1;
    verts(3,iind) = uu;
    tri(itind,:) = fliplr(tri0) + istart-1;

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

    a = v1s(1,:).*dxdu0 + v1s(2,:).*dydu0 + v1s(3,:).*dzdu0;
    b = dxdu0.*dxdv0 + dydu0.*dydv0 + dzdu0.*dzdv0;
    c = dxdu0.*dxdu0 + dydu0.*dydu0 + dzdu0.*dzdu0;




    e = v1s(1,:).*dxdv0 + v1s(2,:).*dydv0 + v1s(3,:).*dzdv0;
    f = b;
    g = dxdv0.*dxdv0 + dydv0.*dydv0 + dzdv0.*dzdv0;

    r = sqrt(x.^2 + y.^2 + z.^2);

    drdu = (repmat(a,[np,1]) + rv*b + ru*c)./r;
    drdv = (repmat(e,[np,1]) + ru*f + rv*g)./r;

    npts = npatches*length(ru);
    srcvals = zeros(12,npts);

    xuse = x./r;
    yuse = y./r;
    zuse = z./r;

    srcvals(1,:) = xuse(:)';
    srcvals(2,:) = yuse(:)';
    srcvals(3,:) = zuse(:)';

    dxduuse = (r.*dxdu - x.*drdu)./r./r;
    dyduuse = (r.*dydu - y.*drdu)./r./r;
    dzduuse = (r.*dzdu - z.*drdu)./r./r;

    srcvals(4,:) = dxduuse(:)';
    srcvals(5,:) = dyduuse(:)';
    srcvals(6,:) = dzduuse(:)';

    dxdvuse = (r.*dxdv - x.*drdv)./r./r;
    dydvuse = (r.*dydv - y.*drdv)./r./r;
    dzdvuse = (r.*dzdv - z.*drdv)./r./r;

    srcvals(7,:) = dxdvuse(:)';
    srcvals(8,:) = dydvuse(:)';
    srcvals(9,:) = dzdvuse(:)';

    dn = cross(srcvals(4:6,:),srcvals(7:9,:));
    srcvals(10:12,:) = dn./repmat(vecnorm(dn,2,1),[3,1]);

    obj = surfer(npatches,norder,srcvals);  


end
