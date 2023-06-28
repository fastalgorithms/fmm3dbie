function [obj] = sphere_quad(norder,nu,nref,iptype)
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
    if(nargin < 4)
        iptype = 11;
    end
    
    nuuse = nu*2^nref;
    nverts = 6*(nuuse+1)*(nuuse+1);
    verts = zeros(3,nverts);

    [vv,uu] = meshgrid(0:nuuse);
    uu = uu/(nuuse)*2-1;
    vv = vv/(nuuse)*2-1;
    uu = uu(:);
    vv = vv(:);
    nv0 = length(uu);
    
    np0 = nuuse*nuuse;
    quad0p = zeros(np0,3);
    quad0m = zeros(np0,3);
    iii = 0;
    for i=1:nuuse
        for j=1:nuuse
            iii = iii + 1;
            
            iv1 = (i-1)*(nuuse+1) + j;
            iv2 = (i-1)*(nuuse+1) + j+1;
            iv3 = i*(nuuse+1) + j;
            quad0p(iii,1) = iv1;
            quad0p(iii,2) = iv2;
            quad0p(iii,3) = iv3;
            
            quad0m(iii,1) = iv1;
            quad0m(iii,2) = iv3;
            quad0m(iii,3) = iv2;
            
        end
    end


    [nt,~] = size(quad0p);

    npatches = 6*nt;
    quads = zeros(npatches,3);

    % +z face
    verts(1,1:nv0) = uu/sqrt(3);
    verts(2,1:nv0) = vv/sqrt(3);
    verts(3,1:nv0) = 1/sqrt(3);
    quads(1:nt,:) = quad0p;

    % -z face
    istart = nv0+1;
    iind = (istart):(istart+nv0-1);
    itind = (nt+1):(2*nt);
    verts(1,iind) = uu/sqrt(3);
    verts(2,iind) = vv/sqrt(3);
    verts(3,iind) = -1/sqrt(3);
    quads(itind,:) = quad0m + istart-1;

    % + x face
    istart = 2*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (2*nt+1):(3*nt);
    verts(1,iind) = 1/sqrt(3);
    verts(2,iind) = uu/sqrt(3);
    verts(3,iind) = vv/sqrt(3);
    quads(itind,:) = quad0p + istart-1;


    % -x face
    istart = 3*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (3*nt+1):(4*nt);
    verts(1,iind) = -1/sqrt(3);
    verts(2,iind) = uu/sqrt(3);
    verts(3,iind) = vv/sqrt(3);
    quads(itind,:) = quad0m + istart-1;


    % + y face
    istart = 4*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (4*nt+1):(5*nt);
    verts(1,iind) = vv/sqrt(3);
    verts(2,iind) = 1/sqrt(3);
    verts(3,iind) = uu/sqrt(3);
    quads(itind,:) = quad0p + istart-1;


    % -y face
    istart = 5*nv0 + 1;
    iind = (istart):(istart+nv0-1);
    itind = (5*nt+1):(6*nt);
    verts(1,iind) = vv/sqrt(3);
    verts(2,iind) = -1/sqrt(3);
    verts(3,iind) = uu/sqrt(3);
    quads(itind,:) = quad0m + istart-1;

    % collect all list of vertices
    v1s = verts(1:3,quads(:,1));
    v2s = verts(1:3,quads(:,2));
    v3s = verts(1:3,quads(:,3));

    % Generate rv nodes on a single patch
    if(iptype == 11)
        rnodes = polytens.lege_nodes(norder);
    else
        rnodes = polytens.cheb_nodes(norder);
    end
 
    % compute the projection along with all derivatives
    ru  = rnodes(1,:)';
    rv = rnodes(2,:)';
    np = length(ru);

    x = repmat(v1s(1,:),[np,1]) + (1+ru)/2*(v2s(1,:)-v1s(1,:)) + ...
         (1+rv)/2*(v3s(1,:)-v1s(1,:));

    y = repmat(v1s(2,:),[np,1]) + (1+ru)/2*(v2s(2,:)-v1s(2,:)) + ...
         (1+rv)/2*(v3s(2,:)-v1s(2,:));

    z = repmat(v1s(3,:),[np,1]) + (1+ru)/2*(v2s(3,:)-v1s(3,:)) + ...
         (1+rv)/2*(v3s(3,:)-v1s(3,:));


    dxdu0 = (v2s(1,:)-v1s(1,:))/2;
    dydu0 = (v2s(2,:)-v1s(2,:))/2;
    dzdu0 = (v2s(3,:)-v1s(3,:))/2;

    dxdv0 = (v3s(1,:)-v1s(1,:))/2;
    dydv0 = (v3s(2,:)-v1s(2,:))/2;
    dzdv0 = (v3s(3,:)-v1s(3,:))/2;


    dxdu = repmat(dxdu0,[np,1]);
    dydu = repmat(dydu0,[np,1]);
    dzdu = repmat(dzdu0,[np,1]);

    dxdv = repmat(dxdv0,[np,1]);
    dydv = repmat(dydv0,[np,1]);
    dzdv = repmat(dzdv0,[np,1]);

    
    r = sqrt(x.^2 + y.^2 + z.^2);
    drdu = (x.*dxdu + y.*dydu + z.*dzdu)./r;
    drdv = (x.*dxdv + y.*dydv + z.*dzdv)./r;

    
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

    obj = surfer(npatches,norder,srcvals,iptype);  


end
