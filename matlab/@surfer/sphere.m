function [obj] = sphere(norder, nu, nref, iptype)
% Generate the triangulation/quadrilaterization of a sphere
% through projecting out a cube discretized with
% nu \times nu patches in each direction, refined
% by nref levels of refinement.
%
%
% Input
%  - norder: order of discretization on each patch
%  - nu: number of patches in u,v direction on each face of cube
%  - nref: additional refinement of patches in each direction
%  - iptype: type of surface
%            * iptype = 1, triangulation with RV nodes
%            * iptype = 11, quadrilateralization with GL nodes
%            * iptype = 12, quadrilaterization with Cheb nodes
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
        iptype = 1;
    end

    if iptype == 1
        obj = surfer.sphere_tri(norder, nu, nref);
    else if iptype == 11 || iptype == 12
        obj = surfer.sphere_quad(norder, nu, nref, iptype);
    else
        error('SURFER:sphere: Invalid patch type, patch type must be 1, 11, or 12\n')
    end
    
end
