function [uvst,xinterp,pols,dumat,dvmat] = patch_sub(norder, iptype)
%
%  PATCH_SUB
%
%  Syntax:
%    [uvst,xinterp,pols,dumat,dvmat] = patch_sub(norder, iptype)
%
%  This subroutine computes the nodes and interpolation matrices
%  needed to split a single patch into four equal sub-patches.
%  The original patch nodes are mapped into each quadrant and the
%  returned matrices interpolate function values and derivatives
%  from the original nodes to the combined set of sub-patch nodes.
%
%  Input arguments:
%    * norder: integer, polynomial order of the patch
%    * iptype: integer, patch type
%        iptype == 1,  Koornwinder-based triangular patch
%        iptype == 11, Legendre-based quadrilateral patch
%        iptype == 12, Chebyshev-based quadrilateral patch
%
%  Output arguments:
%    * uvst: (2, 4*npols) parameter coordinates of the sub-patch
%        nodes on the original patch, where npols is the number of
%        nodes per patch at order norder
%    * xinterp: (4*npols, npols) interpolation matrix mapping
%        function values at the original nodes to values at uvst
%    * pols: (4*npols, npols) basis polynomial values at uvst
%    * dumat: (4*npols, npols) interpolation matrix mapping parent
%        dr/du tangent vectors at the original nodes to dr/dt tangent
%        vectors at the sub-patch nodes, where t is the sub-patch's
%        own local parameter. Accounts for the chain rule ds/dt = +-1/2
%        of the affine sub-patch mapping (sign -1/2 for iptype==1 quad 2,
%        +1/2 for all other quads).
%    * dvmat: (4*npols, npols) same as dumat but for the v component
%

    if iptype == 1
        [uvs] = koorn.rv_nodes(norder);
        uvs1 = uvs/2;
        uvs2 = -uvs/2 + [0.5;0.5];        
        uvs3 = uvs/2+[0.5;0];
        uvs4 = uvs/2+[0;0.5];
        [uvst] = [uvs1,uvs2,uvs3,uvs4];
        amat = koorn.vals2coefs(norder,uvs);
        pols = koorn.pols(norder,uvst);
        xinterp = pols.'*amat;
    elseif iptype == 11
        [uvs] = polytens.lege.nodes(norder);
        uvs1 = uvs/2+[-0.5;-0.5];
        uvs2 = uvs/2+[0.5;-0.5];
        uvs3 = uvs/2+[-0.5;0.5];
        uvs4 = uvs/2+[0.5;0.5];
        [uvst] = [uvs1,uvs2,uvs3,uvs4];
        amat = polytens.lege.vals2coefs(norder,uvs);
        pols = polytens.lege.pols(norder,uvst);
        xinterp = pols.'*amat;
    elseif iptype == 12
        [uvs] = polytens.cheb.nodes(norder);
        uvs1 = uvs/2+[-0.5;-0.5];
        uvs2 = uvs/2+[0.5;-0.5];
        uvs3 = uvs/2+[-0.5;0.5];
        uvs4 = uvs/2 + [0.5;0.5];
        [uvst] = [uvs1,uvs2,uvs3,uvs4];
        amat = polytens.cheb.vals2coefs(norder,uvs);
        pols = polytens.cheb.pols(norder,uvst);
        xinterp = pols.'*amat;
    end
    pols = pols.';
    npols = size(xinterp, 2);
    if iptype == 1
        % Quad 2 uses uvs2 = -uvs/2 + [0.5;0.5], so ds/dt = -1/2 for
        % both components; quads 1,3,4 use uvs/2+offset, so ds/dt = +1/2.
        signs = [ones(npols,1); -ones(npols,1); ones(2*npols,1)];
        dumat = (signs .* xinterp) / 2;
        dvmat = dumat;
    else
        % iptype 11,12: all quads use uvs/2+offset, so ds/dt = +1/2.
        dumat = xinterp / 2;
        dvmat = xinterp / 2;
    end

end