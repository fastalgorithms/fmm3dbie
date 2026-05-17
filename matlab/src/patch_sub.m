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
%    * dumat: (4*npols, npols) interpolation matrix for the
%        u-derivative, mapping original node values to u-derivatives
%        at the sub-patch nodes
%    * dvmat: (4*npols, npols) interpolation matrix for the
%        v-derivative, mapping original node values to v-derivatives
%        at the sub-patch nodes
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

        [~,deru,derv] = koorn.ders(norder,uvs);
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

        [~,deru,derv] = polytens.lege.ders(norder,uvs);
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

        [~,deru,derv] = polytens.cheb.ders(norder,uvs);
    end
    pols = pols.';
dumat = kron(eye(4),deru.'*amat);
dvmat = kron(eye(4),derv.'*amat);    

end