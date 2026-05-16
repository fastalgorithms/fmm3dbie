function [Sout,islice,Sshift] = vertex_surfer(S,vert,idvertpatch,tol)
%
%  VERTEX_SURFER
%
%  Syntax:
%    [Sout,islice] = vertex_surfer(S,vert)
%    [Sout,islice] = vertex_surfer(S,vert,idvertpatch)
%    [Sout,islice] = vertex_surfer(S,vert,idvertpatch,tol)
%    [Sout,islice,Sshift] = vertex_surfer(S,vert)
%
%  This subroutine extracts all patches of S that share the vertex
%  vert and returns them as a new surfer object Sout with positions
%  shifted so that vert is at the origin. The shifted positions on
%  vertex patches are recomputed accurately using an integration-based
%  approach rather than simple translation.
%
%  Input arguments:
%    * S: surfer object describing the source surface
%    * vert: (3,1) coordinates of the target vertex
%    * idvertpatch: (1,m) integer array of patch indices known to
%        share vert (optional). If not provided, the patches are
%        identified automatically by comparing S.end_pt_verts to vert
%        within tol.
%    * tol: tolerance for identifying vertex patches (optional,
%        default 1e-5)
%
%  Output arguments:
%    * Sout: surfer object containing only the m vertex patches,
%        with positions shifted so that vert maps to the origin
%    * islice: indices into S point arrays corresponding to the
%        vertex patches
%    * Sshift: (optional) surfer object with the same connectivity
%        as S but with all positions shifted by -vert, and with the
%        vertex patch positions replaced by the accurately recomputed
%        shifted values
%
    assert(size(vert,2)==1)

    if nargin < 4
        tol = 1e-5;
    end

    if nargin < 3
        idvertpatch = [];
        for i = 1:S.npatches
            if min(vecnorm(S.end_pt_verts{i} - vert)) < tol
                idvertpatch = [idvertpatch,i];
            end
        end
    end

    islice = [];
    srcvals = [];
    for i = idvertpatch
        iinds = S.ixyzs(i):S.ixyzs(i+1)-1;
        islice = [islice,iinds];
        iptype = S.iptype(i);
        nleg = 16;
        verts_i = S.end_pt_verts{i};
        dists = vecnorm(verts_i - vert);

        if iptype == 1
            norder = S.norders(i);
            if dists(1) < tol, pt = [0;0]; end
            if dists(2) < tol, pt = [0;1]; end
            if dists(3) < tol, pt = [1;0]; end
            uvs = koorn.rv_nodes(norder);
            [xintmatu,xintmatv,~] = koorn.koorn_int_mat(uvs,norder,nleg,pt);
        elseif iptype == 11
            norder = S.norders(i);
            epts = [-1,1,1,-1;-1,-1,1,1];
            [~,ivert] = min(dists);
            pt = epts(:,ivert);
            uvs = polytens.lege_nodes(norder);
            [xintmatu,xintmatv,~] = polytens.lege.int_mat(uvs,norder,nleg,pt);
        elseif iptype == 12
            norder = S.norders(i);
            epts = [-1,1,1,-1;-1,-1,1,1];
            [~,ivert] = min(dists);
            pt = epts(:,ivert);
            uvs = polytens.cheb_nodes(norder);
            [xintmatu,xintmatv,~] = polytens.cheb.int_mat(uvs,norder,nleg,pt);
        end
        du = S.du(:,iinds);
        dv = S.dv(:,iinds);

        r = du*xintmatu.' + dv*xintmatv.';
        n = cross(du,dv);
        n = n./vecnorm(n);

        srcvals = [srcvals, [r;du;dv;n]];
    end
    
    if isempty(idvertpatch)
        Sout = [];
    else
        Sout = surfer(numel(idvertpatch), S.norders(idvertpatch), srcvals, S.iptype(idvertpatch));
    end

    if nargout > 2
        srcvalsshift = [S.r-vert;S.du;S.dv;S.n];
        srcvalsshift(1:3,islice) = srcvals(1:3,:);
        Sshift = surfer(S.npatches,S.norders,srcvalsshift,S.iptype);

    end
end