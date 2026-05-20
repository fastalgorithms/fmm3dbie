function [S, chnkr] = triangulate_chunker_interior(chnkr, nord, Hmax, curvtol)
% TRIANGULATE_CHUNKER  Build a high-order surfer mesh for the region
%  enclosed by a chunker boundary curve.
%
% [S, chnkr] = triangulate_chunker_interior(chnkr, nord, Hmax)
% [S, chnkr] = triangulate_chunker_interior(chnkr, nord, Hmax, curvtol)
%
% Uses the MATLAB PDE Toolbox to generate a triangular mesh of the domain
% bounded by chnkr, then fits high-order Koornwinder polynomial patches on
% each triangle to produce a surfer object.
%
% Boundary patches are rotated so that the boundary corresponds to the edge
% v = 0 (for edge-on patches), or the vertex (0,1) (for vertex-on patches).
%
% Inputs:
%   chnkr   - chunker object describing the boundary curve
%   nord    - polynomial order for the surfer patches
%   Hmax    - maximum triangle edge length for the PDE Toolbox mesh
%   curvtol - (optional) tolerance for adaptive chunk refinement based on
%             curvature. If omitted, no refinement is performed.
%
% Outputs:
%   S       - surfer object representing the interior domain
%   chnkr   - chunker object (possibly refined if curvtol was supplied)
%
% Requirements:
%   MATLAB PDE Toolbox (uses generateMesh)
%
% See also SURFER

%% Check that PDE Toolbox is installed

addons = matlab.addons.installedAddons();
if ~any(strcmp(addons.Name, 'Partial Differential Equation Toolbox'))
    error('triangulate_chunker_interior:missingToolbox', ...
        ['The Partial Differential Equation Toolbox is required. ' ...
         'Please install it via the Add-On Explorer.']);
end

%% Adaptively refine chunks where curvature is not resolved

if nargin >= 4
    k = 16;
    [xl, ~, vals2coefs, ~] = lege.exps(k);

    kappa = signed_curvature(chnkr);

    coefs = vals2coefs * reshape(kappa, k, []);
    tails = max(abs(coefs(end-2:end, :)));
    iref  = find(tails > curvtol);

    opts = []; opts.splitchunks = iref;
    chnkr = refine(chnkr, opts);
    chnkr = sort(chnkr);
else
    k = 16;
    [xl, ~, vals2coefs, ~] = lege.exps(k);
end

%% Build high-order polygon from chunker nodes

rt = chnkr.r;
x  = reshape(rt(1,:), k, []);
y  = reshape(rt(2,:), k, []);

ps = lege.pols(-1, k-1);
vx = ps.' * vals2coefs * x;
vy = ps.' * vals2coefs * y;

[rs, ~] = chunkends(chnkr);
rv = squeeze(rs(:,1,:));

%% Generate triangular mesh

pgon = polyshape(vx, vy);
tr   = triangulation(pgon);
gm   = fegeometry(tr);
gm   = addVertex(gm, 'Coordinates', rv.');
gm   = generateMesh(gm, 'GeometricOrder', 'linear', 'Hmax', Hmax);

msh = gm.Mesh;
ee  = msh.Elements;
nn  = msh.Nodes;

%% Identify boundary edges and map to chunker segments

tr_tmp = triangulation(ee.', nn.');
ff = freeBoundary(tr_tmp).';

node_vert = findNodes(msh, 'nearest', rv);
nv = numel(node_vert);

n1   = node_vert(1);
ind1 = find(ff(1,:) == n1);
ff   = circshift(ff, [0, -ind1+1]);

[I, J]       = meshgrid(node_vert, ff(1,:));
[inode, iff] = find(I.' == J.');
[~, iord]    = sort(inode);
iff = iff(iord);

node_vert = [node_vert, node_vert(1)];

xs     = zeros(k, size(ff,2));
ys     = zeros(k, size(ff,2));
iverts = zeros(2, size(ff,2));
strts  = zeros(size(ff,2), 1);
finis  = zeros(size(ff,2), 1);
iis    = zeros(size(ff,2), 1);

icount = 1;
for ii = 1:nv
    i0 = iff(ii);
    if ii ~= nv
        i1   = iff(ii+1);
        inds = ff(1, i0:(i1-1));
    else
        i1   = iff(1);
        inds = ff(1, i0:end);
    end

    v1 = nn(:, ff(1,i0));
    v2 = nn(:, ff(1,i1));

    istrts = nn(:, inds);
    iends  = nn(:, [inds(2:end), ff(1,i1)]);
    sdist1 = vecnorm(istrts - v1);
    sdist2 = vecnorm(istrts - v2);
    spar   = sdist1 ./ (sdist1 + sdist2);

    fdist1 = vecnorm(iends - v1);
    fdist2 = vecnorm(iends - v2);
    fpar   = fdist1 ./ (fdist1 + fdist2);

    xloc = x(:,ii);
    yloc = y(:,ii);

    % Snap interior boundary nodes onto the chunker curve
    for iii = 2:numel(inds)
        tt = spar(iii)*2 - 1;
        nn(:, inds(iii)) = [(lege.pols(tt,k-1)).' * vals2coefs * xloc; ...
                            (lege.pols(tt,k-1)).' * vals2coefs * yloc];
    end

    inds = [inds, ff(1,i1)];
    for iii = 1:numel(inds)-1
        strt = 2*spar(iii) - 1;
        fini = 2*fpar(iii) - 1;
        ts   = (xl+1)/2*(fini-strt) + strt;
        strts(icount)    = strt;
        finis(icount)    = fini;
        iis(icount)      = ii;
        xs(:,icount)     = (lege.pols(ts,k-1)).' * vals2coefs * xloc;
        ys(:,icount)     = (lege.pols(ts,k-1)).' * vals2coefs * yloc;
        iverts(1,icount) = inds(iii);
        iverts(2,icount) = inds(iii+1);
        icount = icount + 1;
    end
end

%% Build high-order surfer patches

norder = 20;
uvs   = koorn.rv_nodes(norder);
nuvs  = size(uvs, 2);
nelem = size(ee, 2);

xx = zeros(nuvs, nelem);
yy = zeros(nuvs, nelem);

bnd_nodes = unique(iverts(:));

for ii = 1:nelem
    nds = ee(:,ii);
    n1 = nds(1); n2 = nds(2); n3 = nds(3);

    [~, icl] = find((iverts == n1) + (iverts == n2) + (iverts == n3));
    [i1, ~]  = find((icl == icl.') - eye(numel(icl)));

    if numel(i1) > 0
        % This triangle has an edge on the boundary.
        % Use Gordon-Hall transfinite interpolation so that the curved
        % boundary edge lies along v = 0 of the reference triangle.
        enum = icl(i1(1));
        inot = setdiff([n1,n2,n3], iverts(:,enum));
        v1   = nn(:, inot);

        xt = xs(:, enum);
        yt = ys(:, enum);

        gam0 = [(lege.pols(-1,k-1)).' * vals2coefs * xt; ...
                (lege.pols(-1,k-1)).' * vals2coefs * yt];
        gamL = [(lege.pols( 1,k-1)).' * vals2coefs * xt; ...
                (lege.pols( 1,k-1)).' * vals2coefs * yt];
        gamx = [(lege.pols((1-uvs(1,:).')*2-1, k-1)).' * vals2coefs * xt, ...
                (lege.pols((1-uvs(1,:).')*2-1, k-1)).' * vals2coefs * yt].';

        ntrm = 1 - uvs(1,:) - uvs(2,:);
        xi   = uvs(1,:);
        eta  = uvs(2,:);
        dtrm = 1 - uvs(1,:);
        tran = ntrm.*gamL + xi.*gam0 + eta.*v1 + ...
               ntrm./dtrm.*(gamx - dtrm.*gamL - xi.*gam0);
        xx(:,ii) = tran(1,:);
        yy(:,ii) = tran(2,:);
    else
        % Linear patch.  Permute vertices so that any boundary vertex
        % (a node that appears in iverts) is at (u,v) = (0,1) in the
        % reference triangle
        [~, ibnd] = intersect([n1, n2, n3], bnd_nodes);
        if ~isempty(ibnd)
            % Rotate the triple so the boundary vertex lands at index 2.
            rot = mod(ibnd(1) - 2, 3);
            tri = circshift([n1, n2, n3], -rot);
            n1 = tri(1); n2 = tri(2); n3 = tri(3);
        end
        pa = nn(:,n1); pb = nn(:,n2); pc = nn(:,n3);
        uvec = pc - pa;
        vvec = pb - pa;
        xx(:,ii) = pa(1) + uvec(1)*uvs(1,:) + vvec(1)*uvs(2,:);
        yy(:,ii) = pa(2) + uvec(2)*uvs(1,:) + vvec(2)*uvs(2,:);
    end
end

%% Fit polynomial coefficients and build surfer

amat = koorn.vals2coefs(norder, uvs);
cfx  = amat * xx;
cfy  = amat * yy;

uv = koorn.rv_nodes(nord);
wt = koorn.rv_weights(nord);
[pols, dersu, dersv] = koorn.ders(norder, uv);

srcx   = pols.' * cfx;
srcy   = pols.' * cfy;
srcdxu = dersu.' * cfx;
srcdxv = dersv.' * cfx;
srcdyu = dersu.' * cfy;
srcdyv = dersv.' * cfy;

srcvals = zeros(12, numel(wt)*nelem);
srcvals(1,:)  = srcx(:);
srcvals(2,:)  = srcy(:);
srcvals(4,:)  = srcdxu(:);
srcvals(5,:)  = srcdyu(:);
srcvals(7,:)  = srcdxv(:);
srcvals(8,:)  = srcdyv(:);
srcvals(12,:) = -1;

S    = surfer(nelem, nord, srcvals, 1);
errs = vecnorm(surf_fun_error(S, srcvals), inf);
fprintf('triangulate_chunker_interior: max resolution error = %.2e\n', max(errs))

end
