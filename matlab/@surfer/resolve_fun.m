function obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, p, opts)
%  RESOLVE_FUN  Refine a surfer until the function f(obj) is resolved.
%
%  Syntax:
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax)
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, p)
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, p, opts)
%
%  Input arguments:
%    * obj:       surfer object to be resolved
%    * tol:       resolution tolerance
%    * f:         function handle evaluated as f(obj), returning
%                 (nfuns, npts) values at the surfer nodes; [] only balances
%    * nordcheck: number of leading orders treated as "resolved";
%                 the tail from nordcheck onwards is checked
%    * nrefmax:   maximum refinement passes in each phase; 0 returns unchanged
%    * p:         (optional, default 2) norm used on the tails
%    * opts:      (optional) options structure:
%        lvlr:         'a' (default) enables spatial level restriction;
%                      'n' disables it and preserves the previous behavior
%        lvlrfac:      maximum nearby patch-radius ratio (default 2.5)
%        lvlrnear:     projection interval multiplier (default 1)
%        npatchmax:    maximum total patch count (default Inf)
%
%  Example (disable balancing):
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, 2, struct('lvlr','n'));
%

if nargin < 6
    p = 2;
end

if nargin < 7
    opts = [];
end

lvlr = 'a';
lvlrfac = 2.5;
lvlrnear = 1;
npatchmax = Inf;
if isfield(opts, 'lvlr'); lvlr = opts.lvlr; end
if isfield(opts, 'lvlrfac'); lvlrfac = opts.lvlrfac; end
if isfield(opts, 'lvlrnear'); lvlrnear = opts.lvlrnear; end
if isfield(opts, 'npatchmax'); npatchmax = opts.npatchmax; end

lvlr = validatestring(lvlr, {'a', 'n'});
validateattributes(lvlrfac, {'numeric'}, {'scalar', 'real', 'finite', '>', 1});
validateattributes(lvlrnear, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(nrefmax, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
validateattributes(npatchmax, {'numeric'}, {'scalar', 'real', 'positive', 'nonnan'});
if isfinite(npatchmax) && npatchmax ~= fix(npatchmax)
    error('surfer:resolve_fun:InvalidPatchLimit', ...
        'opts.npatchmax must be a positive integer or Inf.');
end

if obj.npatches > npatchmax
    error('surfer:resolve_fun:PatchLimit', ...
        'Input surface exceeds opts.npatchmax = %g.', npatchmax);
end
if nrefmax == 0
    return
end
nordcheck = round(mean(nordcheck));

% Resolve the function first; an empty function only balances the mesh.
if ~isempty(f)
    srcvals_g = [];
    norders_g = [];
    iptypes_g = [];

    for jj = 1:nrefmax
        norder_iptype = [obj.norders(:), obj.iptype(:)];
        no_ip_uni = unique(norder_iptype, 'rows');
        nuni = size(no_ip_uni, 1);

        srcvals_b = [];
        norders_b = [];
        iptypes_b = [];

        npatchtotal = numel(iptypes_g) + obj.npatches;
        for i = 1:nuni
            ip0 = no_ip_uni(i,2);
            n0  = no_ip_uni(i,1);
            ipatchuse = find((obj.norders(:).' == n0) & (obj.iptype(:).' == ip0));

            if nuni == 1
                Stmp = obj;
            else
                idpts = [];
                for j = ipatchuse(:).'
                    idpts = [idpts, obj.ixyzs(j):(obj.ixyzs(j+1)-1)];
                end
                srcvalstmp = [obj.r(:,idpts); obj.du(:,idpts); obj.dv(:,idpts); obj.n(:,idpts)];

                Stmp = surfer(length(ipatchuse), n0, srcvalstmp, ip0);
            end

            [srcvals_g0, npatches_g, srcvals_b0, npatches_b] = ...
                resolve_fun_unif(Stmp, tol, f, nordcheck, p, npatchmax-npatchtotal);
            npatchtotal = npatchtotal + 3*(numel(ipatchuse)-npatches_g);

            srcvals_g = [srcvals_g, srcvals_g0];
            iptypes_g = [iptypes_g, ip0*ones(1,npatches_g)];
            norders_g = [norders_g, n0*ones(1,npatches_g)];

            srcvals_b = [srcvals_b, srcvals_b0];
            iptypes_b = [iptypes_b, ip0*ones(1,npatches_b)];
            norders_b = [norders_b, n0*ones(1,npatches_b)];
        end

        if isempty(iptypes_b)
            break
        end
        obj = surfer(length(iptypes_b), norders_b(:), srcvals_b, iptypes_b);
    end

    srcvals_g = [srcvals_g, srcvals_b];
    iptypes_g = [iptypes_g, iptypes_b];
    norders_g = [norders_g, norders_b];
    obj = surfer(length(iptypes_g), norders_g(:), srcvals_g, iptypes_g);
end

% Find nearby centers using the existing spatial search, not all patch pairs.
if strcmpi(lvlr, 'n') || obj.npatches == 0
    return
end

% Sample patch boundaries once per order/type, in addition to existing nodes.
[no_ip_uni, ~] = unique([obj.norders(:), obj.iptype(:)], 'rows');
edgepols = cell(size(no_ip_uni,1), 1);
for i = 1:size(no_ip_uni,1)
    norder = no_ip_uni(i,1);
    iptype = no_ip_uni(i,2);
    if iptype == 1
        t = linspace(0, 1, norder+1);
        uv = [t, 1-t, 0*t; 0*t, t, 1-t];
        edgepols{i} = koorn.pols(norder, uv);
    else
        t = linspace(-1, 1, norder+1);
        uv = [t, 1+0*t, -t, -1+0*t; -1+0*t, t, 1+0*t, -t];
        if iptype == 11
            edgepols{i} = polytens.lege.pols(norder, uv);
        else
            edgepols{i} = polytens.cheb.pols(norder, uv);
        end
    end
end

for iter = 0:nrefmax
    radii = obj.rads(:);
    if any(~isfinite(radii) | radii <= 0)
        error('surfer:resolve_fun:InvalidPatchRadius', ...
            'Spatial level restriction requires finite, positive patch radii.');
    end
    if max(radii) <= lvlrfac*min(radii)
        return
    end
    [~, igroup] = ismember([obj.norders(:), obj.iptype(:)], no_ip_uni, 'rows');
    points = cell(obj.npatches, 1);
    directions = cell(obj.npatches, 1);
    sample_radii = zeros(obj.npatches, 1);
    for i = 1:obj.npatches
        ids = obj.ixyzs(i):obj.ixyzs(i+1)-1;
        coefs = obj.srccoefs{i};
        edges = coefs(1:3,:)*edgepols{igroup(i)};
        points{i} = [obj.r(:,ids), edges] - obj.cms(:,i);
        sample_radii(i) = max(vecnorm(points{i}, 2, 1));

        uvcenter = [0;0];
        if obj.iptype(i) == 1; uvcenter = [1/3;1/3]; end
        [~, icenter] = min(sum((obj.uvs_targ(:,ids) - uvcenter).^2, 1));
        t1 = obj.du(:,ids(icenter));
        t1 = t1/norm(t1);
        t3 = cross(t1, obj.dv(:,ids(icenter)));
        t3 = t3/norm(t3);
        t2 = cross(t3, t1);

        verts = obj.end_pt_verts{i};
        edgevecs = verts(:,[2:end,1]) - verts;
        sideaxes = cross(repmat(t3, 1, size(verts,2)), edgevecs, 1);
        lengths = vecnorm(sideaxes, 2, 1);
        keep = lengths > 64*eps(max(1,radii(i)));
        sideaxes = sideaxes(:,keep)./lengths(keep);
        directions{i} = [t1, t2, t3, sideaxes];
    end

    % Radii find candidates; oriented projection bounds decide proximity.
    searchfac = max(1,lvlrnear)*max(sample_radii./radii)*(1 + 1/lvlrfac);
    near = getnear(obj, obj.cms, searchfac);
    isrc = near.col_ind;
    itarg = repelem((1:obj.npatches).', diff(near.row_ptr));
    bad = radii(isrc) > lvlrfac*radii(itarg);
    isrc = isrc(bad);
    itarg = itarg(bad);
    splitflag = false(obj.npatches, 1);
    geotol = 1e-10*max(sample_radii) + 64*eps(max(1,max(abs(obj.r(:)))));
    for k = 1:numel(isrc)
        i = isrc(k);
        j = itarg(k);
        if splitflag(i); continue; end
        dirs = [directions{i}, directions{j}];
        proji = dirs.'*points{i};
        projj = dirs.'*points{j};
        lo_i = min(proji, [], 2); hi_i = max(proji, [], 2);
        lo_j = min(projj, [], 2); hi_j = max(projj, [], 2);
        delta = dirs.'*(obj.cms(:,j) - obj.cms(:,i));
        centers = (lo_j + hi_j - lo_i - hi_i)/2 + delta;
        widths = (hi_i - lo_i + hi_j - lo_j)/2;
        splitflag(i) = all(abs(centers) <= lvlrnear*widths + geotol);
    end
    isplit = find(splitflag);
    if isempty(isplit)
        return
    end
    if iter == nrefmax
        error('surfer:resolve_fun:LevelRestrictionLimit', ...
            'Spatial level restriction exceeded nrefmax = %d.', nrefmax);
    end
    if obj.npatches + 3*numel(isplit) > npatchmax
        error('surfer:resolve_fun:PatchLimit', ...
            'Spatial level restriction would exceed opts.npatchmax = %g.', ...
            npatchmax);
    end
    obj = split_patches(obj, isplit);
end
end


function [srcvals_good, npatches_good, srcvals_bad, npatches_bad] = ...
        resolve_fun_unif(S, tol, f, nordcheck, p, npatchroom)
% Refine a surfer where every patch has the same type and order.

norders = S.norders(1);
iptype  = S.iptype(1);
if iptype == 1
    npoly = (norders+1)*(norders+2)/2;
else
    npoly = (norders+1)^2;
end

% check resolution
errps = surf_fun_error(S, f(S), p, nordcheck);
errps = vecnorm(errps, inf, 1);
ibad  = find((errps > tol) | isnan(errps));
if 3*numel(ibad) > npatchroom
    error('surfer:resolve_fun:PatchLimit', ...
        'Function refinement would exceed opts.npatchmax.');
end
if isempty(ibad)
    srcvals_good = [S.r; S.du; S.dv; S.n];
    npatches_good = S.npatches;
    srcvals_bad = zeros(12,0);
    npatches_bad = 0;
    return
end

r  = reshape(S.r,  3, npoly, []);
du = reshape(S.du, 3, npoly, []);
dv = reshape(S.dv, 3, npoly, []);
n  = reshape(S.n,  3, npoly, []);

% separate good and bad patches
rbad = reshape(r(:,:,ibad),  3, []);
ubad = reshape(du(:,:,ibad), 3, []);
vbad = reshape(dv(:,:,ibad), 3, []);
nbad = reshape(n(:,:,ibad),  3, []);

r(:,:,ibad)  = [];
du(:,:,ibad) = [];
dv(:,:,ibad) = [];
n(:,:,ibad)  = [];

srcvals_good  = [reshape(r,3,[]); reshape(du,3,[]); reshape(dv,3,[]); reshape(n,3,[])];
npatches_good = size(srcvals_good, 2) / npoly;

% split bad patches: interpolate position, tangents, and normal to sub-nodes
[~, xinterp, ~, dumat, dvmat] = patch_sub(norders, iptype);

nbad_patches = size(rbad, 2) / npoly;

% Each of rbad, ubad, vbad, nbad is (3, npoly*nbad_patches).
% Reshape to (3*npoly, nbad_patches), apply the interpolation matrix
% blockwise across the 3 components by stacking them, then reshape back.
% Equivalent to: for each xyz component, apply xinterp (or dumat/dvmat)
% to the (npoly, nbad_patches) slice.

rbad = reshape(rbad, 3, npoly, nbad_patches);
ubad = reshape(ubad, 3, npoly, nbad_patches);
vbad = reshape(vbad, 3, npoly, nbad_patches);
nbad = reshape(nbad, 3, npoly, nbad_patches);

% Apply interpolation: (4*npoly, npoly) * (npoly, nbad_patches) per component
rbad = reshape(permute(rbad, [2,3,1]), npoly, []);   % (npoly, nbad_patches*3)
ubad = reshape(permute(ubad, [2,3,1]), npoly, []);
vbad = reshape(permute(vbad, [2,3,1]), npoly, []);
nbad = reshape(permute(nbad, [2,3,1]), npoly, []);

rbad = reshape(xinterp * rbad, 4*npoly, nbad_patches, 3);
ubad = reshape(dumat   * ubad, 4*npoly, nbad_patches, 3);
vbad = reshape(dvmat   * vbad, 4*npoly, nbad_patches, 3);
nbad = reshape(xinterp * nbad, 4*npoly, nbad_patches, 3);

% Permute back to (3, 4*npoly*nbad_patches)
rbad = reshape(permute(rbad, [3,1,2]), 3, []);
ubad = reshape(permute(ubad, [3,1,2]), 3, []);
vbad = reshape(permute(vbad, [3,1,2]), 3, []);
nbad = reshape(permute(nbad, [3,1,2]), 3, []);

nbad = nbad ./ vecnorm(nbad);

srcvals_bad  = [rbad; ubad; vbad; nbad];
npatches_bad = size(srcvals_bad, 2) / npoly;

end
