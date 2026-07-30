function obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, p)
%  RESOLVE_FUN  Refine a surfer until the function f(obj) is resolved.
%
%  Syntax:
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax)
%    obj = resolve_fun(obj, tol, f, nordcheck, nrefmax, p)
%
%  Input arguments:
%    * obj:       surfer object to be resolved
%    * tol:       resolution tolerance
%    * f:         function handle evaluated as f(obj), returning
%                 (nfuns, npts) values at the surfer nodes
%    * nordcheck: number of leading orders treated as "resolved";
%                 the tail from nordcheck onwards is checked
%    * nrefmax:   maximum number of refinement passes
%    * p:         (optional, default 2) norm used on the tails
%

if nargin < 6
    p = 2;
end

nordcheck = round(mean(nordcheck));
if nrefmax == 0
    return
end

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

    for i = 1:nuni
        ip0 = no_ip_uni(i,2);
        n0  = no_ip_uni(i,1);
        ipatchuse = find((obj.norders(:).' == n0) & (obj.iptype(:).' == ip0));

        idpts = [];
        for j = ipatchuse(:).'
            idpts = [idpts, obj.ixyzs(j):(obj.ixyzs(j+1)-1)];
        end
        srcvalstmp = [obj.r(:,idpts); obj.du(:,idpts); obj.dv(:,idpts); obj.n(:,idpts)];

        Stmp = surfer(length(ipatchuse), n0, srcvalstmp, ip0);

        [srcvals_g0, npatches_g, srcvals_b0, npatches_b] = resolve_fun_unif(Stmp, tol, f, nordcheck, p);

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


function [srcvals_good, npatches_good, srcvals_bad, npatches_bad] = ...
        resolve_fun_unif(S, tol, f, nordcheck, p)
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
