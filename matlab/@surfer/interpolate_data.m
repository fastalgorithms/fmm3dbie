function dens_int = interpolate_data(obj, dens, ipatch_ids, uvs_targ)
% INTERPOLATE_DATA interpolates a density on surfer object 
%
% dens_int = interpolate_data(obj,dens,ipatch_id,uvs_targ) where dens_int is a 
% vector with values of dens interpolated to the given ipatch_ids and uvs_targ
%
% Input: obj    - a surfer object
%        dens   - a density defined on obj
%        ipatch_id - patch_ids of targs
%        uvs_targ - uv_coordinates
%
% Output: values of the interpolated density

    if size(dens,1) == obj.npts
       fl = false;
    elseif size(dens,2) == obj.npts
       dens = dens.';
       fl = true;
    else
        error("INTERPOLATE_DATA:error in size");
    end

    upids = unique(ipatch_ids);
    npatch = length(upids);
    dens_int = zeros([length(ipatch_ids), size(dens,2)]);

for ii = 1:npatch

    patch_id = upids(ii);
    norder = obj.norders(patch_id);
    iptype = obj.iptype(patch_id);

    vals = dens(obj.patch_id == patch_id,:);
    uvs = obj.uvs_targ(:,obj.patch_id == patch_id);
    uvs_int = uvs_targ(:,ipatch_ids == patch_id);

    if iptype == 1
     %   Rokhlin-Vioreanu nodes

     v2c = koorn.vals2coefs(norder, uvs);
     c2v = koorn.coefs2vals(norder, uvs_int);
     dens_int(ipatch_ids == patch_id,:) = c2v*(v2c*vals);

    elseif iptype == 11
     %   Gauss-Legendre nodes

    v2c = polytens.lege.vals2coefs(norder,uvs);
    c2v = polytens.lege.coefs2vals(norder,uvs_int);
    dens_int(ipatch_ids == patch_id,:) = c2v*(v2c*vals);

    elseif iptype == 12
     %   Chebyshev nodes

    v2c = polytens.cheb.vals2coefs(norder,uvs);
    c2v = polytens.cheb.coefs2vals(norder,uvs_int);
    dens_int(ipatch_ids == patch_id,:) = c2v*(v2c*vals);

    end
    
end

if fl
    dens_int = dens_int .';
end

end