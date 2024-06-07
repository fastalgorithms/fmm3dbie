function [fcoefs] = vals_to_coefs_surface_fun(obj, f)
%
%  Oversamples the given geometry given a current discretization
%  and a patchwise oversampling parameter, or a common oversampling
%  order for all patches
%
    npatches = obj.npatches;
    npts = obj.npts;
    [m1, m2] = size(f);
    if m1 == npts
       fuse = f.';
    elseif m2 == npts
      fuse = f;
    else
      error('vals_to_coefs_surface_fun: invalid size of input array');
    end
    [m1, m2] = size(fuse);

    fcoefs = zeros(m1, m2);

    
    ntmp = zeros(obj.npatches,2);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = obj.iptype;
    [ntmp_uni, ~, intmp] = unique(ntmp,'rows');
    nuni = size(ntmp_uni,1);
    
    umats = cell(nuni,1);
    for i=1:nuni
        
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);
        if(iptype == 1)
            rnodes = koorn.rv_nodes(norder);
            umats{i} = koorn.vals2coefs(norder,rnodes);  
        elseif(iptype == 11)
            rnodes = polytens.lege_nodes(norder);
            umats{i} = polytens.lege_vals2coefs(norder,rnodes);
        elseif(iptype == 12)
            rnodes = polytens.cheb_nodes(norder);
            umats{i} = polytens.cheb_vals2coefs(norder,rnodes);
        end
    end

    for i=1:npatches
        istart = obj.ixyzs(i);
        iend = obj.ixyzs(i+1)-1;
        iind = istart:iend;
        fcoefs(:, iind) = fuse(:, iind)*umats{intmp(i)}';
    end   
    fcoefs = reshape(fcoefs, size(f));
end
