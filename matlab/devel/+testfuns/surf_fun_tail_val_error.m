function [errps] = surf_fun_tail_val_error(obj, fun, lptype, norderhead, ...
    functional_integrand)
%  SURF_FUN_TAIL_VAL_ERROR  Per-patch L^p norm of the high-degree tail,
%  measured in VALUE space rather than coefficient space.
%
%  Same tail selection as SURF_FUN_ERROR, but the tail coefficients are
%  mapped back to the nodes and the norm is taken over values. The result
%  is a genuine L^p(P_j) norm.  No (|P_j|/|S|)^(1/p) scaling is applied.
%
%    p = 0 or Inf : max |tail| over the patch nodes
%    p >= 1       : ( sum |tail|^p .* wts )^(1/p)

    if nargin < 3 || isempty(lptype), lptype = 0; end
    if nargin < 4, norderhead = []; end
    if nargin < 5, functional_integrand = []; end

    npatches = obj.npatches;

    fcoefs = obj.vals2coefs(fun);
    if size(fcoefs,2) == obj.npts
        fcoefs = fcoefs.';
    end
    nfuns = size(fcoefs,2);

    errps = zeros(npatches, nfuns);

    ntmp = [obj.norders(:), obj.iptype(:)];
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');

 %   patcharea = full(sparse(obj.patch_id(:),1,obj.wts(:),npatches,1));

    for i = 1:size(ntmp_uni,1)
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);
        inds   = find(intmp == i);
        indv   = obj.ixyzs(inds);

        if isempty(norderhead)
            head = norder - 1;
        else
            head = norderhead;
        end

        if iptype == 1
            npols = (norder+1)*(norder+2)/2;
            ihead = 1:((head+1)*(head+2)/2);
            bmat  = koorn.coefs2vals(norder, koorn.rv_nodes(norder));
        elseif iptype == 11 || iptype == 12
            npols = (norder+1)^2;
            [N1,N2] = ndgrid(0:norder, 0:norder);
            ihead   = find(N1(:) + N2(:) <= head);
            if iptype == 11
                bmat = polytens.lege.coefs2vals(norder, polytens.lege.nodes(norder));
            else
                bmat = polytens.cheb.coefs2vals(norder, polytens.cheb.nodes(norder));
            end
        else
            error('SURF_FUN_TAIL_VAL_ERROR: unsupported iptype %d', iptype);
        end

        [I,V] = meshgrid(indv, 1:npols);
        itot  = I(:) + V(:) - 1;

        cf = reshape(fcoefs(itot,:), npols, numel(inds), nfuns);
        cf(ihead,:,:) = 0;                       % keep tail only

        v = bmat * reshape(cf, npols, []);       % tail values at nodes
        v = reshape(v, npols, numel(inds), nfuns);

        if ~isempty(functional_integrand)
            x = reshape(obj.r(1,itot),npols,numel(inds));
            v = functional_integrand(v,x);
        end
        if lptype == 0 || isinf(lptype)
            e = reshape(max(abs(v),[],1), numel(inds), nfuns);
        else
            w = reshape(obj.wts(itot), npols, numel(inds));
            e = reshape(sum(abs(v).^lptype .* w,1), ...
            numel(inds),nfuns).^(1/lptype);
    %        e = e ./ patcharea(inds).^(1/lptype);
        end

        errps(inds,:) = e;
    end

    errps = errps.';
end
