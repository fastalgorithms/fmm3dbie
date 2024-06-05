function [errps] = surf_fun_error(obj,fun,p)
% SURF_FUN_ERROR  estimate pointwise approximation error of func on surface.
%
% [errps] = surf_fun_error(obj,fun,p)
%  
%       This routine computes an estimate of the approximation error
%       in a collection of functions defined on a surfer object.
%
%       Inputs:
%       obj     --- the surfer object
%       fun     --- (nfuns*npts) array containing function values
%                   at the points in obj
%       p       --- optional: the norm to use to measure the coeffs
%                   (default is infinity)
%       
%       iptype = 1 : Koornwinder (look at last norder)
%       iptype = 11: Gauss-Legendere (look at upper triangle > norder)
%       iptype = 12: Chebyshev (look at upper triangle > norder)

    npatches = obj.npatches;

%%%%%%%%%%%%%%

    ntmp = zeros(npatches,2);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = obj.iptype;

    [fcoefs] = vals_to_coefs_surface_fun(obj, fun);
    fcoefs = fcoefs.';
    nfuns = size(fcoefs,2);

    errps = zeros(npatches,nfuns);
    parea = zeros(npatches,1);
    [ntmp_uni,~,intmp] = unique(ntmp,'rows');
    nuni = size(ntmp_uni,1);
    umats = cell(nuni,1);
    for i=1:nuni
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);

        inds = find(intmp == i);
        indv = obj.ixyzs(inds);

        if (iptype == 1)
            norderhead = norder-1;
            i1    = (norderhead+1)*(norderhead+2)/2;
            npols = (norder+1)*(norder+2)/2;
            ipols = (i1+1):npols;
        end
        if (iptype == 11)
            nords1  = 0:norder;
            nords2  = 0:norder;
            [N1,N2] = ndgrid(nords1,nords2);
            degs = N1(:) + N2(:);
            ipols = find(degs > norder);
        end
        if (iptype == 12)
            nords1  = 0:norder;
            nords2  = 0:norder;
            [N1,N2] = ndgrid(nords1,nords2);
            degs = N1(:) + N2(:);
            ipols = find(degs > norder);
        end

        [I,V] = meshgrid(indv,ipols);
        itot = I(:)+V(:)-1;
        cftmp = reshape(fcoefs(itot,:),numel(ipols),[]);
        if (nargin >2)
            nrms = vecnorm(cftmp,p,1);
        else
            nrms = vecnorm(cftmp,Inf,1);
        end
        errs = reshape(nrms,[],nfuns);
        errps(inds,:) = errs;

        areas = sum([obj.weights{inds}],1);
        parea(inds) = areas;
    end
    
    if (nargin >2)
        a     = area(obj);
        errps = (parea/a).^(1/p).*errps;
    end

    errps = errps.';


end