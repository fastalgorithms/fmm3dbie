function [errps] = surf_fun_error(obj, fun, p, norderhead)
%  SURF_FUN_ERROR  Estimate pointwise approximation error of a function
%  on a surfer object.
%
%  Syntax:
%    errps = surf_fun_error(obj, fun)
%    errps = surf_fun_error(obj, fun, p)
%    errps = surf_fun_error(obj, fun, p, norderhead)
%
%  Input arguments:
%    * obj:        surfer object
%    * fun:        (nfuns, npts) array of function values at the nodes of obj
%    * p:          (optional) norm used to measure the tail coefficients
%                  (default: Inf)
%    * norderhead: (optional) if provided, only the basis coefficients
%                  corresponding to degrees strictly above norderhead are
%                  treated as the "tail". When omitted the full high-degree
%                  tail (degrees > norder-1 for iptype 1, > norder for
%                  iptype 11/12) is used, matching the default behaviour.
%
%  Output arguments:
%    * errps: (nfuns, npatches) estimated per-patch approximation errors

    npatches = obj.npatches;

    ntmp = zeros(npatches, 2);
    ntmp(:,1) = obj.norders;
    ntmp(:,2) = obj.iptype;

    [fcoefs] = obj.vals2coefs(fun);

    [~, m2] = size(fcoefs);
    if m2 == obj.npts
        fcoefs = fcoefs.';
    end
    nfuns = size(fcoefs, 2);

    errps = zeros(npatches, nfuns);
    parea = zeros(npatches, 1);
    [ntmp_uni, ~, intmp] = unique(ntmp, 'rows');
    nuni = size(ntmp_uni, 1);

    for i = 1:nuni
        norder = ntmp_uni(i,1);
        iptype = ntmp_uni(i,2);

        inds = find(intmp == i);
        indv = obj.ixyzs(inds);

        % Determine which polynomial indices form the "tail".
        % If norderhead is supplied, use it as the cutoff; otherwise
        % default to checking only the very last order (norder-1 for
        % iptype 1, or total degree > norder for iptype 11/12).
        if nargin < 4
            head = norder - 1;
        else
            head = norderhead;
        end

        if iptype == 1
            i1    = (head+1)*(head+2)/2;
            npols = (norder+1)*(norder+2)/2;
            ipols = (i1+1):npols;
        elseif iptype == 11 || iptype == 12
            nords1 = 0:norder;
            nords2 = 0:norder;
            [N1,N2] = ndgrid(nords1, nords2);
            degs  = N1(:) + N2(:);
            ipols = find(degs > head);
        end

        [I, V] = meshgrid(indv, ipols);
        itot   = I(:) + V(:) - 1;
        cftmp  = reshape(fcoefs(itot,:), numel(ipols), []);

        if nargin > 2 && ~isempty(p)
            nrms = vecnorm(cftmp, p, 1);
        else
            nrms = vecnorm(cftmp, Inf, 1);
        end

        errs = reshape(nrms, [], nfuns);
        errps(inds,:) = errs;

        areas = sum([obj.weights{inds}], 1);
        parea(inds) = areas;
    end

    if nargin > 2 && ~isempty(p)
        a     = area(obj);
        errps = (parea/a).^(1/p) .* errps;
    end

    errps = errps.';

end
