function [errps] = surf_fun_koorn_shell_error(obj, fun, lptype, ...
    norderhead, functional_integrand)
%  SURF_FUN_KOORN_SHELL_ERROR  Per-patch Koornwinder shell estimator.
%
%  If norderhead = q, use the Koornwinder coefficients of total degree
%  q+1. For each patch and function, return
%
%    p = 0 or Inf : maximum absolute coefficient on the shell
%    p = 1        : mean absolute coefficient on the shell
%    p = 2        : RMS coefficient on the shell
%
%  With norderhead empty, q = norder-1, so the highest-degree shell is
%  used. No coefficient-to-value reconstruction or quadrature weights are
%  used.

    if nargin < 3 || isempty(lptype), lptype = 0; end
    if nargin < 4, norderhead = []; end
    if nargin < 5, functional_integrand = []; end

    if ~isempty(functional_integrand)
        if size(fun,2) == obj.npts
            fun = fun.';
        end
        x = obj.r(1,:).';
        fun = functional_integrand(fun,x);
    end

    patcharea = full(sparse(obj.patch_id(:),1,obj.wts(:), ...
    obj.npatches,1));

    fcoefs = obj.vals2coefs(fun);
    if size(fcoefs,2) == obj.npts
        fcoefs = fcoefs.';
    end

    nfuns = size(fcoefs,2);
    errps = zeros(obj.npatches,nfuns);

    [norders,~,iorder] = unique(obj.norders(:));

    for i = 1:numel(norders)
        norder = norders(i);
        inds = find(iorder == i);
        indv = obj.ixyzs(inds);
        npols = (norder+1)*(norder+2)/2;

        if isempty(norderhead)
            head = norder-1;
        else
            head = norderhead;
        end

     
        % Keep all degrees above norderhead, matching the value-tail estimator.
        ishell = ((head+1)*(head+2)/2 + 1):npols;

        [I,V] = meshgrid(indv,1:npols);
        itot = I(:)+V(:)-1;

        cf = reshape(fcoefs(itot,:),npols,numel(inds),nfuns);
        cshell = cf(ishell,:,:);

        if lptype == 0 || isinf(lptype)
            e = max(abs(cshell),[],1);
        elseif lptype == 1
            e = mean(abs(cshell),1);
        elseif lptype == 2
            e = sqrt(mean(abs(cshell).^2,1));
        end


        e = reshape(e,numel(inds),nfuns);

        if lptype == 1 || lptype == 2
            e = e .* patcharea(inds).^(1/lptype);
        end

        errps(inds,:) = e;
    end

    errps = errps.';
end
