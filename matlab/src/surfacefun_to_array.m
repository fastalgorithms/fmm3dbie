function vals = surfacefun_to_array(fvals,dom,S)
    
    iptype = S.iptype;
    [n,~] = size(dom.x{1});
    
    if(iptype == 11)
        xint = legpts(n);
    elseif(iptype == 12)
        xint = chebpts(n,[-1,1],1);
    else
        fprintf('Invalid iptype, iptype must be 11 or 12\n');
        fprintf('Exiting\n');
        S = [];
        return
    end
    xcheb = chebpts(n,[-1,1]);
    
    eval = barymat(xint,xcheb);
    npatches = S.npatches;
    npts = S.npts;
    if(isa(fvals,'surfacefun'))
        m = 1;
    elseif(isa(fvals,'surfacefunv'))
        m = 3;
    else
        fprtinf('In SURFACE_FUN_TO_ARRAY: invalid input type for fvals\n');
        fprintf('returning empty array\n');
        vals = [];
    end
    vals = zeros(npts,m);
    if(m == 1)
        for i=1:npatches
            istart = S.ixyzs(i);
            iend = S.ixyzs(i+1)-1;
           
            vtmp = eval * fvals.vals{i} * eval.';
            vtmp = vtmp.';
            vals(istart:iend,1) = vtmp(:);
        end  
    else
        for m=1:3
            fx = fvals.components{m};
            vals(:,m) = surfacefun_to_array(fx,dom,S);
        end
        
    end
    
    
end