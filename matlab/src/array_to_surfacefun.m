function fvals = array_to_surfacefun(vals,dom,S)
% vals: double (npts) or double(npts,3);
% dom: domain defined as surfacemesh
% S: domain defined as surfer
%
% The values are assumed to be given on tensor product chebyshev 
% or legendre nodes depending on iptype of S. These are the first
% kind chebyshev nodes or Legendre nodes
%
% Note that the subroutine currently assumes that the order of
% discretization for each patch is the same

    
       
    iptype = S.iptype;
    npatches = S.npatches;
    
    [np,m] = size(vals);
    if(np~=S.npts) 
        fprintf('In ARRAY_TO_SURFACE_FUN: Invalid leading dimension of vals array\n');
        fprintf('returning empty array\n');
        if(m==3)
            fvals = surfacefunv(dom);
        elseif(m==1)
            fvals = surfacefun(dom);
        else
            fprintf('In ARRAY_TO_SURFACE_FUN: Invalid trailing dimension of vals array\n');
            fprintf('Returning empty surfacefun\n');
            fvals = surfacefun(dom);
        end
        return
    end
    
    if(m~=3 && m~=1)
        fprintf('In ARRAY_TO_SURFACE_FUN: Invalid trailing dimension of vals\n');
        fprintf('returning empty surfacefun\n')
        fvals = surfacefun(dom);
        return
    end
    
    
    
    [n,~] = size(dom.x{1});
    
    if(iptype == 11)
        [xint, ~, w] = polytens.lege.pts(n);
    elseif(iptype == 12)
        [xint, ~, w] = polytens.cheb.pts(n);
    else
        fprintf('Invalid iptype, iptype must be 11 or 12\n');
        fprintf('returning empty surface fun\n');
        if(issurfacefunv)
            fvals = surfacefunv(dom);
        else
            fvals = surfacefun(dom);
        end
        return
    end
    opts_use.kind = 2;
    xcheb = polytens.cheb.pts(n,opts_use);
    
    eval = barymat(xcheb,xint,w);
    
    if(m==1)
        valsuse = cell(npatches,1);
    
        for i=1:npatches
           istart = S.ixyzs(i);
           iend = S.ixyzs(i+1)-1;
           valstmp = vals(istart:iend,:);
           
           valstmp = reshape(valstmp(:),[n,n]);
           valstmp = valstmp.';
           valsuse{i} = eval * valstmp * eval.';
           
        end
        
        fvals = surfacefun(valsuse,dom);
    else
        fvalsx = array_to_surfacefun(vals(:,1),dom,S);
        fvalsy = array_to_surfacefun(vals(:,2),dom,S);
        fvalsz = array_to_surfacefun(vals(:,3),dom,S);
        fvals = surfacefunv(fvalsx,fvalsy,fvalsz);
    end


    
end
