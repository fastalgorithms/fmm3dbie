function [S,varargout] = surfacemesh_to_surfer(dom,opts)
    if(nargin==1)
        opts = [];
    end
    iptype = 11;
    if(isfield(opts,'iptype'))
        iptype = opts.iptype;
    end
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
    [npatches,~] = size(dom.x);
    
    npts = npatches*n*n;
    srcvals = zeros(12,npts);
    npols = n*n;
    
%
%  reorient domain 
% 
    rndom = normal(dom);
    x = dom.x;
    y = dom.y;
    z = dom.z;
%  Estimate normals via xu \times xv    
    xyzu  = surfacefunv(surfacefun(dom.xu,dom), ...
                surfacefun(dom.yu,dom),surfacefun(dom.zu,dom));
    xyzv  = surfacefunv(surfacefun(dom.xv,dom), ... 
                surfacefun(dom.yv,dom),surfacefun(dom.zv,dom));

    rndom2 = cross(xyzu,xyzv);
    rndom2 = rndom2./norm(rndom2);
    normm = norm(rndom-rndom2);
    normp = norm(rndom+rndom2);
    
    for i=1:npatches
        if(norm(normp.vals{i})<norm(normm.vals{i}))
            x{i} = x{i}';
            y{i} = y{i}';
            z{i} = z{i}';
        end
    end
    
    dom_out = surfacemesh(x,y,z);
    
    rndom3 = normal(dom_out);
    xyzu  = surfacefunv(surfacefun(dom_out.xu,dom_out), ...
                surfacefun(dom_out.yu,dom_out),surfacefun(dom_out.zu,dom_out));
    xyzv  = surfacefunv(surfacefun(dom_out.xv,dom_out), ... 
                surfacefun(dom_out.yv,dom_out),surfacefun(dom_out.zv,dom_out));

    rndom4 = cross(xyzu,xyzv);
    rndom4 = rndom4./norm(rndom4);
    
    err = norm(norm(rndom3-rndom4));
    assert(err<1e-13,'Incorrect normals post orientation fix\n');
    
    
    

    for i=1:npatches
        x2 = eval * dom_out.x{i} * eval.';
        y2 = eval * dom_out.y{i} * eval.';
        z2 = eval * dom_out.z{i} * eval.';
        
        dxdu = eval * dom_out.xu{i} * eval.';
        dydu = eval * dom_out.yu{i} * eval.';
        dzdu = eval * dom_out.zu{i} * eval.';
        
        dxdv = eval * dom_out.xv{i} * eval.';
        dydv = eval * dom_out.yv{i} * eval.';
        dzdv = eval * dom_out.zv{i} * eval.';
        
        x2 = x2';
        y2 = y2';
        z2 = z2';
        
        dxdu = dxdu';
        dydu = dydu';
        dzdu = dzdu';
        
        
        dxdv = dxdv';
        dydv = dydv';
        dzdv = dzdv';
        
        istart = (i-1)*npols+1;
        iend = i*npols;
        
        srcvals(1,istart:iend) = x2(:);
        srcvals(2,istart:iend) = y2(:);
        srcvals(3,istart:iend) = z2(:);
        
        srcvals(4,istart:iend) = dxdu(:);
        srcvals(5,istart:iend) = dydu(:);
        srcvals(6,istart:iend) = dzdu(:);
        
        srcvals(7,istart:iend) = dxdv(:);
        srcvals(8,istart:iend) = dydv(:);
        srcvals(9,istart:iend) = dzdv(:);
       
    end
    
    dn = cross(srcvals(4:6,:),srcvals(7:9,:));
    srcvals(10:12,:) = dn./repmat(vecnorm(dn,2,1),[3,1]);
    
    S = surfer(npatches,n-1,srcvals,iptype);
    varargout{1} = dom_out;

end