function val = apply_v2v(S,zk,mu,Av2v_cor,nover,eps)
    
    if zk >= 1e-8
        error('No support for flexural problem')
    end

    [S_over,xinterp] = oversample(S,nover);
    test_fn_over = xinterp*mu;
    
    % lap2d_s = kernel('l','s');
    % A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
    % 
    % sol2 = A_native*test_fn_over + v2v_cor*test_fn;
    
    targ = S.r(1:2,:);

    x2 = S.r(1,:).^2+S.r(2,:).^2;
    y2 = S_over.r(1,:).^2+S_over.r(2,:).^2;

    srcinfo = []; 
    srcinfo.sources = S_over.r(1:2,:);
    srcinfo.charges = test_fn_over(:).'.*S_over.wts(:).';
    Ux = lfmm2d(eps,srcinfo,0,targ,1);
    Ux = (x2.*Ux.pottarg).';

    srcinfo.charges = S_over.r(1,:).*test_fn_over(:).'.*S_over.wts(:).';
    Uxy1 = lfmm2d(eps,srcinfo,0,targ,1);
    Uxy1 = (S.r(1,:).*Uxy1.pottarg).';

    srcinfo.charges = S_over.r(2,:).*test_fn_over(:).'.*S_over.wts(:).';
    Uxy2 = lfmm2d(eps,srcinfo,0,targ,1);
    Uxy2 = (S.r(2,:).*Uxy2.pottarg).';    

    srcinfo.charges = y2.*test_fn_over(:).'.*S_over.wts(:).';
    Uy = lfmm2d(eps,srcinfo,0,targ,1);
    Uy = (Uy.pottarg).';  

    val = 1/(8*pi)*(Ux-2*Uxy1-2*Uxy2+Uy);
    val = val + Av2v_cor*mu(:);

end