function val = apply_v2v(S,mu,Av2v_cor,nover,eps)
    
    [S_over,xinterp] = oversample(S,nover);
    test_fn_over = xinterp*mu;
    
    % lap2d_s = kernel('l','s');
    % A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
    % 
    % sol2 = A_native*test_fn_over + v2v_cor*test_fn;
    
    srcinfo = []; 
    srcinfo.sources = S_over.r(1:2,:);
    srcinfo.charges = test_fn_over(:).'.*S_over.wts(:).';
    targ = S.r(1:2,:);
    U = lfmm2d(eps,srcinfo,0,targ,1);
    val = -1/(2*pi)*U.pottarg.';
    val = val + Av2v_cor*mu(:);

end