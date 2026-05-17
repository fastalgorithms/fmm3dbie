function val = apply_b2v_neu(S,chnkr,rho,Ab2v_cor,eps)
    
    % lap2d_s = kernel('l','s');
    % A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
    % 
    % sol2 = A_native*test_fn_over + v2v_cor*test_fn;
    
    srcinfo = []; 
    srcinfo.sources = chnkr.r(:,:);
    srcinfo.dipstr = rho(:).'.*chnkr.wts(:).';
    srcinfo.dipvec = 
    targ = S.r(1:2,:);
    U = lfmm2d(eps,srcinfo,0,targ,1);
    val = -1/(2*pi)*U.pottarg.';
    val = val + Ab2v_cor*rho(:);

end