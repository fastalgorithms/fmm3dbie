function val = apply_b2v_neu(S,zk,chnkr,rho,Ab2v_cor,eps)
    
    % lap2d_s = kernel('l','s');
    % A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
    % 
    % sol2 = A_native*test_fn_over + v2v_cor*test_fn;
    
    srcinfo = []; 
    srcinfo.sources = chnkr.r(:,:);
    srcinfo.charges = rho(:).'.*chnkr.wts(:).';
    targ = S.r(1:2,:);
    U = hfmm2d(eps,zk,srcinfo,0,targ,1);
    val = U.pottarg.';
    val = val + Ab2v_cor*rho(:);
    
end