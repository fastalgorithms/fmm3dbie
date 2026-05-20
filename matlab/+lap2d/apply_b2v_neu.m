function val = apply_b2v_neu(S,chnkr,rho,Ab2v_cor,eps)
    
    srcinfo = []; 
    srcinfo.sources = chnkr.r(:,:);
    srcinfo.charges = rho(:).'.*chnkr.wts(:).';
    targ = S.r(1:2,:);
    U = lfmm2d(eps,srcinfo,0,targ,1);
    val = -1/(2*pi)*U.pottarg.';
    val = val + Ab2v_cor*rho(:);

end