function val = apply_b2v_dir(S,zk,chnkr,rho,Ab2v_cor,eps)
    
    srcinfo = []; 
    srcinfo.sources = chnkr.r(:,:);
    srcinfo.dipstr = rho(:).'.*chnkr.wts(:).';
    srcinfo.dipvec = chnkr.n(:,:);
    targ = S.r(1:2,:);
    U = hfmm2d(eps,zk,srcinfo,0,targ,1);
    val = U.pottarg.';
    val = val + Ab2v_cor*rho(:);
    
end