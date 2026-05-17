function val = apply_v2b_dir(S,zk,targinfo,mu,Av2b_cor,nover,eps)
    
    [S_over,xinterp] = oversample(S,nover);
    mu_over = xinterp*mu;
    
    % lap2d_s = kernel('l','s');
    % A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
    % 
    % sol2 = A_native*test_fn_over + v2v_cor*test_fn;
    
    srcinfo = []; 
    srcinfo.sources = S_over.r(1:2,:);
    srcinfo.charges = mu_over(:).'.*S_over.wts(:).';
    targ = targinfo.r(1:2,:);
    U = hfmm2d(eps,zk,srcinfo,0,targ,1);
    val = U.pottarg.';
    val = val + Av2b_cor*mu(:);

end