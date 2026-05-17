function val = apply_v2b_neu(S,targinfo,mu,Av2b_cor,nover,eps)
    
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
    U = lfmm2d(eps,srcinfo,0,targ,2);
    val = -1/(2*pi)*(targinfo.n(1,:).*U.gradtarg(1,:) + ...
                targinfo.n(2,:).*U.gradtarg(2,:)).';
    val = val + Av2b_cor*mu(:);

end