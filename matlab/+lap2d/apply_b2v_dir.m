function val = apply_b2v_dir(S, chnkr, rho, Ab2v_cor, eps)
%
%  lap2d.apply_b2v_dir
%
%  FMM-accelerated evaluation of the Laplace double-layer potential from
%  boundary to volume targets:
%
%    val(r) = \int (d/dn') G(r,r') rho(r') ds(r'),   r in S
%
%  Inputs:
%    S        - surfer object (target domain)
%    chnkr    - chunker object (source boundary)
%    rho      - boundary density (npt x 1)
%    Ab2v_cor - near-quadrature correction sparse matrix
%    eps      - FMM precision
%
%  Output:
%    val      - potential values at S nodes (npts x 1)

    srcinfo = [];
    srcinfo.sources = chnkr.r(:,:);
    srcinfo.dipstr  = rho(:).' .* chnkr.wts(:).';
    srcinfo.dipvec  = chnkr.n(:,:);
    targ = S.r(1:2,:);
    U = lfmm2d(eps, srcinfo, 0, targ, 1);
    val = -1/(2*pi) * U.pottarg.';
    val = val + Ab2v_cor * rho(:);

end