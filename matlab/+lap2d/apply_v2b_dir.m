function val = apply_v2b_dir(S, targinfo, mu, Av2b_cor, nover, eps)
%
%  lap2d.apply_v2b_dir
%
%  FMM-accelerated evaluation of the Laplace single-layer volume potential
%  at boundary (or off-surface) targets:
%
%    val(r) = \int G(r,r') mu(r') dA(r'),   r in targinfo
%
%  where G(r,r') = -1/(2*pi) log|r-r'| is the 2D Laplace SLP kernel.
%
%  Inputs:
%    S        - surfer object (source domain)
%    targinfo - target struct with field r (2 or 3 x ntarg)
%    mu       - density values at S nodes (npts x 1)
%    Av2b_cor - near-quadrature correction sparse matrix
%    nover    - oversampling order for smooth part
%    eps      - FMM precision
%
%  Output:
%    val      - potential values at targets (ntarg x 1)

    [S_over, xinterp] = oversample(S, nover);
    mu_over = xinterp * mu;

    srcinfo = [];
    srcinfo.sources = S_over.r(1:2,:);
    srcinfo.charges = mu_over(:).' .* S_over.wts(:).';
    targ = targinfo.r(1:2,:);
    U = lfmm2d(eps, srcinfo, 0, targ, 1);
    val = -1/(2*pi) * U.pottarg.';
    val = val + Av2b_cor * mu(:);

end
