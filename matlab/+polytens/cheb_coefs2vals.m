function bmat = cheb_coefs2vals(norder,uvs)
    bmat = polytens.cheb_pols(norder,uvs);
    bmat = bmat.';
end