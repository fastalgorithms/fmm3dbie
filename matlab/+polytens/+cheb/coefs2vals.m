function bmat = coefs2vals(norder,uvs)
    bmat = polytens.cheb.pols(norder,uvs);
    bmat = bmat.';
end