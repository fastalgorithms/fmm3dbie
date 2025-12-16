function bmat = coefs2vals(norder,uvs)
    bmat = polytens.lege.pols(norder,uvs);
    bmat = bmat.';
end