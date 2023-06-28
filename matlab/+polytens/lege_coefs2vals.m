function bmat = lege_coefs2vals(norder,uvs)
    bmat = polytens.lege_pols(norder,uvs);
    bmat = bmat.';
end