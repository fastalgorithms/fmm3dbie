function amat = lege_vals2coefs(norder,uvs)
    bmat = polytens.lege_coefs2vals(norder,uvs);
    amat = inv(bmat);
end