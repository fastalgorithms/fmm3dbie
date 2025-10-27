function amat = vals2coefs(norder,uvs)
    bmat = polytens.lege.coefs2vals(norder,uvs);
    amat = inv(bmat);
end