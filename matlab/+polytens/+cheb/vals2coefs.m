function amat = vals2coefs(norder,uvs)
    bmat = polytens.cheb.coefs2vals(norder,uvs);
    amat = inv(bmat);
end