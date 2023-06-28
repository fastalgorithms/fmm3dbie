function amat = cheb_vals2coefs(norder,uvs)
    bmat = polytens.cheb_coefs2vals(norder,uvs);
    amat = inv(bmat);
end