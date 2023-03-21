function amat = vals2coefs(norder,uvs)
%
%  This function returns the values to 
%  coefficients matrix for a given
%  set of nodes uvs. Note that
%  there must be (norder+1)*(norder+2)/2 
%  nodes in uvs
%
  bmat = koorn.coefs2vals(norder,uvs);
  amat = inv(bmat);

end
