function bmat = coefs2vals(norder,uvs)
%
%  This function returns the values to 
%  coefficients matrix for a given
%  set of nodes uvs. Note that
%  there must be (norder+1)*(norder+2)/2 
%  nodes in uvs
%
     bmat = koorn.pols(norder,uvs);
     bmat = bmat';

end
