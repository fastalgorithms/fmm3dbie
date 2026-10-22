function val = gsgrav(rho, ejs, src, targ)
%SURFWAVE.GRAVITY.GSGRAV   gravity Green's function,
%
%   G_S(x,y) = ej*rho/(2 pi |x-y|)
%              + (ej*rho^2/4) [ -K0(rho|x-y|) + 2i H0^(1)(rho|x-y|) ]
%
%   G_phi = G_S/g with g = 2*rho, i.e. gsgrav(rho, 1/g, ...). rho may be a
%   vector (sum over roots) with ejs elementwise. See
%   SURFWAVE.GRAVITY.GSGRAV_DER for derivatives.
%
% See also SURFWAVE.GRAVITY.GREEN

val = surfwave.gravity.gsgrav_der(rho, ejs, src, targ);
end
