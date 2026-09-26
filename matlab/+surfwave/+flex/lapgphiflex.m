function val = lapgphiflex(rts,ejs,src,targ)
%LAPGPHIFLEX  Laplacian of the flexural Green's function G_phi.
%
%   Same partial-fraction loop as GPHIFLEX, with the weight ej*rhoj
%   replaced by -ej*rhoj^3.
%
%   Derivation.  With Q(z) = alpha z^5 + gamma z - 1, rts its roots and
%   ejs(j) = 1/Q'(rts(j)),
%
%       Ghat_phi(xi) = sum_j ej/(|xi| - rho_j) = 1/Q(|xi|) ,
%
%   so   -|xi|^2 Ghat_phi = -sum_j ej (|xi| + rho_j + rho_j^2/(|xi|-rho_j))
%                         = -sum_j ej rho_j^2/(|xi| - rho_j) ,
%   the first two terms dropping by the moment conditions
%   sum_j ej = sum_j ej rho_j = 0.  Multiplying the weight by rho_j^2 in
%   the spectral variable multiplies it by rho_j^2 in real space as well,
%   because the Struve function K_0 = H_0 - Y_0
%   satisfies the inhomogeneous Bessel equation
%
%       Lap K_0(a r) = -a^2 K_0(a r) + 2a/(pi r) ,
%
%   whose 1/r remainder cancels under sum_j ej rho_j^2 = 0.  Hence
%
%       Lap G_phi = -1/4 sum_j ej rho_j^3 [ bracket_j ] ,
%
%   with the same branch brackets used by GPHIFLEX.
%
% See also SURFWAVE.FLEX.GPHIFLEX, SURFWAVE.FLEX.S3DGPHIFLEX.

if isstruct(src)
    src = src.r;
end
if isstruct(targ)
    targ = targ.r;
end

val = 0;

for i = 1:5

    rhoj = rts(i);
    ej   = ejs(i);

    if abs(angle(rhoj)) < 1e-8

       sk0 = surfwave.flex.struveKdiffgreen(rhoj,src,targ);
       h0  = surfwave.flex.helmdiffgreen(rhoj,src,targ);

       h0 = -4i*h0;

       val = val + ej*rhoj^3*(-sk0 + 2i*h0);

    else

       sk0 = surfwave.flex.struveKdiffgreen(-rhoj,src,targ);

       val = val + ej*rhoj^3*sk0;

    end

end

val = -1/4*val;

end
