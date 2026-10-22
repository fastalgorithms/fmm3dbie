function val = lapgphiflex(rts,ejs,src,targ)
%SURFWAVE.FLEX.LAPGPHIFLEX   Laplacian of the flexural G_phi, (nt,ns).
%
%   With Q(z) = alpha z^5 + gamma z - 1 and ejs = 1/Q'(rts),
%   -|xi|^2 Ghat_phi = -sum_j ej rho_j^2/(|xi| - rho_j) by the moment
%   conditions sum_j ej = sum_j ej rho_j = 0, so the Laplacian is the root
%   sum of G_phi with weights ej rho_j replaced by -ej rho_j^3.
%
% See also SURFWAVE.FLEX.GREEN, SURFWAVE.ROOTSUM

val = surfwave.rootsum(rts(:), -ejs(:).*rts(:).^3/4, src, targ, 1e-8);
end
