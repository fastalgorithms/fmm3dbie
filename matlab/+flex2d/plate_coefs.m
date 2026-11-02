function cf = plate_coefs(p, nu, zk)
%FLEX2D.PLATE_COEFS  unpack the plate coefficient struct p (fields alpha,
% dalpha, d2alpha, beta) into the (7,nt) target coefficients of
%
%   L u = alpha L0 u + c1 d_x Lap u + c2 d_y Lap u + c3 Lap u
%         + c4 u_yy + c5 u_xx + c6 u_xy + c7 u
%
% applied to G, with L0 = Delta^2 + b0 Delta + c0 the operator of
% FLEX2D.GREEN for the wavenumbers zk, and
%
%   c1 = 2 alpha_x,  c2 = 2 alpha_y,  c3 = Lap alpha - b0 alpha,
%   c4 = -(1-nu) alpha_xx,  c5 = -(1-nu) alpha_yy,
%   c6 = 2 (1-nu) alpha_xy,  c7 = -c0 alpha - beta.

zk = zk(:).';
if isscalar(zk)
    zk = [zk, 1i*zk];
end
zk(abs(zk) < 1e-6) = 0;
b0 = zk(1)^2 + zk(2)^2;
c0 = zk(1)^2*zk(2)^2;

alpha = p.alpha(:).';
nt    = numel(alpha);
da    = reshape(p.dalpha, 2, nt);
d2a   = reshape(p.d2alpha, 3, nt);
beta  = p.beta(:).';

cf = [2*da(1,:); ...
      2*da(2,:); ...
      d2a(1,:) + d2a(3,:) - b0*alpha; ...
      -(1-nu)*d2a(1,:); ...
      -(1-nu)*d2a(3,:); ...
      2*(1-nu)*d2a(2,:); ...
      -c0*alpha - beta];

end
