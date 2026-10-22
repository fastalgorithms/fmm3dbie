function varargout = gsgrav_der(rts, ejs, src, targ)
%SURFWAVE.GRAVITY.GSGRAV_DER   gravity Green's function and its target
%   derivatives,
%
%      G(x,y) = sum_j ej rho_j/(2 pi r) + ej rho_j^2/4 K(rho_j r),
%
%   K the radial profile of SURFWAVE.ROOTSUM (singularities removed in
%   closed form, so accurate as r -> 0). rts = rho (= g/2) and ejs = 1
%   give G_S, ejs = 1/g gives G_phi.
%
%   [val,grad,hess,third,fourth]: val (nt,ns), grad(:,:,1:2),
%   hess(:,:,1:3), third(:,:,1:4), fourth(:,:,1:5), packed as
%   SURFWAVE.HELMGREEN.
%
% See also SURFWAVE.GRAVITY.GREEN, SURFWAVE.GRAVITY.GSGRAV

rts = rts(:);  ejs = ejs(:);
nout = max(nargout, 1);
varargout = cell(1, nout);
[varargout{:}] = surfwave.rootsum(rts, ejs.*rts.^2/4, src, targ);
L = cell(1, nout);
[L{:}] = surfwave.lap3dgreen(src, targ);
c = 2*sum(ejs.*rts);                  % ej rho/(2 pi r) = 2 ej rho/(4 pi r)
for k = 1:nout, varargout{k} = varargout{k} + c*L{k}; end
end
