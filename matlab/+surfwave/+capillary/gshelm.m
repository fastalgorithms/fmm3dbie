function varargout = gshelm(rts,ejs,src,targ)
%SURFWAVE.CAPILLARY.GSHELM   capillary Green's function G_S and its target
%   derivatives: val (nt,ns), grad(:,:,1:2) = G_x, G_y and hess(:,:,1:3) =
%   G_xx, G_xy, G_yy, for the dispersion roots rts and residues ejs
%   (SURFWAVE.CAPILLARY.FIND_ROOTS_CAPILLARY).
%
% See also SURFWAVE.CAPILLARY.GREEN, SURFWAVE.ROOTSUM

rts = rts(:);  ejs = ejs(:);
[varargout{1:max(nargout,1)}] = surfwave.rootsum(rts, ejs.*rts.^2/2, src, targ);
end
