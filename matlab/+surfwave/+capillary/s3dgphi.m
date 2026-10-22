function val = s3dgphi(rts,ejs,src,targ)
%SURFWAVE.CAPILLARY.S3DGPHI   the 3D single layer applied to the capillary
%   G_phi (without the Laplace part 1/(4 pi r)), (nt,ns).
%
% See also SURFWAVE.CAPILLARY.GREEN, SURFWAVE.ROOTSUM

val = surfwave.rootsum(rts(:), ejs(:)/8, src, targ);
end
