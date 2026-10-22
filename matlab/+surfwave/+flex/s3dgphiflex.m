function val = s3dgphiflex(rts,ejs,src,targ)
%SURFWAVE.FLEX.S3DGPHIFLEX   the 3D single layer applied to the flexural
%   G_phi (without the Laplace part 1/(4 pi r)), (nt,ns).
%
% See also SURFWAVE.FLEX.GREEN, SURFWAVE.ROOTSUM

val = surfwave.rootsum(rts(:), ejs(:)/8, src, targ, 1e-8);
end
