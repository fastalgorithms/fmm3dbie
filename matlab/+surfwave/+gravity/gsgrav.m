function [val] = gsgrav(rts,src,targ)
%
% Evaluates the gravity wave free-surface Green's function G_s, built
% from the single dispersion root of the linear dispersion relation
%   |z| - 1 = 0  (gravity waves, non-dimensionalized).
%
% Input:
%   rts  - dispersion root (scalar complex)
%   src  - source positions, (2,ns) array or struct with field .r
%   targ - target positions, (2,nt) array or struct with field .r
%
% Output:
%   val  - (nt,ns) array of Green's function values G_s(targ,src)
%

% src = src.r;
% targ = targ.r;

rslf = 1e-14;

if isstruct(src)
    src = src.r;
end

if isstruct(targ)
    targ = targ.r;
end

eulergamma = 0.57721566490153286060651209008240243;

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

r = sqrt((xt-xs).^2 + (yt-ys).^2);

% clear xs ys xt yt

val = 0;

rhoj = rts;

if (abs(angle(rhoj)) < rslf) && (abs(rhoj) > rslf)

   zt = r*rhoj;
   [cr0] = surfwave.capillary.struveR(zt);

   h0 = chnk.helm2d.green(rhoj,src,targ);
   % h0(r < rslf) = 1/(2*pi)*(1i*pi/2  - eulergamma + log(2/rhoj));
   h0 = -4i*h0;

   sk0 = -1i*cr0+1i*h0;

   val = val + rhoj^2*(-sk0 + 2i*h0);

end

val = 1/2*val;
val = val + rts./(r.*pi);
val(r < rslf) = 0;

end