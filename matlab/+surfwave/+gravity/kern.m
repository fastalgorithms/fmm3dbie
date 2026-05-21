function submat= kern(rts,ejs,srcinfo,targinfo,type)
%SURFWAVE.GRAVITY.KERN gravity wave layer potential kernels
%
% Syntax: submat = surfwave.gravity.kern(rts,ejs,srcinfo,targinfo,type)
%
% Let x be targets and y be sources.
%
% Input:
%   rts     - dispersion root (scalar complex)
%   ejs     - residue corresponding to rts (unused for gravity, kept for
%             interface consistency with the capillary counterpart)
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%   targinfo - description of targets in ptinfo struct format,
%                ptinfo.r - positions (2,:) array
%   type - string, determines kernel type:
%                type == 'gs_s'    single layer G_s
%                type == 'gphi_s'  single layer G_phi
%
% Output:
%   submat - (ntarg,nsrc) matrix of kernel evaluations
%
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

% single layer Gs
if strcmpi(type,'gs_s')
  submat = surfwave.gravity.gsgrav(rts,src,targ);
end

% single layer Gphi
if strcmpi(type,'gphi_s')
  submat = surfwave.gravity.gphigrav(rts,src,targ);
end




