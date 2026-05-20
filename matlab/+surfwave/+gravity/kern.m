function submat= kern(rts,ejs,srcinfo,targinfo,type)
% 
% Syntax: submat = surfwave.gravity.kern(rts,srcinfo,targinfo,type)
%
% Let x be targets and y be sources for these formulas,
%  
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'g_s'
%                type == 'g_phi'
%
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
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




