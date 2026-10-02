function submat= kern(srcinfo,targinfo,type,varargin)
%LAP2D.KERN standard Laplace volume potential kernels in 2D
% 
% Syntax: submat = lap2d.kern(srcinfo,targinfo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x the unit normal at the target (if defined). The sources are
% points of a flat surface in the z = 0 plane, and only the first two
% coordinates of the sources and targets are used.
%  
% Kernels based on G(x,y) = -log|x-y|/(2*pi)
%
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
%
% Input:
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (3,:) array
%   targinfo - description of targets in ptinfo struct format,
%                sprime requires normal info in targinfo.n.
%   type - string, determines kernel type
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'sgrad', gradient of single layer, returns
%                      (2*nt,ns) with d/dx and d/dy interleaved
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also LAP2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

switch lower(type)

    case {'s', 'single'}
        submat = lap2d.green(src,targ);

    case {'sp', 'sprime'}
        [~,grad] = lap2d.green(src,targ,true);
        submat = grad(:,:,1).*targinfo.n(1,:).' + ...
                 grad(:,:,2).*targinfo.n(2,:).';

    case {'sg', 'sgrad'}
        [~,grad] = lap2d.green(src,targ,true);
        submat = reshape(permute(grad,[3,1,2]),2*nt,ns);

    otherwise
        error('LAP2D.KERN: unknown kernel type ''%s''.', type);

end

end
