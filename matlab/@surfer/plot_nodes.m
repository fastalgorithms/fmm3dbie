function h = plot_nodes(S, v, varargin)
% PLOT_NODES  show nodes of a surfer object, or function on nodes
%
% plot_nodes(S) shows all nodes of surfer object S as a point cloud
%  on the current figure. Returns handle to a scatter3 object.
%
% plot_nodes(S,v) uses the values in v, a vector of length S.npts, to color
%  the nodes.
%
% plot_nodes(S,[],...) or plot_nodes(S,v,...) passes any other arguments to
%  scatter3.
%
% Some view and properties of the current axes are set for convenience of
%  3D viewing.
%
% Example
%   S = geometries.sphere(1.0);
%   plot_nodes(S)
%   h = plot_nodes(S, S.r(3,:))            % color by z-coord
%   h.SizeData = 1000;                     % make larger blobs
%
% See also PLOT, SCATTER3

% Barnett 6/5/24
  
if nargin==1 || isempty(v)
  h = scatter3(S.r(1,:),S.r(2,:), S.r(3,:), '.', varargin{:});
else
  if length(v)~=S.npts
    error('v vector must have length = S.npts!');
  end
  h = scatter3(S.r(1,:),S.r(2,:), S.r(3,:), 100, v, '.', varargin{:});
end
view(3)
axis equal

