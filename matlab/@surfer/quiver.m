function varargout = quiver(S, v, noder, varargin)
% QUIVER  plot vector field on a surfer object using quiver3
%
% quiver(S, v, noder) oversamples the geometry S and the vector field v to
%   order noder, then calls quiver3 to plot arrows at each node.
%
%   S      - surfer object
%   v      - (3 x S.npts) array of vectors at the nodes of S, or
%            (3*S.npts x 1) vector in the same ordering
%   noder  - oversampling order (integer). Pass S.norders (or 0) for no
%            oversampling.
%
% quiver(S, v, noder, ...) passes any extra arguments through to quiver3.
%
% h = quiver(...) returns the handle to the quiver3 object.
%
% Example
%   S = geometries.sphere(1.0);
%   quiver(S, S.n, 3)                      % plot normals, oversampled to order 3
%   h = quiver(S, S.n, S.norders(1));      % no oversampling
%
% See also PLOT_NODES, OVERSAMPLE

% Barnett 5/21/26

v = reshape(v, 3, S.npts);

% Oversample geometry and interpolate v to oversampled nodes
[Sover, Pmat] = oversample(S, noder);
vover = v * Pmat';   % (3 x nptso)

x = Sover.r(1,:);
y = Sover.r(2,:);
z = Sover.r(3,:);
u = vover(1,:);
vv = vover(2,:);
w = vover(3,:);

h = quiver3(x, y, z, u, vv, w, varargin{:});

if norm(Sover.r(3,:)) == 0
    view(0,90)
else
    view(3)
end
axis equal

if nargout >= 1
    varargout{1} = h;
end
