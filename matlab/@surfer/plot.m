function plot(obj, varargin)
%PLOT   Plot the surface.

parser = inputParser;
parser.KeepUnmatched = true;
parser.addParameter('surface', 'auto', @(s) contains(lower(s), {'auto', 'on', 'off'}));
parse(parser, varargin{:});
showSurface = parser.Results.surface;
varargin = namedargs2cell(parser.Unmatched);

defaultStyle = {'Color', 'k', 'LineStyle', '-', 'LineWidth', 1};

holdState = ishold();

if ( (~holdState && strcmpi(showSurface, 'auto')) || strcmpi(showSurface, 'on') )
    % Plot the surface, making it slightly smaller so lines show up
    % more clearly.
    hold on
    for k = 1:obj.npatches
        scl = 0.01;
        norder = obj.norders(1);
        rv_nodes = koorn.rv_nodes(norder);
        [~,n] = size(rv_nodes);
        rv_nodes_all = zeros(2,n+3);
        rv_nodes_all(:,1:n) = rv_nodes;
        rv_nodes_all(:,n+2) = [1,0];
        rv_nodes_all(:,n+3) = [0,1];
        T = delaunay(rv_nodes_all(1,:),rv_nodes_all(2,:));
        x = zeros(n+3,1);
        y = zeros(n+3,1);
        z = zeros(n+3,1);
        x(1:n) = obj.srcvals{k}(1,:);
        y(1:n) = obj.srcvals{k}(2,:);
        z(1:n) = obj.srcvals{k}(3,:);
        
        nend = (n+1):(n+3);
        pols = koorn.pols(norder,rv_nodes_all(:,nend));
        x(nend) = obj.srccoefs{k}(1,:)*pols;
        y(nend) = obj.srccoefs{k}(2,:)*pols;
        z(nend) = obj.srccoefs{k}(3,:)*pols;
        
        trisurf(T,x,y,z,z,...
                'EdgeColor', 'None', ...
                'AmbientStrength', 0.6, 'DiffuseStrength', 0.4, 'SpecularStrength', 0.3);
    end
end

if ( ~holdState )
    axis equal
end

if ( ~holdState )
    hold off
end

end
