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
    [norder_uni,~,inuni] = unique(obj.norders);
    nuni = length(norder_uni);
    T = cell(nuni,1);
    pols = cell(nuni,1);
    ns = cell(nuni,1);
    for i=1:nuni
        norder = norder_uni(i);
        rv_nodes = koorn.rv_nodes(norder);
        [~,n] = size(rv_nodes);
        rv_nodes_all = zeros(2,n+3);
        ns{i} = n;
        rv_nodes_all(:,1:n) = rv_nodes;
        rv_nodes_all(:,n+2) = [1,0];
        rv_nodes_all(:,n+3) = [0,1];
        T{i} = delaunay(rv_nodes_all(1,:),rv_nodes_all(2,:));
        nend = (n+1):(n+3);
        pols{i} = koorn.pols(norder,rv_nodes_all(:,nend));
        
    end
    
    for k = 1:obj.npatches
        scl = 0.01;
        
        
        n = ns{inuni(k)};
        x = zeros(n+3,1);
        y = zeros(n+3,1);
        z = zeros(n+3,1);
        x(1:n) = obj.srcvals{k}(1,:);
        y(1:n) = obj.srcvals{k}(2,:);
        z(1:n) = obj.srcvals{k}(3,:);
        
        nend = (n+1):(n+3);
        x(nend) = obj.srccoefs{k}(1,:)*pols{inuni(k)};
        y(nend) = obj.srccoefs{k}(2,:)*pols{inuni(k)};
        z(nend) = obj.srccoefs{k}(3,:)*pols{inuni(k)};
        
        trisurf(T{inuni(k)},x,y,z,z,...
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
