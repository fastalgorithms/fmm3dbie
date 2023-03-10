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
    
    
    
    ntot = sum([ns{inuni(1:obj.npatches)}]) + 3*obj.npatches;
    Ttot = cat(1,T{inuni(1:obj.npatches)});
    xall = zeros(ntot,1);
    yall = zeros(ntot,1);
    zall = zeros(ntot,1);
    
    
    istart = 0;
    itstart = 0;
    
    for k = 1:obj.npatches
        scl = 0.01;
        
        
        n = ns{inuni(k)};
        
        xall(istart+(1:n)) = obj.srcvals{k}(1,:);
        yall(istart+(1:n)) = obj.srcvals{k}(2,:);
        zall(istart+(1:n)) = obj.srcvals{k}(3,:);
        
        nend = istart + ((n+1):(n+3));
        xall(nend) = obj.srccoefs{k}(1,:)*pols{inuni(k)};
        yall(nend) = obj.srccoefs{k}(2,:)*pols{inuni(k)};
        zall(nend) = obj.srccoefs{k}(3,:)*pols{inuni(k)};
        
        [ntt,~] = size(T{inuni(k)});
        Ttot(itstart+(1:ntt),:) = Ttot(itstart+(1:ntt),:) + istart;
        istart = istart + n+3;
        itstart = itstart + ntt;
        
    end
    
    trisurf(Ttot,xall,yall,zall,zall,...
            'EdgeColor', 'None', ...
            'AmbientStrength', 0.6, 'DiffuseStrength', 0.4, 'SpecularStrength', 0.3);
end

if ( ~holdState )
    axis equal
end

if ( ~holdState )
    hold off
end

end
