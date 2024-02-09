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
    norder_iptype = zeros(obj.npatches,2);
    norder_iptype(:,1) = obj.norders;
    norder_iptype(:,2) = obj.iptype;
    [no_ip_uni,~,inuni] = unique(norder_iptype,'rows');
    nuni = size(no_ip_uni,1);
    T = cell(nuni,1);
    pols = cell(nuni,1);
    ns = cell(nuni,1);
    nstot = cell(nuni,1);
    for i=1:nuni
        ip0 = no_ip_uni(i,2);
        norder = no_ip_uni(i,1);
        if ip0 == 1
            rnodes = koorn.rv_nodes(norder);
            [~,n] = size(rnodes);
            rnodes_all = zeros(2,n+3);
            
            rnodes_all(:,1:n) = rnodes;
            rnodes_all(:,n+2) = [1,0];
            rnodes_all(:,n+3) = [0,1];
            nend = (n+1):(n+3);
            pols{i} = koorn.pols(norder,rnodes_all(:,nend));
            nstot{i} = n+3;
        elseif ip0 == 11
            rnodes = polytens.lege_nodes(norder);
            [~,n] = size(rnodes);
            rnodes_all = zeros(2,n+4);
            rnodes_all(:,1:n) = rnodes;
            rnodes_all(:,n+1) = [-1,-1];
            rnodes_all(:,n+2) = [-1,1];
            rnodes_all(:,n+3) = [1,1];
            rnodes_all(:,n+4) = [1,-1];
            nend = (n+1):(n+4);
            pols{i} = polytens.lege_pols(norder,rnodes_all(:,nend));
            nstot{i} = n+4;
        elseif ip0 == 12
            rnodes = polytens.cheb_nodes(norder);
            [~,n] = size(rnodes);
            rnodes_all = zeros(2,n+4);
            rnodes_all(:,1:n) = rnodes;
            rnodes_all(:,n+1) = [-1,-1];
            rnodes_all(:,n+2) = [-1,1];
            rnodes_all(:,n+3) = [1,1];
            rnodes_all(:,n+4) = [1,-1];
            nend = (n+1):(n+4);
            pols{i} = polytens.cheb_pols(norder,rnodes_all(:,nend));
            nstot{i} = n+4;
        end
        ns{i} = n;
        T{i} = delaunay(rnodes_all(1,:),rnodes_all(2,:));     
    end
    
    
    
    ntot = sum([nstot{inuni(1:obj.npatches)}]);
    Ttot = cat(1,T{inuni(1:obj.npatches)});
    xall = zeros(ntot,1);
    yall = zeros(ntot,1);
    zall = zeros(ntot,1);
    
    
    istart = 0;
    itstart = 0;
    
    for k = 1:obj.npatches
        
        n = ns{inuni(k)};
        xall(istart+(1:n)) = obj.srcvals{k}(1,:);
        yall(istart+(1:n)) = obj.srcvals{k}(2,:);
        zall(istart+(1:n)) = obj.srcvals{k}(3,:);
        nextra = nstot{inuni(k)} - n;
        nend = istart + ((n+1):(n+nextra));
        xall(nend) = obj.srccoefs{k}(1,:)*pols{inuni(k)};
        yall(nend) = obj.srccoefs{k}(2,:)*pols{inuni(k)};
        zall(nend) = obj.srccoefs{k}(3,:)*pols{inuni(k)};
        
        [ntt,~] = size(T{inuni(k)});
        Ttot(itstart+(1:ntt),:) = Ttot(itstart+(1:ntt),:) + istart;
        istart = istart + n+nextra;
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
