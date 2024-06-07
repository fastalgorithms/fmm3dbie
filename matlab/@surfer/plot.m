function varargout = plot(obj, varargin)
% PLOT  Plot surface patches of a surfer object in current figure
%
% plot(S) shows all patches, plotting in color the curvature function on nodes.
%
% plot(S,v) where v is vector of length S.npatches plots each patch
%  in a color given by the elements of v. If instead v has length S.npts,
%  plots as a function on the nodes. Both use current colormap.
%
% plot(S,v,...) passes other arguments (...) to trisurf.
%
% Some view and properties of the current axes are set for convenience of
%  3D viewing.
%
% Example
%   S = geometries.sphere(1.0);
%   plot(S,rand(1,S.npatches))    % shows patches in random colors
%
% See also TRISURF
  
[srcvals,srccoefs,~,~,~,~] = extract_arrays(obj);
f = obj.mean_curv.';
fcoefs = vals_to_coefs_surface_fun(obj, f);
istart = 1;
if nargin > 1
    if isnumeric(varargin{1})
        if length(varargin{1}(:)) == obj.npts
            f = reshape(varargin{1}, [1, obj.npts]);
            fcoefs = vals_to_coefs_surface_fun(obj, f);
            istart = 2;
        elseif length(varargin{1}(:)) == obj.npatches
            fuse = reshape(varargin{1}, [1, obj.npatches]);
            f = zeros(1, obj.npts);
            for i = 1:obj.npatches
                i1 = obj.ixyzs(i);
                i2 = obj.ixyzs(i+1)-1;
                f(i1:i2) = fuse(i);
            end
            fcoefs = vals_to_coefs_surface_fun(obj, f);
            istart = 2;
        else
            error('plot array does not have length = npts or npatches');
        end
    end
end
parser = inputParser;
parser.KeepUnmatched = true;
parser.parse(varargin{istart:end});
argnames = fieldnames(parser.Unmatched);
argvals = struct2cell(parser.Unmatched);
args = [argnames(:).'; argvals(:).'];
varargin_use = args(:).';




% Plot the surface, making it slightly smaller so lines show up
% more clearly.

norder_iptype = zeros(obj.npatches,2);
norder_iptype(:,1) = obj.norders;
norder_iptype(:,2) = obj.iptype;
[no_ip_uni,~,iuni] = unique(norder_iptype,'rows');
nuni = size(no_ip_uni,1);
T = cell(nuni,1);
pols = cell(nuni,1);
nstot = cell(nuni,1);

nplotpts = 10;
x1 = 0:1/nplotpts:1;
[x, y] = meshgrid(x1);
tri_quad = delaunay(x,y);


xx = x(x+y<=1);
yy = y(x+y<=1);
tri_tri = delaunay(xx,yy);

x = 2*x - 1;
y = 2*y - 1;

for i=1:nuni
    ip0 = no_ip_uni(i,2);
    norder = no_ip_uni(i,1);
    if ip0 == 1
        nuse = length(xx(:));
        rnodes = zeros(2,nuse);
        rnodes(1,:) = xx(:);
        rnodes(2,:) = yy(:);
        pols{i} = koorn.pols(norder, rnodes);
        nstot{i} = nuse;
        T{i} = tri_tri;
    elseif ip0 == 11
        nuse = length(x(:));
        rnodes = zeros(2,nuse);
        rnodes(1,:) = x(:);
        rnodes(2,:) = y(:);
        pols{i} = polytens.lege_pols(norder, rnodes);
        nstot{i} = nuse;
        T{i} = tri_quad;
    elseif ip0 == 12
        nuse = length(x(:));
        rnodes = zeros(2,nuse);
        rnodes(1,:) = x(:);
        rnodes(2,:) = y(:);
        pols{i} = polytens.cheb_pols(norder, rnodes);
        nstot{i} = nuse;
        T{i} = tri_quad;
    end
end

   
ntot = sum([nstot{iuni(1:obj.npatches)}]);
Ttot = cat(1,T{iuni(1:obj.npatches)});
xall = zeros(ntot,1);
yall = zeros(ntot,1);
zall = zeros(ntot,1);
fall = zeros(ntot,1);


istart = 0;
itstart = 0;

for k = 1:obj.npatches
    i1 = obj.ixyzs(k);
    i2 = obj.ixyzs(k+1)-1;

    fcuse = fcoefs(i1:i2);
    fvals = pols{iuni(k)}.'*fcuse(:);

    n = nstot{iuni(k)};
    xall(istart+(1:n)) = obj.srccoefs{k}(1,:)*pols{iuni(k)};
    yall(istart+(1:n)) = obj.srccoefs{k}(2,:)*pols{iuni(k)};
    zall(istart+(1:n)) = obj.srccoefs{k}(3,:)*pols{iuni(k)};
    fall(istart+(1:n)) = fvals(:);
    
           
    [ntt,~] = size(T{iuni(k)});
    Ttot(itstart+(1:ntt),:) = Ttot(itstart+(1:ntt),:) + istart;
    istart = istart + n;
    itstart = itstart + ntt;
    
end


h = trisurf(Ttot, xall, yall, zall, fall,...
        'EdgeColor', 'black', 'EdgeAlpha', 0.1, ...
        'AmbientStrength', 0.6, 'DiffuseStrength', 0.4, ...
        'SpecularStrength', 0.3, varargin_use{:});


view(3)
shading interp
axis equal
% axis vis3d    % might be useful
grid on

if nargout == 1
    varargout{1} = h;
end

end
