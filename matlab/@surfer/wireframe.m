function wireframe(S,opts)
% WIREFRAME  Plot surface edges of patches of a surfer object in current figure
%
% wireframe(S) shows all patches, plotting in color the curvature function on nodes.
%
% wireframe(S,opts) where opts is a struct with options
%
%   opts.nsign - plot edges at location r + 5e-3 * nsign * n. (default 1)
%       to avoid z fighting. 
%   opts.wfill - if true, also plot surface in white
%       (default false)
%
% Some view and properties of the current axes are set for convenience of
%  3D viewing.
%
% Example
%   S = geometries.sphere(1.0);
%   wireframe(S)
%
% See also PLOT

if nargin < 2
    opts = [];
end

norder_iptype = zeros(S.npatches,2);
norder_iptype(:,1) = S.norders;
norder_iptype(:,2) = S.iptype;
[no_ip_uni,~,iuni] = unique(norder_iptype,'rows');
nuni = size(no_ip_uni,1);
xinterp = cell(nuni,1);

for i = 1:nuni
    iptype = no_ip_uni(i,2);
    norder = no_ip_uni(i,1);
if iptype == 1
    ts = linspace(0,1,norder+1);
    uvs1 = [ts;0*ts];
    uvs2 = [flip(ts);1-flip(ts)];
    uvs3 = [0*ts;flip(ts)];
    [uvst] = [uvs1,uvs2,uvs3];
    [uvs] = koorn.rv_nodes(norder);
    
    amat = koorn.vals2coefs(norder,uvs);
    pols = koorn.pols(norder,uvst);
    xinterp{i} = pols.'*amat;
elseif iptype == 11
    ts = linspace(-1,1,norder+1);
    [uvs] = polytens.lege.nodes(norder);
    uvs1 = [ts; 0*ts-1];
    uvs2 = [1+0*ts;ts];
    uvs3 = [flip(ts); 0*ts+1];
    uvs4 = [-1+0*ts;flip(ts)];
    [uvst] = [uvs1,uvs2,uvs3,uvs4];
    amat = polytens.lege.vals2coefs(norder,uvs);
    pols = polytens.lege.pols(norder,uvst);
    xinterp{i} = pols.'*amat;
elseif iptype == 12
    ts = linspace(-1,1,norder+1);
    [uvs] = polytens.cheb.nodes(norder);
    uvs1 = [ts; 0*ts-1];
    uvs2 = [1+0*ts;ts];
    uvs3 = [flip(ts); 0*ts+1];
    uvs4 = [-1+0*ts;flip(ts)];
    [uvst] = [uvs1,uvs2,uvs3,uvs4];
    amat = polytens.cheb.vals2coefs(norder,uvs);
    pols = polytens.cheb.pols(norder,uvst);
    xinterp{i} = pols.'*amat;
end

end

nsign = 1;
if isfield(opts, 'nsign')
    nsign = opts.nsign;
end

rwire_plot = [];
for i = 1:S.npatches
    iinds = S.ixyzs(i):S.ixyzs(i+1)-1;
rwire = xinterp{iuni(i)}*S.r(:,iinds).';
nwire = xinterp{iuni(i)}*S.n(:,iinds).';
rwire = rwire + nsign*5e-3*nwire;
rwire = [rwire;NaN*ones(3,size(rwire,2))];
rwire = reshape(rwire,[],3).';
rwire_plot = [rwire_plot,rwire];
end

plot3(rwire_plot(1,:),rwire_plot(2,:),rwire_plot(3,:),'k-',LineWidth=1)

wfill = 0;
if isfield(opts,'wfill')
    wfill = opts.wfill;
end

if wfill
hold on
h = plot(S,NaN*S.r(1,:));
h.FaceColor = 'w';
h.AmbientStrength= 0.6;
h.DiffuseStrength= 0.4;
h.SpecularStrength= 0.3;

hold off
end
end