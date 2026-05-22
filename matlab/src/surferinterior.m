function in = surferinterior(S,targs, opts)
% Identify which targets are inside a surfer. Assumes normals are pointed
% outwards. 
%
% If the surface is interesting, use opts.quadrature = true (slow, but accurate)

if nargin < 3
    opts = [];
end

quadrature = false;
if isfield(opts,'quadrature')
    quadrature = opts.quadrature;
end

targ = targs.r;

eps = 1e-4;
% try
%     if quadrature, error('use quadrature');end
%     pg = 0;
%     pgt = 1;
% 
%     srcinfo.sources = S.r;
%     srcinfo.dipoles = S.n.*S.wts(:).';
% 
%     U = lfmm3d(eps,srcinfo,pg,targ,pgt);
% 
%     us = U.pottarg / 4/pi;
% 
% catch
%     warning('fmm3d not detected, surferinterior will be slow')
    dens = ones(S.npts,1);
    us = lap3d.dirichlet.eval(S,dens,targs,eps,[0,1]);
% end


in = -us;
ifix = (abs(in - round(in))>0.1) | abs(in)>1.5;
if any(ifix)
    targfix = [];
    targfix.r = targ(:,ifix);
    [sxyz, patch_inds, uvsloc, dists, flags] = get_closest_pts(S, targfix);
    % if any(flags>0) || any(~(dists<1e2)), error('failed to identify point region'), end
    if any(~(dists<1e2)), warning('failed to identify point region'), end
    ibad = ~(dists<1e2);


    ns = 0*targfix.r;
    for i = 1:length(targfix)
        ipatch = patch_inds(i);
        iptype = S.iptype(ipatch);
        if iptype == 1
        pols = koorn.pols(S.norders(ipatch),uvsloc(:,i));
        elseif iptype == 11
        pols = polytens.lege.pols(S.norders(ipatch), uvsloc(:,i));
        elseif iptype == 12
        pols = polytens.cheb.pols(S.norders(ipatch), uvsloc(:,i));
        else
        error('unsupported patch type')
        end
        du = S.srccoefs{ipatch}(4:6,:)*pols;
        dv = S.srccoefs{ipatch}(7:9,:)*pols;
        rtmp = cross(du,dv);
        ns(:,i) = rtmp./norm(rtmp);
    end

 d = targfix.r - sxyz;
 signs = sign(dot(d,ns));
 in(ifix) = (1-signs)/2;
 in(ibad) = -1;
end
in = logical(round(in));
end
