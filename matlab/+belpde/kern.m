function submat = kern(zk2,srcinfo,targinfo,type)
%BELPDE.KERN parametrix and remainder kernels for PDEs on surfaces
%
% Syntax: submat = belpde.kern(zk2,srcinfo,targinfo,type)
%
% Let x be targets and y be sources for these formulas, with
% n_x the unit normal at the target and H_x the mean curvature there.
%
% Kernels based on K(x,y) = log(|x-y|)/(2*pi), or the 1/4i*H_0^(1)(|x-y),
% and the corresponding remainders
%
% R(x,y) = \Delta_\Gamma K(x,y)
%
% Input:
%   zk2 - complex number, Helmholtz wave number squared, ignored for
%               Laplace-Beltrami
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (3,:) array
%                ptinfo.du - first derivative with respect to u in
%                     underlying parameterization (3,:)
%                ptinfo.dv - first derivative with respect to v in
%                     underlying parameterization (3,:)
%                ptinfo.n - normals (3,:)
%   targinfo - description of targets in ptinfo struct format,
%                the remainder kernels require normal info in
%                targinfo.n and mean curvature in targinfo.mean_curv
%   type - string, determines kernel type
%                type == 'klb', Laplace-Beltrami parametrix
%                type == 'rlb', Laplace-Beltrami remainder
%                type == 'khb', Helmholtz-Beltrami parametrix
%                type == 'rhb', Helmholtz-Beltrami remainder
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets.

src = srcinfo.r;
targ = targinfo.r;

dx = targ(1,:)' - src(1,:);
dy = targ(2,:)' - src(2,:);
dz = targ(3,:)' - src(3,:);

switch lower(type)

  case {'klb'}
    dr2 = dx.^2 + dy.^2 + dz.^2;
    submat = log(dr2)/2/(2*pi);

  case {'rlb'}
    dr2 = dx.^2 + dy.^2 + dz.^2;
    targnorm = targinfo.n;
    mean_curv = targinfo.mean_curv(:);

    tmp = (targnorm(1,:)'.*dx+targnorm(2,:)'.*dy+targnorm(3,:)'.*dz)./dr2;
    submat = (2*tmp.*tmp - 2.*mean_curv.*tmp)/(2*pi);

  case {'khb'}
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    cdr = sqrt(zk2).*dr;

    submat = -1i*besselh(0,cdr)/4;

  case {'rhb'}
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    cdr = sqrt(zk2).*dr;
    targnorm = targinfo.n;
    mean_curv = targinfo.mean_curv(:);

    tmp = (targnorm(1,:)'.*dx+targnorm(2,:)'.*dy+targnorm(3,:)'.*dz)./dr;

    tmp2 = -zk2.*besselh(2,cdr).*tmp.*tmp;
    submat = (tmp2+sqrt(zk2).*besselh(1,cdr).*2.*mean_curv.*tmp)/4i;

  otherwise
    error('BELPDE.KERN: unknown kernel type ''%s''.', type);

end

end
