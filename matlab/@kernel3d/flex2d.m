function obj = flex2d(type, zk, nu, pfun)
%KERNEL3D.FLEX2D   Construct a flexural (thin-plate) kernel in 2D.
%
%   KERNEL3D.FLEX2D(type,zk), KERNEL3D.FLEX2D(type,zk,nu) or
%   KERNEL3D.FLEX2D(type,zk,nu,pfun), for a flat surfer in z = 0, where
%   zk = [zk1, zk2] (scalar zk: [zk, 1i*zk]; 0: biharmonic), nu is the
%   Poisson ratio and type is one of:
%      's'                   - G,                         [1 1]
%      'sp'                  - d/dn_x G,                  [1 1]
%      'sg'                  - [d/dx G; d/dy G],          [2 1]
%      'clamped_plate_bcs'   - [G; d/dn_x G],             [2 1]
%      'supported_plate_bcs' - [G; M_nn + nu*M_tt],       [2 1]
%      'free_plate_bcs'      - [M_nn + nu*M_tt; shear],   [2 1]
%      'flex_varalpha'       - L G, see below,            [1 1]
%                              (alias 'varcoef')
%
%   flex_varalpha is L applied to G at the target,
%      L G = Delta(alpha Delta G) - beta G - (1-nu)(alpha_xx G_yy
%            - 2 alpha_xy G_xy + alpha_yy G_xx) = alpha delta(x-y) + K,
%   with alpha delta(x-y) in obj.diag. pfun(t) returns alpha (1,nt),
%   dalpha (2,nt), d2alpha (3,nt: xx, xy, yy) and beta (1,nt), see
%   FLEX2D.PLATE_COEFS.
%
%   Returns a KERNEL3D.RADGREEN kernel on FLEX2D.GREEN (flex_varalpha: on
%   G and on Lap G), with an fmm if hfmm2d and rfmm2d are on the path
%   (none for free).
%
% See also FLEX2D.GREEN, KERNEL3D.RADGREEN

if ( nargin < 2 )
    error('KERNEL3D.FLEX2D: requires type and zk arguments.');
end
if ( nargin < 3 ), nu = []; end
if ( nargin < 4 ), pfun = []; end

zk = zk(:).';
if ( isscalar(zk) && abs(zk) >= 1e-6 )
    zk = [zk, 1i*zk];
end

switch lower(type)
    case {'s', 'single', 'slp'}
        rtype = 's';
    case {'sp', 'sprime'}
        rtype = 'sp';
    case {'sg', 'sgrad'}
        rtype = 'grad';
    case {'clamped_plate_bcs', 'clamped'}
        rtype = 'clamped_bc';
    case {'supported_plate_bcs', 'supported'}
        rtype = 'supported_bc';
    case {'free_plate_bcs', 'free'}
        rtype = 'free_bc';
    case {'flex_varalpha', 'varalpha', 'varcoef', 'var'}
        rtype = 'flex_varalpha';
        if ( isempty(nu) || isempty(pfun) )
            error('KERNEL3D.FLEX2D: ''flex_varalpha'' requires nu and pfun.');
        end
    otherwise
        error('KERNEL3D.FLEX2D: unknown plate kernel type ''%s''.', type);
end
if ( any(strcmp(rtype, {'supported_bc', 'free_bc'})) && isempty(nu) )
    error('KERNEL3D.FLEX2D: ''%s'' requires nu.', type);
end

hasfmm = exist('hfmm2d', 'file') && exist('rfmm2d', 'file');
fmmg = {};  fmmh = {};
if ( hasfmm )
    fg = @(eps,s,pg,t,pgt) flex2d.fmm(eps, zk, s, pg, t, pgt);
    fh = @(eps,s,pg,t,pgt) flex2d.fmm(eps, zk, s, pg, t, pgt, true);
    fmmg = {'fmm', @(varargin) fmm2dto3d(fg, varargin{:})};
    fmmh = {'fmm', @(varargin) fmm2dto3d(fh, varargin{:})};
end
g = @(s,t) flex2d.green(zk, s, t);

if ( strcmp(rtype, 'flex_varalpha') )
    % G up to second derivatives, plus grad of H = Lap G
    args = {'dim', 2, 'zk', zk};
    KG = kernel3d(g, @(t) wgreen(t, nu, pfun, zk), 'nderiv', 2, args{:}, fmmg{:});
    KH = kernel3d(@(s,t) lapgreen(zk, s, t), @(t) wlap(t, nu, pfun, zk), ...
                  'nderiv', 1, args{:}, fmmh{:});
    obj = KG + KH;
    obj.name = 'flex2d';
    obj.type = 'flex_varalpha';
    obj.diag = @(t) reshape(pfun(t).alpha, [], 1);
    return
end

args = {'dim', 2, 'zk', zk};
if ( ~isempty(nu) ), args = [args, {'nu', nu}]; end
if ( ~strcmp(rtype, 'free_bc') ), args = [args, fmmg]; end
obj = kernel3d(g, rtype, args{:});

end


function W = wgreen(t, nu, pfun, zk)
% weights on [G | G_x G_y | G_xx G_xy G_yy], see FLEX2D.PLATE_COEFS
c = flex2d.plate_coefs(pfun(t), nu, zk).';
z = zeros(size(c,1), 1);
W = [c(:,7), z, z, c(:,3) + c(:,5), c(:,6), c(:,3) + c(:,4)];
end


function W = wlap(t, nu, pfun, zk)
% weights on [H | H_x H_y], H = Lap G
c = flex2d.plate_coefs(pfun(t), nu, zk).';
W = [zeros(size(c,1), 1), c(:,1), c(:,2)];
end


function varargout = lapgreen(zk, s, t)
% H = Lap G and its gradient
out = cell(1, 4);
[out{:}] = flex2d.green(zk, s, t);
varargout{1} = out{3}(:,:,1) + out{3}(:,:,3);
varargout{2} = cat(3, out{4}(:,:,1) + out{4}(:,:,3), out{4}(:,:,2) + out{4}(:,:,4));
end
