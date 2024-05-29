function [E, H] = incoming_sources(zk, srcinfo, targinfo, type, varargin)
%EM3D.INCOMING_FIELDS generates incoming E and H fields
% for electric, and magnetic dipoles
%
% Syntax: [E,H] = em3d.incoming_sources(zk, srcinfo, targinfo, 
%                    'electric dipole')
%         [E,H] = em3d.incoming_sources(zk, srcinfo, targinfo, 
%                    'magnetic dipole')
%
% Electric dipole, due to a point source at the origin
% E_{ed} = - \nabla \times \nabla \times G_{k}(r) edip
% H_{ed} =  i k \nabla \times G_{k}(r) edip
%
% Magnetic dipole, due to a point source at the origin
% E_{md} = i k \nabla \times G_{k}(r) hdip
% H_{md} = \nabla \times \nabla \times G_{k}(r) hdip
%
% where G_{k}(r) is the Helmholtz Green's function = exp(ikr)/(4\pi r)
%
% Input:
%   zk       - complex number, Maxwell wavenumber
%   srcinfo  - description of sources in ptinfo struct format
%               ptinfo.r - positions (3,:) array
%               ptinfo.edips - electric dipole strengths (3,:) 
%                 complex array (must be provided if type is ed or ehd)
%               ptinfo.hdips - magnetic dipole strengths (3,:)
%                 complex array (must be provided if type is hd or ehd)
%   targinfo - description of targets in ptinfo struct format,
%               only ptinfo.r needed
%   type     - string, determine source type
%                 type == 'ed' or 'electric dipole', electric dipole
%                 type == 'hd' or 'magnetic dipole', magnetic dipole
%                 type == 'ehd' or 'hed' or 
%                 'electric and magnetic dipole', or
%                 'magnetic and electric dipole', generate fields
%                 due to both electric and magneetic dipoles
%
%   varargin - currently unused
%
% Output:
%   E,H - electric and magnetic fields at the target locations
%

targs = targinfo.r;
src = srcinfo.r;

[~, nt] = size(targs);
[~, ns] = size(src);

xs = repmat(src(1,:), nt, 1);
ys = repmat(src(2,:), nt, 1);
zs = repmat(src(3,:), nt, 1);

xt = repmat(targs(1,:).', 1, ns);
yt = repmat(targs(2,:).', 1, ns);
zt = repmat(targs(3,:).', 1, ns);

E = complex(zeros(3,nt));
H = complex(zeros(3,nt));

rx = xt - xs;
ry = yt - ys;
rz = zt - zs;

r = sqrt(rx.^2 + ry.^2 + rz.^2);
rinv = 1./r;
zexp = exp(1i*zk*r)/4/pi;
ztmp1 = (1j*zk*rinv.^2 - rinv.^3);
ztmp2 = (-zk^2*rinv.^3 - 3*1j*zk*rinv.^4 + 3*rinv.^5);
ztmp3 = (zk^2.*rinv + ztmp1).*zexp;

ztmp1_use = 1j*zk*zexp.*ztmp1;

if strcmpi(type, 'ed') || strcmpi(type, 'electric dipole') || ...
        strcmpi(type, 'ehd') || strcmpi(type, 'hed') || ...
        strcmpi(type, 'electric and magnetic dipole') || ...
        strcmpi(type, 'magnetic and electric dipole')
    edips = srcinfo.edips;
    edipsx = edips(1,:);
    edipsy = edips(2,:);
    edipsz = edips(3,:);

    ztmp = ztmp2.*zexp;
    zxx = ztmp.*rx.*rx;
    zxy = ztmp.*rx.*ry;
    zxz = ztmp.*rx.*rz;

    zyy = ztmp.*ry.*ry;
    zyz = ztmp.*ry.*rz;

    zzz = ztmp.*rz.*rz;
    
    Ex = -(ztmp3*edipsx(:) + zxx*edipsx(:) + zxy*edipsy(:) + ...
                       zxz*edipsz(:));
    Ey = -(ztmp3*edipsy(:) + zxy*edipsx(:) + zyy*edipsy(:) + ...
                       zyz*edipsz(:));
    Ez = -(ztmp3*edipsz(:) + zxz*edipsx(:) + zyz*edipsy(:) + ...
                       zzz*edipsz(:));
    E(1,:) = E(1,:) + Ex.'; 
    E(2,:) = E(2,:) + Ey.';
    E(3,:) = E(3,:) + Ez.';

    zx = ztmp1_use.*rx;
    zy = ztmp1_use.*ry;
    zz = ztmp1_use.*rz;

    H(1,:) = H(1,:) + (zy*edipsz(:) - zz*edipsy(:)).';
    H(2,:) = H(2,:) + (zz*edipsx(:) - zx*edipsz(:)).';
    H(3,:) = H(3,:) + (zx*edipsy(:) - zy*edipsx(:)).';
end



if strcmpi(type, 'hd') || strcmpi(type, 'magnetic dipole') || ...
        strcmpi(type, 'ehd') || strcmpi(type, 'hed') || ...
        strcmpi(type, 'electric and magnetic dipole') || ...
        strcmpi(type, 'magnetic and electric dipole')
    hdips = srcinfo.hdips;
    hdipsx = hdips(1,:);
    hdipsy = hdips(2,:);
    hdipsz = hdips(3,:);

    ztmp = ztmp2.*zexp;
    zxx = ztmp.*rx.*rx;
    zxy = ztmp.*rx.*ry;
    zxz = ztmp.*rx.*rz;

    zyy = ztmp.*ry.*ry;
    zyz = ztmp.*ry.*rz;

    zzz = ztmp.*rz.*rz;
    
    Hx = (ztmp3*hdipsx(:) + zxx*hdipsx(:) + zxy*hdipsy(:) + ...
                       zxz*hdipsz(:));
    Hy = (ztmp3*hdipsy(:) + zxy*hdipsx(:) + zyy*hdipsy(:) + ...
                       zyz*hdipsz(:));
    Hz = (ztmp3*hdipsz(:) + zxz*hdipsx(:) + zyz*hdipsy(:) + ...
                       zzz*hdipsz(:));
    H(1,:) = H(1,:) + Hx.';
    H(2,:) = H(2,:) + Hy.';
    H(3,:) = H(3,:) + Hz.';

    zx = ztmp1_use.*rx;
    zy = ztmp1_use.*ry;
    zz = ztmp1_use.*rz;


    E(1,:) = E(1,:) + (zy*hdipsz(:) - zz*hdipsy(:)).';
    E(2,:) = E(2,:) + (zz*hdipsx(:) - zx*hdipsz(:)).';
    E(3,:) = E(3,:) + (zx*hdipsy(:) - zy*hdipsx(:)).';
end
