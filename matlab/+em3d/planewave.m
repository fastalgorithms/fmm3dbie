function [E, H] = planewave(zk, dir, pol, targinfo)
%EM3D.planewave returns the fields due to plane waves
%
% Input arguments:
%   zk - wave number
%   dir(2) - incident planewave direction in polar coordinates
%            (\theta_d, \phi_d)
%   pol(2) - polarization vector 
%   targinfo - either (3,n) array or a struct
%      with field targinfo.r which is a (3,n) array
%
%   E = (\that pol(1) + \phat pol(2)) exp(ik x \rhat)
%   H = (-that pol(2) + \phat pol(1)) exp(ik x \rhat)
%   
    if isa(targinfo, 'double')
      targs = targinfo; 
    elseif isa(targinfo, 'struct')  || isa(targinfo, 'surfer')
      targs = targinfo.r; 
    end

    [n1, n2] = size(targs);
    if n1~=3
      error('HELM3D.planewave: invalid size of target array');
    end

    thet = dir(1);
    phi = dir(2);

    rhat = [sin(thet)*cos(phi); sin(thet)*sin(phi); cos(thet)];
    that = [cos(thet)*cos(phi); cos(thet)*sin(phi); -sin(thet)];
    phat = [-sin(phi); cos(phi); 0];
    uexp = exp(1i*zk*sum(rhat(:).*targs));
    
    E = (pol(1)*that + pol(2)*phat)*uexp;
    H = (-pol(2)*that + pol(1)*phat)*uexp;

end
