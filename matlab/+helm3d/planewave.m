function [u, gradu] = planewave(zk, dir, targinfo)
%HELM3D.PLANEWAVE  Potential and gradient due to plane wave.
%
%  [u, gradu] = planewave(zk, dir, targinfo)
%
% Input arguments:
%   zk - wave number
%   dir(2) or dir(3) - incident plane wave direction
%     if it is a two vector then dir = (theta, phi), 
%     d = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
%
%     if ir is a three vector then
%     d = (dir(1), dir(2), dir(3))/norm(dir) 
%   targinfo - either (3,n) array or a struct
%      with field targinfo.r which is a (3,n) array
% Outputs:
%  u - potential (N,1) at each point in targinfo.r
%  gradu - gradient of potential at each point in targinfo.r
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

    if length(dir(:)) == 2
       thet = dir(1);
       phi = dir(2);
       d = [sin(thet)*cos(phi),; sin(thet)*sin(phi); cos(thet)];
    elseif length(dir(:)) == 3
       rn = sqrt(dir(1)^2 + dir(2)^2 + dir(3)^3);
       d = dir(:);
       d = d/rn;
    else
      error('HELM3D.planewave: invalid size of direction vector');
    end
    uexp = exp(1i*zk*sum(d.*targs));
    
    u = uexp.';
    if nargout >= 2
        gradu = complex(zeros(3,n2));
        gradu(1,:) = 1i*zk*d(1)*uexp;
        gradu(2,:) = 1i*zk*d(2)*uexp;
        gradu(3,:) = 1i*zk*d(3)*uexp;
    end
    

end
