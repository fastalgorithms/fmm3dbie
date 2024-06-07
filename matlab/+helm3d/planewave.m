function [u, gradu] = planewave(zk, dir, targinfo)
%HELM3D.PLANEWAVE  Potential and gradient due to plane wave.
%
%  [u, gradu] = planewave(zk, dir, targinfo)
%
% Input arguments:
%   zk - wave number
%   dir(3) - incident planewave direction (should be a unit vector)
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
    uexp = exp(1i*zk*sum(dir(:).*targs));
    
    u = uexp.';
    if nargout >= 2
        gradu = complex(zeros(3,n2));
        gradu(1,:) = 1i*zk*dir(1)*uexp;
        gradu(2,:) = 1i*zk*dir(2)*uexp;
        gradu(3,:) = 1i*zk*dir(3)*uexp;
    end
    

end
