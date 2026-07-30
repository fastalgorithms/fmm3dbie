function v = volume(obj)
% VOLUME Signed volume enclosed by surfer object (divergence theorem).
%
% v = volume(S)
%
% Uses the divergence theorem with vector field F = (x,y,z)/3,
% giving v = integral over S of (nx*x + ny*y + nz*z)/3 dA.

    vfrm_dn = (obj.n(1,:).*obj.r(1,:) + obj.n(2,:).*obj.r(2,:) + obj.n(3,:).*obj.r(3,:)).' / 3;
    v = sum(vfrm_dn .* obj.wts);
end
