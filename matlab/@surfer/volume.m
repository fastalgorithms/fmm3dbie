function v = volume(obj)
% VOLUME Signed volume enclosed by surfer object (divergence theorem).
%
% v = volume(S)
%
% Uses the divergence theorem with vector field F = (x,y,0)/2,
% giving v = integral over S of (nx*x + ny*y)/2 dA.

    vfrm_dn = (obj.n(1,:).*obj.r(1,:) + obj.n(2,:).*obj.r(2,:)).' / 2;
    v = sum(vfrm_dn .* obj.wts);
end
