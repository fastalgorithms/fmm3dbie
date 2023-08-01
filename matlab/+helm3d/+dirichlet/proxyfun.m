function [Kpxy,lst] = proxyfun(x,slf,lst,proxy_dict,l,ctr,zpars,rn,wts)
    proxy = proxy_dict.proxy;
    weigt = l*proxy_dict.weigt;
    norms = proxy_dict.norms;
    % Shift and scale precomputed proxy surface
    pxy = bsxfun(@plus, proxy*l, ctr');
    
    % set up sourceinfo and targetinfo
    srcinfo = [];
    srcinfo.r = x(:,slf);
    srcinfo.n = rn(:,slf);
    w_sqrt = sqrt(wts(slf).');
    
    targinfo = [];
    targinfo.r = pxy;
    targinfo.n = norms;
    
    Kpxy1 = bsxfun(@times,helm3d.kern(zpars(1),srcinfo, ...
      targinfo,'s'),w_sqrt);
    Kpxy1 = bsxfun(@times, weigt.', Kpxy1);
    Kpxy2 = bsxfun(@times,helm3d.kern(zpars(1),srcinfo, ...
      targinfo,'c',zpars(2),zpars(3)),w_sqrt);
    Kpxy2 = bsxfun(@times, weigt.', Kpxy2);
    
    Kpxy = [Kpxy1;Kpxy2];
    dxyz = abs(x(1:3,lst)-ctr(1:3)')/l;
    lst = lst(max(dxyz) < 2.5);

end