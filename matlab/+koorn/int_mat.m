function [xintmatu,xintmatv,xinterp] = int_mat(uvs,norder,nleg,pt)

    xintmatu = zeros(size(uvs,2),(norder+1)*(norder+2)/2);
    xintmatv = zeros(size(uvs,2),(norder+1)*(norder+2)/2);
    xinterp = zeros(size(uvs,2),(norder+1)*(norder+2)/2);
    uvskoorn     = koorn.rv_nodes(norder);
    amat = koorn.vals2coefs(norder,uvskoorn);
    np = size(uvs,2);

    [xlege,wlege] = lege.exps(nleg);

    for ii=1:np

    up = (xlege+1)/2*(uvs(1,ii)-pt(1))+pt(1);
    vp = (xlege+1)/2*(uvs(2,ii)-pt(2))+pt(2);
    wl = wlege/2;
    pp =wl.'*koorn.pols(norder,[up.';vp.']).'*amat;
    xintmatu(ii,:) = pp*(uvs(1,ii)-pt(1));
    xintmatv(ii,:) = pp*(uvs(2,ii)-pt(2));
    
    pol = koorn.pols(norder,uvs(:,ii)).'*amat;
    xinterp(ii,:) = pol;
    end

end
