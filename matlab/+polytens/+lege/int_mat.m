function [xintmatu,xintmatv,xinterp] = int_mat(uvs,norder,nleg,pt)

    xintmatu = zeros(size(uvs,2),(norder+1)*(norder+1));
    xintmatv = zeros(size(uvs,2),(norder+1)*(norder+1));
    xinterp = zeros(size(uvs,2),(norder+1)*(norder+1));
    uvs_lege     = polytens.lege.nodes(norder);
    amat = polytens.lege.vals2coefs(norder,uvs_lege);
    np = size(uvs,2);

    [xlege,wlege] = polytens.lege.rts(nleg);

    for ii=1:np

    up = (xlege+1)/2*(uvs(1,ii)-pt(1))+pt(1);
    vp = (xlege+1)/2*(uvs(2,ii)-pt(2))+pt(2);
    wl = wlege/2;
    pp =wl.'*polytens.lege.pols(norder,[up.';vp.']).'*amat;
    xintmatu(ii,:) = pp*(uvs(1,ii)-pt(1));
    xintmatv(ii,:) = pp*(uvs(2,ii)-pt(2));

    pol = polytens.lege.pols(norder,uvs(:,ii)).'*amat;
    xinterp(ii,:) = pol;
    end

end
