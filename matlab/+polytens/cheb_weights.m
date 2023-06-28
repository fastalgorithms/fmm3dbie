function wts = cheb_weights(norder)
    [~,w] = chebpts(norder+1,[-1,1],1);
    wts = w'*w;
    wts = wts(:).'; 
end