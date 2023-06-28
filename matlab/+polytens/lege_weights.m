function wts = lege_weights(norder)
    [~,w] = legpts(norder+1);
    wts = w'*w;
    wts = wts(:).'; 
end