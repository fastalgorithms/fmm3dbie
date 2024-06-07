function wts = lege_weights(norder)
    [~,w] = polytens.lege.pts(norder+1);
    wts = w'*w;
    wts = wts(:).'; 
end
