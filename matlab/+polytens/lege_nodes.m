function xvals = lege_nodes(norder)
    [x] = legpts(norder+1);
    [yy,xx] = meshgrid(x);
    xx = xx(:);
    yy = yy(:);
    xvals = zeros(2,(norder+1)*(norder+1));
    xvals(1,:) = xx;
    xvals(2,:) = yy;
end