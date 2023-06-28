function xvals = cheb_nodes(norder,kind)
    kind0 = 1;
    if(nargin == 2)
        kind0 = kind;
    end
    [x] = chebpts(norder+1,[-1,1],kind0);
    [yy,xx] = meshgrid(x);
    xx = xx(:);
    yy = yy(:);
    xvals = zeros(2,(norder+1)*(norder+1));
    xvals(1,:) = xx;
    xvals(2,:) = yy;
end