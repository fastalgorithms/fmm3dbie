function xvals = nodes(norder,kind)
    kind0 = 1;
    if(nargin == 2)
        kind0 = kind;
    end
    opts = [];
    opts.kind = kind0;
    [x] = polytens.cheb.pts(norder+1,opts);
    [yy,xx] = meshgrid(x);
    xx = xx(:);
    yy = yy(:);
    xvals = zeros(2,(norder+1)*(norder+1));
    xvals(1,:) = xx;
    xvals(2,:) = yy;
end
