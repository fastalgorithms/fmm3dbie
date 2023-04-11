function Q = init_empty_quadrature_correction(targinfo,opts)
    Q = [];
    [~,ntarg] = size(targinfo.xyzs);
    Q.targinfo = targinfo;
    Q.row_ptr = ones(ntarg+1,1);
    Q.nnz = 0;
    Q.col_ind = 0;
    ifcomplex = 0
    Q.wnear = 0;
    if(isfield(opts,'type'))
       if(strcmpi(opts.type,'complex'))
          Q.wnear = 0+0j;
       end
    end
    Q.iquad = 1;
    Q.nquad = 1;
    if(isfield(opts,'rfac'))
        Q.rfac = opts.rfac;
    end
    Q.kernel_order = -1;
    Q.wavenumber = 0;
end
