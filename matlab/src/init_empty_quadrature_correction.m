function Q = init_empty_quadrature_correction(targinfo,opts)
    if nargin < 2
        opts = [];
    end
    Q = [];
    [~,ntarg] = size(targinfo.r);
    Q.targinfo = targinfo;
    Q.row_ptr = ones(ntarg+1,1);
    Q.nnz = 0;
    Q.col_ind = 0;
    ifcomplex = 0;

    nker = 1;
    if (isfield(opts,'nker'))
      nker = opts.nker;
    end

    Q.wnear = zeros(nker,1);
    if(isfield(opts,'type'))
       if(strcmpi(opts.type,'complex'))
          Q.wnear = complex(zeros(nker,1));
       end
    end
    Q.iquad = [1,1];
    Q.nquad = 1;
    if(isfield(opts,'rfac'))
        Q.rfac = opts.rfac;
    end
    Q.kernel_order = -1;
    Q.wavenumber = 0;
    Q.rfac = 3;
    Q.format = 'rsc';
end
