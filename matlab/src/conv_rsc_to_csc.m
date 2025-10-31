function [col_ptr, row_ind, iper] = conv_rsc_to_csc(npatches,row_ptr,col_ind)
% reorder rsc format to csc format
% row_ind(col_ptr(i)):(row_ind(col_ptr(i+1))-1) are the indices of targets near patch i

    ntarg = size(row_ptr,1)-1;
    nnz = length(col_ind);
    [jsort, iper] = sort(col_ind);

    row_ind_exp = zeros(nnz,1);
    for i = 1:ntarg
        row_ind_exp(row_ptr(i):(row_ptr(i+1)-1)) = i;
    end
    row_ind = row_ind_exp(iper);

    col_ptr = zeros(npatches+1,1);
    
    icur = 1;
    col_ptr(1) = 1;
    for i = 1:nnz
        if (jsort(i) > icur) 
            icur0 = icur;
            icur = jsort(i);
            col_ptr((icur0+1):icur) = i;
        end
    end
    col_ptr((icur+1):end) = nnz+1;
end