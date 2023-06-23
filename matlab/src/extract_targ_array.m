function [targs] = extract_targ_array(targinfo)
%
%  extract_targ_array
%     subroutine to extract flattened target array
%     from its struct format
%
%  Syntax
%     targs = extract_targ_arrays(targinfo)
%  
%  Input arguments
%     targinfo - targinfo struct
%
%  Output arguments
%     targs (ndtarg,ntarg) - flattened target array with
%         ndtarg = 3/12, 3 if only targinfo.r is present
%         and 12 if r,du,dv,and n are present in the struct
%
    [~,ntarg] = size(targinfo.r);
    if(isfield(targinfo,'du') && isfield(targinfo,'dv') && isfield(targinfo,'n'))
       ndtarg = 12;
       targs = zeros(ndtarg,ntarg);
       targs(1:3,:) = targinfo.r;
       targs(4:6,:) = targinfo.du;
       targs(7:9,:) = targinfo.dv;
       targs(10:12,:) = targinfo.n;
    else
       ndtarg = 3;
       targs = targinfo.r;
    end

end
