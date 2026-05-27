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
    has_n  = isfield(targinfo,'n') || isprop(targinfo,'n') ;
    has_du = isfield(targinfo,'du') || isprop(targinfo,'du') ;
    has_dv = isfield(targinfo,'dv') || isprop(targinfo,'dv') ;
    targs(1:3,:) = targinfo.r;
   if has_du,  targs(4:6,:)   = targinfo.du;  end
   if has_dv,  targs(7:9,:)   = targinfo.dv;  end
   if has_n,   targs(10:12,:) = targinfo.n;   end

    if size(targs,1)> 3 && size(targs,1) < 12
        targs = [targs;zeros(12-size(targs,1),size(targs,2))];
    end

end
