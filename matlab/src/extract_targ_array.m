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
%         ndtarg = 3/12/13/14, 3 if only targinfo.r is present
%          12 if r,du,dv,and n are present in the struct
%          13 if kappa (curvature) is present
%          14+ if extra data is present
%
    [~,ntarg] = size(targinfo.r);
    targs = targinfo.r;    

    if isfield(targinfo,'du')
    targs(4:6,:) = targinfo.du;
    end  
    if isfield(targinfo,'dv')
    targs(7:9,:) = targinfo.dv;
    end
    if isfield(targinfo,'n')
    targs(10:12,:) = targinfo.n;
    end
    if isfield(targinfo,'tau')
    targs(4:6,:) = targinfo.tau;
    end            
    if isfield(targinfo,'kappa')
    targs(13,:) = targinfo.kappa;
    end 
    if isfield(targinfo,'data')
    targs(14:13+size(targinfo.data,1),:) = targinfo.data;
    end 

    ndtarg = size(targs,1);

end