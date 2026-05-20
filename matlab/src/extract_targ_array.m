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
    
    targs = targinfo.r(:,:);    
    

    if isfield(targinfo,'du') || isprop(targinfo,'du')
        nddu = size(targinfo.du(:,:),1);
        targs(4+(0:nddu-1),:) = targinfo.du(:,:);
    end  
    if isfield(targinfo,'dv') || isprop(targinfo,'dv')
        nddv = size(targinfo.dv(:,:),1);
        targs(7+(0:nddv-1),:) = targinfo.dv(:,:);
    end
    if isfield(targinfo,'n') || isprop(targinfo,'n')
        ndn = size(targinfo.n(:,:),1);
        targs(10+(0:ndn-1),:) = targinfo.n(:,:);
    end
    if isfield(targinfo,'tau') || isprop(targinfo,'tau')
        ndtau = size(targinfo.tau(:,:),1);
        targs(4+(0:ndtau-1),:) = targinfo.tau(:,:);
    end            
    if isfield(targinfo,'d') || isprop(targinfo,'d')
        ndd = size(targinfo.d(:,:),1);
        targs(4+(0:ndd-1),:) = targinfo.d(:,:);
    end            
    if mod(size(targs,1),3) == 2
        targs = [targs;zeros(1,size(targs,2))]; 
    end

    if isfield(targinfo,'kappa')
    targs(13,:) = targinfo.kappa;
    end 
    if isfield(targinfo,'data')
    targs(14:13+size(targinfo.data,1),:) = targinfo.data;
    end 

    ndtarg = size(targs,1);

end