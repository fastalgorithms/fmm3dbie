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
    
    has_geom = any([has_targ_field(targinfo,'du'), ...
        has_targ_field(targinfo,'dv'), has_targ_field(targinfo,'n')]);
    has_tau = has_targ_field(targinfo,'tau');
    has_d = has_targ_field(targinfo,'d');
    if sum([has_geom, has_tau, has_d]) > 1
        error('FMM3DBIE:extract_targ_array:overlap', ...
            'Target fields du/dv/n, tau, and d use overlapping rows.');
    end

    if has_targ_field(targinfo,'du')
        nddu = size(targinfo.du(:,:),1);
        targs(4+(0:nddu-1),:) = targinfo.du(:,:);
    end  
    if has_targ_field(targinfo,'dv')
        nddv = size(targinfo.dv(:,:),1);
        targs(7+(0:nddv-1),:) = targinfo.dv(:,:);
    end
    if has_targ_field(targinfo,'n')
        ndn = size(targinfo.n(:,:),1);
        targs(10+(0:ndn-1),:) = targinfo.n(:,:);
    end
    if has_targ_field(targinfo,'tau')
        ndtau = size(targinfo.tau(:,:),1);
        targs(4+(0:ndtau-1),:) = targinfo.tau(:,:);
    end            
    if has_targ_field(targinfo,'d')
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

function tf = has_targ_field(targinfo, name)
tf = (isstruct(targinfo) && isfield(targinfo,name)) || isprop(targinfo,name);
end
