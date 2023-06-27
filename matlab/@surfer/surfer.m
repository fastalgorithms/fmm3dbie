classdef surfer
%SURFER class which describes a surface divided into triangles/quads (patches). 
% On each patch the surface is represented by the values of its position, 
% its two tangential derivatives, and the corresponding basis function
% expansions on either an RV grid on triangles or a tensor product GL/Chebyshev
% grid on quads
%


    
    properties
        iptype
        weights
        norders
        curv
        ffform
        ffforminv
        npatches
        npts
        srccoefs
    end
    
   	properties(SetAccess=private)
        r
        du
        dv
        n
    end
    properties (Access = private)
        srcvals
        
    end    
    properties (Hidden)
        ixyzs
        ifcurv
        ifffform
        
       
    end
    
    methods
        function obj = surfer(npatches,norders,srcvals,iptype)
% This needs fixing particularly handling the varargin part, 
% and also dealing with quad patches in the future
            obj.npatches = npatches;
            obj.srcvals = cell(npatches,1);
            obj.srccoefs = cell(npatches,1);
            obj.weights = cell(npatches,1);
            if (nargin == 3)
                obj.iptype = ones(npatches,1);
            else
                if(length(iptype)==1)
                    obj.iptype = iptype*ones(npatches,1);
                elseif(length(iptype)==npatches)
                    obj.iptype = iptype;
                else
                    fprintf('Incompatible size of iptype, returning\n');
                    return;     
                end
                
            end
            
            if(length(norders)==1)
                
                obj.norders = ones(npatches,1)*norders;
            elseif(length(norders) == npatches)
                obj.norders = norders;
            else
                fprintf('Incompatible size of norders, returning\n');
                return;     
            end
                
% extract unique set of interpolation matrices
            norder_iptype = zeros(npatches,2);
            norder_iptype(:,1) = obj.norders;
            norder_iptype(:,2) = obj.iptype;
            [no_ip_uni,~,iuse] = unique(norder_iptype,'rows');
            nuni = size(no_ip_uni,1);
            rwts = cell(nuni,1);
            umats = cell(nuni,1);
            npols = cell(nuni,1);
            for i=1:nuni
                ip0 = no_ip_uni(i,2);
                if(ip0 == 1)
                    rnodes = koorn.rv_nodes(no_ip_uni(i,1));
                    rwts{i} = koorn.rv_weights(no_ip_uni(i,1));
                    umats{i} = koorn.vals2coefs(no_ip_uni(i,1),rnodes); 
                    npols{i} = size(rnodes,2);
                elseif(ip0==11)
                    rnodes = polytens.lege_nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.lege_weights(no_ip_uni(i,1));
                    umats{i} = polytens.lege_vals2coefs(no_ip_uni(i,1),rnodes);
                    npols{i} = size(rnodes,2);
                elseif(ip0==12)
                    rnodes = polytens.cheb_nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.cheb_weights(no_ip_uni(i,1));
                    umats{i} = polytens.cheb_vals2coefs(no_ip_uni(i,1),rnodes);
                    npols{i} = size(rnodes,2);
                else
                    fprintf('Invalid type of patch, exiting\n');
                    return
                end
            end
            
            npts_per_patch = zeros(npatches,1);
            for i=1:npatches
                npts_per_patch(i) = npols{iuse(i)};
            end
            npts_per_patch = [1;npts_per_patch];
            obj.ixyzs = cumsum(npts_per_patch);
            
            obj.npts = obj.ixyzs(end)-1;
            
            for i=1:npatches
                iind = obj.ixyzs(i):(obj.ixyzs(i+1)-1);
                obj.srcvals{i} = srcvals(:,iind);
                obj.srccoefs{i} = obj.srcvals{i}(1:9,:)*umats{iuse(i)}';
                du = obj.srcvals{i}(4:6,:);
                dv = obj.srcvals{i}(7:9,:);
                da = vecnorm(cross(du,dv),2).*rwts{iuse(i)};
                obj.weights{i} = da(:);         
            end
            obj.ifcurv = 0;
            obj.ffform = 0;
            obj.ffform = 0;
            obj.curv = cell(npatches,1);
            obj.ffforminv = cell(npatches,1);
            obj.ffform = cell(npatches,1);
            
            sv = [obj.srcvals{:}];
            obj.r  = sv(1:3,:);
            obj.du = sv(4:6,:);
            obj.dv = sv(7:9,:);
            obj.n  = sv(10:12,:);
        end
        
        
         [varargout] = plot(obj,varargin);
         a = area(obj);
         [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(obj);
         [objout,varargout] = oversample(obj,novers);
         [objout,varargout] = affine_transf(obj,mat,shift);
         [varargout] = scatter(obj,varargin);
         [spmat] = conv_rsc_to_spmat(obj,row_ptr,col_ind,wnear);
         
        
    end

    methods(Static)
        obj = load_from_file(fname,varargin);
        obj = sphere(varargin);
        obj = sphere_quad(varargin);
        obj = ellipsoid(varargin);
        obj = axissym(fcurve,cparams,varargin);
        obj = surface_hps_mesh_to_surfer(dom);
    end
end
