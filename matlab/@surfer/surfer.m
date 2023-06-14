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
        function obj = surfer(npatches,norders,srcvals)
% This needs fixing particularly handling the varargin part, 
% and also dealing with quad patches in the future
            obj.npatches = npatches;
            obj.srcvals = cell(npatches,1);
            obj.srccoefs = cell(npatches,1);
            obj.weights = cell(npatches,1);
            obj.iptype = ones(npatches,1);
            
            if(length(norders)==1)
                
                obj.norders = ones(npatches,1)*norders;
                
                rwts = cell(1);
                umats = cell(1);
                rnodes = koorn.rv_nodes(norders);
                rwts{1} = koorn.rv_weights(norders);
                umats{1} = koorn.vals2coefs(norders,rnodes);
                
                iuse = ones(npatches,1);
                
                
                
            elseif(length(norders)==npatches)
                obj.norders = norders;
                [norders_uni,~,iuse] = unique(norders);
                nuni = length(norders_uni);
                rwts = cell(nuni,1);
                umats = cell(nuni,1);
                for i=1:nuni
                    rnodes = koorn.rv_nodes(norders_uni(i));
                    rwts{i} = koorn.rv_weights(norders_uni(i));
                    umats{i} = koorn.vals2coefs(norders_uni(i),rnodes);    
                end
            else
                fprintf('Incompatible size of norders, returning\n');
                return;     
            end
            npts_per_patch = (obj.norders+1).*(obj.norders+2)/2;
            npts_per_patch = [1;npts_per_patch];
            obj.ixyzs = cumsum(npts_per_patch);
            
            obj.npts = obj.ixyzs(end)-1;
            
            for i=1:npatches
                iind = obj.ixyzs(i):(obj.ixyzs(i+1)-1);
                obj.srcvals{i} = srcvals(:,iind);
                obj.srccoefs{i} = obj.srcvals{i}(1:9,:)*umats{iuse(i)}';
                du = obj.srcvals{i}(4:6,:);
                dv = obj.srcvals{i}(7:9,:);
                %du = obj.du{i}(:,:);
                %dv = obj.dv{i}(:,:);
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
        obj = ellipsoid(varargin);
        obj = axissym(fcurve,cparams,varargin);
    end
end
