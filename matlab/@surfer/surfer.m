classdef surfer
%SURFER class which describes a surface divided into triangles/quads (patches). 
% On each patch the surface is represented by the values of its position, 
% its two tangential derivatives, and the corresponding basis function
% expansions on either an RV grid on triangles or a tensor product GL/Chebyshev
% grid on quads
%
% To see fields and available methods: doc surfer

    
    properties
        iptype        % type of each patch (1,11,12,..)
        weights       % cell array of quadrature weight arrays per patch
        norders       % expansion order of each patch
        
        npatches      % number of patches
        npts          % total number of quadrature nodes
        srccoefs      % cell array of orthog poly coeffs of [r;du;dv] per patch
    end
    
   	properties(SetAccess=private)
        r             % quadrature node locations (3*npts)
        du            % dr/du tangential vectors at nodes (3*npts)
        dv            % dr/dv tangential vectors at nodes (3*npts)
        dru           % normalized dr/du vectors at nodes (3*npts)
        drv           % normalized dr/dv vectors at nodes (3*npts)
        n             % unit outward normals at nodes (3*npts)
        wts           % quadrature weights (surface element) (npts*1)
        patch_id      % which patch each node belongs to (npts*1)
        uvs_targ      % (u,v) param coords of all nodes (2*npts)
        curv
        ffform
        ffforminv
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
    % create a surfer object
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
            rnodes = cell(nuni,1);
            rwts = cell(nuni,1);
            umats = cell(nuni,1);
            npols = cell(nuni,1);
            dumats = cell(nuni,1);
            dvmats = cell(nuni,1);
            for i=1:nuni
                ip0 = no_ip_uni(i,2);
                if(ip0 == 1)
                    rnodes{i} = koorn.rv_nodes(no_ip_uni(i,1));
                    rwts{i} = koorn.rv_weights(no_ip_uni(i,1));
                    umats{i} = koorn.vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = koorn.ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
                elseif(ip0==11)
                    rnodes{i} = polytens.lege_nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.lege_weights(no_ip_uni(i,1));
                    umats{i} = polytens.lege_vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = polytens.lege_ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
                elseif(ip0==12)
                    rnodes{i} = polytens.cheb_nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.cheb_weights(no_ip_uni(i,1));
                    umats{i} = polytens.cheb_vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = polytens.cheb_ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
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
            obj.patch_id = zeros(obj.npts,1);
            obj.uvs_targ = zeros(2,obj.npts);

            ffform = zeros(2,2,obj.npts);
            sfform = zeros(2,2,obj.npts);
            ffforminv = zeros(2,2,obj.npts);
            obj.curv = zeros(obj.npts,1);
            obj.ffforminv = cell(npatches,1);
            obj.ffform = cell(npatches,1);
            
            for i=1:npatches
                iind = obj.ixyzs(i):(obj.ixyzs(i+1)-1);
                obj.srcvals{i} = srcvals(:,iind);
                obj.srccoefs{i} = obj.srcvals{i}(1:9,:)*umats{iuse(i)}';

                du = obj.srcvals{i}(4:6,:);
                dv = obj.srcvals{i}(7:9,:);
                dn = obj.srcvals{i}(10:12,:);

                dxucoefs = obj.srccoefs{i}(4:6,:);
                dxvcoefs = obj.srccoefs{i}(7:9,:);

                dxuu = dxucoefs*dumats{iuse(i)};
                dxuv = 0.5*(dxucoefs*dvmats{iuse(i)} + dxvcoefs*dumats{iuse(i)});
                dxvv = dxvcoefs*dvmats{iuse(i)};

                L = dxuu(1,:).*dn(1,:) + dxuu(2,:).*dn(2,:) + dxuu(3,:).*dn(3,:);
                M = dxuv(1,:).*dn(1,:) + dxuv(2,:).*dn(2,:) + dxuv(3,:).*dn(3,:);
                N = dxvv(1,:).*dn(1,:) + dxvv(2,:).*dn(2,:) + dxvv(3,:).*dn(3,:);

                dunormsq = sum(du.*du,1);
                dvnormsq = sum(dv.*dv,1);
                duv = sum(du.*dv,1);

                ddinv = 1.0./(dunormsq.*dvnormsq - duv.*duv);

                ffform(1,1,iind) = dunormsq;
                ffform(1,2,iind) = duv;
                ffform(2,1,iind) = duv;
                ffform(2,2,iind) = dvnormsq;

                sfform(1,1,iind) = L;
                sfform(1,2,iind) = M;
                sfform(2,1,iind) = M;
                sfform(2,2,iind) = N;
                
                ffforminv(1,1,iind) = dvnormsq.*ddinv;
                ffforminv(1,2,iind) = -duv.*ddinv;
                ffforminv(2,1,iind) = -duv.*ddinv;
                ffforminv(2,2,iind) = dunormsq.*ddinv;

                obj.curv(iind) = -0.5*(L.*dvnormsq - ...
                   2*M.*duv + dunormsq.*N).*ddinv;

                obj.ffform{i} = ffform(:,:,iind);
                obj.ffforminv{i} = ffforminv(:,:,iind);
                
                da = vecnorm(cross(du,dv),2).*rwts{iuse(i)};
                obj.weights{i} = da(:);      
                obj.patch_id(iind) = i;
                obj.uvs_targ(1:2,iind) = rnodes{iuse(i)}; 
            end
            
            sv = [obj.srcvals{:}];
            obj.wts = cat(1,obj.weights{:});
            obj.r  = sv(1:3,:);
            obj.du = sv(4:6,:);
            obj.dv = sv(7:9,:);
            obj.n  = sv(10:12,:);
            obj.dru = obj.du./repmat(vecnorm(obj.du,2,1),[3,1]);
            drv = cross(obj.n,obj.du);
            obj.drv = drv./repmat(vecnorm(drv,2,1),[3,1]);
        end
        
        
         [varargout] = plot(obj,varargin);
         a = area(obj);
         [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(obj);
         [objout,varargout] = oversample(obj,novers);
         [objout,varargout] = affine_transf(obj,mat,shift);
         [varargout] = scatter(obj,varargin);
         [spmat] = conv_rsc_to_spmat(obj,row_ptr,col_ind,wnear);
         [objout,varargout] = rotate(obj, eul);
         
        
    end

    methods(Static)
        obj = load_from_file(fname,varargin);
        [obj,varargout] = surfacemesh_to_surfer(dom);
    end
end
