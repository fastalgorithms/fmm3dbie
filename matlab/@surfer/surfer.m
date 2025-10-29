classdef surfer
%SURFER class which describes a surface divided into triangles/quads (patches).
%
% On each patch the surface is represented by the values of its position, 
% its two tangential derivatives, and the corresponding basis function
% expansions on either an RV grid on triangles or a tensor product GL/Chebyshev
% grid on quads.
%
% Notes:
%
% Cell arrays have one element per patch. For scalar quantities per patch,
% a plain array is used.
% 
% The types for each patch are:
%        * iptype = 1,  triangular patch discretized using
%                       Rokhlin-Vioreanu nodes
%        * iptype = 11, quadrangular patch discretized using tensor
%                       product Gauss-Legendre nodes
%        * iptype = 12, quadrangular patch discretized using tensor
%                       product Chebyshev
% For order p=norder, a triangle patch uses (p+1)(p+2)/2 nodes, whereas
%  either quad patch type uses p^2 nodes. The number of coeffs equals this
%  number of nodes.
%
% Surfer properties:
%  obj.iptype      - type of each patch
%  obj.weights     - cell array of quadrature weight arrays per patch
%  obj.norders     - expansion order p of each patch
%  obj.npatches    - numbers of patches
%  obj.npts        - total number of discretization points
%  obj.srcccoefs   - cell array of orthogonal polynomial coefs of
%                       [r; du; dv] per patch
%  obj.r           - discretization node locations (3,npts)
%  obj.du          - dr/du tangent vectors at nodes (3,npts)
%  obj.dv          - dr/dv tangent vectors at nodes (3,npts)
%  obj.n           - unit outward normals at nodes (3,npts)
%  obj.wts         - quadrature weights for integrating smooth functions
%                    on surface (surface element) (npts,1)
%  obj.patch_id    - which patch each node belongs to (npts,1)
%  obj.uvs_targ    - (u,v) param coords of nodes within own patch (2,npts)
%  obj.mean_curv   - mean curvatures at nodes (npts*1)
%  obj.ffform      - cell array of first fundamental forms at nodes (2,2,n)
%  obj.ffforminv   - cell array of inverses of ffforms at nodes (2,2,n)
%
%
% Surfer methods
%   plot(obj, varargin)      - plot the surface, by default plots the mean
%                              curvature
%   scatter(obj,s,c, ...)    - scatter plot of surface
%   plot_nodes(obj,v, ...)   - plot discretization nodes describing the
%                              surface
%   extract_arrays(obj)      - extract flattened arrays for srcvals,
%                              srccoefs, norders, ixyzs, iptype, and wts
%   oversample(obj, novers)  - oversample patch(i) on the surface
%                              to expansion order novers(i)
%   affine_transf(obj,mat,s) - apply affine transformation, and translation
%   conv_rsc_spmat(obj, ...) - convert row sparse representation of 
%                              near quadrature to sparse matrix
%                              representation
%   rotate(obj, eul)         - rotate surface based on euler angles
%   scale(obj, sf)           - scale object
%   interpolate_data(obj,...)- interpolate density to uv_targs and
%                              patch_ids
%   translate(obj, r)        - translate object
%   merge([array of objs])   - merge an array of surface objects
%   area(obj)                - compute the surface area of object
%   surf_fun_error(obj,f,p)  - estimate error in function approximation
%                              on surface via basis expansion tails
%   vals2coefs(obj, vals)    - construct basis function expansions
%                              of function on surface
% author:
    
    properties
        iptype        % type of each patch (integer, 1,11,12,...)
        weights       % cell array of quadrature weight arrays per patch
        norders       % expansion order p of each patch (length npatches)
        npatches      % number of patches
        npts          % total number of discretization points
        srccoefs      % cell array of orthog poly coeffs of [r;du;dv] per patch
    end
    
   	properties(SetAccess=private)
        r             % discretization node locations (3,npts)
        du            % dr/du tangent vectors at nodes (3,npts)
        dv            % dr/dv tangent vectors at nodes (3,npts)
        dru           % normalized dr/du vectors at nodes (3,npts)
        drv           % normalized dru \times n at nodes (3,npts)
        n             % unit outward normals at nodes (3,npts)
        wts           % quadrature weights (surface element) (npts*1)
        patch_id      % which patch each node belongs to (npts*1)
        uvs_targ      % (u,v) param coords of nodes within own patch (2,npts)
        mean_curv     % mean curvatures at nodes (npts*1)
        ffform        % cell array of first fundamental forms at nodes (2,2,npts)
        ffforminv     % cell array of inverses of ffforms at nodes (2,2,npts)
        aspect_ratio  % measure of local aspect ratio of patch, (npts,1);
        patch_distortion % integral of aspect ratio
        cms           % centroid location of patches (3,npatches)
        rads          % radii of bounding (npatches,1)
        end_pt_verts  % vertices of patch 
    end
    properties (Access = private)
        srcvals      % cell array of nodes vals of [r;du;dv;n] per patch
        
    end    
    properties (Hidden)
        ixyzs
        ifcurv
        ifffform
        
       
    end
    
    methods
    function obj = surfer(npatches, norders, srcvals, iptype)
    % SURFER Create a surfer object.
    %
    % S = surfer(npatches,norders,srcvals,iptype) creates a surfer object
    %  with npatches patches, of types iptype (integer vector of length
    %  npatches) and orders norders (integer vector of length npatches),
    %  using the 12*npts array srcvals, where npts must match the expected
    %  total number of nodes in all the patches (as determined from their
    %  types and orders). If norders or iptype have length=1, the same value
    %  is applied to all patches.
    %
    % Note: under the hood this recreates all surface node and patch info
    %  from srcvals, duplicating code in the fortran library. Its format is
    %  srcvals = [r; du; dv; n] where
    %             r is (3,npts) node coords
    %             du is (3,npts) tangent dr/du coords
    %             dv is (3,npts) tangent dr/dv coords
    %             n is (3,npts) unit outward normal coords
      
% This needs fixing particularly handling the varargin part, 
% and also dealing with quad patches in the future
            obj.npatches = npatches;
            obj.srcvals = cell(npatches,1);
            obj.srccoefs = cell(npatches,1);
            obj.weights = cell(npatches,1);
            obj.end_pt_verts = cell(npatches,1);
            if (nargin == 3)
                obj.iptype = ones(npatches,1);
            else
                if isscalar(iptype)
                    obj.iptype = iptype*ones(npatches,1);
                elseif(length(iptype)==npatches)
                    obj.iptype = iptype;
                else
                    fprintf('Incompatible size of iptype, returning\n');
                    return;     
                end
                
            end
            
            if isscalar(norders)
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
            % matrices for interpolating to patch vertices
            p_vert_mats = cell(nuni,1); 
            for i=1:nuni
                ip0 = no_ip_uni(i,2);
                if(ip0 == 1)
                    rnodes{i} = koorn.rv_nodes(no_ip_uni(i,1));
                    rwts{i} = koorn.rv_weights(no_ip_uni(i,1));
                    umats{i} = koorn.vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = koorn.ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
                    epts = [0,0,1;0,1,0];
                    p_vert_mats{i} = koorn.pols(no_ip_uni(i,1), epts);
                elseif(ip0 == 11)
                    rnodes{i} = polytens.lege.nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.lege.weights(no_ip_uni(i,1));
                    umats{i} = polytens.lege.vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = polytens.lege.ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
                    epts = [-1,1,1,-1;-1,-1,1,1];
                    p_vert_mats{i} = polytens.lege.pols(no_ip_uni(i,1), epts);
                elseif(ip0 == 12)
                    rnodes{i} = polytens.cheb.nodes(no_ip_uni(i,1));
                    rwts{i} = polytens.cheb.weights(no_ip_uni(i,1));
                    umats{i} = polytens.cheb.vals2coefs(no_ip_uni(i,1),rnodes{i});
                    [~, dumats{i}, dvmats{i}] = polytens.cheb.ders(no_ip_uni(i,1), rnodes{i});
                    npols{i} = size(rnodes{i},2);
                    epts = [-1,1,1,-1;-1,-1,1,1];
                    p_vert_mats{i} = polytens.cheb.pols(no_ip_uni(i,1), epts);
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
            obj.mean_curv = zeros(obj.npts,1);
            obj.ffforminv = cell(npatches,1);
            obj.ffform = cell(npatches,1);
            obj.aspect_ratio = zeros(obj.npts,1);
            obj.patch_distortion = zeros(obj.npatches,1);
            obj.cms = zeros(3,obj.npatches);
            obj.rads = zeros(obj.npatches,1);
            
            
            for i=1:npatches
                iind = obj.ixyzs(i):(obj.ixyzs(i+1)-1);
                obj.srcvals{i} = srcvals(:,iind);
                obj.srccoefs{i} = obj.srcvals{i}(1:9,:)*umats{iuse(i)}';

                du = obj.srcvals{i}(4:6,:);
                dv = obj.srcvals{i}(7:9,:);
                dn = obj.srcvals{i}(10:12,:);


                dunormsq = sum(du.*du,1);
                dvnormsq = sum(dv.*dv,1);
                duv = sum(du.*dv,1);
         
                ddinv = 1.0./(dunormsq.*dvnormsq - duv.*duv);
                tt = sqrt((dunormsq - dvnormsq).^2 + 4.*(duv.^2));
                
                % compute aspect ratio
                 
                asp_rat = (dunormsq + dvnormsq + tt).'./ ...
                          (dunormsq + dvnormsq - tt).';  
                obj.aspect_ratio(iind) = asp_rat;

                % compute first fundamental form and its inverse
                ffform(1,1,iind) = dunormsq;
                ffform(1,2,iind) = duv;
                ffform(2,1,iind) = duv;
                ffform(2,2,iind) = dvnormsq;
                obj.ffform{i} = ffform(:,:,iind);

                ffforminv(1,1,iind) = dvnormsq.*ddinv;
                ffforminv(1,2,iind) = -duv.*ddinv;
                ffforminv(2,1,iind) = -duv.*ddinv;
                ffforminv(2,2,iind) = dunormsq.*ddinv;
                obj.ffforminv{i} = ffforminv(:,:,iind);

                % compute second fundamental form
                dxucoefs = obj.srccoefs{i}(4:6,:);
                dxvcoefs = obj.srccoefs{i}(7:9,:);

                dxuu = dxucoefs*dumats{iuse(i)};
                dxuv = 0.5*(dxucoefs*dvmats{iuse(i)} + dxvcoefs*dumats{iuse(i)});
                dxvv = dxvcoefs*dvmats{iuse(i)};

                L = dxuu(1,:).*dn(1,:) + dxuu(2,:).*dn(2,:) + dxuu(3,:).*dn(3,:);
                M = dxuv(1,:).*dn(1,:) + dxuv(2,:).*dn(2,:) + dxuv(3,:).*dn(3,:);
                N = dxvv(1,:).*dn(1,:) + dxvv(2,:).*dn(2,:) + dxvv(3,:).*dn(3,:);

                sfform(1,1,iind) = L;
                sfform(1,2,iind) = M;
                sfform(2,1,iind) = M;
                sfform(2,2,iind) = N;

                % Mean curvature (LG-2MF+NE)/2(EG-F^2):
                obj.mean_curv(iind) = -0.5*(L.*dvnormsq - ...
                                       2*M.*duv + dunormsq.*N).*ddinv;

                % compute centroid and bounding sphere radii
                everts = obj.srccoefs{i}(1:3,:)*p_vert_mats{iuse(i)};
                obj.end_pt_verts{i} = everts;
                obj.cms(1:3,i) = mean(everts, 2);

                rr = sqrt((obj.cms(1,i) - everts(1,:)).^2 + ... 
                          (obj.cms(2,i) - everts(2,:)).^2 + ...
                          (obj.cms(3,i) - everts(3,:)).^2);
                obj.rads(i) = max(rr);

                
                
                
                
                da = vecnorm(cross(du,dv),2).*rwts{iuse(i)};
                obj.weights{i} = da(:);      
                obj.patch_id(iind) = i;
                obj.uvs_targ(1:2,iind) = rnodes{iuse(i)}; 

                rsum = sum(da(:));
                rp = sum(da(:).*asp_rat);
                obj.patch_distortion(i) = rp/rsum;
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
         dens_int = interpolate_data(obj, dens, ipatch_ids, uvs_targ);
         [objout,varargout] = affine_transf(obj,mat,shift);
         [varargout] = scatter(obj,varargin);
         [spmat] = conv_rsc_to_spmat(obj,row_ptr,col_ind,wnear);
         [objout,varargout] = rotate(obj, eul);
         [varargout] = plot_nodes(S, v, varargin);
         [objout,varargout] = scale(obj,sf);
         [objout,varargout] = translate(obj, r);
         [objout] = merge(Sarray);
         [coefs] = vals2coefs(obj,vals);
         [errps] = surf_fun_error(obj,fun,p);    
        
    end

    methods(Static)
        obj = load_from_file(fname,varargin);
        [obj,varargout] = surfacemesh_to_surfer(dom,varargin);
    end
end
