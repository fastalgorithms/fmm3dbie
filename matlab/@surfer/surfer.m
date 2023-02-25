classdef surfer
%SURFER class which describes a surface divided into triangles/quads (patches). 
% On each patch the surface is represented by the values of its position, 
% its two tangential derivatives, and the corresponding basis function
% expansions on either an RV grid on triangles or a tensor product GL/Chebyshev
% grid on quads
%


    
    properties 
        srcvals
        srccoefs
        iptype
        weights
        norders
        curv
        ffform
        ffforminv
        npatches
        npts
    end
    
    properties (Hidden)
        ixyzs
        ifcurv
        ifffform
        
       
    end
    
    methods
        function obj = surfer(npatches,norder,srcvals)
% This needs fixing particularly handling the varargin part            
            obj.npatches = npatches;
            obj.norders = ones(npatches,1)*norder;
            npols = (norder+1)*(norder+2)/2;
            
            obj.npts = npols*npatches;
            obj.srcvals = cell(npatches,1);
            obj.srccoefs = cell(npatches,1);
            obj.weights = cell(npatches,1);
            rnodes = koorn.rv_nodes(norder);
            rwts = koorn.rv_weights(norder);
            umat = koorn.vals2coefs(norder,rnodes);
            obj.iptype = ones(npatches,1);
            for i=1:npatches
                
                obj.srcvals{i} = srcvals(:,((i-1)*npols+1):i*npols);
                obj.srccoefs{i} = obj.srcvals{i}(1:9,:)*umat';
                ru = obj.srcvals{i}(4:6,:);
                rv = obj.srcvals{i}(7:9,:);
                
                da = vecnorm(cross(ru,rv),2).*rwts;
                obj.weights{i} = da(:);
                
                
            end
            ixyzs = 1:npols:obj.npts+1;
            obj.ixyzs = ixyzs(:);
            obj.ifcurv = 0;
            obj.ffform = 0;
            obj.ffform = 0;
            obj.curv = cell(npatches,1);
            obj.ffforminv = cell(npatches,1);
            obj.ffform = cell(npatches,1);
            
        end
        
        
         [varargout] = plot(obj,varargin);
         a = area(obj);
        
    end
    methods(Static)
        obj = load_from_file(fname,varargin)
        
    end
end
