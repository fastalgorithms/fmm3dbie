function S = multiscale_mesher(fnamein, norder, opts)
% 
%  MULTISCALE_MESHER creates a smooth high order triangulated surface 
%    based on an input mesh file containing first/second order
%    triangles/quads. If the mesh comprises of quads, then each quad
%    is split into two triangles. 
%
%  Supported mesh formats include, .gidmsh, .tri, .msh, gmshv2, 
%  and gmshv4
%  
%  Syntax
%    S = multiscale_mesher(fnamein, norder)
%    S = multiscale_mesher(fnamein, norder, opts)
%
%  Input arguments:
%    * fnamein: input mesh file name
%    * norder: order of discretization for smooth surface
%    * opts: options struct (optional)
%        opts.nquad (12), quadrature order for computing level set
%        opts.rlam (10), smoothing parameter (should be between 2.5 and 10)
%        opts.adapt_sigma (1), adapt_sigma = 0, uses uniform sigma for mollifier
%                              adapt_sigma = 1, uses adaptive sigma
%        opts.nrefine (0), number of refinements
%        opts.filetype, type of file
%          filetype = 1, for .msh from gidmsh
%          filetype = 2, for .tri
%          filetype = 3, for .gidmsh
%          filetype = 4, for .msh gmsh v2
%          filetype = 5, for .msh gmsh v2
%                               
%      
%  
    d = dir(fnamein);
    if isempty(d)
        error('MULTISCALE_MESHER: invalid file\n');
    end
    fnameuse = fullfile(d.folder, d.name);

    if nargin < 3
        opts = [];
    end


    fnameoutuse = fullfile(d.folder, '/tmp');
    
    norder_skel = 12;
    if isfield(opts, 'nquad')
        norder_skel = opts.nquad;
    end

    norder_smooth = norder;
    nrefine = 0;
    if isfield(opts, 'nrefine')
        nrefine = opts.nrefine;
    end


    adapt_sigma = 1;
    if isfield(opts, 'adapt_sigma')
        adapt_sigma = opts.adapt_sigma;
    end
    rlam = 10;
    if isfield(opts, 'rlam')
        rlam = opts.rlam;
    end

    if isfield(opts, 'filetype')
        ifiletype = opts.filetype;
    else
        ier = 0;
        ifiletype = 0;
        mex_id_ = 'MWF77_get_filetype(i cstring[x], io int[x], io int[x])';
[ifiletype, ier] = fmm3dbie_routs(mex_id_, fnameuse, ifiletype, ier, 1000, 1, 1);
        if ier > 0
            error('MULTISCALE_MESHER: error determining file type\n');
        end
    end

    if norder > 20 || norder < 1
        error('MULTISCALE_MESHER: norder too high, must be less than 20');
    end
    
    if norder_skel > 20 || norder < 1
        error('MULTISCALE_MESHER: opts.nquad too high, must be less than 20');
    end
    ier = 0;
    mex_id_ = 'MWF77_multiscale_mesher(i cstring[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x], i cstring[x], io int[x])';
[ier] = fmm3dbie_routs(mex_id_, fnameuse, ifiletype, norder_skel, norder_smooth, nrefine, adapt_sigma, rlam, fnameoutuse, ier, 1000, 1, 1, 1, 1, 1, 1, 1000, 1);

    if ier == 3 || ier == 4
        error_message = ['MULTISCALE_MESHER: error in main smoothing routine', ...
                'try a larger value of rlam (if not set, use opts.rlam and set it to ',...
                'a value greater than 10)\n' ...
                'if that does not work, then mesh cannot be smooth with the surface smoother\n'];
        error(error_message);
    end

    S = cell(1,nrefine+1); 
    for i=0:nrefine
        fname = fullfile(d.folder, ['/tmp_o' num2str(norder,'%02.f') '_r' num2str(i,'%02.f') '.go3']);
        S{i+1} = surfer.load_from_file(fname);
        delete(fname);
    end
end
%
%%   Common routines
%%
%
%-------------------------------------------------
