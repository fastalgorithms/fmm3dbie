% @SURFER
%
% Files
%   affine_transf         - Affine transformation of surfer object.
%   area                  - Surface area of surfer object (sum of quadrature weights).
%   conv_rsc_to_spmat     - Make sparse matrix from patchwise sparse format.
%   extract_arrays        - extract and concatenate panel arrays from surfer object.
%   load_from_file        - Create surfer object by loading a .go3 file
%   merge                 - Combine array of surfer objects into a single such object.
%   oversample            - Creates surfer object with patches of higher order than given.
%   plot                  - Plot surface patches of a surfer object in current figure
%   plot_nodes            - show nodes of a surfer object, or function on nodes
%   rotate                - Return a surfer object rotated by Euler angles
%   scale                 - Create a rescaled copy of a surfer object.
%   scatter               - Plots values of a surfer object, with colors specified by data.
%   surf_fun_error        - estimate pointwise approximation error of func on surface.
%   surfacemesh_to_surfer - Create a surfer object from structured quad array
%   surfer                - class which describes a surface divided into triangles/quads (patches).
%   translate             - translate a surfer object by given vector(s)
%   vals2coefs            - convert values to coeffs on given surfer object
