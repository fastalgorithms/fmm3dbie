function lapn = lap_n_op(obj)
%  lapn = shape_op(S,p)
%     This subroutine evaluates the surface Laplacian of the surface normal
%     vector
%
%  Input arguments:
%    * obj: surfer object
%
%  Output arguments:
%    * lapn: double (3,S.npts) 
%
    n = obj.n;
    [lapnx] = get_surface_laplacian(obj,n(1,:));
    [lapny] = get_surface_laplacian(obj,n(2,:));
    [lapnz] = get_surface_laplacian(obj,n(3,:));

    lapn = [lapnx;lapny;lapnz];
end
