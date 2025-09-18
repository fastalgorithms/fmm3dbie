function lapn = lap_n_op(obj)

    n = obj.n;
    [lapnx] = get_surface_laplacian(obj,n(1,:));
    [lapny] = get_surface_laplacian(obj,n(2,:));
    [lapnz] = get_surface_laplacian(obj,n(3,:));

    lapn = [lapnx;lapny;lapnz];

end

function c = pmt(a,b)
    c = pagemtimes(a,b);
end