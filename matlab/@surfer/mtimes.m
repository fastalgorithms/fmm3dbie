function objout = mtimes(A, obj)
% MTIMES Apply a 3x3 matrix transformation to a surfer object.
%
% S2 = A * S applies the linear map r -> A*r to the surfer S, updating
%  node positions, tangent vectors, and normals. A must be a 3x3 matrix.
%  A scalar A is interpreted as A*eye(3).
%
% See also AFFINE_TRANSF, PLUS, MINUS, SCALE, ROTATE

if isa(A, 'numeric') && isa(obj, 'surfer')
    if isscalar(A)
        objout = scale(obj, A);
        return
    end
    [m, n] = size(A);
    if m == 3 && n == 3
        objout = affine_transf(obj, A);
    else
        error('SURFER:mtimes:invalid', ...
            'Matrix must be 3x3 to transform a surfer object.');
    end
elseif isa(obj, 'numeric') && isa(A, 'surfer')
    % Allow S * A by delegating to A' * S
    objout = mtimes(obj.', A.');
else
    error('SURFER:mtimes:invalid', ...
        'mtimes requires a 3x3 numeric matrix and a surfer object.');
end

end
