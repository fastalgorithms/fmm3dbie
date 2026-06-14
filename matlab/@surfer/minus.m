function objout = minus(obj, v)
% MINUS Translate a surfer object by subtracting a constant vector.
%
% S2 = S - v translates the surfer S by -v, i.e. r -> r - v.
%
% See also AFFINE_TRANSF, PLUS, TRANSLATE, MTIMES

if ~isa(obj, 'surfer')
    % v - S: treat as -(S - v)
    objout = affine_transf(v, -eye(3), obj(:));
    return
end

if ~(isa(v, 'numeric') && numel(v) == 3)
    error('SURFER:minus:invalid', ...
        'minus only supported for a surfer object minus a 3-element vector.');
end

objout = affine_transf(obj, eye(3), -v(:));

end
