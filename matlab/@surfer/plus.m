function objout = plus(obj, v)
% PLUS Translate a surfer object by a constant vector.
%
% S2 = S + v translates the surfer S by the 3x1 (or 3-element) vector v,
%  i.e. r -> r + v.
%  Also supports v + S.
%
% See also AFFINE_TRANSF, MINUS, TRANSLATE, MTIMES

if ~isa(obj, 'surfer')
    % v + S: swap and recurse
    objout = plus(v, obj);
    return
end

if ~(isa(v, 'numeric') && numel(v) == 3)
    error('SURFER:plus:invalid', ...
        'plus only supported for a surfer object and a 3-element vector.');
end

objout = affine_transf(obj, eye(3), v(:));

end
