function [kernels, coefs, diagonal] = radpp_sum_parts(kernels, coefs)
%RADPP_SUM_PARTS Separate the kernel and local terms of a radial sum.
if nargin < 2, coefs = ones(1,numel(kernels)); end
kernels = reshape(kernels,1,[]);
coefs = reshape(coefs,1,[]);
assert(numel(kernels) == numel(coefs), 'Need one coefficient per kernel.');
assert(all(cellfun(@(K) isequal(K.opdims,kernels{1}.opdims),kernels)), ...
       'kernel3d dimensions must agree to add');

diagonal = [];
if nargout > 2
    hasdiag = cellfun(@(K) ~isempty(K.diag),kernels) & coefs ~= 0;
    if any(hasdiag)
        terms = cellfun(@(K) K.diag,kernels(hasdiag),'UniformOutput',false);
        diagonal = radpp_sum(coefs(hasdiag),terms{:});
    end
end

neutral = cellfun(@(K) strcmp(K.name,'eye') || K.iszero,kernels);
kernels = kernels(~neutral);
coefs = coefs(~neutral);
for k = 1:numel(kernels)
    kernels{k}.diag = [];
end
end
