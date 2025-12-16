function varargout = dealm(x)

n = size(x, 2);
if ( n == 1 )
    n = size(x, 1);
    x = x(:).';
end

if ( nargout > n )
    error('Too many outputs.');
end

varargout = cell(1, nargout);
for k = 1:nargout
    varargout{k} = x(:,k);
end

end
