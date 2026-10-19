function Kg = grad(K)
%GRAD   Target gradient of a kernel from KERNEL3D.RADPP, KERNEL3D.RADGREEN,
% KERNEL3D.LAP3D or KERNEL3D.HELM3D, the kernel of the gradient type, e.g.
% 'grad' (alias 'sg') for 's' and 'cgrad' for 'c'.
%
% See also KERNEL3D.NORMALTRACE.

if ( ~any(strcmp(K.type, {'s', 'd', 'c', 'trans_rep'})) )
    error('KERNEL3D.GRAD: Currently only supported for the types s, d, c and trans_rep.');
end
type = [strrep(K.type, 'trans_rep', 'trep') 'g'];

p = K.params;
switch K.name
    case {'radpp', 'radgreen'}
        % the new type with the same green, fit range (radpp) and options
        opts = {'dim', p.dim, 'opts', p.fitopts, 'zk', K.zk};
        if ( any(strcmp(K.type, {'c', 'trans_rep'})) )
            opts = [opts, {'coefs', p.kcoefs}];
        end
        if ( ~isempty(p.fmm) ), opts = [opts, {'fmm', p.fmm}]; end
        if ( strcmp(K.name, 'radpp') )
            Kg = kernel3d.radpp(p.green, [p.rmin p.rmax], type, opts{:});
        else
            Kg = kernel3d.radpp(p.green, [], type, opts{:});
            q = Kg.params;
            Kg.params.key = {q.spec.name, q.kcoefs, q.nu, q.dim};
        end
    case 'laplace'
        coefs = {};
        if ( isfield(p, 'coefs') ), coefs = {p.coefs}; end
        Kg = kernel3d.lap3d(type, coefs{:});
    case 'helmholtz'
        coefs = {};
        if ( isfield(p, 'coefs') ), coefs = {p.coefs}; end
        Kg = kernel3d.helm3d(type, K.zk, coefs{:});
    otherwise
        error('KERNEL3D.GRAD: Currently only supported for kernel3d.radpp, kernel3d.radgreen, kernel3d.lap3d and kernel3d.helm3d.');
end

end
