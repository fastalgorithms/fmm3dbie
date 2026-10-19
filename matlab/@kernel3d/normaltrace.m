function Kn = normaltrace(K)
%NORMALTRACE   Target normal derivative of a kernel from KERNEL3D.RADPP,
% KERNEL3D.RADGREEN, KERNEL3D.LAP3D or KERNEL3D.HELM3D, the kernel of the
% prime type, e.g. 'sp' for 's' and 'sp_diff' for 's_diff'.
%
% See also KERNEL3D.GRAD, KERNEL3D.NTRACE.

p = K.params;

% the prime type, e.g. 'sp' for 's' and 'sp_diff' for 's_diff'
type = regexprep(K.type, '_diff$', '');
if ( ~any(strcmp(type, {'s', 'd', 'c', 'trans_rep'})) )
    error('KERNEL3D.NORMALTRACE: Currently only supported for the types s, d, c and trans_rep, and their _diff.');
end
type = [strrep(type, 'trans_rep', 'trep') 'p'];

% the difference kernels of KERNEL3D.HELM3D
if ( endsWith(K.type, '_diff') )
    Kn = kernel3d.helm3d([type '_diff'], p.zks, p.coefs);
    return
end

switch K.name
    case {'radpp', 'radgreen'}
        % the new type with the same green, fit range (radpp) and options
        opts = {'dim', p.dim, 'opts', p.fitopts, 'zk', K.zk};
        if ( any(strcmp(K.type, {'c', 'trans_rep'})) )
            opts = [opts, {'coefs', p.kcoefs}];
        end
        if ( ~isempty(p.fmm) ), opts = [opts, {'fmm', p.fmm}]; end
        if ( strcmp(K.name, 'radpp') )
            Kn = kernel3d.radpp(p.green, [p.rmin p.rmax], type, opts{:});
        else
            Kn = kernel3d.radpp(p.green, [], type, opts{:});
            q = Kn.params;
            Kn.params.key = {q.spec.name, q.kcoefs, q.nu, q.dim};
        end
    case 'laplace'
        coefs = {};
        if ( isfield(p, 'coefs') ), coefs = {p.coefs}; end
        Kn = kernel3d.lap3d(type, coefs{:});
    case 'helmholtz'
        coefs = {};
        if ( isfield(p, 'coefs') ), coefs = {p.coefs}; end
        Kn = kernel3d.helm3d(type, K.zk, coefs{:});
    otherwise
        error('KERNEL3D.NORMALTRACE: Currently only supported for kernel3d.radpp, kernel3d.radgreen, kernel3d.lap3d and kernel3d.helm3d.');
end

end
