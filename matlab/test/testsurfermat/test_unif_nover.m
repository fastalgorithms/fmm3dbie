%TEST_UNIF_NOVER   Verify that opts.unif_nover=1 in surfermat produces
%   oversampling orders that are constant across target surfaces for each
%   source surface j: novers{i,j} must be the same for all i.

run ../../startup.m

S1 = geometries.sphere(1, 1, [0;0;0], 4);
S2 = geometries.sphere(1, 1, [2;0;0], 4);
S3 = geometries.sphere(1, 1, [4;0;0], 4);

kern = kernel3d('laplace', 'c', [1;1]);
eps  = 1e-10;
opts = struct('corrections', 1, 'unif_nover', 1);
[~, novers] = surfermat([S1, S2, S3], kern, eps, opts);

pass = true;
for j = 1:3
    for i = 1:3
        for k = 1:3
        if i== k, continue, end
    if ~isequal(novers{i,j}, novers{k,j})
        fprintf('FAIL: novers{%d,%d} ~= novers{$d,%d}\n', i,j,k, j);
        pass = false;
    end
        end
    end
end

if pass
    fprintf('PASS\n');
end

opts = struct('corrections', 1);
[~, novers_nonuni] = surfermat([S1, S2, S3], kern, eps, opts);
