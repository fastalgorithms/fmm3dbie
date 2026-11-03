function [rts,ejs] = find_roots(a,gamma,n)
%SURFWAVE.FIND_ROOTS   Roots of a z^n + gamma z - 1 and the partial
%   fraction coefficients ejs of its inverse, sorted by abs(angle(rts)).
%   n = 3: capillary (a = beta), n = 5: flexural (a = alpha), n = 1:
%   gravity (a = 2/g, gamma = 0).

    if ( n == 1 )
        rts = 1/(a + gamma);  ejs = rts;
        return
    end

    d1 = gamma ./ a;
    d0 = - 1 ./ a;

    % companion matrix, then Newton
    C = diag(ones(n-1,1),-1);
    C(1:2,end) = [-d0; -d1];

    f = @(z) z.^n + d1*z + d0;
    fp = @(z) n*z.^(n-1) + d1;

    rts = eig(C);
    for j = 1:3
        rts = rts - f(rts)./fp(rts);
    end

    ejs = 1./fp(rts);
    ejs = ejs./a;

    [~, sort_id] = sort(abs(angle(rts)));
    rts = rts(sort_id);
    ejs = ejs(sort_id);

end
