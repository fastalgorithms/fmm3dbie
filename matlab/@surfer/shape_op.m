function shape = shape_op(obj)
%  shape = shape_op(S,p)
%     This subroutine evaluates the shape operator of a surface in
%     cartesian coordinates.
%
%     The shape operator is defined to be surface gradient of the surface
%     normal.
%
%  Input arguments:
%    * obj: surfer object
%
%  Output arguments:
%    * shape: double (3,3,S.npts) 
%
    ffrm  = get_first_fundamental_form(obj);
    ffrm2 = get_second_fundamental_form(obj);
    

    dffrm = ffrm(1,1,:).*ffrm(2,2,:)-ffrm(1,2,:).*ffrm(2,1,:);
    iffrm = ffrm;
    iffrm(1,1,:) =  ffrm(2,2,:)./dffrm;
    iffrm(2,2,:) =  ffrm(1,1,:)./dffrm;
    iffrm(1,2,:) = -ffrm(1,2,:)./dffrm;
    iffrm(2,1,:) = -ffrm(2,1,:)./dffrm;

    shape_red = pagemtimes(iffrm,ffrm2);
    shape_red = pagemtimes(shape_red,iffrm);

    du = obj.du;
    dv = obj.dv;

    X = zeros(3,2,obj.npts);
    XT= zeros(2,3,obj.npts);

    X(:,1,:) = du;
    X(:,2,:) = dv;
    XT(1,:,:) = du;
    XT(2,:,:) = dv;

    shape = pagemtimes(pagemtimes(X,shape_red),XT);
end