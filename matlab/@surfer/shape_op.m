function shape = shape_op(obj)

    ffrm  = get_first_fundamental_form(obj);
    ffrm2 = get_second_fundamental_form(obj);
    

    dffrm = ffrm(1,1,:).*ffrm(2,2,:)-ffrm(1,2,:).*ffrm(2,1,:);
    iffrm = ffrm;
    iffrm(1,1,:) =  ffrm(2,2,:)./dffrm;
    iffrm(2,2,:) =  ffrm(1,1,:)./dffrm;
    iffrm(1,2,:) = -ffrm(1,2,:)./dffrm;
    iffrm(2,1,:) = -ffrm(2,1,:)./dffrm;

    shape_red = pmt(iffrm,ffrm2);

    du = obj.du;
    dv = obj.dv;

    X = zeros(3,2,obj.npts);
    XT= zeros(2,3,obj.npts);

    X(:,1,:) = du;
    X(:,2,:) = dv;
    XT(1,:,:) = du;
    XT(2,:,:) = dv;

    shape = pmt(pmt(X,shape_red),XT);



end

function c = pmt(a,b)
    c = pagemtimes(a,b);
end