function shape = shape_op(obj)
%  shape = shape_op(S)
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
    ffrm2 = cell2mat(obj.sfform);
    ffrm2 = reshape(ffrm2, 2, obj.npatches, 2,[]); 
    ffrm2 = permute(ffrm2,[1,3,4,2]); ffrm2= ffrm2(:,:,:);

    iffrm = cell2mat(obj.ffforminv);
    iffrm = reshape(iffrm,  2, obj.npatches, 2,[]);
    iffrm = permute(iffrm,[1,3,4,2]); iffrm= iffrm(:,:,:);

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