function [ier,geo] = readtri(filename)

ier = 0;

try

    fid = fopen(filename);
    frewind(fid);
    inputs = textscan(fid,"%f %f",1);
    fgetl(fid);

    m = inputs{1};
    N = inputs{2};

    inputs = textscan(fid,"%f %f %f",m);
    points = zeros(3,m);
    for ii=1:3
        points(ii,:) = inputs{ii};
    end

    inputs = textscan(fid,"%f %f %f",N);
    elems = zeros(3,N);
    for ii=1:3
        elems(ii,:) = inputs{ii};
    end

    ptarr = zeros(3,m+3*N);
    ptarr(:,1:m) = points;
    triarr= zeros(6,N);

    ipointer=m+1;

    for j=1:N
        triarr(1:3,j) = elems(:,j);
        triarr(4:6,j) = ipointer+(0:2);
        ptarr(:,ipointer)   = (ptarr(:,triarr(1,j))+ptarr(:,triarr(2,j)))/2;
        ptarr(:,ipointer+1) = (ptarr(:,triarr(2,j))+ptarr(:,triarr(3,j)))/2;
        ptarr(:,ipointer+2) = (ptarr(:,triarr(1,j))+ptarr(:,triarr(3,j)))/2;
        ipointer=ipointer+3;
    end

    fclose(fid);

    geo = {};
    geo.npoints = m+3*N;
    geo.ntri = N;
    geo.npe = 6;
    geo.ndim = 3;
    geo.points = ptarr;
    geo.tris   = triarr;
    geo.basepts= points;

catch
    ier = 1;
    geo  = {};

end