function [ier,geo] = readmsh(filename)

ier = 0;

try

    fid = fopen(filename);
    inputs=textscan(fid,'%f',5);
    % dimension
    ndim = inputs{1}(1);
    % number of points per element
    npe = inputs{1}(3);
    % number of points
    m = inputs{1}(4);
    % number of elements
    N = inputs{1}(5);

    %skip the end of the first line
    fgetl(fid);

    ptarr = zeros(ndim,m);
    for ii=1:m
        fgetl(fid);
        inputs=textscan(fid,'%f',ndim);
        fgetl(fid);
        ptarr(:,ii) = inputs{1};
    end

    fgetl(fid);
    triarr = zeros(npe,N);
    for ii=1:N
        fgetl(fid);
        inputs=textscan(fid,'%f',npe);
        triarr(:,ii) = inputs{1};
        fgetl(fid);
    end

    fclose(fid);

    geo = {};
    geo.npoints = m;
    geo.ntri = N;
    geo.npe = npe;
    geo.ndim = ndim;
    geo.points = ptarr;
    geo.tris   = triarr;

catch
    ier = 1;
    geo  = {};

end