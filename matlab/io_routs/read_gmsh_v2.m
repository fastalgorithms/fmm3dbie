function [ier,geo] =read_gmsh_v2(filename)

ier = 0;

try

    fid = fopen(filename);
    frewind(fid)
    ifnodefnd = false;
    ifelemfnd = false;

    itype_tri3 = 2;
    itype_tri6 = 9;
    itype_quad4 = 4;
    itype_quad9 = 10;
    itype_quad8 = 16;

    while ((~ifelemfnd) || (~ifnodefnd))
        fline = fgetl(fid);
        ifnodes = contains(fline,'$Nodes');

        if (ifnodes)
            ifnodefnd = true;
            inputs=textscan(fid,'%f',1);
            nnodes = inputs{1};
            fgetl(fid);

            points = zeros(3,nnodes);
            inputs=textscan(fid,'%f %f %f %f',nnodes);
            points(1,:) = inputs{2};
            points(2,:) = inputs{3};
            points(3,:) = inputs{4};
            fgetl(fid);
        end

        ifelems = contains(fline,'$Elements');

        if(ifelems)
            ntri = 0;
            npts = nnodes;

            ifelemfnd = true;
            inputs=textscan(fid,'%f',1);
            nelems = inputs{1};
            fgetl(fid);
            for ii=1:nelems
                inputs=textscan(fid,'%f',3);
                itri = inputs{1}(2);
                if (itri == itype_tri3)
                    ntri = ntri + 1;
                    npts = npts + 3;
                elseif (itri == itype_tri6)
                    ntri = ntri + 1;
                elseif (itri == itype_quad4)
                    ntri = ntri + 2;
                    npts = npts + 5;
                elseif (itri == itype_quad8)
                    ntri = ntri + 2;
                    npts = npts + 1;
                elseif (itri == itype_quad9)
                    ntri = ntri + 2;
                end
                fgetl(fid);
            end
        end

    end

    m = npts;
    N = ntri;
    ptarr = zeros(3,m);
    triarr = zeros(6,N);

    ptarr(:,1:nnodes) = points;
    npts_use = nnodes;

    frewind(fid);

    ifelemfnd = false;
    while ((~ifelemfnd))
        fline = fgetl(fid);
        ifelems = contains(fline,'$Elements');
        if(ifelems)
            ifelemfnd = true;
        end

    end

    fline = fgetl(fid);

    nel_tri3 = 3;
    nel_tri6 = 6;
    nel_quad4 = 4;
    nel_quad9 = 9;
    nel_quad8 = 8;
    ntri = 0;

    xs = ptarr(1,:);
    ys = ptarr(2,:);
    zs = ptarr(3,:);

    for ii=1:nelems

        inputs=textscan(fid,'%f',3);
        itri = inputs{1}(2);
        ntag = inputs{1}(3);

        if (itri == itype_tri3)

            nel = ntag + nel_tri3;
            inputs=textscan(fid,'%f',nel);
            itmp = inputs{1};

            inode1 = itmp(ntag+1);
            inode2 = itmp(ntag+2);
            inode3 = itmp(ntag+3);

            ntri = ntri + 1;
            triarr(1,ntri) = inode1;
            triarr(2,ntri) = inode2;
            triarr(3,ntri) = inode3;
            triarr(4,ntri) = npts_use + 1;
            triarr(5,ntri) = npts_use + 2;
            triarr(6,ntri) = npts_use + 3;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode1) + xs(inode2))/2;
            ptarr(2,npts_use) = (ys(inode1) + ys(inode2))/2;
            ptarr(3,npts_use) = (zs(inode1) + zs(inode2))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode2) + xs(inode3))/2;
            ptarr(2,npts_use) = (ys(inode2) + ys(inode3))/2;
            ptarr(3,npts_use) = (zs(inode2) + zs(inode3))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode3) + xs(inode1))/2;
            ptarr(2,npts_use) = (ys(inode3) + ys(inode1))/2;
            ptarr(3,npts_use) = (zs(inode3) + zs(inode1))/2;

        elseif (itri == itype_tri6)

            nel = ntag + nel_tri6;
            inputs=textscan(fid,'%f',nel);
            itmp = inputs{1};

            ntri = ntri + 1;
            triarr(:,ntri) = itmp(ntag + (1:6));

        elseif (itri == itype_quad4)

            nel = ntag + nel_quad4;
            inputs=textscan(fid,'%f',nel);
            itmp = inputs{1};

            inode1 = itmp(ntag+1);
            inode2 = itmp(ntag+2);
            inode3 = itmp(ntag+3);
            inode4 = itmp(ntag+4);

            ntri = ntri + 1;
            triarr(1,ntri) = inode1;
            triarr(2,ntri) = inode2;
            triarr(3,ntri) = inode4;
            triarr(4,ntri) = npts_use + 1;
            triarr(5,ntri) = npts_use + 5;
            triarr(6,ntri) = npts_use + 4;

            ntri = ntri + 1;
            triarr(1,ntri) = inode2;
            triarr(2,ntri) = inode4;
            triarr(3,ntri) = inode3;
            triarr(4,ntri) = npts_use + 5;
            triarr(5,ntri) = npts_use + 3;
            triarr(6,ntri) = npts_use + 2;


            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode1) + xs(inode2))/2;
            ptarr(2,npts_use) = (ys(inode1) + ys(inode2))/2;
            ptarr(3,npts_use) = (zs(inode1) + zs(inode2))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode2) + xs(inode3))/2;
            ptarr(2,npts_use) = (ys(inode2) + ys(inode3))/2;
            ptarr(3,npts_use) = (zs(inode2) + zs(inode3))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode3) + xs(inode4))/2;
            ptarr(2,npts_use) = (ys(inode3) + ys(inode4))/2;
            ptarr(3,npts_use) = (zs(inode3) + zs(inode4))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode4) + xs(inode1))/2;
            ptarr(2,npts_use) = (ys(inode4) + ys(inode1))/2;
            ptarr(3,npts_use) = (zs(inode4) + zs(inode1))/2;

            npts_use = npts_use+1;
            ptarr(1,npts_use) = (xs(inode1) + xs(inode2) + ...
                xs(inode3) + xs(inode4))/4;
            ptarr(2,npts_use) = (ys(inode1) + ys(inode2) + ...
                ys(inode3) + ys(inode4))/4;
            ptarr(3,npts_use) = (zs(inode1) + zs(inode2) + ...
                zs(inode3) + zs(inode4))/4;

        elseif (itri == itype_quad8)

            nel = ntag + nel_quad8;
            inputs=textscan(fid,'%f',nel);
            itmp = inputs{1};

            inode8 = itmp((ntag+1):(ntag+8));

            ntri = ntri + 1;
            triarr(1,ntri) = inode8(1);
            triarr(2,ntri) = inode8(2);
            triarr(3,ntri) = inode8(4);
            triarr(4,ntri) = inode8(5);
            triarr(5,ntri) = npts_use + 1;
            triarr(6,ntri) = inode8(8);


            ntri = ntri + 1;
            triarr(1,ntri) = inode8(3);
            triarr(2,ntri) = inode8(4);
            triarr(3,ntri) = inode8(2);
            triarr(4,ntri) = inode8(7);
            triarr(5,ntri) = npts_use + 1;
            triarr(6,ntri) = inode8(6);

            npts_use = npts_use + 1;
            ptarr(1:3,npts_use) = 0;

            ptarr(:,npts_use) = [sum(xs(inode8));sum(ys(inode8));...
                sum(zs(inode8))]/8;

        elseif (itri == itype_quad9)

            nel = 3 + ntag + nel_quad9;
            inputs=textscan(fid,'%f',nel);
            itmp = inputs{1};

            ntri = ntri + 1;
            triarr(1,ntri) = itmp(3+ntag+1);
            triarr(2,ntri) = itmp(3+ntag+2);
            triarr(3,ntri) = itmp(3+ntag+4);
            triarr(4,ntri) = itmp(3+ntag+5);
            triarr(5,ntri) = itmp(3+ntag+9);
            triarr(6,ntri) = itmp(3+ntag+8);


            ntri = ntri + 1;
            triarr(1,ntri) = itmp(3+ntag+3);
            triarr(2,ntri) = itmp(3+ntag+4);
            triarr(3,ntri) = itmp(3+ntag+2);
            triarr(4,ntri) = itmp(3+ntag+7);
            triarr(5,ntri) = itmp(3+ntag+9);
            triarr(6,ntri) = itmp(3+ntag+6);


        end
        fgetl(fid);
    end


    fclose(fid);

    geo = {};
    geo.nnodes = nnodes;
    geo.npoints = size(ptarr,2);
    geo.ntri = N;
    geo.ndim = 3;
    geo.points = ptarr;
    geo.tris   = triarr;

catch
    ier = 1;
    geo  = {};
end