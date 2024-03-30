function [ier,geo] =read_gmsh_v4(filename)

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
itype_quad8 = NaN;

while ((~ifelemfnd) || (~ifnodefnd))
    fline = fgetl(fid);
    ifnodes = contains(fline,'$Nodes');

    if (ifnodes)
        ifnodefnd = true;
        inputs=textscan(fid,'%f %f %f %f',1);
        nentity = inputs{1};
        nodetot = inputs{2};
        fgetl(fid);

        points = zeros(3,nodetot);
        for ii=1:nentity
            inputs=textscan(fid,'%f %f %f %f',1);
            fgetl(fid);
            nnode = inputs{4};
            inputs=textscan(fid,'%f',nnode);
            iindvec = inputs{1};
            inputs=textscan(fid,'%f %f %f',nnode);
            xs = inputs{1};
            ys = inputs{2};
            zs = inputs{3};
            points(1,iindvec) = xs;
            points(2,iindvec) = ys;
            points(3,iindvec) = zs;
        end


        fgetl(fid);
        nnodes = nodetot;
    end



    ifelems = contains(fline,'$Elements');

    if(ifelems)
        inputs=textscan(fid,'%f %f %f %f',1);
        nentity = inputs{1};
        neletot = inputs{2};
        fgetl(fid);

        ntri = 0;
        npts = 0;

        ifelemfnd = true;

        for ii=1:nentity
            inputs=textscan(fid,'%f %f %f %f',1);
            fgetl(fid);
            ielemtype = inputs{3};
            nelem = inputs{4};
            if(ielemtype==itype_tri3)
                ntri = ntri + nelem;
                npts = npts + 3*nelem;
            end

            if(ielemtype==itype_tri6)
                ntri = ntri + nelem;
            end

            if(ielemtype==itype_quad4)
                ntri = ntri + 2*nelem;
                npts = npts + 5*nelem;
            end

            if(ielemtype==itype_quad8)
                ntri = ntri + 2*nelem;
                npts = npts + nelem;
            end

            if(ielemtype==itype_quad9)
                ntri = ntri + 2*nelem;
            end
            for jj=1:nelem
                fgetl(fid);
            end
        end
    end

end

npts = npts + nnodes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nel_tri3 = 3;
nel_tri6 = 6;
nel_quad4 = 4;
nel_quad9 = 9;

xs = points(1,:);
ys = points(2,:);
zs = points(3,:);

ptarr=zeros(3,npts);
triarr=zeros(6,ntri);
npts_use = nnodes;

ptarr(:,1:nnodes) = points;

frewind(fid);

ifelemfnd = false;
while ((~ifelemfnd) || (~ifnodefnd))
    fline = fgetl(fid);
    ifelems = contains(fline,'$Elements');
    if(ifelems)
        ifelemfnd = true;
    end
end

inputs=textscan(fid,'%f %f %f %f',1);
nentity = inputs{1};
nodetot = inputs{2};
fgetl(fid);

points = zeros(3,npts);
for ii=1:nentity
    inputs=textscan(fid,'%f %f %f %f',1);
    fgetl(fid);
    ielemtype = inputs{3};
    nelem = inputs{4};

    if(ielemtype==itype_tri3)
        nel = 1 + nel_tri3;
        for j=1:nelem
            inputs=textscan(fid,'%f',nel);
            fgetl(fid);
            itmp = inputs{1};
            inode1 = itmp(2);
            inode2 = itmp(3);
            inode3 = itmp(4);

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
        end
    end

    if(ielemtype==itype_tri6)
        nel = 1 + nel_tri6;
        for j=1:nelem
            inputs=textscan(fid,'%f',nel);
            fgetl(fid);
            itmp = inputs{1};
            ntri = ntri + 1;
            for l = 1:6
                triarr(l,ntri) = itmp(l+1);
            end
        end
    end

    if(ielemtype==itype_quad4)
        nel = 1 + nel_quad4;
        for j=1:nelem
            inputs=textscan(fid,'%f',nel);
            fgetl(fid);
            itmp = inputs{1};
            inode1 = itmp(2);
            inode2 = itmp(3);
            inode3 = itmp(4);
            inode4 = itmp(5);

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
        end
    end

    if(ielemtype==itype_quad8)
        nel = 1 + nel_quad8;
        for j=1:nelem
            inputs=textscan(fid,'%f',nel);
            fgetl(fid);
            itmp = inputs{1};
            inode8 = itmp(2:9);

            ntri = ntri + 1;
            triarr(1,ntri) = itmp(2);
            triarr(2,ntri) = itmp(3);
            triarr(3,ntri) = itmp(5);
            triarr(4,ntri) = itmp(6);
            triarr(5,ntri) = npts_use + 1;
            triarr(6,ntri) = itmp(9);


            ntri = ntri + 1;
            triarr(1,ntri) = itmp(4);
            triarr(2,ntri) = itmp(5);
            triarr(3,ntri) = itmp(3);
            triarr(4,ntri) = itmp(8);
            triarr(5,ntri) = npts_use + 1;
            triarr(6,ntri) = itmp(7);

            npts_use = npts_use + 1;
            ptarr(1:3,npts_use) = 0;
            for l=1:8
                ptarr(1,npts_use) = ptarr(1,npts_use) + ...
                    xs(inode8(l));
                ptarr(2,npts_use) = ptarr(1,npts_use) + ...
                    ys(inode8(l));
                ptarr(3,npts_use) = ptarr(1,npts_use) + ...
                    zs(inode8(l));
            end
            ptarr(1:3,npts_use) = ptarr(1:3,npts_use)/8;
        end
    end

    if(ielemtype==itype_quad9)
        nel = 1 + nel_quad9;
        for j=1:nelem
            inputs=textscan(fid,'%f',nel);
            fgetl(fid);
            itmp = inputs{1};
            ntri = ntri + 1;
            triarr(1,ntri) = itmp(2);
            triarr(2,ntri) = itmp(3);
            triarr(3,ntri) = itmp(5);
            triarr(4,ntri) = itmp(6);
            triarr(5,ntri) = itmp(10);
            triarr(6,ntri) = itmp(9);


            ntri = ntri + 1;
            triarr(1,ntri) = itmp(4);
            triarr(2,ntri) = itmp(5);
            triarr(3,ntri) = itmp(3);
            triarr(4,ntri) = itmp(8);
            triarr(5,ntri) = itmp(10);
            triarr(6,ntri) = itmp(7);
        end
    end

end

    fclose(fid);

    geo = {};
    geo.nnodes = nnodes;
    geo.npoints = size(ptarr,2);
    geo.ntri = ntri;
    geo.ndim = 3;
    geo.points = ptarr;
    geo.tris   = triarr;

catch
    ier = 1;
    geo  = {};
end