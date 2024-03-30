function [ier,geo] = readgidmsh(filename)

ier = 0;

try

    fid = fopen(filename);
    frewind(fid)

    inputs=textscan(fid,'%s %s %f %s %s %s %f',1);
    if (strcmpi(inputs{2},"dimension"))
        dim = inputs{3};
        if ((round(dim)~= dim) || dim<1)
            disp("Error in dimension heading");
            ier = 1;
            dim = 3;
        end
    end
    if (strcmpi(inputs{4},"ElemType"))
        if (~strcmpi(inputs{5},"Triangle"))
            disp("unsupported element type");
            ier = 2;
        end
    end
    if (strcmpi(inputs{6},"Nnode"))
        npertri = inputs{7};
        if ((round(npertri)~= npertri) || npertri<1)
            disp("Error in Nnode heading");
            ier = 3;
            npertri = 6;
        end
    end

    ifn1 = false;
    ifn2 = false;

    nnodes = 0;

    while (~(ifn1 & ifn2))
        fline = fgetl(fid);
        ifnodes = contains(fline,'Coordinates');
        if (ifnodes)
            if (ifn1)
                ifn2 = true;
            else
                ifn1 = true;
            end
        else
            if (ifn1)
                nnodes = nnodes + 1;
            end
        end
    end

    ife1 = false;
    ife2 = false;

    frewind(fid)
    nelems = 0;
    while (~(ife1 & ife2))
        fline = fgetl(fid);
        ifelems = contains(fline,'Elements');
        if (ifelems)
            if (ife1)
                ife2 = true;
            else
                ife1 = true;
            end
        else
            if (ife1)
                nelems = nelems + 1;
            end
        end
    end

    ptarr = zeros(dim,nnodes);
    triarr= zeros(npertri,nelems);

    frewind(fid);
    ifnodes = false;
    while (~ifnodes)
        fline = fgetl(fid);
        ifnodes = contains(fline,'Coordinates');
        if (ifnodes)
            inputs = textscan(fid,join(repmat("%f ",1,dim+1)),nnodes);
        end
    end

    for ii=1:dim
        ptarr(ii,:) = inputs{ii+1};
    end

    frewind(fid)
    ifelems = false;
    while (~(ifelems))
        fline = fgetl(fid);
        ifelems = contains(fline,'Elements');
        if (ifelems)
            inputs = textscan(fid,join(repmat("%f ",1,npertri+1)),nelems);
        end
    end

    for ii=1:npertri
        triarr(ii,:) = inputs{ii+1};
    end

    fclose(fid);

    geo = {};
    geo.npoints = nnodes;
    geo.ntri = nelems;
    geo.npe = npertri;
    geo.ndim = dim;
    geo.points = ptarr;
    geo.tris   = triarr;

catch
    ier = -1;
    geo  = {};

end
