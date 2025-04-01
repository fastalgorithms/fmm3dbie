        subroutine struve102(z,cval0,cval1)
        implicit real *8 (a-h,o-z)
        complex *16 z,cval0,cval1

        if (abs(z) .lt. 95) then
        call eval_struve_quad(z,cval0,cval1)
        return
        endif

        if (abs(z) .ge. 95) then
        call eval_struve_asymp(z,cval0,cval1)
        return
        endif

        return
        end
c
c
c
c
c
        subroutine eval_struve_quad(z,cval0,cval1)
        implicit real *8 (a-h,o-z)
        dimension roots(38),weights(38)
        complex *16 cval0,cval1,z,ima
        data ima /(0,1)/

        done = 1
        pi   = atan(done)*4

c        
c        Nodes:
c        
        data roots/
     1  0.23287797554979764574239513845863310D-02,
     2  0.77797210528565602202196751376561240D-02,
     3  0.16271618033986618198702811217156061D-01,
     4  0.27707503723264139261127948179293856D-01,
     5  0.41963012934216623038977890242892569D-01,
     6  0.58891418218743405778560394085441146D-01,
     7  0.78329391689422605486978192014671943D-01,
     8  0.10010263860402900542980504967723438D+00,
     9  0.12403097627111796188202638019389090D+00,
     *  0.14993260442201184905565140071361611D+00,
     1  0.17762745410598842069997492381278595D+00,
     2  0.20693961217903811539183343424801665D+00,
     3  0.23769889224165259793060177883684293D+00,
     4  0.26974166243662189627173897294466337D+00,
     5  0.30291105282421662136522246065862774D+00,
     6  0.33705665857467929257135214342578155D+00,
     7  0.37203383756330568046403009390787936D+00,
     8  0.40770267783662548059749513310209640D+00,
     9  0.44392668553654047387383102654870569D+00,
     *  0.48057121923925020977835009680085893D+00,
     1  0.51750167311684907657185822394561155D+00,
     2  0.55458138900459695607645970207925428D+00,
     3  0.59166925629421969496594497355751397D+00,
     4  0.62861693883643152084181844298409790D+00,
     5  0.66526565104397613812338690028125968D+00,
     6  0.70144239471352053795506975786674479D+00,
     7  0.73695557165015248739431902109333251D+00,
     8  0.77158992088379493410487935842608357D+00,
     9  0.80510082300500561137518182658171708D+00,
     *  0.83720822069456348292247823421660014D+00,
     1  0.86759080958908695616059731434000299D+00,
     2  0.89588187731948532410575618640971064D+00,
     3  0.92166932590546832099280447416876834D+00,
     4  0.94450394921279951711390249056169951D+00,
     5  0.96392127329142721670042025053186995D+00,
     6  0.97948100688506504715664233260675227D+00,
     7  0.99082007648743246982849951429630583D+00,
     8  0.99769783994110755421389441552630474D+00/
c        
c        Weights:
c        
        data weights/
     1  0.39054246308472108537572143018542489D-02,
     2  0.69851156794182489718899661398675367D-02,
     3  0.99824773011900623033011655549735733D-02,
     4  0.12868477385036746802691008127225257D-01,
     5  0.15617988869632857012575714934677317D-01,
     6  0.18211505504860765237462396137799908D-01,
     7  0.20635307068313352195120159032980171D-01,
     8  0.22881084883173884340069747783455705D-01,
     9  0.24945234599334783449814705461144249D-01,
     *  0.26827969985166285752458832247831409D-01,
     1  0.28532384616938756151580122811620477D-01,
     2  0.30063553632169448782541566530842284D-01,
     3  0.31427731341652938389422977827806917D-01,
     4  0.32631669241865390472468304928884732D-01,
     5  0.33682055948710753297697772662680929D-01,
     6  0.34585065860587188374649279101325372D-01,
     7  0.35345995387008851284279099044910354D-01,
     8  0.35968962282864375803451612771640899D-01,
     9  0.36456643112443903367680001996253688D-01,
     *  0.36810024688899325062925443352680342D-01,
     1  0.37028146561546252296714166666237650D-01,
     2  0.37107812810054908793043848494316541D-01,
     3  0.37043252593480939850414399158839917D-01,
     4  0.36825710732858259789984880122231590D-01,
     5  0.36442953627131857006501666929733034D-01,
     6  0.35878685118130530049529049930153823D-01,
     7  0.35111887371695664011105788319219062D-01,
     8  0.34116143950671667563908738183976545D-01,
     9  0.32859083990848582685090102345858933D-01,
     *  0.31302235858241964317952738389628315D-01,
     1  0.29401832538196193841361929053862802D-01,
     2  0.27111495402561458561401180300153319D-01,
     3  0.24388182350537953056835591739052492D-01,
     4  0.21202989623306897961392790227859606D-01,
     5  0.17557342613948635157400627184362894D-01,
     6  0.13500760099128853258969852643313634D-01,
     7  0.91371732847652890752352921528928191D-02,
     8  0.46001586784217374251199154849668279D-02/
        
        cval0 = 0
        cval1 = 0
        npts = 38
        do i=1,npts
            r = roots(i)
            w = weights(i)
            cval0 = cval0 + 1.0d0/(sqrt(1.0d0-r*r))*
     1      (exp(ima*r*z)-(1+r*(exp(ima*z)-1)))*w
            cval1 = cval1 + 1.0d0/(sqrt(1.0d0-r*r))*
     1      (ima*r*exp(ima*r*z)-(ima*r*(exp(ima*z))))*w
        enddo

        cval0 = cval0 + pi/2 + (exp(ima*z)-1)
        cval1 = cval1 + ima*exp(ima*z)

        cval0 = 2/pi*cval0
        cval1 = 2/pi*cval1
        cval1 = -cval1 + ima*2/pi

        return
        end
c
c
c
c
c
c

        subroutine eval_struve_asymp(zvin,cval0,cval1)
        implicit real *8 (a-h,o-z)
        complex *16 ima,zvin,zv,zvv,h0,h1,ct,ct1,cval0,
     1      cval1
        data ima/(0,1)/

        zv = -zvin
        
        eps = 1.0d-15

        done = 1
        pi   = atan(done)*4.0d0

        cval0 = 2.0d0/pi/zv 
        ct    = cval0

        cval1 = 2.0d0/pi
        ct1   = cval1

        r1 = abs(ct)

        do i=1,1000
            di = i
            ct =-ct/zv/zv/4.0d0/di/di*
     1              (2.0d0*di*(2.0d0*di-1.0d0))**2
            ct1= -ct1/zv/zv*4.0d0*
     1              ((di-1.0d0)**2-0.25d0)

            r2 = abs(ct)
ccc            call prin2('r2 = *',r2,1)
            if (r2 .gt. r1) goto 2100
            cval0 = cval0 + ct
            cval1 = cval1 + ct1
            if (r2 .lt. eps) goto 2100
            r1 = r2
        enddo

 2100 continue

ccc        call prinf('i = *',i,1)
        zvv = zvin
        call hank103(zvv,h0,h1,1)
ccc        call prin2('h0 = *',h0,2)

        if (real(zvv) .gt. 0) then
            cval0 = cval0+ima*h0
            cval1 = cval1-ima*h1
        endif
        if (real(zvv) .le. 0) then
            cval0 = cval0-ima*h0
            cval1 = cval1+ima*h1
        endif

        cval0 = -ima*cval0
        cval1 =  ima*cval1


        return
        end
c
c
c
c
c
        subroutine sk0(rhoj,dr,val)
        implicit real *8 (a-h,o-z)
        complex *16 rhoj,val,ima,cr0,cr1
        complex *16 h0,h1,zt
        real *8 dr
        logical ilow
        data ima /(0,1)/

        ilow = .false.

        if (imag(rhoj) .lt. 0) then
            ilow = .true.
            rhoj = dconjg(rhoj)
        endif

        zt = rhoj*dr
    
        call hank103(zt,h0,h1,1)
        call struve102(zt,cr0,cr1)

        val = -ima*cr0+ima*h0

        if (ilow) then
            val = dconjg(val)
            rhoj = dconjg(rhoj)
        endif

        return
        end
c
c
c
c
        subroutine gradsk0(rhoj,dx,dy,sk0x,sk0y)
        implicit real *8 (a-h,o-z)
        complex *16 rhoj,val,ima,cr0,cr1
        complex *16 h0,h1,zt,h0x,h0y,cr0x,cr0y
        complex *16 sk0x,sk0y
        real *8 dx,dy,dr,dr2,pi
        logical ilow 
        data ima /(0,1)/

        ilow = .false.
        pi = atan(1.0d0)*4.0d0

        if (imag(rhoj) .lt. 0) then
            ilow = .true.
            rhoj = dconjg(rhoj)
        endif

        dr2 = dx**2 + dy**2
        dr = sqrt(dr2)

        zt = rhoj*dr

        call hank103(zt,h0,h1,1)
        call struve102(zt,cr0,cr1)

        cr0x = 2/pi*ima*rhoj*dx/dr - rhoj*dx/dr*cr1
        cr0y = 2/pi*ima*rhoj*dy/dr - rhoj*dy/dr*cr1

        h0x = -rhoj*dx/dr*h1
        h0y = -rhoj*dy/dr*h1

        sk0x = -ima*cr0x+ima*h0x
        sk0y = -ima*cr0y+ima*h0y

        if (ilow) then
            sk0x = dconjg(sk0x)
            sk0y = dconjg(sk0y)
            rhoj = dconjg(rhoj)
        endif

        return
        end
c
c
c
c
c
        subroutine lapsk0(rhoj,dr,lap)
        implicit real *8 (a-h,o-z)
        complex *16 rhoj,val,ima,cr0,cr1
        complex *16 h0,h1,zt
        complex *16 lapcr0,laph0,lap
        real *8 dr,pi
        logical ilow 
        data ima /(0,1)/

        ilow = .false.
        pi = atan(1.0d0)*4.0d0

        if (imag(rhoj) .lt. 0) then
            ilow = .true.
            rhoj = dconjg(rhoj)
        endif

        zt = rhoj*dr

        call hank103(zt,h0,h1,1)
        call struve102(zt,cr0,cr1)

        laph0 = -rhoj*rhoj*h0
        lapcr0 = 2/pi*ima*rhoj/dr - rhoj**2*cr0

        lap = -ima*lapcr0+ima*laph0

        if (ilow) then
            lap = dconjg(lap)
            rhoj = dconjg(rhoj)
        endif

        return
        end
c
c
c
c
