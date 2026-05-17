        
        subroutine hankdiff(zs,nzs,nder,nsub,
     1      cvals)
        complex *16 zs(nzs),cvals(nder+1,nzs),z

        rcut = 1
        do ii=1,nzs
        z = zs(ii)
        if (abs(z) .le. rcut) then
        call hankdiff_loc(z,1,nder,nsub,
     1      cvals(1,ii))
        else
        call hankdiff_far(z,1,nder,nsub,cvals(1,ii))
        endif

        enddo




        return
        end
c
c
c
c
c
        subroutine hankdifftest(zs,nzs,nder,nsub,
     1      cvals)
        complex *16 zs(nzs),cvals(nder+1,nzs),z,
     1      ctmp(nder+1)

        rcut = 1

        do ii=1,nzs
        
        z = zs(ii)

        call hankdiff_loc(z,1,nder,nsub,
     1      cvals(1,ii))
        call prin2('cvals = *',cvals(1,ii),2*nder+2)
        call hankdiff_far(z,1,nder,nsub,ctmp)
        call prin2('cvals = *',ctmp,2*nder+2)

        do jj=1,nder + 1
        cvals(jj,ii) = cvals(jj,ii)-ctmp(jj)
        enddo

        enddo



        return
        end
c
c
c
c
c
        subroutine hankdiff_loc(zs,nzs,nder,nsub,
     1      cvals)
        complex *16 zs(*),cfs(11),cfslog(11),
     1      cfd(11,nder),cfdlog(11,nder+1),
     1      cvals(nzs,nder+1),
     1      zt,clog,z,cvt

        data cfs/
     1    (  1.00000000000000000d0, -0.07380429510868722527d0       ),
     1    ( -0.25000000000000000d0,  0.17760601686906714209d0       ),
     1    (  0.01562500000000000d0, -0.016073968025938425623d0      ),
     1    ( -0.43402777777777777778d-3, 0.0005386026668616549570d0  ),
     1    (  6.781684027777777778d-6, -9.495005205221546473d-6      ),
     1    ( -6.7816840277777777778d-8, 1.03584760336280966885d-7    ),
     1    (  4.7095027970679012346d-10, -7.6930799009029318530d-10  ),
     1    ( -2.4028075495244394054d-12, 4.1435657365127097585d-12   ),
     1    (  9.385966990329841427d-15, -1.6932715179356949520d-14   ), 
     1    ( -2.8969033920771115516d-17, 5.4310606578547997903d-17   ),  
     1    (  7.242258480192778879d-20, -1.4038708139145750726d-19   )/ 

        data cfslog/
     1     ( 0,  0.6366197723675813431d0        ),        
     1     ( 0, -0.15915494309189533577d0       ), 
     1     ( 0,  0.009947183943243458486d0      ),
     1     ( 0, -0.00027631066509009606904d0    ),
     1     ( 0,  4.3173541420327510788d-6       ),
     1     ( 0, -4.3173541420327510788d-8       ),
     1     ( 0,  2.9981625986338549158d-10      ),
     1     ( 0, -1.5296747952213545489d-12      ),
     1     ( 0,  5.9752921688334162066d-15      ),
     1     ( 0, -1.8442259780350050020d-17      ),
     1     ( 0,  4.6105649450875125051d-20      )/


        nord   = 11
        nterms = 11

        do ii=1,nord
            cfd(ii,1)    = 2*(ii-1)*cfs(ii) 
            if (ii .gt. nsub) then
                cfd(ii,1) = cfd(ii,1) + cfslog(ii)
            endif
            cfdlog(ii,1) = 0
            if (ii .gt. nsub) then
                cfdlog(ii,1) = 2*(ii-1)*cfslog(ii)
            endif
        enddo

        do jj=2,nder
            do ii=1,nord
                cfd(ii,jj)    = (2*(ii-1)-jj+1)*cfd(ii,jj-1) 
                cfd(ii,jj) = cfd(ii,jj) + cfdlog(ii,jj-1)
                cfdlog(ii,jj) = (2*(ii-1)-jj+1)*cfdlog(ii,jj-1)
            enddo

        enddo


        do jj=1,nzs

        z = zs(jj)
        zt   = 1
        clog = log(z) 
        

        cvt = 0

        do ii=1,nterms
        
        cvt = cvt + zt*cfs(ii)

        if (ii > nsub) then
            cvt = cvt + zt*clog*cfslog(ii)
        endif

        zt  = zt*(z*z)

        enddo
        
        cvals(jj,1) = cvt

        do kk=1,nder

        zt   = z**(-kk)

        cvt = 0

        do jjj=1,nterms
        
        cvt = cvt + zt*cfd(jjj,kk)

        cvt = cvt + zt*clog*cfdlog(jjj,kk)

        zt  = zt*(z*z)

        enddo
        
        cvals(jj,kk+1) = cvt


c           .   .   .   end of nder
        enddo

c           .   .   .   end of zs loop
        enddo

        return
        end
c
c
c
c
c
        subroutine hankdiff_far(zs,nzs,nder,nsub,cvals)
        implicit real *8 (a-h,o-z)
        complex *16 zs(*),cfs(11),cfslog(11),
     1      cfd(11,nder),cfdlog(11,nder+1),
     1      cvals(nzs,nder+1),
     1      zt,clog,z,cvt,h0,h1,h2,hd,hanks(10),
     1      hders(10)

        cvals = 0
         
        data cfslog/
     1     ( 0,  0.6366197723675813431d0        ),        
     1     ( 0, -0.15915494309189533577d0       ), 
     1     ( 0,  0.009947183943243458486d0      ),
     1     ( 0, -0.00027631066509009606904d0    ),
     1     ( 0,  4.3173541420327510788d-6       ),
     1     ( 0, -4.3173541420327510788d-8       ),
     1     ( 0,  2.9981625986338549158d-10      ),
     1     ( 0, -1.5296747952213545489d-12      ),
     1     ( 0,  5.9752921688334162066d-15      ),
     1     ( 0, -1.8442259780350050020d-17      ),
     1     ( 0,  4.6105649450875125051d-20      )/


        nord   = 11
        nterms = 11
        do ii=(nsub+1),nterms
            cfslog(ii) = 0
        enddo

        do ii=1,nord
            cfd(ii,1)    = 0
            cfd(ii,1) = cfd(ii,1) + cfslog(ii)
            cfdlog(ii,1) = 0
            cfdlog(ii,1) = 2*(ii-1)*cfslog(ii)
        enddo

        do jj=2,nder
            do ii=1,nord
                cfd(ii,jj)    = (2*(ii-1)-jj+1)*cfd(ii,jj-1) 
                cfd(ii,jj) = cfd(ii,jj) + cfdlog(ii,jj-1)
                cfdlog(ii,jj) = (2*(ii-1)-jj+1)*cfdlog(ii,jj-1)
            enddo

        enddo

        do jj=1,nzs

        z = zs(jj)
        zt   = 1
        clog = log(z) 
        
        
        call hank101(z,h0,h1)
        cvals(jj,1) = h0
        cvals(jj,2) =-h1

        hders(1) = -h0+h1/z
        hders(2) =  h0/z - (2-z*z)*h1/(z*z)
        hders(3) =  (-3+z*z)*h0/(z*z) - 2*(z*z-3)*h1/(z*z*z)
        hders(4) =  -2*(-6+z*z)*h0/(z*z*z) - 
     1          (24-7*z*z+z*z*z*z)*h1/(z*z*z*z)

        do kk=2,nder
            cvals(jj,kk+1) = hders(kk-1)
        enddo

        cvt = 0

        do ii=1,nterms
        
        if (ii .le. nsub) then
            cvt = cvt + zt*clog*cfslog(ii)
        endif

        zt  = zt*(z*z)

        enddo
        
        cvals(jj,1) = cvals(jj,1)-cvt

        do kk=1,nder

        zt   = z**(-kk)

        cvt = 0

        do jjj=1,nterms
        
        cvt = cvt + zt*cfd(jjj,kk)

        cvt = cvt + zt*clog*cfdlog(jjj,kk)

        zt  = zt*(z*z)

        enddo
        
        cvals(jj,kk+1) = cvals(jj,kk+1)-cvt


c           .   .   .   end of nder
        enddo


c           .   .   .   end of zs loop
        enddo

        return
        end
