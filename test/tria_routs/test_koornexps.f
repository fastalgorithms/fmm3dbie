
      subroutine test_koornf_ders(ipass)
      implicit real *8 (a-h,o-z)
      implicit integer *8 (i-n)
      integer *8 ntest
      real *8, allocatable :: pols(:),ders(:,:)
      real *8, allocatable :: polsex(:),dersex(:,:),ders2ex(:,:)
      real *8, allocatable :: rat1(:,:),rat2(:,:,:),rsc1(:,:)
      real *8 uv(2)
      integer *8 nordert(5)

      call prini(6,13)

      inorder = 5
      nordert(1) = 3
      nordert(2) = 4
      nordert(3) = 10
      nordert(4) = 17
      nordert(5) = 22
      
      ntest = 100
      thresh = 1.0d-12

      ipass = 1

      do i=1,inorder
        norder = nordert(i)
        npols = (norder+1)*(norder+2)/2
        allocate(pols(npols),ders(2,npols),polsex(npols),
     1     dersex(2,npols),rat1(2,0:norder),rat2(3,0:norder,0:norder),
     2     rsc1(0:norder,0:norder))
        call koornf_init(norder,rat1,rat2,rsc1)

        do j=1,ntest
          uv(1) = hkrand(0)
          uv(2) = hkrand(0)*(1-uv(1))

          call koorn_ders(uv,norder,npols,polsex,dersex)
          call koornf_ders(uv,norder,npols,pols,ders,rat1,rat2,rsc1)
          do l=1,npols
            err1 = abs(pols(l)-polsex(l))
            r1 = abs(polsex(l))
            if(r1.gt.1) err1=err1/r1

            err2 = abs(ders(1,l)-dersex(1,l))+
     1         abs(ders(2,l)-dersex(2,l))
            r2 = abs(dersex(1,l))+abs(dersex(2,l))
            if(r2.gt.1) err2 = err2/r2
            errm = max(err1,err2)
            if(errm.gt.thresh) then
              ipass = 0
              call prinf('failed test for*',i,0)
              call prin2('uv=*',uv,2)
              call prinf('pol number=*',l,1)
              call prinf('norder=*',norder,1)
              print *, ders(1,l),dersex(1,l)
              print *, ders(2,l),dersex(2,l)
              print *, pols(l),polsex(l)
              goto 1000
            endif
          enddo
        enddo
        deallocate(pols,ders,polsex,dersex,rat1,rat2,rsc1)
      enddo

      call prinf('successfully cleared fast koorn ders test*',i,0)
 1000 continue     

      return
      end
