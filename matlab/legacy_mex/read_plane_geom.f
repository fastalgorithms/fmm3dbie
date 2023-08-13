      subroutine get_plane_geom(iref,npatches,norder,npts,norders,ixyzs,
     1   iptype,srcvals,srccoefs,wts)
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: iref,npatches,norder,npts
      integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
      integer, intent(out) :: iptype(npatches)
      real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)
      real *8, intent(out) :: wts(npts)

      character *300 fname
      
      write(fname,'(a,i1,a)') 
     1    '../../fmm3dbie/geometries/A380_Final_o03_r0',iref,'.go3'
      
      call open_gov3_geometry(trim(fname),npatches,norders,ixyzs,
     1   iptype,npts,
     1   srcvals,srccoefs,wts)
      

      
      return
      end
