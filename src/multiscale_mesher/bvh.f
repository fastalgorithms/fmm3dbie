ccccccccccccccccccccccccccccccccccccccccc
      integer function delta(n,mid,leaves_sorted,idx)
      implicit none
      integer n,idx
      integer mid(n)
      integer leaves_sorted(n)

      delta = xor(mid(leaves_sorted(idx+1)),mid(leaves_sorted(idx)))

      return
      end
ccccccccccccccccccccccccccccccccccccccccc

      module bvh
      contains
ccccccccccccccccccccccccccccccccccccccccc
      subroutine normalize(n, p_aabb, p, aabb)
      implicit none
      integer n,i,j
      double precision p_aabb(3,n), p(3,n), aabb(3,2)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j)
      do i = 1,n
        do j = 1,3
          p_aabb(j,i) = (p(j,i) - aabb(j,1)) / (aabb(j,2) - aabb(j,1))
        enddo
      enddo
C$OMP END PARALLEL DO

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine overlap_aabb(aabb1, aabb2,if_overlap)
      implicit none
      logical if_overlap
      double precision aabb1(3,2), aabb2(3,2)
      integer i

      do i = 1,3
        if(aabb1(i,1).gt.aabb2(i,2) .or. aabb1(i,2).lt.aabb2(i,1)) then
          if_overlap = .false.
          return
        endif
      enddo

      if_overlap = .true.
      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine expand_bits(v)
      implicit none
      integer n1,n2,n3,n4,n5,n6,n7,n8
      integer v
      DATA n1 /Z'00010001'/, n2 /Z'FF0000FF'/
      DATA n3 /Z'00000101'/, n4 /Z'0F00F00F'/
      DATA n5 /Z'00000011'/, n6 /Z'C30C30C3'/
      DATA n7 /Z'00000005'/, n8 /Z'49249249'/

c      print *, "n1"
c      call print_bits_32(n1)
c      print *, "n2"
c      call print_bits_32(n2)
c      print *, "n3"
c      call print_bits_32(n3)
c      print *, "n4"
c      call print_bits_32(n4)
c      print *, "n5"
c      call print_bits_32(n5)
c      print *, "n6"
c      call print_bits_32(n6)
c      print *, "n7"
c      call print_bits_32(n7)
c      print *, "n8"
c      call print_bits_32(n8)
c
c      print *, "v0"
c      call print_bits_32(v)

      v = and(v * n1, n2)
c      print *, "v1"
c      call print_bits_32(v)

      v = and(v * n3, n4)
c      print *, "v2"
c      call print_bits_32(v)

      v = and(v * n5, n6)
c      print *, "v3"
c      call print_bits_32(v)

      v = and(v * n7, n8)
c      print *, "v4"
c      call print_bits_32(v)

c      v = and(v * X'00010001',X'FF0000FF')
c      v = and(v * X'00000101',X'0F00F00F')
c      v = and(v * X'00000011',X'C30C30C3')
c      v = and(v * X'00000005',X'49249249')

c      v = (v * 0x00010001u) & 0xFF0000FFu;
c      v = (v * 0x00000101u) & 0x0F00F00Fu;
c      v = (v * 0x00000011u) & 0xC30C30C3u;
c      v = (v * 0x00000005u) & 0x49249249u;
      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine interleave_morton_3d(ix,iy,iz,m)
      implicit none
      integer ix,iy,iz
      integer m
       
      call expand_bits(ix)
      call expand_bits(iy)
      call expand_bits(iz)
     
      m = or(or(lshift(ix,2),lshift(iy,1)),iz)

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine morton_3d(n,p,m,aabb)
      implicit none
      integer n
      double precision p(3,n)
      double precision aabb(3,2)
      double precision bsize(3)
      integer m(n)
      integer istart, iend
      integer ix,iy,iz,i
      integer nthreads
      integer, external :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
      integer, external :: OMP_GET_MAX_THREADS
      

      nthreads = 1
C$    nthreads = OMP_GET_MAX_THREADS()
      print *, nthreads

c     compute size of the global bounding box
      do i = 1,3
        bsize(i) = aabb(i,2) - aabb(i,1)
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,ix,iy,iz)
      do i = 1,n 
        ix = min(max(int((p(1,i)-aabb(1,1))/bsize(1)*1024.0),0),1023)
        iy = min(max(int((p(2,i)-aabb(2,1))/bsize(2)*1024.0),0),1023)
        iz = min(max(int((p(3,i)-aabb(3,1))/bsize(3)*1024.0),0),1023)
        call interleave_morton_3d(ix,iy,iz,m(i))
      enddo
C$OMP END PARALLEL DO

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine sort_morton_3d(n,m,idx)
      implicit none
      integer n,m(n),idx(n)

      call sorti_para(n,m,idx)

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine is_leaf(idx, is_leaf_node, child_idx)
      implicit none
      integer idx, child_idx
      logical is_leaf_node

      if(idx .lt. 0) then
        is_leaf_node = .true.
        child_idx = -idx
      else
        is_leaf_node = .false.
        child_idx = idx
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine bvh_overlap_aabb(n,leaf_box,root,bbox,nodes_box,
     1           l_child,r_child,aabb,col_list,ncol)
      implicit none
      integer n,root,ncol
      double precision aabb(3,2),leaf_box(3,2*n)
      double precision bbox(3,2)
      double precision nodes_box(3,2*(n-1))
      integer l_child(n-1)
      integer r_child(n-1)
      integer col_list(*)
      integer stack(64),stack_ptr,invalid,node
      logical overlap_l, overlap_r
      logical is_l_leaf, is_r_leaf
      integer col_cnt, l_child_idx, r_child_idx
      logical traverse_l, traverse_r

      col_cnt = 0
      stack_ptr = 1
      stack(stack_ptr) = invalid
      stack_ptr = stack_ptr + 1

c     node must be internal node(not leaf) for now, to fix
      node = root
      do while (node .ne. invalid)
        call is_leaf(l_child(node), is_l_leaf, l_child_idx)
        call is_leaf(r_child(node), is_r_leaf, r_child_idx)

        call overlap_aabb(aabb,nodes_box(1,2*l_child_idx-1),overlap_l)
        call overlap_aabb(aabb,nodes_box(1,2*r_child_idx-1),overlap_r)
        

        if(overlap_l .and. is_l_leaf) then
          col_cnt = col_cnt + 1
          col_list(col_cnt) = l_child_idx
        endif
        if(overlap_r .and. is_r_leaf) then
          col_cnt = col_cnt + 1
          col_list(col_cnt) = r_child_idx
        endif

        traverse_l = (overlap_l .and. .not.is_l_leaf)
        traverse_r = (overlap_r .and. .not.is_r_leaf)
        if(.not.traverse_l .and. .not.traverse_r) then
          stack_ptr = stack_ptr - 1
          node = stack(stack_ptr)
        else
          if(traverse_l) then
            node = l_child_idx
          else
            node = r_child_idx
          endif
          if(traverse_l .and. traverse_r) then
            stack(stack_ptr) = r_child_idx
            stack_ptr = stack_ptr + 1
          endif
        endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine bvh_overlap_aabbs(n,leaf_box,root,bbox,nodes_box,
     1           l_child,r_child,m,aabbs,col_pair,npair)
      implicit none
      integer n, m, root, npair, i, j, cnt
      double precision leaf_box(3,2*n), aabbs(3,2*m)
      double precision bbox(3,2)
      double precision nodes_box(3,2*(n-1))
      integer l_child(n-1)
      integer r_child(n-1)
      integer, allocatable :: i_npair(:)
      integer, allocatable :: col_pair(:,:), col_list(:,:)
      
      allocate(col_list(20,m))
      allocate(i_npair(m))
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,i_npair)
      do i = 1,m
        call bvh_overlap_aabb(n,leaf_box,root,bbox,nodes_box,
     1       l_child,r_child,aabbs(1,2*i-1),col_list(1,i),i_npair(i))

      enddo
C$OMP END PARALLEL DO

      npair = 0

      do i = 1,m
        npair = npair + i_npair(i)
      enddo
      
      allocate(col_pair(2,npair))

      cnt = 0
      do i = 1,m
        do j = 1,i_npair(i)
          col_pair(1,cnt+j) = i
          col_pair(2,cnt+j) = col_list(j,i)
        enddo
        cnt = cnt + i_npair(i)
      enddo


      deallocate(i_npair)
      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine bvh_init(n,aabbs)
      implicit none
      integer n
      double precision aabbs(3,2*n)

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine init_bbox(bbox)
      implicit none
      double precision bbox(3,2)
      integer i

      do i = 1,3
c        bbox(i,1) =  1.797693D+308
c        bbox(i,2) = -1.797693D+308
        bbox(i,1) =  huge(bbox(1,1))
        bbox(i,2) = -huge(bbox(1,1))
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
c     omp thread level atomic
      subroutine update_min_atomic(min_val,val)
      implicit none
      double precision min_val, val

C$OMP ATOMIC
      min_val = max(min_val,val)
C$OMP END ATOMIC

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
c     omp thread level atomic
      subroutine update_max_atomic(max_val,val)
      implicit none
      double precision max_val, val

C$OMP ATOMIC
      max_val = max(max_val,val)
C$OMP END ATOMIC

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine update_min(min_val,val)
      implicit none
      double precision min_val, val

      min_val = min(min_val, val)

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine update_max(max_val,val)
      implicit none
      double precision  max_val, val

      max_val = max(max_val, val)

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine atomic_exchange(x,v,prev)
      implicit none
      integer x,v,prev

C$OMP ATOMIC CAPTURE
      prev = x
      x = v
C$OMP END ATOMIC

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
c     expand bbox with box
      subroutine expand_bbox_atomic(bbox,box)
      implicit none
      double precision bbox(3,2), box(3,2)
      integer i

      do i = 1,3
        call update_min_atomic(bbox(i,1), box(i,1))
        call update_max_atomic(bbox(i,2), box(i,2))
      enddo

      return
      end


ccccccccccccccccccccccccccccccccccccccccc
      subroutine copy_box(box_src,box_trg)
      implicit none
      integer i
      double precision box_src(6), box_trg(6)

      do i = 1,6
        box_trg(i) = box_src(i)
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine compute_centers(n,aabbs,centers)
      implicit none
      integer n,i,j
      double precision aabbs(3,2*n),centers(3,n)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j)
      do i = 1,n
        do j = 1,3
          centers(j,i) = (aabbs(j,2*i-1)+aabbs(j,2*i))/2.0
        enddo
      enddo
C$OMP END PARALLEL DO

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine bvh_build(n,aabbs,root,bbox,nodes_box,l_child,r_child)
      implicit none
      integer, external :: delta
      integer n,i,j,root,leaf_id,left,right,current,idx,nnodes
      integer previous,invalid,boxid,d1,d2,parent
      double precision aabbs(3,2*n)
      double precision bbox(3,2)
      double precision box(3,2)
      integer, allocatable :: mid(:),leaves_sorted(:)
c      integer, allocatable :: leaves_id(:)
      integer, allocatable :: l_child(:),r_child(:)
      integer, allocatable :: other_bounds(:)
      double precision, allocatable :: centers(:,:), nodes_box(:,:)
      logical is_leaf

      invalid = huge(invalid)
      nnodes = n-1

      call init_bbox(bbox)
c     calculate the global bounding box
c     todo: omp it
      do i = 1,n
        do j = 1,3
          call update_min(bbox(j,1), aabbs(j,2*i-1))
          call update_max(bbox(j,2), aabbs(j,2*i))
        enddo
      enddo

      allocate(centers(3,n))
      allocate(mid(n), leaves_sorted(n))
c      allocate(mid(n), leaves_sorted(n), leaves_id(n))
      allocate(l_child(nnodes), r_child(nnodes))
      allocate(nodes_box(3,2*(nnodes)))
      allocate(other_bounds(nnodes))

      do i = 1,nnodes
        other_bounds(i) = invalid
        call init_bbox(nodes_box(1,2*i-1))
      enddo
      call compute_centers(n,aabbs,centers)
      call morton_3d(n,centers,mid,bbox)
      call sort_morton_3d(n,mid,leaves_sorted)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(leaf_id,current,boxid,box,idx,previous,left,right,parent)
      do leaf_id = 1,n
        current = leaf_id

        boxid = leaves_sorted(leaf_id)
        call copy_box(aabbs(1,2*boxid-1),box)
        is_leaf = .true.
        do
          if(left.eq.1 .and. right.eq.n) then
            root = current
            exit
          endif

c         negative for leaf index and positive for internal node index
          if(is_leaf) then
            idx = -leaves_sorted(current)
          else
            idx =  current
          endif

c         choose parent
          if((left.eq.1) .or. ((right.ne.n) .and.
     1       (delta(n,mid,leaves_sorted,right) .lt.
     2        delta(n,mid,leaves_sorted,left-1)))) then
            parent =  right
            call atomic_exchange(other_bounds(parent),left,previous)
            if(invalid .ne. previous) then
              right = previous
            endif
            l_child(parent) = idx
          else
            parent = left - 1
            call atomic_exchange(other_bounds(parent),right,previous)
            if(invalid .ne. previous) then
              left = previous
            endif
            r_child(parent) = idx
          endif

c         expand current node box
          call expand_bbox_atomic(nodes_box(1,parent*2-1),box)

c         terminate this thread
          if(invalid .eq. previous) then
            exit
          endif
          
          current =  parent
          call copy_box(nodes_box(1,parent*2-1),box)
          is_leaf = .false.

        enddo

      enddo
C$OMP END PARALLEL DO

      return
      end

ccccccccccccccccccccccccccccccccccccccccc
      subroutine print_bits_32(val)
      implicit none
      integer val
      character(len=32) :: bin

      write(bin, '(B32.32)') val
      print *, bin(1:32)


      return
      end
ccccccccccccccccccccccccccccccccccccccccc
      end module

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
      module alloc
      contains
      subroutine alloc_n_ints(n,arr)
        implicit none
        integer n,i
        integer, allocatable :: arr(:)
        allocate(arr(n))
        do i = 1,n
          arr(i) = 2*i
        enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
      subroutine dealloc_n_ints(arr)
        implicit none
        integer, allocatable :: arr(:)
        deallocate(arr)
      return
      end
ccccccccccccccccccccccccccccccccccccccccc
      end module

ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc
c
c
c     main driver
c     testing
c
c
c
      use alloc
      use bvh
      implicit none
      character(4) b
      character c
      integer a,i
      integer  M, ICOUNT/1/, JCOUNT
      integer x,y,z,mid
      integer, allocatable :: arr(:),indx(:)
      integer n
      double precision d


      REAL*4  TEMP
      character(len=32) :: bin

      print *, "parallel bvh"

c      JCOUNT = B'01000000000000000000000000000000'
c      JCOUNT = X'0000001'
c      JCOUNT = or(JCOUNT,X'0000001')
      JCOUNT = 1023
      print *, "jcount 1"
      call print_bits_32(JCOUNT)

      call expand_bits(JCOUNT)
      print *, "jcount 2"
      call print_bits_32(JCOUNT)
      print *, "jcount 3"
      JCOUNT = lshift(JCOUNT,2)
      call print_bits_32(JCOUNT)

      d = 0.5

c      x = 1023
c      y = 1023
c      z = 1023
c      call interleave_morton_3d(x,y,z,mid)
c      call print_bits_32(mid)
c
c      n = 100
c      allocate(arr(n),indx(n))
c      arr = 0
c      arr(1) = 100
c      arr(2) = 20
c      call sorti_para(n,arr,indx)
c
c      print *, indx(1:n)
c
c
c      print *, sizeof(n)

      x = min(int(0.9999*1024),1023)
      d = 1.797693D+308
      n = huge(n)
      d = huge(d)
      print *, n
      print *, d
      n = 10
      
      call alloc_n_ints(n, arr)
      do i = 1,n
        print *, arr(i)
      enddo
      call dealloc_n_ints(arr)

      return
      end
