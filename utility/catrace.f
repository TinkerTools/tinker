c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine catrace  --  build protein alpha-carbon trace  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "catrace" constructs a coordinates file containing only the
c     alpha-carbon atoms of a protein structure
c
c
      program catrace
      use sizes
      use atomid
      use atoms
      use couple
      use files
      use inform
      implicit none
      integer i,j,k,m
      integer jj,kk,frame
      integer ixyz,icat
      integer freeunit
      logical nfind,cfind
      character*240 xyzfile
      character*240 catfile
c
c
c     get the initial coordinates file to be converted
c
      call initial
      call getxyz
      call attach
      call field
      call katom
c
c     reopen coordinates file and open alpha carbon output
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      icat = freeunit ()
      catfile = filename
      call suffix (catfile,'xyz','new')
      open (unit=icat,file=catfile,status ='new')
c
c     extract the alpha carbons from the coordinates file
c
      do while (.not. abort)
         m = 0
         do i = 1, n
            if (atomic(i) .eq. 6) then
               nfind = .false.
               cfind = .false.
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (atomic(k) .eq. 7)  nfind = .true.
                  if (atomic(k) .eq. 6) then
                     do jj = 1, n12(k)
                        kk = i12(jj,k)
                        if (atomic(kk) .eq. 8)  cfind = .true.
                     end do
                  end if
               end do
               if (nfind .and. cfind) then
                  m = m + 1
                  name(m) = 'CA'
                  x(m) = x(i)
                  y(m) = y(i)
                  z(m) = z(i)
                  type(m) = type(i)
               end if
            end if
         end do
c
c     assign connectivities between consecutive alpha carbons
c
         n = m
         do i = 1, n
            n12(i) = 0
            if (i .ne. 1) then
               n12(i) = n12(i) + 1
               i12(n12(i),i) = i - 1
            end if
            if (i .ne. n) then
               n12(i) = n12(i) + 1
               i12(n12(i),i) = i + 1
            end if
         end do
c
c     print alpha-carbon trace and read next coordinates file
c
         call prtxyz (icat)
         call readxyz (ixyz)
      end do
      close (unit=ixyz)
      close (unit=icat)
c
c     perform any final tasks before program exit
c
      call final
      end
