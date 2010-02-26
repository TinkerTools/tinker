c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtxyz  --  output of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtxyz" writes out a set of Cartesian coordinates
c     to an external disk file
c
c
      subroutine prtxyz (ixyz)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'titles.i'
      integer i,k,ixyz
      logical opened
      character*120 xyzfile
c
c
c     open output unit if not already done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         write (ixyz,10)  n
   10    format (i6)
      else
         write (ixyz,20)  n,title(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     finally, write the coordinates for each atom
c
      if (digits .le. 6) then
         do i = 1, n
            write (ixyz,30)  i,name(i),x(i),y(i),z(i),type(i),
     &                       (i12(k,i),k=1,n12(i))
   30       format (i6,2x,a3,3f12.6,9i6)
         end do
      else if (digits .le. 8) then
         do i = 1, n
            write (ixyz,40)  i,name(i),x(i),y(i),z(i),type(i),
     &                       (i12(k,i),k=1,n12(i))
   40       format (i6,2x,a3,3f14.8,9i6)
         end do
      else
         do i = 1, n
            write (ixyz,50)  i,name(i),x(i),y(i),z(i),type(i),
     &                       (i12(k,i),k=1,n12(i))
   50       format (i6,2x,a3,3f16.10,9i6)
         end do
      end if
      if (.not. opened)  close (unit=ixyz)
      return
      end
