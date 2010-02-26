c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prtint  --  output of internal coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtint" writes out a set of Z-matrix internal
c     coordinates to an external disk file
c
c
      subroutine prtint (izmt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'titles.i'
      include 'zclose.i'
      include 'zcoord.i'
      integer i,k,izmt
      logical opened
      character*120 zmtfile
c
c
c     open output unit if not already done
c
      inquire (unit=izmt,opened=opened)
      if (.not. opened) then
         zmtfile = filename(1:leng)//'.int'
         call version (zmtfile,'new')
         open (unit=izmt,file=zmtfile,status='new')
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         write (izmt,10)  n
   10    format (i6)
      else
         write (izmt,20)  n,title(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     output of first three atoms is handled separately
c
      if (n .ge. 1) then
         write (izmt,30)  1,name(1),type(1)
   30    format (i6,2x,a3,i6)
      end if
      if (n .ge. 2) then
         if (digits .le. 6) then
            write (izmt,40)  2,name(2),type(2),iz(1,2),zbond(2)
   40       format (i6,2x,a3,2i6,f10.5)
         else if (digits .le. 8) then
            write (izmt,50)  2,name(2),type(2),iz(1,2),zbond(2)
   50       format (i6,2x,a3,2i6,f12.7)
         else
            write (izmt,60)  2,name(2),type(2),iz(1,2),zbond(2)
   60       format (i6,2x,a3,2i6,f14.9)
         end if
      end if
      if (n .ge. 3) then
         if (digits .le. 6) then
            write (izmt,70)  3,name(3),type(3),iz(1,3),zbond(3),
     &                       iz(2,3),zang(3)
   70       format (i6,2x,a3,2i6,f10.5,i6,f10.4)
         else if (digits .le. 8) then
            write (izmt,80)  3,name(3),type(3),iz(1,3),zbond(3),
     &                       iz(2,3),zang(3)
   80       format (i6,2x,a3,2i6,f12.7,i6,f12.6)
         else
            write (izmt,90)  3,name(3),type(3),iz(1,3),zbond(3),
     &                       iz(2,3),zang(3)
   90       format (i6,2x,a3,2i6,f14.9,i6,f14.8)
         end if
      end if
c
c     convert the torsional angles to lie in standard range
c
      do i = 4, n
         if (iz(4,i) .eq. 0) then
            dowhile (ztors(i) .lt. -180.0d0)
               ztors(i) = ztors(i) + 360.0d0
            end do
            dowhile (ztors(i) .gt. 180.0d0)
               ztors(i) = ztors(i) - 360.0d0
            end do
         end if
      end do
c
c     now, output the fourth through final atoms
c
      if (digits .le. 6) then
         do i = 4, n
            write (izmt,100)  i,name(i),type(i),iz(1,i),zbond(i),
     &                        iz(2,i),zang(i),iz(3,i),ztors(i),iz(4,i)
  100       format (i6,2x,a3,2i6,f10.5,i6,f10.4,i6,f10.4,i6)
         end do
      else if (digits .le. 8) then
         do i = 4, n
            write (izmt,110)  i,name(i),type(i),iz(1,i),zbond(i),
     &                        iz(2,i),zang(i),iz(3,i),ztors(i),iz(4,i)
  110       format (i6,2x,a3,2i6,f12.7,i6,f12.6,i6,f12.6,i6)
         end do
      else
         do i = 4, n
            write (izmt,120)  i,name(i),type(i),iz(1,i),zbond(i),
     &                        iz(2,i),zang(i),iz(3,i),ztors(i),iz(4,i)
  120       format (i6,2x,a3,2i6,f14.9,i6,f14.8,i6,f14.8,i6)
         end do
      end if
c
c     finally, add or delete bonds as required
c
      if (nadd.ne.0 .or. ndel.ne.0)  write (izmt,130)
  130 format ()
      do i = 1, nadd
         write (izmt,140)  (iadd(k,i),k=1,2)
  140    format (2i6)
      end do
      if (ndel .ne. 0)  write (izmt,150)
  150 format ()
      do i = 1, ndel
         write (izmt,160)  (idel(k,i),k=1,2)
  160    format (2i6)
      end do
      if (.not. opened)  close (unit=izmt)
      return
      end
