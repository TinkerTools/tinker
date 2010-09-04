c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine readxyz  --  input of Cartesian coordinates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "readxyz" gets a set of Cartesian coordinates from
c     an external disk file
c
c
      subroutine readxyz (ixyz)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,j,k,m,ixyz
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer list(maxatm)
      logical exist,opened
      logical quit,reorder
      logical clash
      character*120 xyzfile
      character*120 record
      character*120 string
c
c
c     initialize the total number of atoms in the system
c
      n = 0
c
c     open the input file if it has not already been done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            open (unit=ixyz,file=xyzfile,status='old')
            rewind (unit=ixyz)
         else
            write (iout,10)
   10       format (/,' READXYZ  --  Unable to Find the Cartesian',
     &                 ' Coordinates File')
            call fatal
         end if
      end if
c
c     read first line and return if already at end of file
c
      quit = .false.
      abort = .true.
      size = 0
      do while (size .eq. 0)
         read (ixyz,20,err=60,end=60)  record
   20    format (a120)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
c
c     parse the title line to get the number of atoms
c
      i = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=60,end=60)  n
c
c     extract the title and determine its length
c
      string = record(next:120)
      first = nexttext (string)
      last = trimtext (string)
      if (last .eq. 0) then
         title = ' '
         ltitle = 0
      else
         title = string(first:last)
         ltitle = trimtext (title)
      end if
c
c     check for too many total atoms in the file
c
      if (n .gt. maxatm) then
         write (iout,30)  maxatm
   30    format (/,' READXYZ  --  The Maximum of',i8,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         type(i) = 0
         do j = 1, maxval
            i12(j,i) = 0
         end do
      end do
c
c     read the coordinates and connectivities for each atom
c
      do i = 1, n
         next = 1
         size = 0
         do while (size .eq. 0)
            read (ixyz,40,err=60,end=60)  record
   40       format (a120)
            size = trimtext (record)
         end do
         read (record,*,err=60,end=60)  tag(i)
         call getword (record,name(i),next)
         string = record(next:120)
         read (string,*,err=50,end=50)  x(i),y(i),z(i),type(i),
     &                                  (i12(j,i),j=1,maxval)
   50    continue
      end do
      quit = .false.
   60 continue
      if (.not. opened)  close (unit=ixyz)
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,70)  i
   70    format (/,' READXYZ  --  Error in Coordinate File at Atom',i6)
         call fatal
      end if
c
c     for each atom, count and sort its attached atoms
c
      do i = 1, n
         n12(i) = 0
         do j = maxval, 1, -1
            if (i12(j,i) .ne. 0) then
               n12(i) = j
               goto 80
            end if
         end do
   80    continue
         call sort (n12(i),i12(1,i))
      end do
c
c     check for scrambled atom order and attempt to renumber
c
      reorder = .false.
      do i = 1, n
         list(tag(i)) = i
         if (tag(i) .ne. i)  reorder = .true.
      end do
      if (reorder) then
         write (iout,90)
   90    format (/,' READXYZ  --  Atom Labels not Sequential,',
     &              ' Attempting to Renumber')
         do i = 1, n
            tag(i) = i
            do j = 1, n12(i)
               i12(j,i) = list(i12(j,i))
            end do
            call sort (n12(i),i12(1,i))
         end do
      end if
c
c     check for atom pairs with identical coordinates
c
      clash = .false.
      call chkxyz (clash)
c
c     make sure that all connectivities are bidirectional
c
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            do m = 1, n12(k)
               if (i12(m,k) .eq. i)  goto 110
            end do
            write (iout,100)  k,i
  100       format (/,' READXYZ  --  Check Connection of Atom',
     &                 i6,' to Atom',i6)
            call fatal
  110       continue
         end do
      end do
      return
      end
