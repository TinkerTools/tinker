c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readxyz  --  input of XYZ-format coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readxyz" gets a set of Cartesian coordinates from an external
c     file in either Tinker XYZ format or simple generic XYZ format
c
c
      subroutine readxyz (ixyz)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use iounit
      use titles
      implicit none
      integer i,j,k,m
      integer ixyz,nmax
      integer next,size
      integer first,last
      integer nexttext
      integer trimtext
      integer, allocatable :: list(:)
      real*8 xtmp,ytmp,ztmp
      real*8 atmp,btmp,gtmp
      logical exist,opened
      logical done,quit
      logical simple,clash
      logical reorder
      character*240 xyzfile
      character*240 record
      character*240 string
c
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
      simple = .false.
      size = 0
      do while (size .eq. 0)
         read (ixyz,20,err=130,end=130)  record
   20    format (a240)
         size = trimtext (record)
      end do
      abort = .false.
      quit = .true.
c
c     parse the first line to get the total number of atoms
c
      n = 0
      next = 1
      call gettext (record,string,next)
      read (string,*,err=130,end=130)  n
c
c     extract the title and determine its length
c
      string = record(next:240)
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
c     check for too few or too many total atoms in the file
c
      if (n .le. 0) then
         write (iout,30)
   30    format (/,' READXYZ  --  The Coordinate File Does Not',
     &              ' Contain Any Atoms')
         call fatal
      else if (n .gt. maxatm) then
         write (iout,40)  maxatm
   40    format (/,' READXYZ  --  The Maximum of',i9,' Atoms',
     &              ' has been Exceeded')
         call fatal
      end if
c
c     initialize coordinates and connectivities for each atom
c
      do i = 1, n
         tag(i) = 0
         name(i) = '   '
         x(i) = 0.0d0
         y(i) = 0.0d0
         z(i) = 0.0d0
         type(i) = 0
         n12(i) = 0
         do j = 1, maxval
            i12(j,i) = 0
         end do
      end do
c
c     read second and following lines until atom lines begin
c
      call unitcell
      done = .false.
      dowhile (.not. done)
         read (ixyz,50,err=130,end=130)  record
   50    format (a240)
         size = trimtext (record)
         if (size .ne. 0) then
c
c     check for an initial atom in the Tinker XYZ format
c
            if (.not. done) then
               read (record,*,err=70,end=70)  tag(1)
               next = 1
               call getword (record,name(1),next)
               if (name(1) .eq. '   ')  goto 70
               string = record(next:240)
               read (string,*,err=60,end=60)  x(1),y(1),z(1),type(1),
     &                                        (i12(j,1),j=1,maxval)
   60          continue
               done = .true.
   70          continue
            end if
c
c     check for an initial atom in the simple XYZ format
c
            if (.not. done) then
               tag(1) = 1
               next = 1
               call getword (record,name(1),next)
               if (next .eq. 1)  goto 80
               string = record(next:240)
               read (string,*,err=80,end=80)  x(1),y(1),z(1)
               done = .true.
               simple = .true.
   80          continue
            end if
c
c     check for optional dimensions of the periodic box
c
            if (.not. done) then
               xtmp = 0.0d0
               ytmp = 0.0d0
               ztmp = 0.0d0
               atmp = 0.0d0
               btmp = 0.0d0
               gtmp = 0.0d0
               read (record,*,err=90,end=90)  xtmp,ytmp,ztmp,
     &                                        atmp,btmp,gtmp
   90          continue
               if (xtmp .ne. 0.0d0) then
                  use_bounds = .true.
                  xbox = xtmp
                  ybox = ytmp
                  zbox = ztmp
                  alpha = atmp
                  beta = btmp
                  gamma = gtmp
                  if (ytmp .eq. 0.0d0)  ybox = xbox
                  if (ztmp .eq. 0.0d0)  zbox = xbox
                  if (atmp .eq. 0.0d0)  alpha = 90.0d0
                  if (btmp .eq. 0.0d0)  beta = 90.0d0
                  if (gtmp .eq. 0.0d0)  gamma = 90.0d0
                  call lattice
               end if
            end if
         end if
      end do
c
c     read second and following atom lines from input file
c
      if (simple) then
         do i = 2, n
            read (ixyz,100,err=130,end=130)  record
  100       format (a240)
            tag(i) = i
            next = 1
            call getword (record,name(i),next)
            string = record(next:240)
            read (string,*,err=130,end=130)  x(i),y(i),z(i)
         end do
         quit = .false.
      else
         do i = 2, n
            read (ixyz,110,err=130,end=130)  record
  110       format (a240)
            read (record,*,err=130,end=130)  tag(i)
            next = 1
            call getword (record,name(i),next)
            string = record(next:240)
            read (string,*,err=120,end=120)  x(i),y(i),z(i),type(i),
     &                                       (i12(j,i),j=1,maxval)
  120       continue
         end do
         quit = .false.
      end if
  130 continue
      if (.not. opened)  close (unit=ixyz)
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,140)  i
  140    format (/,' READXYZ  --  Error in Coordinate File at Atom',i9)
         call fatal
      end if
c
c     for each atom, count and sort its attached atoms
c
      if (.not. abort) then
         do i = 1, n
            do j = maxval, 1, -1
               if (i12(j,i) .ne. 0) then
                  n12(i) = j
                  goto 150
               end if
            end do
  150       continue
            call sort (n12(i),i12(1,i))
         end do
c
c     perform dynamic allocation of some local arrays
c
         nmax = 0
         do i = 1, n
            nmax = max(tag(i),nmax)
            do j = 1, n12(i)
               nmax = max(i12(j,i),nmax)
            end do
         end do
         allocate (list(nmax))
c
c     check for scrambled atom order and attempt to renumber
c
         reorder = .false.
         do i = 1, n
            list(tag(i)) = i
            if (tag(i) .ne. i)  reorder = .true.
         end do
         if (reorder) then
            write (iout,160)
  160       format (/,' READXYZ  --  Atom Labels not Sequential,',
     &                 ' Attempting to Renumber')
            do i = 1, n
               tag(i) = i
               do j = 1, n12(i)
                  i12(j,i) = list(i12(j,i))
               end do
               call sort (n12(i),i12(1,i))
            end do
         end if
c
c     perform deallocation of some local arrays
c
         deallocate (list)
c
c     check for atom pairs with identical coordinates
c
         clash = .false.
         if (n .le. 10000)  call chkxyz (clash)
c
c     make sure all atom connectivities are bidirectional
c
         do i = 1, n
            do j = 1, n12(i)
               k = i12(j,i)
               do m = 1, n12(k)
                  if (i12(m,k) .eq. i)  goto 180
               end do
               write (iout,170)  k,i
  170          format (/,' READXYZ  --  Check Connection of Atoms',
     &                    i9,' and',i9)
               call fatal
  180          continue
            end do
         end do
      end if
      return
      end
