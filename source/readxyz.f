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
      real*8 hx,hy,hz
      logical exist,opened
      logical done,quit
      logical simple,clash
      logical reorder,second
      logical held,comtitle
      character*3 hname
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
         read (ixyz,20,err=190,end=190)  record
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
      read (string,*,err=190,end=190)  n
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
c     read second and following lines until atom lines begin; the
c     line after the atom count is the comment of a simple XYZ file
c     unless it is a Tinker atom or box line, and a comment that can
c     be read as a simple atom is held in case the comment is missing
c
      call unitcell
      done = .false.
      second = .true.
      held = .false.
      comtitle = .false.
      dowhile (.not. done)
         read (ixyz,50,err=190,end=120)  record
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
               read (string,*,err=70,end=70)  x(1),y(1),z(1)
               read (string,*,err=60,end=60)  x(1),y(1),z(1),type(1),
     &                                        (i12(j,1),j=1,maxval)
   60          continue
               done = .true.
   70          continue
            end if
c
c     check for an initial atom in the simple XYZ format
c
            if (.not.done .and. .not.second) then
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
c     a held comment followed by a non-atom line is the first atom
c
            if (.not.done .and. held) then
               backspace (unit=ixyz)
               goto 120
            end if
c
c     check for optional dimensions of the periodic box, which need
c     three values unless the line ends after fewer of them
c
            if (.not. done) then
               xtmp = 0.0d0
               ytmp = 0.0d0
               ztmp = 0.0d0
               atmp = 0.0d0
               btmp = 0.0d0
               gtmp = 0.0d0
               read (record,*,err=90,end=100)  xtmp,ytmp,ztmp,
     &                                        atmp,btmp,gtmp
               goto 100
   90          continue
               if (ztmp .eq. 0.0d0)  xtmp = 0.0d0
  100          continue
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
c
c     otherwise keep the line after the atom count as a comment
c
            if (.not. done) then
               if (second .and. xtmp.eq.0.0d0) then
                  if (ltitle .eq. 0) then
                     first = nexttext (record)
                     last = trimtext (record)
                     title = record(first:last)
                     ltitle = trimtext (title)
                     comtitle = .true.
                  end if
                  next = 1
                  call getword (record,hname,next)
                  if (next .eq. 1)  goto 110
                  string = record(next:240)
                  read (string,*,err=110,end=110)  hx,hy,hz
                  held = .true.
  110             continue
               end if
            end if
         end if
         second = .false.
      end do
      goto 130
c
c     no atom line follows a held comment, so it was the first atom
c
  120 continue
      if (.not. held)  goto 190
      tag(1) = 1
      name(1) = hname
      x(1) = hx
      y(1) = hy
      z(1) = hz
      if (comtitle) then
         title = ' '
         ltitle = 0
      end if
      held = .false.
      simple = .true.
  130 continue
c
c     read second and following atom lines from input file; if the
c     last line is missing and a comment was held, the file has no
c     comment line and the held comment holds the first atom
c
      if (simple) then
         do i = 2, n
            read (ixyz,140,err=190,end=160)  record
  140       format (a240)
            tag(i) = i
            next = 1
            call getword (record,name(i),next)
            string = record(next:240)
            read (string,*,err=150,end=150)  x(i),y(i),z(i)
         end do
         quit = .false.
         goto 190
  150    continue
         backspace (unit=ixyz)
  160    continue
         if (.not.held .or. i.ne.n)  goto 190
         do k = n, 2, -1
            name(k) = name(k-1)
            x(k) = x(k-1)
            y(k) = y(k-1)
            z(k) = z(k-1)
         end do
         tag(n) = n
         name(1) = hname
         x(1) = hx
         y(1) = hy
         z(1) = hz
         if (comtitle) then
            title = ' '
            ltitle = 0
         end if
         quit = .false.
      else
         do i = 2, n
            read (ixyz,170,err=190,end=190)  record
  170       format (a240)
            read (record,*,err=190,end=190)  tag(i)
            next = 1
            call getword (record,name(i),next)
            string = record(next:240)
            read (string,*,err=180,end=180)  x(i),y(i),z(i),type(i),
     &                                       (i12(j,i),j=1,maxval)
  180       continue
         end do
         quit = .false.
      end if
  190 continue
      if (.not. opened)  close (unit=ixyz)
c
c     an error occurred in reading the coordinate file
c
      if (quit) then
         write (iout,200)  i
  200    format (/,' READXYZ  --  Error in Coordinate File at Atom',i9)
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
                  goto 210
               end if
            end do
  210       continue
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
            write (iout,220)
  220       format (/,' READXYZ  --  Atom Labels not Sequential,',
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
                  if (i12(m,k) .eq. i)  goto 240
               end do
               write (iout,230)  k,i
  230          format (/,' READXYZ  --  Check Connection of Atoms',
     &                    i9,' and',i9)
               call fatal
  240          continue
            end do
         end do
      end if
      return
      end
