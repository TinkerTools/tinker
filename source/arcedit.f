c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program arcedit  --  create or extract from an archive  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "arcedit" is a utility program for coordinate files which
c     concatenates multiple coordinate sets into a new archive or
c     performs any of several manipulations on an existing archive
c
c
      program arcedit
      use atoms
      use bound
      use files
      use inform
      use iounit
      use output
      use usage
      implicit none
      integer i,j,k,nask
      integer iarc,ixyz,idcd
      integer start,stop
      integer step,size
      integer nmode,mode
      integer lext,lengb
      integer leng1,leng2
      integer now,freeunit
      integer, allocatable :: list(:)
      real*8 xr,yr,zr
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      logical exist,query
      logical first,opened
      character*1 letter
      character*7 ext,modtyp
      character*240 arcfile
      character*240 basename
      character*240 dcdfile
      character*240 xyzfile
      character*240 record
      character*240 string
c
c
c     initialization and set number of archive modifications
c
      call initial
      nmode = 8
c
c     try to get a filename from the command line arguments
c
      call nextarg (arcfile,exist)
      if (exist) then
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end if
c
c     ask for the user specified input archive filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter the Coordinate Archive File Name :  ',$)
         read (input,20)  arcfile
   20    format (a240)
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     get file format type by inspection of first character
c
      coordtype = 'CARTESIAN'
      iarc = freeunit ()
      open (unit=iarc,file=arcfile,status='old')
      rewind (unit=iarc)
      read (iarc,30)  letter
   30 format (a1)
      archive = .false.
      if (letter .eq. ' ')  archive = .true.
      if (letter.ge.'0' .and. letter.le.'9')  archive = .true.
      binary = (.not. archive)
      close (unit=iarc)
c
c     find out which archive modification to perform
c
      mode = -1
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=40,end=40)  mode
         if (mode.ge.0 .and. mode.le.nmode)  query = .false.
      end if
   40 continue
      if (query) then
         write (iout,50)
   50    format (/,' The Tinker Archive File Utility Can :',
     &           //,4x,'(1) Create an Archive from Individual Frames',
     &           /,4x,'(2) Extract Individual Frames from an Archive',
     &           /,4x,'(3) Trim an Archive to Remove Atoms or Frames',
     &           /,4x,'(4) Enforce Periodic Boundaries for a Trajectory'
     &           /,4x,'(5) Unfold Periodic Boundaries for a Trajectory',
     &           /,4x,'(6) Remove Periodic Box Size from a Trajectory',
     &           /,4x,'(7) Convert Tinker Archive to Binary DCD File',
     &           /,4x,'(8) Convert Binary DCD File to Tinker Archive')
         do while (mode.lt.0 .or. mode.gt.nmode)
            mode = -1
            write (iout,60)
   60       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,70,err=80,end=80)  mode
   70       format (i10)
   80       continue
         end do
      end if
c
c     set code for the type of procedure to be performed
c
      if (mode .eq. 0)  modtyp = 'EXIT'
      if (mode .eq. 1)  modtyp = 'CREATE'
      if (mode .eq. 2)  modtyp = 'EXTRACT'
      if (mode .eq. 3)  modtyp = 'TRIM'
      if (mode .eq. 4)  modtyp = 'FOLD'
      if (mode .eq. 5)  modtyp = 'UNFOLD'
      if (mode .eq. 6)  modtyp = 'UNBOUND'
      if (mode .eq. 7) then
         modtyp = 'ARCDCD'
         archive = .true.
         binary = .false.
      end if
      if (mode .eq. 8) then
         modtyp = 'DCDARC'
         archive = .false.
         binary = .true.
      end if
c
c     create and open a new Tinker formatted archive file
c
      if (modtyp .eq. 'CREATE') then
         iarc = freeunit ()
         call basefile (arcfile)
         basename = arcfile
         lengb = leng
         call suffix (arcfile,'arc','new')
         open (unit=iarc,file=arcfile,status='new')
c
c     open an existing Tinker archive file for processing
c
      else if (archive) then
         iarc = freeunit ()
         call basefile (arcfile)
         basename = arcfile
         lengb = leng
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
         do while (.not. exist)
            write (iout,90)
   90       format (/,' Enter the Coordinate Archive File Name :  ',$)
            read (input,100)  arcfile
  100       format (a240)
            call basefile (arcfile)
            basename = arcfile
            lengb = leng
            call suffix (arcfile,'arc','old')
            inquire (file=arcfile,exist=exist)
         end do
         open (unit=iarc,file=arcfile,status='old')
         rewind (unit=iarc)
         call readxyz (iarc)
         rewind (unit=iarc)
         call active
c
c     open an existing binary DCD trajectory file for processing
c
      else if (binary) then
         idcd = freeunit ()
         dcdfile = arcfile
         call basefile (dcdfile)
         basename = dcdfile
         lengb = leng
         call suffix (dcdfile,'dcd','old')
         inquire (file=dcdfile,exist=exist)
         do while (.not. exist)
            write (iout,110)
  110       format (/,' Enter the DCD Binary Archive File Name :  ',$)
            read (input,120)  dcdfile
  120       format (a240)
            call basefile (dcdfile)
            basename = dcdfile
            lengb = leng
            call suffix (dcdfile,'dcd','old')
            inquire (file=dcdfile,exist=exist)
         end do
         call nextarg (xyzfile,exist)
         if (exist) then
            call basefile (xyzfile)
            call suffix (xyzfile,'xyz','old')
            inquire (file=xyzfile,exist=exist)
         end if
         nask = 0
         do while (.not.exist .and. nask.lt.maxask)
            nask = nask + 1
            write (iout,130)
  130       format (/,' Enter Formatted Coordinate File Name :  ',$)
            read (input,140)  xyzfile
  140       format (a240)
            call basefile (xyzfile)
            call suffix (xyzfile,'xyz','old')
            inquire (file=xyzfile,exist=exist)
         end do
         if (.not. exist)  call fatal
         ixyz = freeunit ()
         open (unit=ixyz,file=xyzfile,status='old')
         rewind (unit=ixyz)
         call readxyz (ixyz)
         close (unit=ixyz)
         open (unit=idcd,file=dcdfile,form='unformatted',status='old')
         rewind (unit=idcd)
         first = .true.
         call readdcd (idcd,first)
         rewind (unit=idcd)
         first = .true.
         call active
      end if
c
c     combine individual files into a single archive file
c
      if (modtyp .eq. 'CREATE') then
         modtyp = 'EXIT'
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=150,end=150)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=150,end=150)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=150,end=150)  step
  150    continue
         if (query) then
            write (iout,160)
  160       format (/,' Numbers of First & Last File and Step',
     &                 ' Increment :  ',$)
            read (input,170)  record
  170       format (a240)
            read (record,*,err=180,end=180)  start,stop,step
  180       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
c
c     cycle over the user specified coordinate files
c
         i = start
         do while (i.ge.start .and. i.le.stop)
            ixyz = freeunit ()
            lext = 3
            call numeral (i,ext,lext)
            xyzfile = basename(1:lengb)//'.'//ext(1:lext)
            call version (xyzfile,'old')
            inquire (file=xyzfile,exist=exist)
            if (.not.exist .and. i.lt.100) then
               lext = 2
               call numeral (i,ext,lext)
               xyzfile = basename(1:lengb)//'.'//ext(1:lext)
               call version (xyzfile,'old')
               inquire (file=xyzfile,exist=exist)
            end if
            if (.not.exist .and. i.lt.10) then
               lext = 1
               call numeral (i,ext,lext)
               xyzfile = basename(1:lengb)//'.'//ext(1:lext)
               call version (xyzfile,'old')
               inquire (file=xyzfile,exist=exist)
            end if
            if (exist) then
               open (unit=ixyz,file=xyzfile,status='old')
               rewind (unit=ixyz)
               call readxyz (ixyz)
               close (unit=ixyz)
               if (i .eq. start)  call active
               nuse = n
               do j = 1, n
                  use(j) = .true.
               end do
               call prtarc (iarc)
            end if
            i = i + step
         end do
      end if
c
c     perform dynamic allocation of some local arrays
c
      if (modtyp .eq. 'TRIM') then
         size = 40
         allocate (list(size))
      end if
c
c     decide whether atoms are to be removed from each frame
c
      if (modtyp .eq. 'TRIM') then
         call active
         if (nuse .eq. n) then
            do i = 1, size
               list(i) = 0
            end do
            i = 1
            query = .true.
            call nextarg (string,exist)
            if (exist) then
               do while (i .le. size)
                  read (string,*,err=190,end=190)  list(i)
                  if (list(i) .eq. 0)  goto 190
                  i = i + 1
                  call nextarg (string,exist)
               end do
  190          continue
               query = .false.
            end if
            if (query) then
               write (iout,200)
  200          format (/,' Numbers of the Atoms to be Removed :  ',$)
               read (input,210)  record
  210          format (a240)
               read (record,*,err=220,end=220)  (list(i),i=1,size)
  220          continue
            end if
            i = 1
            do while (list(i) .ne. 0)
               list(i) = max(-n,min(n,list(i)))
               if (list(i) .gt. 0) then
                  k = list(i)
                  if (use(k)) then
                     use(k) = .false.
                     nuse = nuse - 1
                  end if
                  i = i + 1
               else
                  list(i+1) = max(-n,min(n,list(i+1)))
                  do k = abs(list(i)), abs(list(i+1))
                     if (use(k)) then
                        use(k) = .false.
                        nuse = nuse - 1
                     end if
                  end do
                  i = i + 2
               end if
            end do
         end if
c
c     perform deallocation of some local arrays
c
         deallocate (list)
      end if
c
c     store index to use in renumbering the untrimmed atoms
c
      k = 0
      do i = 1, n
         iuse(i) = 0
         if (use(i)) then
            k = k + 1
            iuse(i) = k
         end if
      end do
c
c     convert Tinker archive to binary DCD trajectory file
c
      if (modtyp .eq. 'ARCDCD') then
         modtyp = 'EXIT'
         first = .true.
         idcd = freeunit ()
         dcdfile = filename(1:leng)//'.dcd'
         call version (dcdfile,'new')
         open (unit=idcd,file=dcdfile,form='unformatted',status='new')
         do while (.true.)
            call readxyz (iarc)
            if (abort)  goto 230
            call prtdcd (idcd,first)
         end do
  230    continue
         close (unit=idcd)
      end if
c
c     convert binary DCD trajectory file to Tinker archive
c
      if (modtyp .eq. 'DCDARC') then
         modtyp = 'EXIT'
         first = .true.
         iarc = freeunit ()
         arcfile = filename(1:leng)//'.arc'
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
         do while (.true.)
            call readdcd (idcd,first)
            if (abort)  goto 240
            call prtarc (iarc)
         end do
  240    continue
         close (unit=iarc)
      end if
c
c     perform dynamic allocation of some local arrays
c
      if (modtyp .ne. 'EXIT') then
         allocate (xold(n))
         allocate (yold(n))
         allocate (zold(n))
      end if
c
c     get the initial and final coordinate frames to process
c
      if (modtyp .ne. 'EXIT') then
         now = 1
         leng1 = 1
         leng2 = leng
         do i = 1, leng
            if (filename(i:i) .eq. '/')  leng1 = i+1
            if (filename(i:i) .eq. ']')  leng1 = i+1
            if (filename(i:i) .eq. ':')  leng1 = i+1
         end do
         do i = leng, leng1, -1
            if (filename(i:i) .eq. '.')  leng2 = i-1
         end do
         leng = leng2 - leng1 + 1
         filename(1:leng) = filename(leng1:leng2)
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=250,end=250)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  step
  250    continue
         if (query) then
            write (iout,260)
  260       format (/,' Numbers of First & Last File and Step',
     &                 ' [<Enter>=Exit] :  ',$)
            read (input,270)  record
  270       format (a240)
            read (record,*,err=280,end=280)  start,stop,step
  280       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
c
c     loop over the individual coordinate files to be extracted
c
         do while (start .ne. 0)
            if (start .le. now) then
               now = 1
               first = .true.
               if (archive)  rewind (unit=iarc)
               if (binary)  rewind (unit=idcd)
            end if
            do k = 1, start-now
               call readcart (iarc,first)
            end do
            i = start
            if (modtyp .eq. 'EXTRACT') then
               do while (i.ge.start .and. i.le.stop)
                  lext = 3
                  call numeral (i,ext,lext)
                  call readcart (iarc,first)
                  if (abort)  goto 290
                  nuse = n
                  do j = 1, n
                     use(j) = .true.
                  end do
                  ixyz = freeunit ()
                  xyzfile = filename(1:leng)//'.'//ext(1:lext)
                  call version (xyzfile,'new')
                  if (archive) then
                     open (unit=ixyz,file=xyzfile,status='new')
                  else if (binary) then
                     open (unit=ixyz,file=xyzfile,form='unformatted',
     &                        status='new')
                  end if
                  first = .true.
                  if (archive)  call prtarc (ixyz)
                  if (binary)  call prtdcd (ixyz,first)
                  close (unit=ixyz)
                  i = i + step
                  do k = 1, step-1
                     call readcart (iarc,first)
                  end do
               end do
            else
               ixyz = freeunit ()
               xyzfile = basename
               call suffix (xyzfile,'arc','new')
               if (archive) then
                  open (unit=ixyz,file=xyzfile,status='new')
               else if (binary) then
                  open (unit=ixyz,file=xyzfile,form='unformatted',
     &                     status='new')
               end if
               do while (i.ge.start .and. i.le.stop)
                  if (modtyp .eq. 'UNBOUND')  use_bounds = .true.
                  call readcart (iarc,first)
                  if (abort)  goto 290
                  if (modtyp .eq. 'FOLD') then
                     call unitcell
                     if (use_bounds) then
                        call lattice
                        call molecule
                        call bounds
                     end if
                  else if (modtyp .eq. 'UNFOLD') then
                     nuse = n
                     do j = 1, n
                        use(j) = .true.
                     end do
                     if (i .eq. start) then
                        call unitcell
                        do j = 1, n
                           xold(j) = x(j)
                           yold(j) = y(j)
                           zold(j) = z(j)
                        end do
                     end if
                     call lattice
                     do j = 1, n
                        xr = x(j) - xold(j)
                        yr = y(j) - yold(j)
                        zr = z(j) - zold(j)
                        if (use_bounds)  call image (xr,yr,zr)
                        x(j) = xold(j) + xr
                        y(j) = yold(j) + yr
                        z(j) = zold(j) + zr
                        xold(j) = x(j)
                        yold(j) = y(j)
                        zold(j) = z(j)
                     end do
                  end if
                  if (i .eq. start)  first = .true.
                  if (modtyp .eq. 'UNBOUND')  use_bounds = .false.
                  if (archive)  call prtarc (ixyz)
                  if (binary)  call prtdcd (ixyz,first)
                  i = i + step
                  do k = 1, step-1
                     call readcart (iarc,first)
                  end do
               end do
               close (unit=ixyz)
            end if
  290       continue
            now = stop
            start = 0
            stop = 0
            step = 0
            query = .true.
            call nextarg (string,exist)
            if (exist) then
               read (string,*,err=300,end=300)  start
               query = .false.
            end if
            call nextarg (string,exist)
            if (exist)  read (string,*,err=300,end=300)  stop
            call nextarg (string,exist)
            if (exist)  read (string,*,err=300,end=300)  step
  300       continue
            if (query) then
               write (iout,310)
  310          format (/,' Numbers of First & Last File and Step',
     &                    ' [<Enter>=Exit] :  ',$)
               read (input,320)  record
  320          format (a240)
               read (record,*,err=330,end=330)  start,stop,step
  330          continue
            end if
            if (stop .eq. 0)  stop = start
            if (step .eq. 0)  step = 1
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (xold)
         deallocate (yold)
         deallocate (zold)
      end if
c
c     perform any final tasks before program exit
c
      inquire (unit=iarc,opened=opened)
      if (opened)  close (unit=iarc)
      inquire (unit=idcd,opened=opened)
      if (opened)  close (unit=idcd)
      call final
      end
