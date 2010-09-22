c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program archive  --  create or extract from an archive  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "archive" is a utility program for coordinate files which
c     concatenates multiple coordinate sets into a single archive
c     file, or extracts individual coordinate sets from an archive
c
c
      program archive
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'usage.i'
      integer i,k,iarc,ixyz
      integer start,stop
      integer step,now
      integer lext,next
      integer lengb
      integer leng1,leng2
      integer freeunit
      integer list(20)
      logical exist,query
      character*1 answer
      character*7 ext,mode
      character*120 arcfile
      character*120 basename
      character*120 xyzfile
      character*120 record
      character*120 string
c
c
c     get the name to use for the coordinate archive file
c
      call initial
      call nextarg (arcfile,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter Name of the Coordinate Archive File :  ',$)
         read (input,20)  arcfile
   20    format (a120)
      end if
c
c     decide whether to create or extract from an archive file
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Create (C), Extract (E) from or Trim (T)',
     &              ' an Archive [C] :  ',$)
         read (input,40)  record
   40    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'E') then
         mode = 'EXTRACT'
      else if (answer .eq. 'T') then
         mode = 'TRIM'
      else
         mode = 'CREATE'
      end if
c
c     open the archive file to be processed
c
      iarc = freeunit ()
      call basefile (arcfile)
      basename = arcfile
      lengb = leng
      call suffix (arcfile,'arc')
      if (mode .eq. 'CREATE') then
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
      else
         call version (arcfile,'old')
      end if
c
c     concatenate individual files into a single archive file
c
      if (mode .eq. 'CREATE') then
         start = 0
         stop = 0
         step = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=50,end=50)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=50,end=50)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=50,end=50)  step
   50    continue
         if (query) then
            write (iout,60)
   60       format (/,' Numbers of First & Last File and Step',
     &                 ' Increment :  ',$)
            read (input,70)  record
   70       format (a120)
            read (record,*,err=80,end=80)  start,stop,step
   80       continue
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
               if (i .eq. start) then
                  call active
                  nuse = n
                  do k = 1, nuse
                     use(k) = .true.
                  end do
               end if
               call prtarc (iarc)
            end if
            i = i + step
         end do
      end if
c
c     extract individual files from a concatenated archive file
c
      if (mode.eq.'EXTRACT' .or. mode.eq.'TRIM') then
         inquire (file=arcfile,exist=exist)
         do while (.not. exist)
            write (iout,90)
   90       format (/,' Enter Name of the Coordinate Archive',
     &                 ' File :  ',$)
            read (input,100)  arcfile
  100       format (a120)
            call basefile (arcfile)
            basename = arcfile
            lengb = leng
            call suffix (arcfile,'arc')
            call version (arcfile,'old')
            inquire (file=arcfile,exist=exist)
         end do
         open (unit=iarc,file=arcfile,status='old')
         rewind (unit=iarc)
         call readxyz (iarc)
         rewind (unit=iarc)
c
c     decide whether atoms are to be removed from each frame
c
         call active
         if (mode .eq. 'EXTRACT') then
            nuse = n
            do i = 1, nuse
               use(i) = .true.
            end do
         else if (mode.eq.'TRIM' .and. nuse.eq.n) then
            do i = 1, 20
               list(i) = 0
            end do
            write (iout,110)
  110       format (/,' Numbers of the Atoms to be Removed :  ',$)
            read (input,120)  record
  120       format (a120)
            read (record,*,err=130,end=130)  (list(i),i=1,20)
  130       continue
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
c     get the initial and final coordinate frames to extract
c
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
            read (string,*,err=140,end=140)  start
            query = .false.
         end if
         call nextarg (string,exist)
         if (exist)  read (string,*,err=140,end=140)  stop
         call nextarg (string,exist)
         if (exist)  read (string,*,err=140,end=140)  step
  140    continue
         if (query) then
            write (iout,150)
  150       format (/,' Numbers of First & Last File and Step',
     &                 ' [<CR>=Exit] :  ',$)
            read (input,160)  record
  160       format (a120)
            read (record,*,err=170,end=170)  start,stop,step
  170       continue
         end if
         if (stop .eq. 0)  stop = start
         if (step .eq. 0)  step = 1
c
c     loop over the individual coordinate files to be extracted
c
         do while (start .ne. 0)
            if (start .le. now) then
               now = 1
               rewind (unit=iarc)
            end if
            do k = 1, start-now
               call readxyz (iarc)
            end do
            i = start
            if (mode .eq. 'EXTRACT') then
               do while (i.ge.start .and. i.le.stop)
                  lext = 3
                  call numeral (i,ext,lext)
                  call readxyz (iarc)
                  if (abort)  goto 180
                  ixyz = freeunit ()
                  xyzfile = filename(1:leng)//'.'//ext(1:lext)
                  call version (xyzfile,'new')
                  open (unit=ixyz,file=xyzfile,status='new')
                  call prtarc (ixyz)
                  close (unit=ixyz)
                  i = i + step
                  do k = 1, step-1
                     call readxyz (iarc)
                  end do
               end do
            else if (mode .eq. 'TRIM') then
               ixyz = freeunit ()
               xyzfile = basename
               call suffix (xyzfile,'arc')
               call version (xyzfile,'new')
               open (unit=ixyz,file=xyzfile,status='new')
               do while (i.ge.start .and. i.le.stop)
                  call readxyz (iarc)
                  if (abort)  goto 180
                  call prtarc (ixyz)
                  i = i + step
                  do k = 1, step-1
                     call readxyz (iarc)
                  end do
               end do
               close (unit=ixyz)
            end if
  180       continue
            now = stop
            start = 0
            stop = 0
            step = 0
            query = .true.
            call nextarg (string,exist)
            if (exist) then
               read (string,*,err=190,end=190)  start
               query = .false.
            end if
            call nextarg (string,exist)
            if (exist)  read (string,*,err=190,end=190)  stop
            call nextarg (string,exist)
            if (exist)  read (string,*,err=190,end=190)  step
  190       continue
            if (query) then
               write (iout,200)
  200          format (/,' Numbers of First & Last File and Step',
     &                    ' [<CR>=Exit] :  ',$)
               read (input,210)  record
  210          format (a120)
               read (record,*,err=220,end=220)  start,stop,step
  220          continue
            end if
            if (stop .eq. 0)  stop = start
            if (step .eq. 0)  step = 1
         end do
      end if
c
c     perform any final tasks before program exit
c
      close (unit=iarc)
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtarc  --  output of a TINKER archive file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtarc" writes out a set of Cartesian coordinates for
c     all active atoms in the TINKER XYZ archive format
c
c
      subroutine prtarc (iarc)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'titles.i'
      include 'usage.i'
      integer i,k,iarc
      logical opened
      character*120 arcfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         arcfile = filename(1:leng)//'.arc'
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         write (iarc,10)  nuse
   10    format (i6)
      else
         write (iarc,20)  nuse,title(1:ltitle)
   20    format (i6,2x,a)
      end if
c
c     finally, write the coordinates for each atom
c
      if (digits .le. 6) then
         do i = 1, n
            if (use(i)) then
               write (iarc,30)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   30          format (i6,2x,a3,3f12.6,9i6)
            end if
         end do
      else if (digits .le. 8) then
         do i = 1, n
            if (use(i)) then
               write (iarc,40)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   40          format (i6,2x,a3,3f14.8,9i6)
            end if
         end do
      else
         do i = 1, n
            if (use(i)) then
               write (iarc,50)  i,name(i),x(i),y(i),z(i),type(i),
     &                          (i12(k,i),k=1,n12(i))
   50          format (i6,2x,a3,3f16.10,9i6)
            end if
         end do
      end if
      if (.not. opened)  close (unit=iarc)
      return
      end
