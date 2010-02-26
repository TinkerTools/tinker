c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program superpose  --  optimal coordinate superposition  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "superpose" takes pairs of structures and superimposes them
c     in the optimal least squares sense; it will attempt to match
c     all atom pairs or only those specified by the user
c
c
      program superpose
      implicit none
      include 'sizes.i'
      include 'align.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'titles.i'
      integer i,ixyz,next
      integer n1,i1,n2,i2
      integer leng1,leng2
      integer ifile1,ifile2
      integer frame1,frame2
      integer start,stop
      integer option,delta
      integer trimtext,freeunit
      integer range(4)
      integer atomic1(maxatm)
      integer atomic2(maxatm)
      real*8 dist,cutoff
      real*8 rmsvalue
      real*8 mass1(maxatm)
      real*8 mass2(maxatm)
      real*8 x1(maxatm),x2(maxatm)
      real*8 y1(maxatm),y2(maxatm)
      real*8 z1(maxatm),z2(maxatm)
      logical header,exist
      logical query,skip,same
      character*1 answer
      character*3 name1(maxatm)
      character*3 name2(maxatm)
      character*120 file1,file2
      character*120 xyzfile
      character*120 record
      character*120 string
c
c
c     get atom names and masses for the first structure type
c
      call initial
      call getxyz
      call field
      call katom
      file1 = filename
      leng1 = trimtext (file1)
      n1 = n
      do i = 1, n1
         name1(i) = name(i)
         atomic1(i) = atomic(i)
         mass1(i) = mass(i)
      end do
c
c     get atom names and masses for the second structure type
c
      call getxyz
      call field
      call katom
      file2 = filename
      leng2 = trimtext (file2)
      n2 = n
      do i = 1, n2
         name2(i) = name(i)
         atomic2(i) = atomic(i)
         mass2(i) = mass(i)
      end do
c
c     get atom pairs to be superimposed from command line
c
      option = 0
      start = 0
      stop = 0
      answer = ' '
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         query = .false.
         read (string,*,err=10,end=10)  option
         if (option .eq. 1) then
            call nextarg (string,exist)
            if (exist) then
               answer = string(1:1)
               read (string,*,err=10,end=10)  start
               answer = ' '
               call nextarg (string,exist)
               if (exist) then
                  answer = string(1:1)
                  read (string,*,err=10,end=10)  stop
                  answer = ' '
               end if
            end if
         end if
      end if
   10 continue
c
c     ask the user which pairs of atoms are to be superimposed
c
      if (query) then
         write (iout,20)
   20    format (/,' Two Options are Available :  (1) Fit atoms',
     &              ' "M" through "N" from structure 1',
     &           /,' to the corresponding atoms of structure 2.',
     &              ' Enter "1,M,N" to use this option.',
     &           /,' If "N" is omitted, the fit uses atoms 1',
     &              ' through "M". If both "M" and "N" are',
     &           /,' omitted, the fit uses all atoms; or (2)',
     &              ' Individual entry of atom range pairs',
     &           /,' to be used in the fitting procedure.')
         write (iout,30)
   30    format (/,' Enter an Option (either 1,M,N or 2',
     &              ' [<CR>=1,0,0]) :  ',$)
         read (input,40)  record
   40    format (a120)
         read (record,*,err=50,end=50)  option,start,stop
   50    continue
         if (option.lt.1 .or. option.gt.2) then
            option = 1
            start = 0
            stop = 0
         end if
      end if
c
c     warning if structures have different numbers of atoms
c
      if (option .eq. 1) then
         if (n1.ne.n2 .and. start.eq.0) then
            write (iout,60)
   60       format (/,' SUPERPOSE  --  The Molecules contain',
     &                 ' Different Numbers of Atoms')
         end if
      end if
c
c     setup automatic superposition with option to omit hydrogens
c
      if (option .eq. 1) then
         if (answer .eq. ' ') then
            call nextarg (answer,exist)
         else
            exist = .true.
         end if
         if (.not. exist) then
            write (iout,70)
   70       format (/,' Include Hydrogen Atoms in the Fitting',
     &                 ' [Y] :  ',$)
            read (input,80)  record
   80       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (start.eq.0 .and. stop.eq.0) then
            start = 1
            stop = min(n1,n2)
         else if (start.ne.0 .and. stop.eq.0) then
            stop = min(n1,n2,start)
            start = 1
         else if (start.ne.0 .and. stop.ne.0) then
            start = max(1,start)
            stop = min(n1,n2,stop)
         end if
         nfit = 0
         do i = start, stop
            skip = .false.
            if (answer .eq. 'N') then
               if (atomic1(i).le.1 .or. atomic2(i).le.1) then
                  skip = .true.
               end if
            end if
            if (.not. skip) then
               nfit = nfit + 1
               ifit(1,nfit) = i
               ifit(2,nfit) = i
            end if
         end do
      end if
c
c     manual input of the pairs of atom ranges to superimpose
c
      if (option .eq. 2) then
         write (iout,90)
   90    format (/,' On successive lines below, enter atom',
     &              ' pairs or pairs of atom ranges to use',
     &           /,' during fitting. Entering "4,7" will fit',
     &              ' atom 4 of structure 1 to atom 7 of',
     &           /,' structure 2, while the entry "4,7,9,12"',
     &              ' will match atoms 4 through 7 from',
     &           /,' structure 1 with atoms 9 through 12 of',
     &              ' structure 2. Hit <RET> to end entry',
     &           /,' of the list of pairs.')
         nfit = 0
         dowhile (.true.)
            do i = 1, 4
               range(i) = 0
            end do
            write (iout,100)
  100       format (/,' Enter a Pair of Atoms or Ranges :  ',$)
            read (input,110)  record
  110       format (a120)
            read (record,*,err=120,end=120)  (range(i),i=1,4)
  120       continue
            if (range(1) .eq. 0) then
               goto 130
            else if (range(2) .eq. 0) then
               nfit = nfit + 1
               ifit(1,nfit) = range(1)
               ifit(2,nfit) = range(1)
            else if (range(3) .eq. 0) then
               nfit = nfit + 1
               ifit(1,nfit) = range(1)
               ifit(2,nfit) = range(2)
            else
               delta = range(3) - range(1)
               do i = range(1), range(2)
                  nfit = nfit + 1
                  ifit(1,nfit) = i
                  ifit(2,nfit) = i + delta
               end do
            end if
         end do
  130    continue
      end if
c
c     decide on the weighting to use for the coordinates
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,140)
  140    format (/,' Use Mass- or Unit-Weighted Coordinates',
     &              ' (M or [U]) :  ',$)
         read (input,150)  record
  150    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'M') then
         do i = 1, nfit
            wfit(i) = 0.5d0 * (mass1(ifit(1,i)) + mass2(ifit(2,i)))
         end do
      else
         do i = 1, nfit
            wfit(i) = 1.0d0
         end do
      end if
c
c     decide whether to write the best fit set of coordinates
c
      call nextarg (answer,exist)
      if (.not. exist) then
         write (iout,160)
  160    format (/,' Write Best-Fit Coordinates of 2nd Molecule',
     &              ' [N] :  ',$)
         read (input,170)  record
  170    format (a120)
         next = 1
         call gettext (record,answer,next)
      end if
      call upcase (answer)
c
c     chose cutoff value for output of atom pair deviations
c
      cutoff = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=180,end=180)  cutoff
  180 continue
      if (cutoff .lt. 0.0d0) then
         write (iout,190)
  190    format (/,' Cutoff Value for Listing RMS Deviations',
     &              ' [0.0] :  ',$)
         read (input,200)  cutoff
  200    format (f20.0)
      end if
c
c     information about structures to be superimposed
c
      write (iout,210)  file1(1:leng1)
  210 format (/,' Structure File 1 :  ',a)
      write (iout,220)  file2(1:leng2)
  220 format (/,' Structure File 2 :  ',a)
c
c     reopen the coordinate files with structures to superimpose
c
      ifile1 = freeunit ()
      call suffix (file1,'xyz')
      call version (file1,'old')
      open (unit=ifile1,file=file1,status ='old')
      rewind (unit=ifile1)
      call suffix (file2,'xyz')
      call version (file2,'old')
      if (file1 .eq. file2) then
         same = .true.
         ifile2 = ifile1
      else
         same = .false.
         ifile2 = freeunit ()
         open (unit=ifile2,file=file2,status ='old')
         rewind (unit=ifile2)
      end if
c
c     read initial structure set from the first coordinate file
c
      frame1 = 1
      call readxyz (ifile1)
      n1 = n
      do i = 1, n1
         x1(i) = x(i)
         y1(i) = y(i)
         z1(i) = z(i)
      end do
c
c     read initial structure set from the second coordinate file
c
      frame2 = 1
      if (same)  frame2 = 2
      call readxyz (ifile2)
      n2 = n
      do i = 1, n2
         x2(i) = x(i)
         y2(i) = y(i)
         z2(i) = z(i)
      end do
      if (abort) then
         abort = .false.
         frame2 = 1
         n2 = n1
         do i = 1, n2
            x2(i) = x1(i)
            y2(i) = y1(i)
            z2(i) = z1(i)
         end do
      end if
c
c     perform the superposition of a structure pair
c
      dowhile (.not. abort)
         write (iout,230)  frame1,frame2
  230    format (/,' File 1 Frame :',i6,13x,'File 2 Frame :',i6)
         write (iout,240)
  240    format (/,' Summary of Results from Structural',
     &              ' Superposition :')
         verbose = .true.
         call impose (n1,x1,y1,z1,n2,x2,y2,z2,rmsvalue)
         write (iout,250)  rmsvalue,frame1,frame2
  250    format (/,' Root Mean Square Distance :',11x,f15.6,2x,2i7)
c
c     write out the results of the superposition
c
         header = .true.
         do i = 1, nfit
            i1 = ifit(1,i)
            i2 = ifit(2,i)
            dist = sqrt((x1(i1)-x2(i2))**2 + (y1(i1)-y2(i2))**2
     &                         + (z1(i1)-z2(i2))**2)
            if (dist .ge. cutoff) then
               if (header) then
                  header = .false.
                  write (iout,260)
  260             format (/,'   Atom in the',9x,'Atom in the',12x,
     &                       'Distance',10x,'Weight'
     &                    /,' First Structure',5x,'Second Structure',
     &                       8x,'Separated',10x,'in Fit'/)
               end if
               write (iout,270)  i1,name1(i1),i2,name2(i2),dist,wfit(i)
  270          format (5x,i5,'-',a3,11x,i5,'-',a3,7x,f13.6,4x,f12.4)
            end if
         end do
         if (.not. header) then
            write (iout,280)  rmsvalue
  280       format (/,' Root Mean Square Distance :',11x,f15.6)
         end if
c
c     create output file for superimposed second structure
c
         if (answer .eq. 'Y') then
            do i = 1, n
               x(i) = x2(i)
               y(i) = y2(i)
               z(i) = z2(i)
            end do
            ixyz = freeunit ()
            xyzfile = file2(1:leng)//'.xyz'
            call version (xyzfile,'new')
            open (unit=ixyz,file=xyzfile,status='new')
            call prtxyz (ixyz)
            close (unit=ixyz)
         end if
c
c     attempt to get next structure pair from coordinate files
c
         frame2 = frame2 + 1
         call readxyz (ifile2)
         n2 = n
         do i = 1, n2
            x2(i) = x(i)
            y2(i) = y(i)
            z2(i) = z(i)
         end do
         if (abort) then
            abort = .false.
            if (same) then
               rewind (unit=ifile1)
               do i = 1, frame1
                  call readxyz (ifile1)
               end do
            end if
            frame1 = frame1 + 1
            call readxyz (ifile1)
            n1 = n
            do i = 1, n1
               x1(i) = x(i)
               y1(i) = y(i)
               z1(i) = z(i)
            end do
            if (.not. abort) then
               frame2 = frame1 + 1
               if (.not. same) then
                  frame2 = 1
                  rewind (unit=ifile2)
               end if
               call readxyz (ifile2)
               n2 = n
               do i = 1, n2
                  x2(i) = x(i)
                  y2(i) = y(i)
                  z2(i) = z(i)
               end do
            end if
         end if
      end do
c
c     perform any final tasks before program exit
c
      close (unit=ifile1)
      if (.not. same)  close (unit=ifile2)
      call final
      end
