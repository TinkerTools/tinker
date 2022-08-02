c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################        
c     ##                                                          ##
c     ##  subroutine getdcd  --  get DCD coordinate archive file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getdcd" asks for a binary DCD trajectory file name and the
c     corresponding Tinker coordinates file, then reads the initial
c     set of DCD coordinates
c
c
      subroutine getdcd (idcd)
      use inform
      use iounit
      use output
      implicit none
      integer idcd,ixyz,nask
      integer freeunit
      logical exist
      character*240 dcdfile
      character*240 xyzfile
c
c
c     try to get a DCD filename from the command line arguments
c
      call nextarg (dcdfile,exist)
      if (exist) then
         call basefile (dcdfile)
         call suffix (dcdfile,'dcd','old')
         inquire (file=dcdfile,exist=exist)
      end if
c
c     ask for the user specified input DCD trajectory filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter DCD Coordinate Archive File Name :  ',$)
         read (input,20)  dcdfile
   20    format (a240)
         call basefile (dcdfile)
         call suffix (dcdfile,'dcd','old')
         inquire (file=dcdfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     try to get a coordinates file from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input coordinates filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,30)
   30    format (/,' Enter Formatted Coordinate File Name :  ',$)
         read (input,40)  xyzfile
   40    format (a240)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     first open and then read the Cartesian coordinates file
c
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      close (unit=ixyz)
c
c     next read the initial set of coordinates from the DCD file
c
      idcd = freeunit ()
      open (unit=idcd,file=dcdfile,form='unformatted',status='old')
      rewind (unit=idcd)
      call readdcd (idcd)
c
c     quit if the DCD trajectory archive file contains no atoms
c
      if (abort) then
         write (iout,50)
   50    format (/,' GETDCD  --  Binary DCD Trajectory File',
     &              ' was not Read Correctly')
         close (unit=idcd)
         call fatal
      end if
      return
      end
