c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine getcart  --  get a Cartesian coordinates file  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "getcart" asks for a Cartesian coordinate file name, then
c     reads the formatted or binary coordinates file
c
c
      subroutine getcart (ixyz)
      use files
      use inform
      use iounit
      use output
      implicit none
      integer ixyz,iaux
      integer nask
      integer freeunit
      logical exist,first
      character*1 letter
      character*240 xyzfile
      character*240 auxfile
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter Cartesian Coordinate File Name :  ',$)
         read (input,20)  xyzfile
   20    format (a240)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     open the file and get format by inspecting first character
c
      filename = xyzfile
      coordtype = 'CARTESIAN'
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      read (ixyz,30)  letter
   30 format (a1)
      archive = .false.
      if (letter .eq. ' ')  archive = .true.
      if (letter.ge.'0' .and. letter.le.'9')  archive = .true.
      binary = (.not. archive)
c
c     read initial Cartesian coordinates from formatted file
c
      if (archive) then
         rewind (unit=ixyz)
         call readxyz (ixyz)
      end if
c
c     get atom types and connectivity from formatted file
c
      if (binary) then
         call nextarg (auxfile,exist)
         if (exist) then
            call basefile (auxfile)
            call suffix (auxfile,'xyz','old')
            inquire (file=auxfile,exist=exist)
         end if
         nask = 0
         do while (.not.exist .and. nask.lt.maxask)
            nask = nask + 1
            write (iout,40)
   40       format (/,' Enter Formatted Coordinate File Name :  ',$)
            read (input,50)  auxfile
   50       format (a240)
            call basefile (auxfile)
            call suffix (auxfile,'xyz','old')
            inquire (file=auxfile,exist=exist)
         end do
         if (.not. exist)  call fatal
         iaux = freeunit ()
         open (unit=iaux,file=auxfile,status='old')
         rewind (unit=iaux)
         call readxyz (iaux)
         close (unit=iaux)
c
c     read initial Cartesian coordinates from binary file
c
         filename = xyzfile
         close (unit=ixyz)
         ixyz = freeunit ()
         open (unit=ixyz,file=xyzfile,form='unformatted',status='old')
         rewind (unit=ixyz)
         first = .true.
         call readdcd (ixyz,first)
      end if
c
c     quit if the Cartesian coordinates file contains no atoms
c
      if (abort) then
         write (iout,60)
   60    format (/,' GETCART  --  Cartesian Coordinate File',
     &              ' was not Read Correctly')
         close (unit=ixyz)
         call fatal
      end if
      return
      end
