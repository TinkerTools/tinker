c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine getarc  --  get coordinate archive or trajectory  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "getarc" asks for a coordinate archive or trajectory file name,
c     then reads in the initial set of coordinates
c
c
      subroutine getarc (iarc)
      use inform
      use iounit
      use output
      implicit none
      integer iarc,nask
      integer freeunit
      logical exist
      character*240 arcfile
c
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
c     ask for the user specified input structure filename
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
c     first open file, then read the initial set of coordinates
c
      coordtype = 'CARTESIAN'
      iarc = freeunit ()
      open (unit=iarc,file=arcfile,status='old')
      rewind (unit=iarc)
      call readxyz (iarc)
c
c     quit if the coordinates archive file contains no atoms
c
      if (abort) then
         write (iout,30)
   30    format (/,' GETARC  --  Coordinates Archive File',
     &              ' does not Contain Any Atoms')
         close (unit=iarc)
         call fatal
      end if
      return
      end
