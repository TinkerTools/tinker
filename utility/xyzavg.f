c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program xyzavg  --  average structure from an archive file  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "xyzavg" compute the coordiantes of the average strurcture from
c     an archive file; assumes the files in the archive have already
c     been superimposed
c
c
      program xyzavg
      use sizes
      use atoms
      use inform
      use iounit
      implicit none
      integer i,ixyz
      integer natom
      integer nframe
      integer freeunit
      real*8 denom
      real*8, allocatable :: xave(:)
      real*8, allocatable :: yave(:)
      real*8, allocatable :: zave(:)
      logical exist
      character*240 xyzfile
c
c
c     try to get a filename from the command line arguments
c
      call initial
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Coordinate Archive File Name :  ',$)
         read (input,20)  xyzfile
   20    format (a240)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do
c
c     open the input file and read the first coordinate set
c
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      call readxyz (ixyz)
      allocate (xave(n))
      allocate (yave(n))
      allocate (zave(n))
      natom = n
      nframe = 1
      do i = 1, n
         xave(i) = x(i)
         yave(i) = y(i)
         zave(i) = z(i)
      end do
c
c     get the remaining coordinate sets from the archive file
c
      dowhile (.true.)
         call readxyz (ixyz)
         if (abort)  goto 30
         natom = n
         nframe = nframe + 1
         do i = 1, n
            xave(i) = xave(i) + x(i)
            yave(i) = yave(i) + y(i)
            zave(i) = zave(i) + z(i)
         end do
      end do
   30 continue
      close (unit=ixyz)
c
c     compute the average coordinates over all structures
c
      n = natom
      denom = dble(nframe)
      do i = 1, n
         x(i) = xave(i) / denom
         y(i) = yave(i) / denom
         z(i) = zave(i) / denom
      end do
c
c     write the average coordinates to a file
c
      call suffix (xyzfile,'xyz','new')
      open (unit=ixyz,file=xyzfile,status='new')
      call prtxyz (ixyz)
      close (unit=ixyz)
      call final
      end
