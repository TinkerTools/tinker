c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getmol2  --  get a Tripos MOL2 format file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getmol2" asks for a Tripos MOL2 molecule file name,
c     then reads the coordinates from the file
c
c
      subroutine getmol2
      use files
      use inform
      use iounit
      implicit none
      integer imol2,nask
      integer freeunit
      logical exist
      character*240 mol2file
c
c
c     try to get a filename from the command line arguments
c
      call nextarg (mol2file,exist)
      if (exist) then
         call basefile (mol2file)
         call suffix (mol2file,'mol2','old')
         inquire (file=mol2file,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter a Tripos MOL2 File Name :  ',$)
         read (input,20)  mol2file
   20    format (a240)
         call basefile (mol2file)
         call suffix (mol2file,'mol2','old')
         inquire (file=mol2file,exist=exist)
      end do
      if (.not. exist)  call fatal
c
c     first open and then read the Tripos MOL2 coordinates file
c
      filename = mol2file
      imol2 = freeunit ()
      open (unit=imol2,file=mol2file,status='old')
      rewind (unit=imol2)
      call readmol2 (imol2)
      close (unit=imol2)
      return
      end
