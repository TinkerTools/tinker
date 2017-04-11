c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2014  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program potplot  --  convert Epot grid to xyz for display  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "potplot" reads a TINKER electrostatic potential grid file
c     and generates an atomic coordinates file for display in FFE
c
c
      program potplot
      use iounit
      implicit none
      integer i,k
      integer ipot,npot
      integer freeunit
      real*8 xpot,ypot,zpot
      logical exist
      character*240 record
      character*240 potfile
c
c
c     find and open the electrostatic potential grid file
c
      call initial
      call nextarg (potfile,exist)
      if (exist) then
         call basefile (potfile)
         call suffix (potfile,'pot','old')
         inquire (file=potfile,exist=exist)
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Electrostatic Potential File Name :  ',$)
         read (input,20)  potfile
   20    format (a240)
         call basefile (potfile)
         call suffix (potfile,'pot','old')
         inquire (file=potfile,exist=exist)
      end do
      ipot = freeunit ()
      open (unit=ipot,file=potfile,status='old')
      rewind (unit=ipot)
c
c     get the grid points and output them as atom coordinates
c
      npot = 0
      read (ipot,30,err=50,end=50)  record
   30 format (a240)
      read (record,*,err=50,end=50)  npot
      write (iout,40)  npot
   40 format (i6,2x,'Epot Grid converted to Atoms for Display')
   50 continue
      do i = 1, npot
         read (ipot,60,err=80,end=80)  record
   60    format (a240)
         read (record,*,err=80,end=80)  k,xpot,ypot,zpot
         write (iout,70)  k,xpot,ypot,zpot
   70    format (i6,2x,'He ',3f12.4,4x,'20')
   80    continue
      end do
      close (unit=ipot)
c
c     perform any final tasks before program exit
c
      call final
      end
