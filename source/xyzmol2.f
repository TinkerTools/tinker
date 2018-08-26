c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1991  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program xyzmol2  --  Cartesian coordinates to Tripos MOL2  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xyzmol2" takes as input a Cartesian coordinates file,
c     converts to and then writes out a Tripos MOL2 file
c
c
      program xyzmol2
      use files
      use iounit
      use titles
      implicit none
      integer imol2,freeunit
      character*240 mol2file
c
c
c     get and read the Cartesian coordinates file
c
      call initial
      call getxyz
      write (iout,10)  title(1:ltitle)
   10 format (' Title :  ',a)
c
c     get a list of the bonds in the structure
c
      call bonds
c
c     open a new version of the Tripos MOL2 file
c
      imol2 = freeunit ()
      mol2file = filename(1:leng)//'.mol2'
      call version (mol2file,'new')
      open (unit=imol2,file=mol2file,status='new')
c
c     output the coordinates into Tripos MOL2 format
c
      call prtmol2 (imol2)
      close (unit=imol2)
c
c     perform any final tasks before program exit
c
      call final
      end
