c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program xyzpdb  --  Cartesian to Protein Data Bank file  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xyzpdb" takes as input a Cartesian coordinates file,
c     then converts to and writes out a Protein Data Bank file
c
c
      program xyzpdb
      implicit none
      include 'files.i'
      include 'inform.i'
      integer ipdb,ixyz
      integer freeunit
      character*120 pdbfile
      character*120 xyzfile
c
c
c     get the Cartesian coordinates file for the system
c
      call initial
      call getxyz
c
c     get atomic number of each atom and count the molecules
c
      call field
      call katom
      call molecule
c
c     open the Protein Data Bank file to be used for output
c
      ipdb = freeunit ()
      pdbfile = filename(1:leng)//'.pdb'
      call version (pdbfile,'new')
      open (unit=ipdb,file=pdbfile,status='new')
c
c     reopen the coordinates file and read the first structure
c
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz')
      call version (xyzfile,'old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     add each successive coordinate frame to the PDB file
c
      dowhile (.not. abort)
         call makepdb
         call prtpdb (ipdb)
         call readxyz (ixyz)
      end do
c
c     perform any final tasks before program exit
c
      close (unit=ipdb)
      call final
      end
