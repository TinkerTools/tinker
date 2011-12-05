c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  pdb.i  --  definition of a Protein Data Bank structure  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     xpdb      x-coordinate of each atom stored in PDB format
c     ypdb      y-coordinate of each atom stored in PDB format
c     zpdb      z-coordinate of each atom stored in PDB format
c     npdb      number of atoms stored in Protein Data Bank format
c     nres      number of residues stored in Protein Data Bank format
c     resnum    number of the residue to which each atom belongs
c     resatm    number of first and last atom in each residue
c     npdb12    number of atoms directly bonded to each CONECT atom
c     ipdb12    atom numbers of atoms connected to each CONECT atom
c     pdblist   list of the Protein Data Bank atom number of each atom
c     pdbtyp    Protein Data Bank record type assigned to each atom
c     atmnam    Protein Data Bank atom name assigned to each atom
c     resnam    Protein Data Bank residue name assigned to each atom
c     chnsym    string with PDB chain identifiers to be included
c     altsym    string with PDB alternate locations to be included
c     instyp    string with PDB insertion records to be included
c
c
      integer npdb,nres
      integer resnum,resatm
      integer npdb12,ipdb12
      integer pdblist
      real*8 xpdb,ypdb,zpdb
      character*1 altsym
      character*3 resnam
      character*4 atmnam
      character*6 pdbtyp
      character*20 chnsym
      character*20 instyp
      common /pdb/ xpdb(maxatm),ypdb(maxatm),zpdb(maxatm),npdb,nres,
     &             resnum(maxatm),resatm(2,maxatm),npdb12(maxatm),
     &             ipdb12(maxval,maxatm),pdblist(maxatm),pdbtyp(maxatm),
     &             atmnam(maxatm),resnam(maxatm),chnsym,altsym,instyp
