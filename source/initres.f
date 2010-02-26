c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine initres  --  setup biopolymer residue names  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "initres" sets names for biopolymer residue types used in
c     PDB file conversion and automated generation of structures
c
c
      subroutine initres
      implicit none
      include 'sizes.i'
      include 'resdue.i'
      integer i
      character*1 acid1(maxamino)
      character*1 base1(maxnuc)
      character*3 acid3(maxamino)
      character*3 base3(maxnuc)
      data acid1  / 'G','A','V','L','I','S','T','C','C','P','F','Y',
     &              'W','H','U','Z','D','N','E','Q','M','K','R','O',
     &              'B','J','f','a','n','m','X' /
      data acid3  / 'GLY','ALA','VAL','LEU','ILE','SER','THR','CYS',
     &              'CYX','PRO','PHE','TYR','TRP','HIS','HID','HIE',
     &              'ASP','ASN','GLU','GLN','MET','LYS','ARG','ORN',
     &              'AIB','PCA','FOR','ACE','NH2','NME','UNK' /
      data base1  / 'A','G','C','U','D','B','I','T','1','2','3','X' /
      data base3  / 'A  ','G  ','C  ','U  ','DA ','DG ','DC ','DT ',
     &              'MP ','DP ','TP ','UNK' /
c
c
c     set values for the 1- and 3-letter amino acid names
c
      do i = 1, maxamino
         amino(i) = acid3(i)
         amino1(i) = acid1(i)
      end do
c
c     set values for the 1- and 3-letter nucleic acid names
c
      do i = 1, maxnuc
         nuclz(i) = base3(i)
         nuclz1(i) = base1(i)
      end do
      return
      end
