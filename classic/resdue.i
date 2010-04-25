c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  resdue.i  --  standard biopolymer residue abbreviations  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     amino    three-letter abbreviations for amino acids types
c     nuclz    three-letter abbreviations for nucleic acids types
c     amino1   one-letter abbreviations for amino acids types
c     nuclz1   one-letter abbreviations for nucleic acids types
c
c
      character*1 amino1,nuclz1
      character*3 amino,nuclz
      common /resdue/ amino(maxamino),nuclz(maxnuc),amino1(maxamino),
     &                nuclz1(maxnuc)
