c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  orbits.i  --  orbital energies for conjugated pisystem  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     q       number of pi-electrons contributed by each atom
c     w       ionization potential of each pisystem atom
c     em      repulsion integral for each pisystem atom
c     nfill   number of filled pisystem molecular orbitals
c
c
      integer nfill
      real*8 q,w,em
      common /orbits/ q(maxpi),w(maxpi),em(maxpi),nfill
