c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  analyz.i  --  energy components partitioned over atoms  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     aesum   total potential energy partitioned over atoms
c     aeb     bond stretch energy partitioned over atoms
c     aea     angle bend energy partitioned over atoms
c     aeba    stretch-bend energy partitioned over atoms
c     aeub    Urey-Bradley energy partitioned over atoms
c     aeaa    angle-angle energy partitioned over atoms
c     aeopb   out-of-plane bend energy partitioned over atoms
c     aeopd   out-of-plane distance energy partitioned over atoms
c     aeid    improper dihedral energy partitioned over atoms
c     aeit    improper torsion energy partitioned over atoms
c     aet     torsional energy partitioned over atoms
c     aept    pi-orbital torsion energy partitioned over atoms
c     aebt    stretch-torsion energy partitioned over atoms
c     aett    torsion-torsion energy partitioned over atoms
c     aev     van der Waals energy partitioned over atoms
c     aec     charge-charge energy partitioned over atoms
c     aecd    charge-dipole energy partitioned over atoms
c     aed     dipole-dipole energy partitioned over atoms
c     aem     multipole energy partitioned over atoms
c     aep     polarization energy partitioned over atoms
c     aer     reaction field energy partitioned over atoms
c     aes     solvation energy partitioned over atoms
c     aelf    metal ligand field energy partitioned over atoms
c     aeg     geometric restraint energy partitioned over atoms
c     aex     extra energy term partitioned over atoms
c
c
      real*8 aesum,aeb,aea,aeba
      real*8 aeub,aeaa,aeopb,aeopd
      real*8 aeid,aeit,aet,aept
      real*8 aebt,aett,aev,aec
      real*8 aecd,aed,aem,aep,aer
      real*8 aes,aelf,aeg,aex
      common /analyz/ aesum(maxatm),aeb(maxatm),aea(maxatm),
     &                aeba(maxatm),aeub(maxatm),aeaa(maxatm),
     &                aeopb(maxatm),aeopd(maxatm),aeid(maxatm),
     &                aeit(maxatm),aet(maxatm),aept(maxatm),
     &                aebt(maxatm),aett(maxatm),aev(maxatm),
     &                aec(maxatm),aecd(maxatm),aed(maxatm),
     &                aem(maxatm),aep(maxatm),aer(maxatm),
     &                aes(maxatm),aelf(maxatm),aeg(maxatm),
     &                aex(maxatm)
