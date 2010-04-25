c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1998  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  khbond.i  --  forcefield parameters for H-bonding terms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnhb   maximum number of hydrogen bonding pair entries
c
c     radhb    radius parameter for hydrogen bonding pairs
c     epshb    well depth parameter for hydrogen bonding pairs
c     khb      string of atom types for hydrogen bonding pairs
c
c
      integer maxnhb
      parameter (maxnhb=500)
      real*8 radhb,epshb
      character*8 khb
      common /khbond/ radhb(maxnhb),epshb(maxnhb),khb(maxnhb)
