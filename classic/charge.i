c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  charge.i  --  partial charges for the current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     pchg      magnitude of the partial charges (e-)
c     nion      total number of partial charges in system
c     iion      number of the atom site for each partial charge
c     jion      neighbor generation site for each partial charge
c     kion      cutoff switching site for each partial charge
c     chglist   partial charge site for each atom (0=no charge)
c
c
      integer nion,iion
      integer jion,kion
      integer chglist
      real*8 pchg
      common /charge/ pchg(maxatm),nion,iion(maxatm),jion(maxatm),
     &                kion(maxatm),chglist(maxatm)
