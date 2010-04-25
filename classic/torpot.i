c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  torpot.i  --  specifics of torsional functional forms  ##
c     ##                                                         ##
c     #############################################################
c
c
c     idihunit  convert improper dihedral energy to kcal/mole
c     itorunit  convert improper torsion amplitudes to kcal/mole
c     torsunit  convert torsional parameter amplitudes to kcal/mole
c     ptorunit  convert pi-orbital torsion energy to kcal/mole
c     storunit  convert stretch-torsion energy to kcal/mole
c     ttorunit  convert stretch-torsion energy to kcal/mole
c
c
      real*8 idihunit,itorunit,torsunit
      real*8 ptorunit,storunit,ttorunit
      common /torpot/ idihunit,itorunit,torsunit,ptorunit,storunit,
     &                ttorunit
