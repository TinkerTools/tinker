c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  imptor.i  --  improper torsions in the current structure  ##
c     ##                                                            ##
c     ################################################################
c
c
c     itors1   1-fold amplitude and phase for each improper torsion
c     itors2   2-fold amplitude and phase for each improper torsion
c     itors3   3-fold amplitude and phase for each improper torsion
c     nitors   total number of improper torsional angles in the system
c     iitors   numbers of the atoms in each improper torsional angle
c
c
      integer nitors,iitors
      real*8 itors1,itors2,itors3
      common /imptor/ itors1(4,maxtors),itors2(4,maxtors),
     &                itors3(4,maxtors),nitors,iitors(4,maxtors)
