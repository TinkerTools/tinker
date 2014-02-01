c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  molcul.i  --  individual molecules within current system  ##
c     ##                                                            ##
c     ################################################################
c
c
c     molmass   molecular weight for each molecule in the system
c     totmass   total weight of all the molecules in the system
c     nmol      total number of separate molecules in the system
c     kmol      contiguous list of the atoms in each molecule
c     imol      first and last atom of each molecule in the list
c     molcule   number of the molecule to which each atom belongs
c
c
      integer nmol,kmol,imol,molcule
      real*8 molmass,totmass
      common /molcul/ molmass(maxatm),totmass,nmol,kmol(maxatm),
     &                imol(2,maxatm),molcule(maxatm)
