c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  align.i  --  information for superposition of structures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     wfit    weights assigned to atom pairs during superposition
c     nfit    number of atoms to use in superimposing two structures
c     ifit    atom numbers of pairs of atoms to be superimposed
c
c
      integer nfit,ifit
      real*8 wfit
      common /align/ wfit(maxatm),nfit,ifit(2,maxatm)
