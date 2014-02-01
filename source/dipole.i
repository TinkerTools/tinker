c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  dipole.i  --  atom & bond dipoles for current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     bdpl      magnitude of each of the dipoles (Debyes)
c     sdpl      position of each dipole between defining atoms
c     ndipole   total number of dipoles in the system
c     idpl      numbers of atoms that define each dipole
c
c
      integer ndipole,idpl
      real*8 bdpl,sdpl
      common /dipole/ bdpl(maxbnd),sdpl(maxbnd),ndipole,idpl(2,maxbnd)
