c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  freeze.i  --  definition of holonomic RATTLE constraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     rateps       convergence tolerance for holonomic constraints
c     krat         ideal distance value for holonomic constraint
c     nrat         number of holonomic distance constraints to apply
c     nratx        number of atom group holonomic constraints to apply
c     irat         atom numbers of atoms in a holonomic constraint
c     iratx        group number of group in a holonomic constraint
c     kratx        spatial constraint type (1=plane, 2=line, 3=point)
c     ratimage     flag to use minimum image for holonomic constraint
c     use_rattle   logical flag to set use of holonomic contraints
c
c
      integer nrat,nratx
      integer irat,iratx
      integer kratx
      real*8 rateps,krat
      logical ratimage,use_rattle
      common /freeze/ rateps,krat(maxfix),nrat,nratx,irat(2,maxfix),
     $                iratx(maxfix),kratx(maxfix),ratimage(maxfix),
     &                use_rattle
