c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  shake.i  --  definition of Shake/Rattle constraints  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     rateps       convergence tolerance for rattle constraints
c     krat         ideal distance value for rattle constraint
c     nrat         number of rattle distance constraints to apply
c     nratx        number of atom group spatial constraints to apply
c     irat         atom numbers of atoms in a rattle constraint
c     iratx        group number of group in a spatial constraint
c     kratx        spatial constraint type (1=plane, 2=line, 3=point)
c     ratimage     flag to use minimum image for rattle constraint
c     use_rattle   logical flag to set use of rattle contraints
c
c
      integer nrat,nratx
      integer irat,iratx
      integer kratx
      real*8 rateps,krat
      logical ratimage,use_rattle
      common /shake/ rateps,krat(maxfix),nrat,nratx,irat(2,maxfix),
     $               iratx(maxfix),kratx(maxfix),ratimage(maxfix),
     &               use_rattle
