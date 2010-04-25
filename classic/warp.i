c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  warp.i  --  parameters for potential surface smoothing  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     m2           second moment of the GDA gaussian for each atom
c     deform       value of the smoothing deformation parameter
c     difft        diffusion coefficient for torsional potential
c     diffv        diffusion coefficient for van der Waals potential
c     diffc        diffusion coefficient for charge-charge potential
c     use_smooth   flag to use a potential energy smoothing method
c     use_dem      flag to use diffusion equation method potential
c     use_gda      flag to use gaussian density annealing potential
c     use_tophat   flag to use analytical tophat smoothed potential
c     use_stophat  flag to use shifted tophat smoothed potential
c
c
      real*8 m2,deform
      real*8 difft,diffv,diffc
      logical use_smooth,use_dem,use_gda
      logical use_tophat,use_stophat
      common /warp/ m2(maxatm),deform,difft,diffv,diffc,use_smooth,
     &              use_dem,use_gda,use_tophat,use_stophat
