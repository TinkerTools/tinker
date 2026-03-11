c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2023  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module extfld  --  applied external electric field vector  ##
c     ##                                                             ##
c     #################################################################
c
c
c     exfreq       frequency of applied external field in gigahertz
c     exfld        components of applied external electric field
c     texfld       components of time dependent applied electric field
c     use_exfld    flag to include applied external electric field
c     use_exfreq   flag to oscillate applied external field
c
c
      module extfld
      implicit none
      real*8 exfreq
      real*8 exfld(3)
      real*8 texfld(3)
      logical use_exfld
      logical use_exfreq
      save
      end
