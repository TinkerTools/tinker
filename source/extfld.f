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
c     exfld       components of applied external electric field
c     use_exfld   flag to include applied external electric field
c
c
      module extfld
      implicit none
      real*8 exfld(3)
      logical use_exfld
      save
      end
