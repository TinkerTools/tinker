c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2020  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module solpot  --  solvation term functional form details  ##
c     ##                                                             ##
c     #################################################################
c
c
c     solvtyp   type of continuum solvation energy model in use
c     borntyp   method to be used for the Born radius computation
c
c
      module solpot
      implicit none
      character*8 solvtyp
      character*8 borntyp
      save
      end
