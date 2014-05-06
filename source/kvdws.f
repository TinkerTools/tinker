c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module kvdws  --  van der Waals term forcefield parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     rad      van der Waals radius parameter for each atom type
c     eps      van der Waals well depth parameter for each atom type
c     rad4     van der Waals radius parameter in 1-4 interactions
c     eps4     van der Waals well depth parameter in 1-4 interactions
c     reduct   van der Waals reduction factor for each atom type
c
c
      module kvdws
      use sizes
      implicit none
      real*8 rad(maxtyp)
      real*8 eps(maxtyp)
      real*8 rad4(maxtyp)
      real*8 eps4(maxtyp)
      real*8 reduct(maxtyp)
      save
      end
