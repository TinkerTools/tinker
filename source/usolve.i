c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  usolve.i  --  inverse preconditioner for induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     mindex   index into preconditioner inverse for CG solver
c     minv     preconditioner inverse for induced dipole CG solver
c
c
      real*8, pointer :: mindex(:)
      real*8, pointer :: minv(:)
      common /usolve/ mindex,minv
