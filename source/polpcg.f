c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2013  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module poltcg  --  induced dipoles via the PCG solver  ##
c     ##                                                         ##
c     #############################################################
c
c
c     mindex    index into preconditioner inverse for PCG solver
c     pcgpeek   value of acceleration factor for PCG peek step
c     minv      preconditioner inverse for induced dipole PCG solver
c     pcgprec   flag to allow use of preconditioner with PCG solver
c     pcgguess  flag to use initial PCG based on direct field
c
c
      module polpcg
      implicit none
      integer, allocatable :: mindex(:)
      real*8 pcgpeek
      real*8, allocatable :: minv(:)
      logical pcgprec
      logical pcgguess
      save
      end
