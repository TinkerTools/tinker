c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2018 by Zhi Wang and Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module poltcg  --  induced dipoles via the TCG solver  ##
c     ##                                                         ##
c     #############################################################
c
c
c     tcgorder   total number of TCG iterations to be used
c     tcgnab     number of mutual induced dipole components
c     tcgpeek    value of acceleration factor for TCG peek step
c     uindt      induced d-dipole components for the TCG method
c     uinpt      induced p-dipole components for the TCG method
c     uad        left-hand side mutual induced d-dipoles
c     uap        left-hand side mutual induced p-dipoles
c     ubd        right-hand side mutual induced d-dipoles
c     ubp        right-hand side mutual induced p-dipoles
c     tcgprec    flag to allow use of a diagonal preconditioner
c     tcgguess   flag to use initial TCG based on direct field
c
c
      module poltcg
      implicit none
      integer tcgorder
      integer tcgnab
      real*8 tcgpeek
      real*8, allocatable :: uindt(:,:)
      real*8, allocatable :: uinpt(:,:)
      real*8, allocatable :: uad(:,:,:)
      real*8, allocatable :: uap(:,:,:)
      real*8, allocatable :: ubd(:,:,:)
      real*8, allocatable :: ubp(:,:,:)
      logical tcgprec
      logical tcgguess
      save
      end
