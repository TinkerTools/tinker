c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2018 by Zhi Wang and Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module poltcg  --  induced dipoles for the TCG method  ##
c     ##                                                         ##
c     #############################################################
c
c
c     tcgorder   total number of TCG iterations to be used
c     tcgnab     number of mutual induced dipole components
c     tcgprec    flag to allow use of a diagonal preconditioner
c     tcgpeek    flag to allow use of a final TCG peek step
c     tcgguess   flag to use initial TCG based on direct field
c     tcgomega   value of acceleration factor for TCG peek step
c     uindt      induced d-dipole components for the TCG method
c     uinpt      induced p-dipole components for the TCG method
c     uad        left-hand side mutual induced d-dipoles
c     uap        left-hand side mutual induced p-dipoles
c     ubd        right-hand side mutual induced d-dipoles
c     ubp        right-hand side mutual induced p-dipoles
c
c
      module poltcg
      implicit none
      integer tcgorder
      integer tcgnab
      logical tcgprec
      logical tcgpeek
      logical tcgguess
      real*8 tcgomega
      real*8, allocatable :: uindt(:,:)
      real*8, allocatable :: uinpt(:,:)
      real*8, allocatable :: uad(:,:,:)
      real*8, allocatable :: uap(:,:,:)
      real*8, allocatable :: ubd(:,:,:)
      real*8, allocatable :: ubp(:,:,:)
      save
      end
