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
c     maxtcg    maximum order for TCG induced dipole iterations
c
c     tcgorder  total number of TCG iterations to be used
c     tcgnab
c     tcgprec   flag to allow conjugate gradient preconditioner
c     tcgpeek   flag to allow use of a final TCG peek step
c     tcgomega  value of acceleration factor for TCG peek step
c     uindt
c     uinpt
c     uad
c     uap
c     ubd
c     ubp
c
c
      module poltcg
      implicit none
      integer maxtcg
      parameter (maxtcg=2)
      integer tcgorder
      integer tcgnab
      logical tcgprec
      logical tcgpeek
      real*8 tcgomega
      real*8, allocatable :: uindt(:,:)
      real*8, allocatable :: uinpt(:,:)
      real*8, allocatable :: uad(:,:,:)
      real*8, allocatable :: uap(:,:,:)
      real*8, allocatable :: ubd(:,:,:)
      real*8, allocatable :: ubp(:,:,:)
      save
      end
