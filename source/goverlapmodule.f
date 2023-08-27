c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  module goverlapmodule  --  gaussian overlap in tree  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     g                     Gaussian representing overlap 
c     level                 level (0=root, 1=atoms, 2=2-body, ...)
c     atom                  the atomic index of the last atom of the 
c                              overlap list (i, j, k, ..., atom) ! =
c                              (Parent, atom)
c     parent_index          index in tree list of parent overlap
c     children_startindex   start index in tree array of children
c     children_count        number of children contains
c     volume                volume of overlap
c     dvv1                  derivative wrt volume of first atom
c     gamma1i               sum gammai for this overlap
c     self_volume           self volume accumulator
c     sfp                   switching function derivatives
c     dv1                   derivative wrt position of first atom
c
c
      module goverlapmodule
      use gaussianvcamodule
      implicit none
c
c
c     gaussian overlap values
c
      type GOverlap
      type(GaussianVca) g
      integer level
      integer atom
      integer parent_index
      integer children_startindex
      integer children_count
      real*8 volume
      real*8 dvv1
      real*8 gamma1i
      real*8 self_volume
      real*8 sfp
      real*8 dv1(3)
      contains
      procedure, pass(this) :: print_overlap 
      end type GOverlap
      contains
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine print_overlap  --  print gaussian overlap  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "print_overlap" prints the guassian overlap in the tree
c
c
      subroutine print_overlap(this)
      class(GOverlap) this
c
c
c     print the overlap details
c
      write(*, "('           ', I6, I7, I7, I8, I8, 11F11.6)")
     &        this%level, this%atom, this%parent_index,
     &        this%children_startindex, this%children_count,
     &        this%volume, this%gamma1i, this%g%a,
     &        this%g%v, this%g%c(1), this%g%c(2), this%g%c(3), 
     &        this%dv1(1), this%dv1(2), this%dv1(3), this%sfp
c      write(*, "('           ', I6, F20.16)") this%level, this%volume
      return
      end
      end module goverlapmodule
