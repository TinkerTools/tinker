c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_bounds  --  periodic wrapping  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "test_bounds" checks the periodic coordinate wrapping case
c     exercised by the tinker-gpu bounds.cpp test
c
c
      subroutine test_bounds
      use atoms
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_bounds')
c
c
      if (skiptest(tname,'amoeba'))  return
      call pushdir ('file/bounds')
      call loadfix ('bounds','bounds.key')
      call bounds
      call assert_real (x(1),-9.0d0,1.0d-5,tname//' x1')
      call assert_real (y(1),-3.0d0,1.0d-5,tname//' y1')
      call assert_real (z(1),4.0d0,1.0d-5,tname//' z1')
      call assert_real (x(2),3.0d0,1.0d-5,tname//' x2')
      call assert_real (y(2),-5.0d0,1.0d-5,tname//' y2')
      call assert_real (z(2),7.0d0,1.0d-5,tname//' z2')
      call popdir
      call final
      return
      end
