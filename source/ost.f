c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  module ost  --  orthogonal space tempering variables  ##
c     ##                                                        ##
c     ############################################################
c
c
c     iost       steps between lambda recursion kernel calculations
c     ostlambda  main lambda value in orthogonal space sampling
c     ostplmda1  sublambda upper bound for polarization
c     ostplmda0  sublambda lower bound for polarization
c     ostelmda1  sublambda upper bound for electrostatics
c     ostelmda0  sublambda lower bound for electrostatics
c     ostvlmda1  sublambda upper bound for van der Waals
c     ostvlmda0  sublambda lower bound for van der Waals
c     use_ost    flag to use orthogonal space tempering
c     use_pol4i  flag to compute polarization lambda deriv for lmda=0
c     use_pol4f  flag to compute polarization lambda deriv for lmda=1
c
c
      module ost
      implicit none
      integer iost
      real*8 ostlambda
      real*8 ostplmda1
      real*8 ostplmda0
      real*8 ostelmda1
      real*8 ostelmda0
      real*8 ostvlmda1
      real*8 ostvlmda0
      logical use_ost
      logical use_pol4i
      logical use_pol4f
      save
      end
