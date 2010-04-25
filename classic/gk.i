c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  gk.i  --  parameters for generalized Kirkwood solvation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     gkr      generalized Kirkwood cavity radii for atom types
c     gkc      tuning parameter exponent in the f(GB) function
c
c
      real*8 gkr,gkc
      common /gk/ gkr(maxtyp),gkc
