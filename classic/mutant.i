c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  mutant.i  --  parameters for free energy perturbation  ##
c     ##                                                         ##
c     #############################################################
c
c
c     lambda    generic weighting of the initial and final states
c     vlambda   weighting of initial and final states for vdw
c     clambda   weighting of initial and final states for charges
c     dlambda   weighting of initial and final states for dipoles
c     mlambda   weighting of initial and final states for multipoles
c     plambda   weighting of initial and final states for polarization
c     nmut      number of atoms mutated from initial to final state
c     imut      atomic sites differing in initial and final state
c     type0     atom type of each atom in the initial state system
c     type1     atom type of each atom in the final state system
c     class0    atom class of each atom in the initial state system
c     class1    atom class of each atom in the final state system
c     mut       true if an atom is to be mutated, false otherwise
c
c
      integer nmut,imut(maxatm)
      integer type0(maxatm),type1(maxatm)
      integer class0(maxatm),class1(maxatm)
      real*8 lambda,vlambda,clambda
      real*8 dlambda,mlambda,plambda
      logical mut(maxatm)
      common /mutant/ lambda,vlambda,clambda,dlambda,mlambda,plambda,
     &                nmut,imut,type0,type1,class0,class1,mut
