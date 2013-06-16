c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mutant.i  --  hybrid atoms for free energy perturbation  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     lambda     generic weighting between initial and final states
c     vlambda    state weighting value for electrostatic potentials
c     elambda    state weighting value for van der Waals potentials
c     scexp
c     scalpha
c     eslvt
c     eslut
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     mut        true if an atom is to be mutated, false otherwise
c     eupdated   flag to mark updated energy value in BAR method
c
c
      integer nmut,imut
      integer type0,class0
      integer type1,class1
      real*8 lambda
      real*8 vlambda,elambda
      real*8 scexp,scalpha
      real*8 eslvt,eslut
      logical mut,eupdated
      common /mutant/ lambda,vlambda,elambda,scexp,scalpha,eslvt,eslut,
     &                nmut,imut(maxatm),type0(maxatm),class0(maxatm),
     &                type1(maxatm),class1(maxatm),mut(maxatm),eupdated
