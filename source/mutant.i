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
c     lambda     weighting of initial state in hybrid Hamiltonian
c     nmut       number of atoms mutated from initial to final state
c     imut       atomic sites differing in initial and final state
c     type0      atom type of each atom in the initial state system
c     class0     atom class of each atom in the initial state system
c     type1      atom type of each atom in the final state system
c     class1     atom class of each atom in the final state system
c     alter      true if an atom is to be mutated, false otherwise
c     eupdated   flag to mark updated energy value in BAR method
c
c
      integer nmut,imut
      integer type0,class0
      integer type1,class1
      real*8 lambda,elmd,vlmd
      real*8 scexp,scalpha
      real*8 eslvt,eslut
      logical alter,eupdated
      common /mutant/ lambda,elmd,vlmd,scexp,scalpha,eslvt,eslut,
     &                nmut,imut(maxatm),type0(maxatm),class0(maxatm),
     &                type1(maxatm),class1(maxatm),alter(maxatm),
     &                eupdated
