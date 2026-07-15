c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module mutant  --  free energy calculation hybrid atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     nmut        number of atoms mutated from initial to final state
c     nmutb       number of atoms in the second ligand group (group B)
c     vcouple     van der Waals lambda type (0=decouple, 1=annihilate)
c     imut        atom sites differing in initial and final state
c     type0       atom type of each atom in the initial state system
c     class0      atom class of each atom in the initial state system
c     type1       atom type of each atom in the final state system
c     class1      atom class of each atom in the final state system
c     mutg        alchemical group of each atom (0=env, 1=ligand A,
c                   2=ligand B) used for relative free energy dual topo
c     lambda      generic weighting between initial and final states
c     vlambda     state weighting value for van der Waals potentials
c     elambda     state weighting value for electrostatic potentials
c     tlambda     state weighting value for torsional potential
c     scexp       scale factor for soft core buffered 14-7 potential
c     scalpha     scale factor for soft core buffered 14-7 potential
c     use_rel     flag to use two-ligand relative dual topology
c     use_subsys  flag that a parameter-zeroed subsystem is active
c     mut         true if an atom is to be mutated, false otherwise
c     subon       true if an atom is active in the subsystem currently
c                   being built for a relative dual topo energy
c
c
      module mutant
      implicit none
      integer nmut
      integer nmutb
      integer vcouple
      integer, allocatable :: imut(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: class0(:)
      integer, allocatable :: type1(:)
      integer, allocatable :: class1(:)
      integer, allocatable :: mutg(:)
      real*8 lambda
      real*8 vlambda
      real*8 elambda
      real*8 tlambda
      real*8 scexp
      real*8 scalpha
      logical use_rel
      logical use_subsys
      logical, allocatable :: mut(:)
      logical, allocatable :: subon(:)
      save
      end
