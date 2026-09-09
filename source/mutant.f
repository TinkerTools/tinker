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
c     class0      atom class of each atom in the initial state system
c     class1      atom class of each atom in the final state system
c     imut        atom sites differing in initial and final state
c     mutg        alchemical group of each atom (0=env, 1=ligand A,
c                   2=ligand B) used for relative free energy dual topo
c     type0       atom type of each atom in the initial state system
c     type1       atom type of each atom in the final state system
c     elambda     state weighting value for electrostatic potentials
c     lambda      generic weighting between initial and final states
c     plambda     state weighting value for polarization potentials
c     scalpham    offset factor for soft core electrostatics
c     scalphav    scale factor for soft core buffered 14-7 potential
c     scexp       scale factor for soft core buffered 14-7 potential
c     tlambda     state weighting value for torsional potential
c     vlambda     state weighting value for van der Waals potentials
c     mutfield    flag restricting dfield to the mpoles by mut atoms
c     setelambda  flag that elambda was set by its own keyword
c     setesoft    flag that the ELE-SOFTCORE keyword was given
c     setplambda  flag that plambda was set by its own keyword
c     setvlambda  flag that vlambda was set by its own keyword
c     use_esoft   flag to soft core real space multipole interactions
c     use_past    flag for absolute single topology polarization
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
      integer, allocatable :: class0(:)
      integer, allocatable :: class1(:)
      integer, allocatable :: imut(:)
      integer, allocatable :: mutg(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: type1(:)
      real*8 elambda
      real*8 lambda
      real*8 plambda
      real*8 scalpham
      real*8 scalphav
      real*8 scexp
      real*8 tlambda
      real*8 vlambda
      logical mutfield
      logical setelambda
      logical setesoft
      logical setplambda
      logical setvlambda
      logical use_esoft
      logical use_past
      logical use_rel
      logical use_subsys
      logical, allocatable :: mut(:)
      logical, allocatable :: subon(:)
      save
      end
