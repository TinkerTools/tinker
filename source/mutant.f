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
c     mutg        alchemical group of each atom (0=not mutated, 1=ligand
c                   or ligand A, 2=ligand B of a relative dual topology)
c     mutcls      class of van der Waals type pairs between two groups
c                   (-1=skipped, 0=full weight, 1 or 2=soft core of that
c                   ligand group)
c     fluxsc      lambda scale applied to the charge flux of each atom
c     pdltflx     unscaled charge flux change in charge of each atom
c     type0       atom type of each atom in the initial state system
c     type1       atom type of each atom in the final state system
c     elambda     state weighting value for electrostatic potentials
c     lambda      generic weighting between initial and final states
c     plambda     state weighting value for polarization potentials
c     scalphar    shift scale for soft core Pauli repulsion in Ang^2
c     scalphalj   shift scale for soft core Lennard-Jones potential
c     scalphav    scale factor for soft core buffered 14-7 potential
c     scexp       scale factor for soft core buffered 14-7 potential
c     scexpm      lambda exponent of soft core buffered 14-7 shift
c     scexpp      one minus lambda exponent of soft core 14-7 shift
c     scexpr      lambda exponent for soft core Pauli repulsion
c     scexprm     lambda exponent of soft core repulsion shift
c     scexprp     one minus lambda exponent of soft core repulsion shift
c     scrdsp      switch radius for soft core damped dispersion in Ang
c     scexpd      lambda exponent for soft core damped dispersion
c     scexpdq     added lambda exponent of soft core dispersion at r=0
c     scexplj     lambda exponent for soft core Lennard-Jones potential
c     scexpljm    lambda exponent of soft core Lennard-Jones shift
c     scexpljp    one minus lambda exponent of soft core LJ shift
c     tlambda     state weighting value for torsional potential
c     vlambda     state weighting value for van der Waals potentials
c     mutfield    flag restricting dfield to the mpoles of mutated atoms
c     mutpsc      flag to collect the mutated partner overlap in alterpol
c     setelambda  flag that elambda was set by its own keyword
c     setplambda  flag that plambda was set by its own keyword
c     setvlambda  flag that vlambda was set by its own keyword
c     use_rel     flag to use the staged two-ligand relative schedule
c                   being built for a relative dual topo energy
c
c
      module mutant
      implicit none
      integer nmut
      integer nmutb
      integer vcouple
      integer mutcls(0:2,0:2)
      integer, allocatable :: class0(:)
      integer, allocatable :: class1(:)
      integer, allocatable :: imut(:)
      integer, allocatable :: mutg(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: type1(:)
      real*8 elambda
      real*8 lambda
      real*8 plambda
      real*8 scalphar
      real*8 scalphalj
      real*8 scalphav
      real*8 scexp
      real*8 scexpm
      real*8 scexpp
      real*8 scexpr
      real*8 scexprm
      real*8 scexprp
      real*8 scrdsp
      real*8 scexpd
      real*8 scexpdq
      real*8 scexplj
      real*8 scexpljm
      real*8 scexpljp
      real*8 tlambda
      real*8 vlambda
      real*8, allocatable :: fluxsc(:)
      real*8, allocatable :: pdltflx(:)
      logical mutfield
      logical mutpsc
      logical setelambda
      logical setplambda
      logical setvlambda
      logical use_rel
      save
      end
