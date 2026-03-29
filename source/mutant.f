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
c     nmut      number of atoms mutated from initial to final state
c     vcouple   van der Waals lambda type (0=decouple, 1=annihilate)
c     emdtexp   multipole exponent for dual topology interpolation
c     epdtexp   polarization exponent for dual topology interpolation
c     imut      atom sites differing in initial and final state
c     type0     atom type of each atom in the initial state system
c     class0    atom class of each atom in the initial state system
c     type1     atom type of each atom in the final state system
c     class1    atom class of each atom in the final state system
c     lambda    generic weighting between initial and final states
c     vlambda   state weighting value for van der Waals potentials
c     elambda   state weighting value for electrostatic potentials
c     tlambda   state weighting value for torsional potential
c     scexp     scale factor for soft core buffered 14-7 potential
c     scalpha   scale factor for soft core buffered 14-7 potential
c     use_emis  logical flag for use of multipole interaction scaling
c     use_epis  logical flag for use of polarization interaction scaling
c     use_emdt  logical flag for use of multipole dual topology
c     use_epdt  logical flag for use of polarization dual topology
c     mut       true if an atom is to be mutated, false otherwise
c
c
      module mutant
      implicit none
      integer nmut
      integer vcouple
      integer emdtexp
      integer epdtexp
      integer, allocatable :: imut(:)
      integer, allocatable :: type0(:)
      integer, allocatable :: class0(:)
      integer, allocatable :: type1(:)
      integer, allocatable :: class1(:)
      real*8 lambda
      real*8 vlambda
      real*8 elambda
      real*8 tlambda
      real*8 scexp
      real*8 scalpha
      logical use_emis
      logical use_epis
      logical use_emdt
      logical use_epdt
      logical, allocatable :: mut(:)
      save
      end
