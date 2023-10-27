c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module xrepel  --  exch repulsion for current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     nxrep      total number of repulsion sites in the system
c     ixrep      number of the atom for each repulsion site
c     xreplist   repulsion multipole site for each atom (0=none)
c     zpxr       nuclear charge parameter value for each atom
c     dmppxr     exchange repulsion alpha damping value for each atom
c     crpxr      ratio of p/s orbital cofficients for each atom
c     cpxr       local pseudo wavefunction coefficient for each atom
c     rcpxr      global pseudo wavefunction coefficient for each atom
c     xrepole    repulsion Cartesian multipoles in the local frame
c
c
      module xrepel
      implicit none
      integer nxrep
      integer, allocatable :: ixrep(:)
      integer, allocatable :: xreplist(:)
      real*8, allocatable :: zpxr(:)
      real*8, allocatable :: dmppxr(:)
      real*8, allocatable :: crpxr(:)
      real*8, allocatable :: cpxr(:,:)
      real*8, allocatable :: rcpxr(:,:)
      real*8, allocatable :: xrepole(:,:)
      save
      end
