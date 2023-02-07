c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module repel  --  Pauli repulsion for current structure  ##        
c     ##                                                           ##
c     ###############################################################
c
c
c     nrep      total number of repulsion sites in the system
c     irep      number of the atom for each repulsion site
c     sizpr     Pauli repulsion size parameter value at each site
c     dmppr     Pauli repulsion alpha damping value at each site
c     elepr     Pauli repulsion valence electrons at each site
c     repole    repulsion Cartesian multipoles in the local frame
c     rrepole   repulsion Cartesian multipoles in the global frame
c
c
      module repel
      implicit none
      integer nrep
      integer, allocatable :: irep(:)
      real*8, allocatable :: sizpr(:)
      real*8, allocatable :: dmppr(:)
      real*8, allocatable :: elepr(:)
      real*8, allocatable :: repole(:,:)
      real*8, allocatable :: rrepole(:,:)
      save
      end
