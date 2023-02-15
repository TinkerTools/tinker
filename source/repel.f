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
c     replist   repulsion multipole site for each atom (0=none)
c     sizpr     Pauli repulsion size parameter value for each atom
c     dmppr     Pauli repulsion alpha damping value for each atom
c     elepr     Pauli repulsion valence electrons for each atom
c     repole    repulsion Cartesian multipoles in the local frame
c     rrepole   repulsion Cartesian multipoles in the global frame
c
c
      module repel
      implicit none
      integer nrep
      integer, allocatable :: irep(:)
      integer, allocatable :: replist(:)
      real*8, allocatable :: sizpr(:)
      real*8, allocatable :: dmppr(:)
      real*8, allocatable :: elepr(:)
      real*8, allocatable :: repole(:,:)
      real*8, allocatable :: rrepole(:,:)
      save
      end
