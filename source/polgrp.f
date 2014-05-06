c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2001  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  module polgrp  --  polarization group connectivity lists   ##
c     ##                                                             ##
c     #################################################################
c
c
c     np11     number of atoms in polarization group of each atom
c     np12     number of atoms in groups 1-2 to each atom
c     np13     number of atoms in groups 1-3 to each atom
c     np14     number of atoms in groups 1-4 to each atom
c     ip11     atom numbers of atoms in same group as each atom
c     ip12     atom numbers of atoms in groups 1-2 to each atom
c     ip13     atom numbers of atoms in groups 1-3 to each atom
c     ip14     atom numbers of atoms in groups 1-4 to each atom
c
c
      module polgrp
      implicit none
      integer, allocatable :: np11(:)
      integer, allocatable :: np12(:)
      integer, allocatable :: np13(:)
      integer, allocatable :: np14(:)
      integer, allocatable :: ip11(:,:)
      integer, allocatable :: ip12(:,:)
      integer, allocatable :: ip13(:,:)
      integer, allocatable :: ip14(:,:)
      save
      end
