c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  module ktorsn  --  torsional angle forcefield parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     maxnt    maximum number of torsional angle parameter entries
c     maxnt5   maximum number of 5-membered ring torsion entries
c     maxnt4   maximum number of 4-membered ring torsion entries
c
c     t1       torsional parameters for standard 1-fold rotation
c     t2       torsional parameters for standard 2-fold rotation
c     t3       torsional parameters for standard 3-fold rotation
c     t4       torsional parameters for standard 4-fold rotation
c     t5       torsional parameters for standard 5-fold rotation
c     t6       torsional parameters for standard 6-fold rotation
c     t15      torsional parameters for 1-fold rotation in 5-ring
c     t25      torsional parameters for 2-fold rotation in 5-ring
c     t35      torsional parameters for 3-fold rotation in 5-ring
c     t45      torsional parameters for 4-fold rotation in 5-ring
c     t55      torsional parameters for 5-fold rotation in 5-ring
c     t65      torsional parameters for 6-fold rotation in 5-ring
c     t14      torsional parameters for 1-fold rotation in 4-ring
c     t24      torsional parameters for 2-fold rotation in 4-ring
c     t34      torsional parameters for 3-fold rotation in 4-ring
c     t44      torsional parameters for 4-fold rotation in 4-ring
c     t54      torsional parameters for 5-fold rotation in 4-ring
c     t64      torsional parameters for 6-fold rotation in 4-ring
c     kt       string of atom classes for torsional angles
c     kt5      string of atom classes for 5-ring torsions
c     kt4      string of atom classes for 4-ring torsions
c
c
      module ktorsn
      implicit none
      integer maxnt
      integer maxnt5
      integer maxnt4
      real*8, allocatable :: t1(:,:)
      real*8, allocatable :: t2(:,:)
      real*8, allocatable :: t3(:,:)
      real*8, allocatable :: t4(:,:)
      real*8, allocatable :: t5(:,:)
      real*8, allocatable :: t6(:,:)
      real*8, allocatable :: t15(:,:)
      real*8, allocatable :: t25(:,:)
      real*8, allocatable :: t35(:,:)
      real*8, allocatable :: t45(:,:)
      real*8, allocatable :: t55(:,:)
      real*8, allocatable :: t65(:,:)
      real*8, allocatable :: t14(:,:)
      real*8, allocatable :: t24(:,:)
      real*8, allocatable :: t34(:,:)
      real*8, allocatable :: t44(:,:)
      real*8, allocatable :: t54(:,:)
      real*8, allocatable :: t64(:,:)
      character*16, allocatable :: kt(:)
      character*16, allocatable :: kt5(:)
      character*16, allocatable :: kt4(:)
      save
      end
