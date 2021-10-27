c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2021  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine spherharm  --  spherical harmonic multipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "spherharm" converts local coordinate frame multipoles to
c     their values in standard spherical harmonics
c
c
      subroutine spherharm
      use atoms
      use math
      use mpole
      integer i
      real*8, allocatable :: spole(:,:)
c
c
c     allocate an array to hold spherical harmonic multipoles
c
      allocate (spole(maxpole,n)
c
c     compute and store the multipoles in spherical harmonics
c     (q -> Q_00, z -> Q_10, x -> Q_11c, y -> Q_11s, zz -> Q_20,
c     xz -> Q_21c, xz -> Q_21c, xx-yy -> Q_22c, xy -> Q_22s)
c
      do i = 1, n
         spole(1,i) = pole(1,i)
         spole(2,i) = pole(4,i)
         spole(3,i) = pole(2,i)
         spole(4,i) = pole(3,i)
         spole(5,i) = pole(13,i)
         spole(6,i) = 2.0d0 * root3 * pole(7,i)
         spole(7,i) = 2.0d0 * root3 * pole(10,i)
         spole(8,i) = root3 * (pole(5,i)-pole(9,i))
         spole(9,i) = 2.0d0 * root3 * pole(6,i)
      end do
      return
      end
