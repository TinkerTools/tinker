c
c
c     ##############################################################
c     ##  COPYRIGHT (C)  2024  by  Moses KJ Chung & Jay W Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine uniquetyp  --  get unique atom type information  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "uniquetyp" determines the number of unique types and map between
c     atom types and unique types
c
c
      subroutine uniquetyp
      use atoms
      use uatom
      implicit none
      integer i,j
      integer at
      real*8 xx
c
c
c     initialize nunique, utype, and utypeinv
c
      nunique = 0
      do i = 1, maxtyp
         utype(i) = 0
         utypeinv(i) = 0
      end do
c
c     loop through atoms to find unique atom types
c
      do i = 1, n
         at = type(i)
c
c     check if the current atom type is already in unique type array
c
         do j = 1, nunique
            if (at .eq. utype(j))  goto 10
         end do
c
c     set unique atom types
c
         nunique = nunique + 1
         utype(nunique) = at
         utypeinv(at) = nunique
10       continue
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine velunique  --  unique atom type velocity  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "velunique" computes the unique atom type velocities
c
c
      subroutine velunique
      use atomid
      use atoms
      use moldyn
      use output
      use uatom
      implicit none
      integer i
      integer ut
c
c
c     zero out velocity
c
      do i = 1, nunique
         utv1(1,i) = 0.0d0
         utv1(2,i) = 0.0d0
         utv1(3,i) = 0.0d0
      end do
!$OMP PARALLEL default(private)
!$OMP& shared(n,utv1,type,utypeinv,v)
!$OMP DO reduction(+:utv1) schedule(guided)
c
c     compute velocity by unique atom type
c
      do i = 1, n
         ut = utypeinv(type(i))
         utv1(1,ut) = utv1(1,ut) + v(1,i)
         utv1(2,ut) = utv1(2,ut) + v(2,i)
         utv1(3,ut) = utv1(3,ut) + v(3,i)
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
