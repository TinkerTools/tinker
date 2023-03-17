c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine exrepel1  --  exch repulsion energy & derivs  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "exrepel1" calculates the exchange repulsion energy and first
c     derivatives with respect to Cartesian coordinates
c
c     literature reference:
c
c     TBD
c
c
      subroutine exrepel1
      use limits
      implicit none
c
c
c     choose the method for summing over pairwise interactions
c
c      if (use_mlist) then
c         call exrepel1b
c      else
c         call exrepel1a
c      end if
      call exrepel1a
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exrepel1a  --  exch repulsion analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exrepel1a" calculates the exchange repulsion energy and first
c     derivatives with respect to Cartesian coordinates using a
c     pairwise double loop
c
c
      subroutine exrepel1a
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use energi
      use group
      use mpole
      use mutant
      use repel
      use reppot
      use shunt
      use units
      use usage
      use virial
      use xrepel
      implicit none
      integer i,j,k
c
c
c     zero out the repulsion energy and derivatives
c
      er = 0.0d0
      do i = 1, n
         der(1,i) = 0.0d0
         der(2,i) = 0.0d0
         der(3,i) = 0.0d0
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     determine pseudo orbital coefficients
c
      call solvcoeff
c
c     rotate the coefficient components into the global frame
c
      call rotcoeff

      if (nrep .eq. 0)  return
      return
      end