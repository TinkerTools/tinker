c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar  --  induced dipole polarization energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar
      use sizes
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use limits
      use math
      use mpole
      use polar
      use potent
      use units
      implicit none
      integer i,j,ii
      real*8 e,epo,f,fi
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
c
c
c     zero out the total polarization energy
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     set the energy conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     initialize local variables for OpenMP calculation
c
      epo = 0.0d0
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,e,fi)
!$OMP DO reduction(+:epo) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            fi = f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = fi * uind(j,i) * udirp(j,i)
               epo = epo + e
            end do
         end if
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      ep = ep + epo
c
c     compute the cell dipole boundary correction term
c
      if (use_ewald) then
         if (boundary .eq. 'VACUUM') then
            f = electric / dielec
            xd = 0.0d0
            yd = 0.0d0
            zd = 0.0d0
            xu = 0.0d0
            yu = 0.0d0
            zu = 0.0d0
            do i = 1, npole
               ii = ipole(i)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               uix = uind(1,i)
               uiy = uind(2,i)
               uiz = uind(3,i)
               xd = xd + dix + rpole(1,i)*x(ii)
               yd = yd + diy + rpole(1,i)*y(ii)
               zd = zd + diz + rpole(1,i)*z(ii)
               xu = xu + uix
               yu = yu + uiy
               zu = zu + uiz
            end do
            e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xu+yd*yu+zd*zu)
            ep = ep + e
         end if
      end if
      return
      end
