c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2001 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program polarize  --  compute the molecular polarizability  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "polarize" computes the molecular polarizability by applying
c     an external field along each axis followed by diagonalization
c     of the resulting polarizability tensor
c
c
      program polarize
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      integer i
      real*8 addu,malpha
      real*8 exfield(3)
      real*8 umol(3),dalpha(3)
      real*8 work1(3),work2(3)
      real*8 alpha(3,3)
      real*8 valpha(3,3)
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call field
      call molecule
      call kpolar
c
c     sum atomic polarizabilities to get additive molecular value
c
      if (.not. use_polar) then
         write (iout,10)
   10    format (/,' POLARIZE  --  Dipole Polarizability',
     &              ' is Not in Use')
         call fatal
      end if
      addu = 0.0d0
      do i = 1, npole
         addu = polarity(i) + addu
      end do
      if (nmol .eq. 1) then
         write (iout,20)  addu
   20    format (/,' Additive Molecular Polarizability :',f15.4)
      else
         write (iout,30)  addu
   30    format (/,' Additive Total Polarizability :',4x,f15.4)
      end if
c
c     compute each column of the polarizability tensor
c
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(1) = 0.001d0
      call moluind (exfield,umol)
      alpha(1,1) = umol(1) / exfield(1)
      alpha(2,1) = umol(2) / exfield(1)
      alpha(3,1) = umol(3) / exfield(1)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(2) = 0.001d0
      call moluind (exfield,umol)
      alpha(1,2) = umol(1) / exfield(2)
      alpha(2,2) = umol(2) / exfield(2)
      alpha(3,2) = umol(3) / exfield(2)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(3) = 0.001d0
      call moluind (exfield,umol)
      alpha(1,3) = umol(1) / exfield(3)
      alpha(2,3) = umol(2) / exfield(3)
      alpha(3,3) = umol(3) / exfield(3)
c
c     print out the full polarizability tensor
c
      if (nmol .eq. 1) then
         write (iout,40)
   40    format (/,' Molecular Polarizability Tensor :',/)
      else
         write (iout,50)
   50    format (/,' Total Polarizability Tensor:',/)
      end if
      write (iout,60)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                 alpha(2,1),alpha(2,2),alpha(2,3),
     &                 alpha(3,1),alpha(3,2),alpha(3,3)
   60 format (15x,3f12.4,/,15x,3f12.4,/,15x,3f12.4)
c
c     diagonalize the tensor and get molecular polarizability
c
      call jacobi (3,3,alpha,dalpha,valpha,work1,work2)
      write (iout,70)
   70 format (/,' Polarizability Tensor Eigenvalues :',/)
      write (iout,80)  dalpha(1),dalpha(2),dalpha(3)
   80 format (15x,3f12.4)
      malpha = (dalpha(1)+dalpha(2)+dalpha(3)) / 3.0d0
      if (nmol .eq. 1) then
         write (iout,90)  malpha
   90    format (/,' Interactive Molecular Polarizability :',f12.4)
      else
         write (iout,100)  malpha
  100    format (/,' Interactive Total Polarizability :',4x,f12.4)
      end if
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine moluind  --  molecular induced dipole in field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "moluind" computes the molecular induced dipole components
c     in the presence of an external electric field
c
c
      subroutine moluind (exfield,umol)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer i,j,iter
      integer maxiter
      real*8 eps,epsold
      real*8 umol(3)
      real*8 exfield(3)
      real*8 udir(3,maxatm)
      real*8 uold(3,maxatm)
      real*8 field(3,maxatm)
      logical done
c
c
c     zero out the atomic and molecular induced dipole components
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
         end do
      end do
      do j = 1, 3
         umol(j) = 0.0d0
      end do
c
c     compute induced dipoles as polarizability times external field
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarity(i) * exfield(j)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 1.0d0
         do i = 1, npole
            do j = 1, 3
               udir(j,i) = uind(j,i)
            end do
         end do
c
c     compute mutual induced dipole moments by an iterative method
c
         dowhile (.not. done)
            do i = 1, npole
               do j = 1, 3
                  field(j,i) = 0.0d0
               end do
            end do
            call ufield (field)
c
c     check to see if the mutual induced dipoles have converged
c
            iter = iter + 1
            epsold = eps
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  uold(j,i) = uind(j,i)
                  uind(j,i) = udir(j,i) + polarity(i)*field(j,i)
                  uind(j,i) = uold(j,i) + polsor*(uind(j,i)-uold(j,i))
                  eps = eps + (uind(j,i)-uold(j,i))**2
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     sum up the total molecular induced dipole components
c
         do i = 1, npole
            umol(1) = umol(1) + uind(1,i)
            umol(2) = umol(2) + uind(2,i)
            umol(3) = umol(3) + uind(3,i)
         end do
c
c     print a warning if induced dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ufield  --  electric field from induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ufield" finds the field at each polarizable site due to the
c     induced dipoles at the other sites using Thole's method to
c     damp the field at close range
c
c
      subroutine ufield (field)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      integer i,j,k
      integer ii,kk
      real*8 r,r2,rr3,rr5
      real*8 xr,yr,zr
      real*8 uir,ukr
      real*8 uix,uiy,uiz
      real*8 ukx,uky,ukz
      real*8 pdi,pti
      real*8 pgamma,damp
      real*8 scale3,scale5
      real*8 fi(3),fk(3)
      real*8 pscale(maxatm)
      real*8 field(3,maxatm)
c
c
c     loop over pairs of sites incrementing the electric field
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         do j = i+1, npole
            pscale(ipole(j)) = 1.0d0
         end do
         do j = 1, np11(ii)
            pscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            pscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            pscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            pscale(ip14(j,ii)) = u4scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            r2 = xr*xr + yr* yr + zr*zr
            ukx = uind(1,k)
            uky = uind(2,k)
            ukz = uind(3,k)
            r = sqrt(r2)
            rr3 = 1.0d0 / (r*r2)
            rr5 = 3.0d0 * rr3 / r2
            uir = xr*uix + yr*uiy + zr*uiz
            ukr = xr*ukx + yr*uky + zr*ukz
c
c     adjust the field to account for polarization damping
c
            scale3 = 1.0d0
            scale5 = 1.0d0
            damp = pdi * pdamp(k)
            if (damp .ne. 0.0d0) then
               pgamma = min(pti,thole(k))
               damp = -pgamma * (r/damp)**3
               if (damp .gt. -50.0d0) then
                  scale3 = 1.0d0 - exp(damp)
                  scale5 = 1.0d0 - (1.0d0-damp)*exp(damp)
               end if
            end if
            fi(1) = -rr3*ukx*scale3 + rr5*ukr*xr*scale5
            fi(2) = -rr3*uky*scale3 + rr5*ukr*yr*scale5
            fi(3) = -rr3*ukz*scale3 + rr5*ukr*zr*scale5
            fk(1) = -rr3*uix*scale3 + rr5*uir*xr*scale5
            fk(2) = -rr3*uiy*scale3 + rr5*uir*yr*scale5
            fk(3) = -rr3*uiz*scale3 + rr5*uir*zr*scale5
            do j = 1, 3
               fi(j) = fi(j) * pscale(kk)
               fk(j) = fk(j) * pscale(kk)
            end do
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               field(j,i) = field(j,i) + fi(j)
               field(j,k) = field(j,k) + fk(j)
            end do
         end do
      end do
      return
      end
