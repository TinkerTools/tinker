c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pressure  --  constant pressure via barostat  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pressure" uses the internal virial to find the pressure
c     in a periodic box and maintains a constant desired pressure
c     via a barostat method
c
c
      subroutine pressure (dt,epot,ekin,temp,pres,stress)
      use sizes
      use bath
      use boxes
      use bound
      use units
      use virial
      implicit none
      integer i,j
      real*8 dt,epot
      real*8 temp,pres
      real*8 factor
      real*8 ekin(3,3)
      real*8 stress(3,3)
c
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     calculate the stress tensor for anisotropic systems
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
c
c     set isotropic pressure to the average of tensor diagonal
c
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     use either the Berendsen or Monte Carlo barostat method
c
      if (isobaric) then
         if (barostat .eq. 'BERENDSEN') then
            call pscale (dt,pres,stress)
         else if (barostat .eq. 'MONTECARLO') then
            call pmonte (epot,temp)
         end if
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine pscale  --  Berendsen barostat via scaling  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "pscale" implements a Berendsen barostat by scaling the
c     coordinates and box dimensions via coupling to an external
c     constant pressure bath
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     S. E. Feller, Y. Zhang, R. W. Pastor, B. R. Brooks, "Constant
c     Pressure Molecular Dynamics Simulation: The Langevin Piston
c     Method", Journal of Chemical Physics, 103, 4613-4621 (1995)
c
c     code for anisotropic pressure coupling was provided by Guido
c     Raos, Dipartimento di Chimica, Politecnico di Milano, Italy
c
c
      subroutine pscale (dt,pres,stress)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use group
      use math
      use mdstuf
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 dt,pres
      real*8 weigh,cosine
      real*8 scale,third
      real*8 xcm,xmove
      real*8 ycm,ymove
      real*8 zcm,zmove
      real*8 stress(3,3)
      real*8 temp(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
c
c
c     find the isotropic scale factor for constant pressure
c
      if (.not. anisotrop) then
         scale = 1.0d0
         third = 1.0d0 / 3.0d0
         scale = (1.0d0 + (dt*compress/taupres)*(pres-atmsph))**third
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
c
c     propagate the new box dimensions to other lattice values
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
         if (integrate .ne. 'RIGIDBODY') then
            do i = 1, n
               if (use(i)) then
                  x(i) = x(i) * scale
                  y(i) = y(i) * scale
                  z(i) = z(i) * scale
               end if
            end do
c
c     couple to pressure bath via center of mass of rigid bodies
c
         else
            scale = scale - 1.0d0
            do i = 1, ngrp
               start = igrp(1,i)
               stop = igrp(2,i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = start, stop
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               xmove = scale * xcm/grpmass(i)
               ymove = scale * ycm/grpmass(i)
               zmove = scale * zcm/grpmass(i)
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
            end do
         end if
c
c     find the anisotropic scale factors for constant pressure
c
      else
         scale = dt*compress / (3.0d0*taupres)
         do i = 1, 3
            do j = 1, 3
               if (j. eq. i) then
                  ascale(j,i) = 1.0d0 + scale*(stress(i,i)-atmsph)
               else
                  ascale(j,i) = scale*stress(j,i)
               end if
            end do
         end do
c
c     modify the current periodic box dimension values
c
         temp(1,1) = xbox
         temp(2,1) = 0.0d0
         temp(3,1) = 0.0d0
         temp(1,2) = ybox * gamma_cos
         temp(2,2) = ybox * gamma_sin
         temp(3,2) = 0.0d0
         temp(1,3) = zbox * beta_cos
         temp(2,3) = zbox * beta_term
         temp(3,3) = zbox * gamma_term
         do i = 1, 3
            do j = 1, 3
               hbox(j,i) = 0.0d0
               do k = 1, 3
                  hbox(j,i) = hbox(j,i) + ascale(j,k)*temp(k,i)
               end do
            end do
         end do
         xbox = sqrt(hbox(1,1)**2 + hbox(2,1)**2 + hbox(3,1)**2)
         ybox = sqrt(hbox(1,2)**2 + hbox(2,2)**2 + hbox(3,2)**2)
         zbox = sqrt(hbox(1,3)**2 + hbox(2,3)**2 + hbox(3,3)**2)
         if (monoclinic) then
            cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                  + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
            beta = radian * acos(cosine)
         else if (triclinic) then
            cosine = (hbox(1,2)*hbox(1,3) + hbox(2,2)*hbox(2,3)
     &                  + hbox(3,2)*hbox(3,3)) / (ybox*zbox)
            alpha = radian * acos(cosine)
            cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                  + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
            beta = radian * acos(cosine)
            cosine = (hbox(1,1)*hbox(1,2) + hbox(2,1)*hbox(2,2)
     &                  + hbox(3,1)*hbox(3,2)) / (xbox*ybox)
            gamma = radian * acos(cosine)
         end if
c
c     propagate the new box dimensions to other lattice values
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
         if (integrate .ne. 'RIGIDBODY') then
            do i = 1, n
               if (use(i)) then
                  x(i) = ascale(1,1)*x(i) + ascale(1,2)*y(i)
     &                      + ascale(1,3)*z(i)
                  y(i) = ascale(2,1)*x(i) + ascale(2,2)*y(i)
     &                      + ascale(2,3)*z(i)
                  z(i) = ascale(3,1)*x(i) + ascale(3,2)*y(i)
     &                      + ascale(3,3)*z(i)
               end if
            end do
c
c     couple to pressure bath via center of mass of rigid bodies
c
         else
            ascale(1,1) = ascale(1,1) - 1.0d0
            ascale(2,2) = ascale(2,2) - 1.0d0
            ascale(3,3) = ascale(3,3) - 1.0d0
            do i = 1, ngrp
               start = igrp(1,i)
               stop = igrp(2,i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = start, stop
                  k = kgrp(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = xcm + y(k)*weigh
                  zcm = xcm + z(k)*weigh
               end do
               xcm = xcm / grpmass(i)
               ycm = ycm / grpmass(i)
               zcm = zcm / grpmass(i)
               xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
     &                    + ascale(1,3)*zcm
               ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
     &                    + ascale(2,3)*zcm
               zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
     &                    + ascale(3,3)*zcm
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
            end do
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pmonte  --  Monte Carlo barostat trial moves  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pmonte" implements a Monte Carlo barostat via random trial
c     changes in the periodic box volume and shape
c
c     literature references:
c
c     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
c     2nd Edition", Academic Press, San Diego, CA, 2002; Section 5.4
c
c     original version written by Alan Grossfield, January 2004;
c     anisotropic modification implemented by Lee-Ping Wang, Stanford
c     University, March 2013
c
c
      subroutine pmonte (epot,temp)
      use sizes
      use atomid
      use atoms
      use bath
      use boxes
      use group
      use math
      use mdstuf
      use molcul
      use units
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 epot,temp,term
      real*8 energy,random
      real*8 third,weigh
      real*8 step,scale
      real*8 enew,rnd6
      real*8 xcm,ycm,zcm
      real*8 vold,vfrac
      real*8 cosine,diff
      real*8 xmove,ymove,zmove
      real*8 xboxold,yboxold,zboxold
      real*8 alphaold,betaold,gammaold
      real*8 kt,de,dv,lnv,expterm
      real*8 temp3(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     save the old lattice parameters, box size and coordinates
c
      if (random() .lt. 1.0d0/dble(voltrial)) then
         xboxold = xbox
         yboxold = ybox
         zboxold = zbox
         alphaold = alpha
         betaold  = beta
         gammaold = gamma
         vold = volbox
         do i = 1, n
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
         end do
         third = 1.0d0 / 3.0d0
c
c     for the isotropic case, change the lattice lengths uniformly
c
         if (.not.anisotrop .or. random().lt.0.5d0) then
            step = volmove * (2.0d0*random()-1.0d0)
            volbox = volbox + step
            scale = (volbox/vold)**third
            xbox = xbox * scale
            ybox = ybox * scale
            zbox = zbox * scale
            call lattice
            if (integrate .eq. 'RIGIDBODY') then
               diff = scale - 1.0d0
               do i = 1, ngrp
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = igrp(1,i)
                  stop = igrp(2,i)
                  do j = start, stop
                     k = kgrp(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xmove = diff * xcm/grpmass(i)
                  ymove = diff * ycm/grpmass(i)
                  zmove = diff * zcm/grpmass(i)
                  do j = start, stop
                     k = kgrp(j)
                     x(k) = x(k) + xmove
                     y(k) = y(k) + ymove
                     z(k) = z(k) + zmove
                  end do
               end do
            else if (volscale .eq. 'MOLECULAR') then
               diff = scale - 1.0d0
               do i = 1, nmol
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = imol(1,i)
                  stop = imol(2,i)
                  do j = start, stop
                     k = kmol(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xmove = diff * xcm/molmass(i)
                  ymove = diff * ycm/molmass(i)
                  zmove = diff * zcm/molmass(i)
                  do j = start, stop
                     k = kmol(j)
                     if (use(k)) then
                        x(k) = x(k) + xmove
                        y(k) = y(k) + ymove
                        z(k) = z(k) + zmove
                     end if
                  end do
               end do
            else
               do i = 1, n
                  if (use(i)) then
                     x(i) = x(i) * scale
                     y(i) = y(i) * scale
                     z(i) = z(i) * scale
                  end if
               end do
            end if
c
c     for anisotropic case alter lattice angles, then scale lengths
c
         else
            rnd6 = 6.0d0*random()
            step  = volmove * (2.0d0*random()-1.0d0)
            scale = (1.0d0+step/vold)**third
            ascale(1,1) = 1.0d0
            ascale(2,2) = 1.0d0
            ascale(3,3) = 1.0d0
            if (monoclinic .or. triclinic) then
               if (rnd6 .lt. 1.0d0) then
                  ascale(1,1) = scale
               else if (rnd6 .lt. 2.0d0) then
                  ascale(2,2) = scale
               else if (rnd6 .lt. 3.0d0) then
                  ascale(3,3) = scale
               else if (rnd6 .lt. 4.0d0) then
                  ascale(1,2) = scale - 1.0d0
                  ascale(2,1) = scale - 1.0d0
               else if (rnd6 .lt. 5.0d0) then
                  ascale(1,3) = scale - 1.0d0
                  ascale(3,1) = scale - 1.0d0
               else
                  ascale(2,3) = scale - 1.0d0
                  ascale(3,2) = scale - 1.0d0
               end if
            else
               if (rnd6 .lt. 2.0d0) then
                  ascale(1,1) = scale
               else if (rnd6 .lt. 4.0d0) then
                  ascale(2,2) = scale
               else
                  ascale(3,3) = scale
               end if
            end if
c
c     modify the current periodic box dimension values
c
            temp3(1,1) = xbox
            temp3(2,1) = 0.0d0
            temp3(3,1) = 0.0d0
            temp3(1,2) = ybox * gamma_cos
            temp3(2,2) = ybox * gamma_sin
            temp3(3,2) = 0.0d0
            temp3(1,3) = zbox * beta_cos
            temp3(2,3) = zbox * beta_term
            temp3(3,3) = zbox * gamma_term
            do i = 1, 3
               do j = 1, 3
                  hbox(j,i) = 0.0d0
                  do k = 1, 3
                     hbox(j,i) = hbox(j,i) + ascale(j,k)*temp3(k,i)
                  end do
               end do
            end do
            xbox = sqrt(hbox(1,1)**2 + hbox(2,1)**2 + hbox(3,1)**2)
            ybox = sqrt(hbox(1,2)**2 + hbox(2,2)**2 + hbox(3,2)**2)
            zbox = sqrt(hbox(1,3)**2 + hbox(2,3)**2 + hbox(3,3)**2)
            if (monoclinic) then
               cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                     + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
               beta = radian * acos(cosine)
            else if (triclinic) then
               cosine = (hbox(1,2)*hbox(1,3) + hbox(2,2)*hbox(2,3)
     &                     + hbox(3,2)*hbox(3,3)) / (ybox*zbox)
               alpha = radian * acos(cosine)
               cosine = (hbox(1,1)*hbox(1,3) + hbox(2,1)*hbox(2,3)
     &                     + hbox(3,1)*hbox(3,3)) / (xbox*zbox)
               beta = radian * acos(cosine)
               cosine = (hbox(1,1)*hbox(1,2) + hbox(2,1)*hbox(2,2)
     &                     + hbox(3,1)*hbox(3,2)) / (xbox*ybox)
               gamma = radian * acos(cosine)
            end if
c
c     find the new box dimensions and other lattice values
c
            call lattice
            vfrac = vold / volbox
            scale = vfrac**third
            xbox = xbox * scale
            ybox = ybox * scale
            zbox = zbox * scale
            call lattice
c
c     scale the coordinates by groups, molecules or atoms
c
            if (integrate .eq. 'RIGIDBODY') then
               ascale(1,1) = ascale(1,1) - 1.0d0
               ascale(2,2) = ascale(2,2) - 1.0d0
               ascale(3,3) = ascale(3,3) - 1.0d0
               do i = 1, ngrp
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = igrp(1,i)
                  stop = igrp(2,i)
                  do j = start, stop
                     k = kgrp(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xcm = xcm / grpmass(i)
                  ycm = ycm / grpmass(i)
                  zcm = zcm / grpmass(i)
                  xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
     &                       + ascale(1,3)*zcm
                  ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
     &                       + ascale(2,3)*zcm
                  zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
     &                       + ascale(3,3)*zcm
                  do j = start, stop
                     k = kgrp(j)
                     x(k) = x(k) + xmove
                     y(k) = y(k) + ymove
                     z(k) = z(k) + zmove
                  end do
               end do
            else if (volscale .eq. 'MOLECULAR') then
               ascale(1,1) = ascale(1,1) - 1.0d0
               ascale(2,2) = ascale(2,2) - 1.0d0
               ascale(3,3) = ascale(3,3) - 1.0d0
               do i = 1, nmol
                  xcm = 0.0d0
                  ycm = 0.0d0
                  zcm = 0.0d0
                  start = imol(1,i)
                  stop = imol(2,i)
                  do j = start, stop
                     k = kmol(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xcm = xcm / molmass(i)
                  ycm = ycm / molmass(i)
                  zcm = zcm / molmass(i)
                  xmove = ascale(1,1)*xcm + ascale(1,2)*ycm
     &                       + ascale(1,3)*zcm
                  ymove = ascale(2,1)*xcm + ascale(2,2)*ycm
     &                       + ascale(2,3)*zcm
                  zmove = ascale(3,1)*xcm + ascale(3,2)*ycm
     &                       + ascale(3,3)*zcm
                  do j = start, stop
                     k = kmol(j)
                     if (use(k)) then
                        x(k) = x(k) + xmove
                        y(k) = y(k) + ymove
                        z(k) = z(k) + zmove
                     end if
                  end do
               end do
            else
               do i = 1, n
                  if (use(i)) then
                     x(i) = ascale(1,1)*x(i) + ascale(1,2)*y(i)
     &                         + ascale(1,3)*z(i)
                     y(i) = ascale(2,1)*x(i) + ascale(2,2)*y(i)
     &                         + ascale(2,3)*z(i)
                     z(i) = ascale(3,1)*x(i) + ascale(3,2)*y(i)
     &                         + ascale(3,3)*z(i)
                  end if
               end do
            end if
         end if
c
c     find the energy change and PV term for the trial move
c
         enew = energy ()
         de = enew - epot
         dv = atmsph * (volbox-vold) / prescon
c
c     set the entropy of mixing term based on system type
c
         if (isothermal) then
            kt = gasconst * kelvin
         else
            kt = gasconst * temp
         end if
         if (integrate .eq. 'RIGIDBODY') then
            lnv = dble(ngrp) * kt * log(volbox/vold)
         else if (volscale .eq. 'MOLECULAR') then
            lnv = dble(nmol) * kt * log(volbox/vold)
         else
            lnv = dble(nuse) * kt * log(volbox/vold)
         end if
c
c     acceptance ratio from energy change, PV and mixing term
c
         term = -(de+dv-lnv) / kt
         expterm = exp(term)
c
c     reject the step, restore old box size and coordinates
c
         if (random() .gt. expterm) then
            xbox = xboxold
            ybox = yboxold
            zbox = zboxold
            call lattice
            do i = 1, n
               x(i) = xold(i)
               y(i) = yold(i)
               z(i) = zold(i)
            end do
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ptest  --  find pressure via finite-difference  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ptest" compares the virial-based value of dE/dV to an estimate
c     from finite-difference volume changes; also finds the isotropic
c     pressure via finite-differences
c
c     original version written by John D. Chodera, University of
c     California, Berkeley, December 2010
c
c
      subroutine ptest
      use sizes
      use atoms
      use bath
      use bound
      use boxes
      use iounit
      use units
      use virial
      implicit none
      integer i
      real*8 energy,third
      real*8 delta,step,scale
      real*8 vold,xboxold
      real*8 yboxold,zboxold
      real*8 epos,eneg
      real*8 dedv_vir,dedv_fd
      real*8 pres_vir,pres_fd
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
c
c
c     set relative volume change for finite-differences
c
      if (.not. use_bounds)  return
      delta = 0.000001d0
      step = volbox * delta
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
c
c     store original box dimensions and coordinate values
c
      xboxold = xbox
      yboxold = ybox
      zboxold = zbox
      vold = volbox
      do i = 1, n
         xold(i) = x(i)
         yold(i) = y(i)
         zold(i) = z(i)
      end do
c
c     get scale factor to reflect a negative volume change
c
      volbox = vold - step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
      xbox = xboxold * scale
      ybox = yboxold * scale
      zbox = zboxold * scale
      call lattice
      do i = 1, n
         x(i) = xold(i) * scale
         y(i) = yold(i) * scale
         z(i) = zold(i) * scale
      end do
c
c     compute potential energy for negative volume change
c
      eneg = energy ()
c
c     get scale factor to reflect a positive volume change
c
      volbox = vold + step
      third = 1.0d0 / 3.0d0
      scale = (volbox/vold)**third
c
c     set new box dimensions and coordinate values
c
      xbox = xboxold * scale
      ybox = yboxold * scale
      zbox = zboxold * scale
      call lattice
      do i = 1, n
         x(i) = xold(i) * scale
         y(i) = yold(i) * scale
         z(i) = zold(i) * scale
      end do
c
c     compute potential energy for positive volume change
c
      epos = energy ()
c
c     restore original box dimensions and coordinate values
c
      xbox = xboxold
      ybox = yboxold
      zbox = zboxold
      call lattice
      do i = 1, n
         x(i) = xold(i)
         y(i) = yold(i)
         z(i) = zold(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
c
c     get virial and finite difference values of dE/dV
c
      dedv_vir = (vir(1,1)+vir(2,2)+vir(3,3)) / (3.0d0*volbox)
      dedv_fd = (epos-eneg) / (2.0d0*delta*volbox)
      write (iout,10)  dedv_vir
   10 format (/,' dE/dV (Virial-based) :',11x,f15.6,' Kcal/mole/A**3')
      write (iout,20)  dedv_fd
   20 format (' dE/dV (Finite Diff) :',12x,f15.6,' Kcal/mole/A**3')
c
c     compute analytical and finite-difference isotropic pressure
c
      pres_vir = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_vir)
      pres_fd = prescon * (dble(n)*gasconst*kelvin/volbox-dedv_fd)
      if (kelvin .eq. 0.0d0) then
         write (iout,30)  pres_vir
         write (iout,40)  pres_fd
   30    format (/,' Pressure (Analytical, 0 K) :',5x,f15.3,
     &              ' Atmospheres')
   40    format (' Pressure (Numerical, 0 K) :',6x,f15.3,
     &              ' Atmospheres')
      else
         write (iout,50)  nint(kelvin),pres_vir
         write (iout,60)  nint(kelvin),pres_fd
   50    format (/,' Pressure (Analytical,',i4,' K) :',3x,f15.3,
     &              ' Atmospheres')
   60    format (' Pressure (Numerical,',i4,' K) :',4x,f15.3,
     &              ' Atmospheres')
      end if
      return
      end
