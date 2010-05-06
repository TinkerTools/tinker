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
      implicit none
      include 'sizes.i'
      include 'bath.i'
      include 'boxes.i'
      include 'bound.i'
      include 'units.i'
      include 'virial.i'
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
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'group.i'
      include 'math.i'
      include 'mdstuf.i'
      include 'usage.i'
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
         xbox = scale * xbox
         ybox = scale * ybox
         zbox = scale * zbox
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
                  x(i) = scale * x(i)
                  y(i) = scale * y(i)
                  z(i) = scale * z(i)
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
c     this version only implements isotropic volume changes
c
c
      subroutine pmonte (epot,temp)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'group.i'
      include 'mdstuf.i'
      include 'molcul.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer start,stop
      real*8 epot,temp
      real*8 energy,random
      real*8 third,weigh
      real*8 step,scale
      real*8 enew,diff
      real*8 xcm,ycm,zcm
      real*8 xmove,ymove,zmove
      real*8 vold,xboxold
      real*8 yboxold,zboxold
      real*8 kt,de,dv,lnv
      real*8 term,expterm
      real*8 xold(maxatm)
      real*8 yold(maxatm)
      real*8 zold(maxatm)
c
c
c     make volume move, save old box size and coordinates
c
      if (random() .lt. 1.0d0/dble(voltrial)) then
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
c     get scale factor that reflects the chosen volume change
c
         step = volmove * (2.0d0*random()-1.0d0)
         volbox = volbox + step
         third = 1.0d0 / 3.0d0
         scale = (volbox/vold)**third
         if (monoclinic) then
            term = 1.0d0 / beta_sin
            scale = 1.0d0 + (scale-1.0d0)*term**third
         else if (triclinic) then
            term = 1.0d0 / (gamma_sin*gamma_term)
            scale = 1.0d0 + (scale-1.0d0)*term**third
         else if (octahedron) then
            term = 2.0d0
            scale = 1.0d0 + (scale-1.0d0)*term**third
         end if
c
c     set the new box dimensions and other lattice values
c
         xbox = xbox * scale
         ybox = ybox * scale
         zbox = zbox * scale
         call lattice
c
c     scale the coordinates by groups, molecules or atoms
c
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
                  x(i) = scale * x(i)
                  y(i) = scale * y(i)
                  z(i) = scale * z(i)
               end if
            end do
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
            kt = gasconst * kelvin0
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
      return
      end