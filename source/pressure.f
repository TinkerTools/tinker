c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine pressure  --  barostat applied at full-step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pressure" uses the internal virial to find the pressure
c     in a periodic box and maintains a constant desired pressure
c     via a barostat method
c
c
      subroutine pressure (dt,epot,ekin,temp,pres,stress)
      use bath
      use boxes
      use bound
      use math
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
      pres = third * (stress(1,1)+stress(2,2)+stress(3,3))
c
c     use the desired barostat to maintain constant pressure
c
      if (isobaric) then
         if (barostat .eq. 'BERENDSEN')  call pscale (dt,pres,stress)
         if (barostat .eq. 'BUSSI')  call pscale (dt,pres,stress)
c        if (barostat .eq. 'MONTECARLO')  call pmonte (epot,temp)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine pressure2  --  barostat applied at half-step  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "pressure2" applies a box size and velocity correction at
c     the half time step as needed for the Monte Carlo barostat
c
c
      subroutine pressure2 (epot,temp)
      use bath
      use bound
      implicit none
      real*8 epot,temp
c
c
c     only necessary if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     use the desired barostat to maintain constant pressure
c
      if (isobaric) then
c        if (barostat .eq. 'BERENDSEN')  call pscale (dt,pres,stress)
c        if (barostat .eq. 'BUSSI')  call pscale (dt,pres,stress)
         if (barostat .eq. 'MONTECARLO')  call pmonte (epot,temp)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine pscale  --  Berendsen & Bussi scaling barostats  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "pscale" implements the Berendsen and Bussi barostats by scaling
c     the box dimensions, coordinates and momenta via coupling to an
c     external constant pressure bath
c
c     literature references:
c
c     H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren,
c     A. DiNola and J. R. Hauk, "Molecular Dynamics with Coupling
c     to an External Bath", Journal of Chemical Physics, 81,
c     3684-3690 (1984)
c
c     M. Bernetti and G. Bussi, "Pressure Control Using Stochastic
c     Cell Rescaling", Journal of Chemical Physics, 153, 114107 (2020)
c
c     V. Del Tatto, "A Fully Anisotropic Formulation of Stochastic
c     Cell Rescaling", arXiv, 2111.06403v1 (2021)
c
c     original code for anisotropic pressure coupling by Guido Raos,
c     Dipartimento di Chimica, Politecnico di Milano, Italy, May 2006
c
c
      subroutine pscale (dt,pres,stress)
      use atomid
      use atoms
      use bath
      use boxes
      use group
      use math
      use mdstuf
      use moldyn
      use rgddyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 dt,pres,weigh
      real*8 eps,deps,term
      real*8 kt,betat,dw
      real*8 scale,scalei
      real*8 scalexy,scalez
      real*8 normal,cosine
      real*8 tension
      real*8 xcm,xmove
      real*8 ycm,ymove
      real*8 zcm,zmove
      real*8 stress(3,3)
      real*8 temp(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
      external normal
c
c
c     find the isotropic scale factor for pressure control
c
      if (prestyp .eq. 'ISOTROPIC') then
         if (barostat .eq. 'BERENDSEN') then
            eps = third * (compress*dt/taupres)
            scale = 1.0d0 + eps*(pres-atmsph)
         else if (barostat .eq. 'BUSSI') then
            kt = gasconst * kelvin
            betat = prescon * compress
            dw = normal ()
            eps = (compress*dt/taupres) * (pres-atmsph)
            deps = sqrt(2.0d0*kt*betat*dt/(volbox*taupres))
            scale = exp(third*(eps+deps*dw))
         end if
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
            do i = 1, nuse
               k = iuse(i)
               x(k) = x(k) * scale
               y(k) = y(k) * scale
               z(k) = z(k) * scale
            end do
            if (barostat .eq. 'BUSSI') then
               do i = 1, nuse
                  k = iuse(i)
                  do j = 1, 3
                     v(j,k) = v(j,k) / scale
                  end do
               end do
            end if
c
c     couple to pressure bath via center of mass of rigid bodies
c
         else
            scalei = scale - 1.0d0
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
               xmove = scalei * xcm/grpmass(i)
               ymove = scalei * ycm/grpmass(i)
               zmove = scalei * zcm/grpmass(i)
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
               if (barostat .eq. 'BUSSI') then
                  do j = 1, 3
                     vcm(j,i) = vcm(j,i) / scale
                     wcm(j,i) = wcm(j,i) / scale
                  end do
               end if
            end do
         end if
c
c     find the semi-isotropic scale factors for pressure control
c
      else if (prestyp .eq. 'SEMIISO') then
         if (barostat .eq. 'BERENDSEN') then
            tension = 0.0d0
            eps = third * (compress*dt/taupres)
            scalexy = 1.0d0 + eps*(0.5d0*(stress(1,1)+stress(2,2))
     &                                +(tension/zbox)-atmsph)
            scalez = 1.0d0 + eps*(stress(3,3)-atmsph)
         else if (barostat .eq. 'BUSSI') then
            tension = 0.0d0
            kt = gasconst * kelvin
            betat = prescon * compress
            eps = third * (compress*dt/taupres)
            deps = sqrt(third2*kt*betat*dt/(volbox*taupres))
            dw = normal ()
            term = 0.5d0*(stress(1,1)+stress(2,2))
     &                + (tension/zbox) - atmsph
            scalexy = 1.0d0 + eps*term + root2*deps*dw
            dw = normal ()
            term = stress(3,3) - atmsph
            scalez = 1.0d0 + eps*term + deps*dw
         end if
c
c     modify the current periodic box dimension values
c
         xbox = xbox * scalexy
         ybox = ybox * scalexy
         zbox = zbox * scalez
c
c     propagate the new box dimensions to other lattice values
c
         call lattice
c
c     couple to pressure bath via atom scaling in Cartesian space
c
         if (integrate .ne. 'RIGIDBODY') then
            do i = 1, nuse
               k = iuse(i)
               x(k) = x(k) * scalexy
               y(k) = y(k) * scalexy
               z(k) = z(k) * scalez
            end do
            if (barostat .eq. 'BUSSI') then
               do i = 1, nuse
                  k = iuse(i)
                  v(1,k) = v(1,k) / scalexy
                  v(2,k) = v(2,k) / scalexy
                  v(3,k) = v(3,k) / scalez
               end do
            end if
c
c     couple to pressure bath via center of mass of rigid bodies
c
         else
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
               xmove = (scalexy-1.0d0) * xcm/grpmass(i)
               ymove = (scalexy-1.0d0) * ycm/grpmass(i)
               zmove = (scalez-1.0d0) * zcm/grpmass(i)
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
               if (barostat .eq. 'BUSSI') then
                  vcm(j,i) = vcm(1,i) / scalexy
                  vcm(2,i) = vcm(2,i) / scalexy
                  vcm(3,i) = vcm(3,i) / scalez
                  wcm(1,i) = wcm(1,i) / scalexy
                  wcm(2,i) = wcm(2,i) / scalexy
                  wcm(3,i) = wcm(3,i) / scalez
               end if
            end do
         end if
c
c     find the anisotropic scale factors for pressure control
c
      else if (prestyp .eq. 'ANISO') then
         if (barostat .eq. 'BERENDSEN') then
            eps = third * (compress*dt/taupres)
            do i = 1, 3
               do j = 1, 3
                  if (j. eq. i) then
                     ascale(j,i) = 1.0d0 + eps*(stress(j,i)-atmsph)
                  else
                     ascale(j,i) = eps * stress(j,i)
                  end if
               end do
            end do
         else if (barostat .eq. 'BUSSI') then
            kt = gasconst * kelvin
            betat = prescon * compress
            eps = third * (compress*dt/taupres)
            deps = sqrt(third2*kt*betat*dt/(volbox*taupres))
            do i = 1, 3
               do j = 1, 3
                  dw = normal ()
                  if (j .eq. i) then
                     term = stress(j,i) - atmsph + prescon*kt/volbox
                     ascale(j,i) = 1.0d0 + eps*term + deps*dw
                  else
c                    ascale(j,i) = eps*stress(j,i) + deps*dw
                     ascale(j,i) = eps * stress(j,i)
                  end if
               end do
            end do
         end if
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
            do i = 1, nuse
               k = iuse(i)
               x(k) = x(k)*ascale(1,1) + y(k)*ascale(1,2)
     &                   + z(k)*ascale(1,3)
               y(k) = x(k)*ascale(2,1) + y(k)*ascale(2,2)
     &                   + z(k)*ascale(2,3)
               z(k) = x(k)*ascale(3,1) + y(k)*ascale(3,2)
     &                   + z(k)*ascale(3,3)
            end do
            if (barostat .eq. 'BUSSI') then
               call invert (3,ascale)
               do i = 1, nuse
                  k = iuse(i)
                  v(1,k) = v(1,k)*ascale(1,1) + v(2,k)*ascale(1,2)
     &                        + v(3,k)*ascale(1,3)
                  v(2,k) = v(1,k)*ascale(2,1) + v(2,k)*ascale(2,2)
     &                        + v(3,k)*ascale(2,3)
                  v(3,k) = v(1,k)*ascale(3,1) + v(2,k)*ascale(3,2)
     &                        + v(3,k)*ascale(3,3)
               end do
            end if
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
               xmove = xcm*ascale(1,1) + ycm*ascale(1,2)
     &                    + zcm*ascale(1,3)
               ymove = xcm*ascale(2,1) + ycm*ascale(2,2)
     &                    + zcm*ascale(2,3)
               zmove = xcm*ascale(3,1) + ycm*ascale(3,2)
     &                    + zcm*ascale(3,3)
               do j = start, stop
                  k = kgrp(j)
                  x(k) = x(k) + xmove
                  y(k) = y(k) + ymove
                  z(k) = z(k) + zmove
               end do
               if (barostat .eq. 'BUSSI') then
                  call invert (3,ascale)
                  vcm(1,i) = vcm(1,i)*ascale(1,1) + vcm(2,i)*ascale(1,2)
     &                          + vcm(3,i)*ascale(1,3)
                  vcm(2,i) = vcm(1,i)*ascale(2,1) + vcm(2,i)*ascale(2,2)
     &                          + vcm(3,i)*ascale(2,3)
                  vcm(3,i) = vcm(1,i)*ascale(3,1) + vcm(2,i)*ascale(3,2)
     &                          + vcm(3,i)*ascale(3,3)
                  wcm(1,i) = wcm(1,i)*ascale(1,1) + wcm(2,i)*ascale(1,2)
     &                          + wcm(3,i)*ascale(1,3)
                  wcm(2,i) = wcm(1,i)*ascale(2,1) + wcm(2,i)*ascale(2,2)
     &                          + wcm(3,i)*ascale(2,3)
                  wcm(3,i) = wcm(1,i)*ascale(3,1) + wcm(2,i)*ascale(3,2)
     &                          + wcm(3,i)*ascale(3,3)
               end if
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
c     changes in the box dimensions and coordinates
c
c     literature references:
c
c     J. Aqvist, P. Wennerstrom, M. Nervall, S. Bjelic, B. O. Brandsal,
c     "Molecular Dynamics Simulations of Water and Biomolecules with
c     a Monte Carlo Constant Pressure Algorithm", Chemical Physics
c     Letters, 384, 288-294 (2004)
c
c     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
c     3rd Edition", Academic Press, San Diego, CA, 2023; Section 6.3
c
c     original version implemented by Alan Grossfield, January 2004;
c     anisotropic modification provided by Lee-Ping Wang, Stanford
c     University, March 2013
c
c
      subroutine pmonte (epot,temp)
      use atomid
      use atoms
      use bath
      use boxes
      use group
      use math
      use mdstuf
      use molcul
      use moldyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer start,stop
      real*8 epot,temp,term
      real*8 energy,random
      real*8 expterm,weigh
      real*8 kt,step,scale
      real*8 eold,rnd6
      real*8 xcm,ycm,zcm
      real*8 volold,cosine
      real*8 dpot,dpv,dkin
      real*8 xmove,ymove,zmove
      real*8 xboxold,yboxold,zboxold
      real*8 alphaold,betaold,gammaold
      real*8 temp3(3,3)
      real*8 hbox(3,3)
      real*8 ascale(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      logical dotrial
      logical isotropic
      logical idealgas
      external random
c
c
c     decide whether to attempt a box size change at this step
c
      dotrial = .false.
      if (random() .lt. 1.0d0/dble(voltrial))  dotrial = .true.
c
c     set constants and decide on type of trial box size change
c
      if (dotrial) then
         kt = gasconst * temp
         if (isothermal)  kt = gasconst * kelvin
         isotropic = .true.
         if (prestyp.eq.'ANISO' .and. random().gt.0.5d0) then
            isotropic = .false.
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (xold(n))
         allocate (yold(n))
         allocate (zold(n))
c
c     save the system state prior to trial box size change
c
         xboxold = xbox
         yboxold = ybox
         zboxold = zbox
         alphaold = alpha
         betaold = beta
         gammaold = gamma
         volold = volbox
         eold = epot
         do i = 1, n
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
         end do
c
c     for the isotropic case, change the lattice lengths uniformly
c
         if (isotropic) then
            step = volmove * (2.0d0*random()-1.0d0)
            volbox = volbox + step
            scale = (volbox/volold)**third
            xbox = xbox * scale
            ybox = ybox * scale
            zbox = zbox * scale
            call lattice
            if (integrate .eq. 'RIGIDBODY') then
               scale = scale - 1.0d0
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
            else if (volscale .eq. 'MOLECULAR') then
               scale = scale - 1.0d0
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
                  xmove = scale * xcm/molmass(i)
                  ymove = scale * ycm/molmass(i)
                  zmove = scale * zcm/molmass(i)
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
               do i = 1, nuse
                  k = iuse(i)
                  x(k) = x(k) * scale
                  y(k) = y(k) * scale
                  z(k) = z(k) * scale
               end do
            end if
c
c     for anisotropic case alter lattice angles, then scale lengths
c
         else
            rnd6 = 6.0d0*random()
            step = volmove * (2.0d0*random()-1.0d0)
            scale = (1.0d0+step/volold)**third
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
c     modify the current periodic box lattice angle values
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
            scale = (volbox/volold)**third
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
                     k = kmol(j)
                     weigh = mass(k)
                     xcm = xcm + x(k)*weigh
                     ycm = ycm + y(k)*weigh
                     zcm = zcm + z(k)*weigh
                  end do
                  xcm = xcm / grpmass(i)
                  ycm = ycm / grpmass(i)
                  zcm = zcm / grpmass(i)
                  xmove = xcm*ascale(1,1) + ycm*ascale(1,2)
     &                       + zcm*ascale(1,3)
                  ymove = xcm*ascale(2,1) + ycm*ascale(2,2)
     &                       + zcm*ascale(2,3)
                  zmove = xcm*ascale(3,1) + ycm*ascale(3,2)
     &                       + zcm*ascale(3,3)
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
                  xmove = xcm*ascale(1,1) + ycm*ascale(1,2)
     &                       + zcm*ascale(1,3)
                  ymove = xcm*ascale(2,1) + ycm*ascale(2,2)
     &                       + zcm*ascale(2,3)
                  zmove = xcm*ascale(3,1) + ycm*ascale(3,2)
     &                       + zcm*ascale(3,3)
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
               do i = 1, nuse
                  k = iuse(i)
                  x(k) = x(k)*ascale(1,1) + y(k)*ascale(1,2)
     &                      + z(k)*ascale(1,3)
                  y(k) = x(k)*ascale(2,1) + y(k)*ascale(2,2)
     &                      + z(k)*ascale(2,3)
                  z(k) = x(k)*ascale(3,1) + y(k)*ascale(3,2)
     &                      + z(k)*ascale(3,3)
               end do
            end if
         end if
c
c     get the potential energy and PV work changes for trial move
c
         epot = energy ()
         dpot = epot - eold
         dpv = atmsph * (volbox-volold) / prescon
c
c     get the kinetic energy contribution for the trial move
c
         idealgas = .true.
c
c     estimate the kinetic energy change as an ideal gas term
c
         if (idealgas) then
            if (integrate .eq. 'RIGIDBODY') then
               dkin = dble(ngrp) * kt * log(volold/volbox)
            else if (volscale .eq. 'MOLECULAR') then
               dkin = dble(nmol) * kt * log(volold/volbox)
            else
               dkin = dble(nmol) * kt * log(volold/volbox)
c              dkin = dble(nuse) * kt * log(volold/volbox)
            end if
c
c     alternatively get the instantaneous kinetic energy change;
c     requires the prior step velocity, which is not available
c
         else
            dkin = 0.0d0
            do i = 1, nuse
               k = iuse(i)
               term = 1.5d0 * mass(k) / ekcal
               do j = 1, 3
c                 dkin = dkin + term*(v(j,k)**2-vold(j,k)**2)
               end do
            end do
            if (integrate .eq. 'RIGIDBODY') then
               dkin = dkin * dble(ngrp)/dble(nuse)
            else if (volscale .eq. 'MOLECULAR') then
               dkin = dkin * dble(nmol)/dble(nuse)
            else
               dkin = dkin * dble(nuse)/dble(nuse)
            end if
         end if
c
c     acceptance ratio from Epot change, Ekin change and PV work
c
         term = -(dpot+dpv+dkin) / kt
         expterm = exp(term)
c
c     reject the step, and restore values prior to trial change
c
         if (random() .gt. expterm) then
            epot = eold
            xbox = xboxold
            ybox = yboxold
            zbox = zboxold
            alpha = alphaold
            beta = betaold
            gamma = gammaold
            call lattice
            do i = 1, n
               x(i) = xold(i)
               y(i) = yold(i)
               z(i) = zold(i)
            end do
         end if
c
c     perform deallocation of some local arrays
c
         deallocate (xold)
         deallocate (yold)
         deallocate (zold)
      end if
      return
      end
