c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine erepel2  --  atomwise Pauli repulsion Hessian  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erepel2" calculates the second derivatives of the Pauli
c     repulsion energy
c
c
      subroutine erepel2 (i)
      use atoms
      use deriv
      use hessn
      use mpole
      use potent
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical twosided
c
c
c     set the default stepsize and accuracy control flags
c
      eps = 1.0d-5
      twosided = .false.
      if (n .le. 300)  twosided = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(npole))
      allocate (d0(3,n))
c
c     perform dynamic allocation of some global arrays
c
      prior = .false.
      if (allocated(der)) then
         prior = .true.
         if (size(der) .lt. 3*n)  deallocate (der)
      end if
      if (.not. allocated(der))  allocate (der(3,n))
c
c     find the multipole definitions involving the current atom;
c     results in a faster but approximate Hessian calculation
c
      nlist = 0
      do k = 1, npole
         if (ipole(k).eq.i .or. zaxis(k).eq.i
     &          .or. xaxis(k).eq.i .or. yaxis(k).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     get multipole first derivatives for the base structure
c
      if (.not. twosided) then
         call erepel2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = der(j,k)
            end do
         end do
      end if
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      if (twosided) then
         x(i) = x(i) - 0.5d0*eps
         call erepel2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = der(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      call erepel2a (nlist,list)
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (der(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      if (twosided) then
         y(i) = y(i) - 0.5d0*eps
         call erepel2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = der(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      call erepel2a (nlist,list)
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (der(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      if (twosided) then
         z(i) = z(i) - 0.5d0*eps
         call erepel2a (nlist,list)
         do k = 1, n
            do j = 1, 3
               d0(j,k) = der(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      call erepel2a (nlist,list)
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (der(j,k)-d0(j,k))/eps
         end do
      end do
c
c     perform deallocation of some global arrays
c
      if (.not. prior)  deallocate (der)
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (d0)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine erepel2a  --  numerical repulsion Hessian  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "erepel2a" computes Pauli repulsion first derivatives for a
c     single atom via a double loop; used to get finite difference
c     second derivatives
c
c
      subroutine erepel2a (nlist,list)
      use atoms
      use bound
      use cell
      use couple
      use deriv
      use group
      use mpole
      use potent
      use repel
      use reppot
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,iii
      integer nlist,jcell
      integer list(*)
      real*8 e,fgrp
      real*8 eterm,de
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,r3,r4,r5
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr,dik,qik
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 diqk,dkqi,qiqk
      real*8 dirx,diry,dirz
      real*8 dkrx,dkry,dkrz
      real*8 dikx,diky,dikz
      real*8 qirx,qiry,qirz
      real*8 qkrx,qkry,qkrz
      real*8 qikx,qiky,qikz
      real*8 qixk,qiyk,qizk
      real*8 qkxi,qkyi,qkzi
      real*8 qikrx,qikry,qikrz
      real*8 qkirx,qkiry,qkirz
      real*8 diqkx,diqky,diqkz
      real*8 dkqix,dkqiy,dkqiz
      real*8 diqkrx,diqkry,diqkrz
      real*8 dkqirx,dkqiry,dkqirz
      real*8 dqikx,dqiky,dqikz
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 sizi,sizk,sizik
      real*8 vali,valk
      real*8 dmpi,dmpk
      real*8 frcx,frcy,frcz
      real*8 taper,dtaper
      real*8 ttri(3),ttrk(3)
      real*8 fix(3),fiy(3),fiz(3)
      real*8 dmpik(11)
      real*8, allocatable :: rscale(:)
      real*8, allocatable :: ter(:,:)
      logical proceed,usei
      character*6 mode
c
c
c     zero out the Pauli repulsion first derivative components
c
      do i = 1, n
         der(1,i) = 0.0d0
         der(2,i) = 0.0d0
         der(3,i) = 0.0d0
      end do
      if (nrep .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     perform dynamic allocation of some local arrays
c
      allocate (rscale(n))
      allocate (ter(3,n))
c
c     initialize connected atom scaling and torque arrays
c
      do i = 1, n
         rscale(i) = 1.0d0
         do j = 1, 3
            ter(j,i) = 0.0d0
         end do
      end do
c
c     set cutoff distances and switching coefficients
c
      mode = 'REPULS'
      call switch (mode)
c
c     compute the Pauli repulsion first derivatives
c
      do iii = 1, nlist
         ii = list(iii)
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         sizi = sizpr(ii)
         dmpi = dmppr(ii)
         vali = elepr(ii)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = r2scale
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = r3scale
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = r4scale
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = r5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            k = ipole(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (ii .eq. kk)  proceed = .false.
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  sizk = sizpr(kk)
                  dmpk = dmppr(kk)
                  valk = elepr(kk)
                  ck = rpole(1,kk)
                  dkx = rpole(2,kk)
                  dky = rpole(3,kk)
                  dkz = rpole(4,kk)
                  qkxx = rpole(5,kk)
                  qkxy = rpole(6,kk)
                  qkxz = rpole(7,kk)
                  qkyy = rpole(9,kk)
                  qkyz = rpole(10,kk)
                  qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  dik = dix*dkx + diy*dky + diz*dkz
                  qik = qix*qkx + qiy*qky + qiz*qkz
                  diqk = dix*qkx + diy*qky + diz*qkz
                  dkqi = dkx*qix + dky*qiy + dkz*qiz
                  qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                      + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     additional intermediates involving moments and distance
c
                  dirx = diy*zr - diz*yr
                  diry = diz*xr - dix*zr
                  dirz = dix*yr - diy*xr
                  dkrx = dky*zr - dkz*yr
                  dkry = dkz*xr - dkx*zr
                  dkrz = dkx*yr - dky*xr
                  dikx = diy*dkz - diz*dky
                  diky = diz*dkx - dix*dkz
                  dikz = dix*dky - diy*dkx
                  qirx = qiz*yr - qiy*zr
                  qiry = qix*zr - qiz*xr
                  qirz = qiy*xr - qix*yr
                  qkrx = qkz*yr - qky*zr
                  qkry = qkx*zr - qkz*xr
                  qkrz = qky*xr - qkx*yr
                  qikx = qky*qiz - qkz*qiy
                  qiky = qkz*qix - qkx*qiz
                  qikz = qkx*qiy - qky*qix
                  qixk = qixx*qkx + qixy*qky + qixz*qkz
                  qiyk = qixy*qkx + qiyy*qky + qiyz*qkz
                  qizk = qixz*qkx + qiyz*qky + qizz*qkz
                  qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz
                  qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz
                  qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz
                  qikrx = qizk*yr - qiyk*zr
                  qikry = qixk*zr - qizk*xr
                  qikrz = qiyk*xr - qixk*yr
                  qkirx = qkzi*yr - qkyi*zr
                  qkiry = qkxi*zr - qkzi*xr
                  qkirz = qkyi*xr - qkxi*yr
                  diqkx = dix*qkxx + diy*qkxy + diz*qkxz
                  diqky = dix*qkxy + diy*qkyy + diz*qkyz
                  diqkz = dix*qkxz + diy*qkyz + diz*qkzz
                  dkqix = dkx*qixx + dky*qixy + dkz*qixz
                  dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
                  dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
                  diqkrx = diqkz*yr - diqky*zr
                  diqkry = diqkx*zr - diqkz*xr
                  diqkrz = diqky*xr - diqkx*yr
                  dkqirx = dkqiz*yr - dkqiy*zr
                  dkqiry = dkqix*zr - dkqiz*xr
                  dkqirz = dkqiy*xr - dkqix*yr
                  dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy
     &                    - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                            -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
                  dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz
     &                    - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                            -qixx*qkxz-qixy*qkyz-qixz*qkzz)
                  dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix
     &                    - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                            -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     get reciprocal distance terms for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  rr11 = 9.0d0 * rr9 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c
                  call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                             11,dmpi,dmpk,dmpik)                  
c
c     calculate intermediate terms needed for the energy
c
                  term1 = vali*valk
                  term2 = valk*dir - vali*dkr + dik
                  term3 = vali*qkr + valk*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  eterm = term1*dmpik(1) + term2*dmpik(3)
     &                       + term3*dmpik(5) + term4*dmpik(7)
     &                       + term5*dmpik(9)
c
c     compute the Pauli repulsion energy for this interaction
c
                  sizik = sizi * sizk * rscale(k)
                  e = sizik * eterm * rr1
c
c     calculate intermediate terms for force and torque
c
                  de = term1*dmpik(3) + term2*dmpik(5)
     &                    + term3*dmpik(7) + term4*dmpik(9)
     &                    + term5*dmpik(11)
                  term1 = -valk*dmpik(3) + dkr*dmpik(5)
     &                       - qkr*dmpik(7)
                  term2 = vali*dmpik(3) + dir*dmpik(5)
     &                       + qir*dmpik(7)
                  term3 = 2.0d0 * dmpik(5)
                  term4 = 2.0d0 * (-valk*dmpik(5) + dkr*dmpik(7)
     &                                - qkr*dmpik(9))
                  term5 = 2.0d0 * (-vali*dmpik(5) - dir*dmpik(7)
     &                                - qir*dmpik(9))
                  term6 = 4.0d0 * dmpik(7)
c     
c     compute the force components for this interaction
c     
                  frcx = de*xr + term1*dix + term2*dkx
     &                      + term3*(diqkx-dkqix) + term4*qix
     &                      + term5*qkx + term6*(qixk+qkxi)
                  frcy = de*yr + term1*diy + term2*dky
     &                      + term3*(diqky-dkqiy) + term4*qiy
     &                      + term5*qky + term6*(qiyk+qkyi)
                  frcz = de*zr + term1*diz + term2*dkz
     &                      + term3*(diqkz-dkqiz) + term4*qiz
     &                      + term5*qkz + term6*(qizk+qkzi)
                  frcx = frcx*rr1 + eterm*rr3*xr
                  frcy = frcy*rr1 + eterm*rr3*yr
                  frcz = frcz*rr1 + eterm*rr3*zr
                  frcx = sizik * frcx
                  frcy = sizik * frcy
                  frcz = sizik * frcz
c
c     compute the torque components for this interaction
c
                  ttri(1) = -dmpik(3)*dikx + term1*dirx
     &                         + term3*(dqikx+dkqirx)
     &                         - term4*qirx - term6*(qikrx+qikx)
                  ttri(2) = -dmpik(3)*diky + term1*diry
     &                         + term3*(dqiky+dkqiry)
     &                         - term4*qiry - term6*(qikry+qiky)
                  ttri(3) = -dmpik(3)*dikz + term1*dirz
     &                         + term3*(dqikz+dkqirz)
     &                         - term4*qirz - term6*(qikrz+qikz)
                  ttrk(1) = dmpik(3)*dikx + term2*dkrx
     &                         - term3*(dqikx+diqkrx)
     &                         - term5*qkrx - term6*(qkirx-qikx)
                  ttrk(2) = dmpik(3)*diky + term2*dkry
     &                         - term3*(dqiky+diqkry)
     &                         - term5*qkry - term6*(qkiry-qiky)
                  ttrk(3) = dmpik(3)*dikz + term2*dkrz
     &                         - term3*(dqikz+diqkrz)
     &                         - term5*qkrz - term6*(qkirz-qikz)
                  ttri(1) = sizik * ttri(1) * rr1
                  ttri(2) = sizik * ttri(2) * rr1
                  ttri(3) = sizik * ttri(3) * rr1
                  ttrk(1) = sizik * ttrk(1) * rr1
                  ttrk(2) = sizik * ttrk(2) * rr1
                  ttrk(3) = sizik * ttrk(3) * rr1
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = fgrp * e
                     frcx = fgrp * frcx
                     frcy = fgrp * frcy
                     frcz = fgrp * frcz
                     do j = 1, 3
                        ttri(j) = fgrp * ttri(j)
                        ttrk(j) = fgrp * ttrk(j)
                     end do
                  end if
c
c     use energy switching if near the cutoff distance
c
                  if (r2 .gt. cut2) then
                     r3 = r2 * r
                     r4 = r2 * r2
                     r5 = r2 * r3
                     taper = c5*r5 + c4*r4 + c3*r3
     &                          + c2*r2 + c1*r + c0
                     dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                           + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                     dtaper = dtaper * e * rr1
                     frcx = frcx*taper - dtaper*xr
                     frcy = frcy*taper - dtaper*yr
                     frcz = frcz*taper - dtaper*zr
                     do j = 1, 3
                        ttri(j) = ttri(j) * taper
                        ttrk(j) = ttrk(j) * taper
                     end do
                  end if
c
c     increment force-based gradient and torque on first site
c
                  der(1,i) = der(1,i) + frcx
                  der(2,i) = der(2,i) + frcy
                  der(3,i) = der(3,i) + frcz
                  ter(1,i) = ter(1,i) + ttri(1)
                  ter(2,i) = ter(2,i) + ttri(2)
                  ter(3,i) = ter(3,i) + ttri(3)
c
c     increment force-based gradient and torque on second site
c
                  der(1,k) = der(1,k) - frcx
                  der(2,k) = der(2,k) - frcy
                  der(3,k) = der(3,k) - frcz
                  ter(1,k) = ter(1,k) + ttrk(1)
                  ter(2,k) = ter(2,k) + ttrk(2)
                  ter(3,k) = ter(3,k) + ttrk(3)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            rscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            rscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            rscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            rscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction components with other unit cells
c
         do iii = 1, nlist
            ii = list(iii)
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            sizi = sizpr(ii)
            dmpi = dmppr(ii)
            vali = elepr(ii)
            ci = rpole(1,ii)
            dix = rpole(2,ii)
            diy = rpole(3,ii)
            diz = rpole(4,ii)
            qixx = rpole(5,ii)
            qixy = rpole(6,ii)
            qixz = rpole(7,ii)
            qiyy = rpole(9,ii)
            qiyz = rpole(10,ii)
            qizz = rpole(13,ii)
            usei = use(i)
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               rscale(i12(j,i)) = r2scale
            end do
            do j = 1, n13(i)
               rscale(i13(j,i)) = r3scale
            end do
            do j = 1, n14(i)
               rscale(i14(j,i)) = r4scale
            end do
            do j = 1, n15(i)
               rscale(i15(j,i)) = r5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = 1, npole
               k = ipole(kk)
               proceed = .true.
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               if (.not. use_intra)  proceed = .true.
               if (proceed)  proceed = (usei .or. use(k))
               if (proceed) then
                  do jcell = 2, ncell
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     call imager (xr,yr,zr,jcell)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (r2 .le. off2) then
                        r = sqrt(r2)
                        sizk = sizpr(kk)
                        dmpk = dmppr(kk)
                        valk = elepr(kk)
                        ck = rpole(1,kk)
                        dkx = rpole(2,kk)
                        dky = rpole(3,kk)
                        dkz = rpole(4,kk)
                        qkxx = rpole(5,kk)
                        qkxy = rpole(6,kk)
                        qkxz = rpole(7,kk)
                        qkyy = rpole(9,kk)
                        qkyz = rpole(10,kk)
                        qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
                        dir = dix*xr + diy*yr + diz*zr
                        qix = qixx*xr + qixy*yr + qixz*zr
                        qiy = qixy*xr + qiyy*yr + qiyz*zr
                        qiz = qixz*xr + qiyz*yr + qizz*zr
                        qir = qix*xr + qiy*yr + qiz*zr
                        dkr = dkx*xr + dky*yr + dkz*zr
                        qkx = qkxx*xr + qkxy*yr + qkxz*zr
                        qky = qkxy*xr + qkyy*yr + qkyz*zr
                        qkz = qkxz*xr + qkyz*yr + qkzz*zr
                        qkr = qkx*xr + qky*yr + qkz*zr
                        dik = dix*dkx + diy*dky + diz*dkz
                        qik = qix*qkx + qiy*qky + qiz*qkz
                        diqk = dix*qkx + diy*qky + diz*qkz
                        dkqi = dkx*qix + dky*qiy + dkz*qiz
                        qiqk = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                            + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     additional intermediates involving moments and distance
c
                        dirx = diy*zr - diz*yr
                        diry = diz*xr - dix*zr
                        dirz = dix*yr - diy*xr
                        dkrx = dky*zr - dkz*yr
                        dkry = dkz*xr - dkx*zr
                        dkrz = dkx*yr - dky*xr
                        dikx = diy*dkz - diz*dky
                        diky = diz*dkx - dix*dkz
                        dikz = dix*dky - diy*dkx
                        qirx = qiz*yr - qiy*zr
                        qiry = qix*zr - qiz*xr
                        qirz = qiy*xr - qix*yr
                        qkrx = qkz*yr - qky*zr
                        qkry = qkx*zr - qkz*xr
                        qkrz = qky*xr - qkx*yr
                        qikx = qky*qiz - qkz*qiy
                        qiky = qkz*qix - qkx*qiz
                        qikz = qkx*qiy - qky*qix
                        qixk = qixx*qkx + qixy*qky + qixz*qkz
                        qiyk = qixy*qkx + qiyy*qky + qiyz*qkz
                        qizk = qixz*qkx + qiyz*qky + qizz*qkz
                        qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz
                        qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz
                        qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz
                        qikrx = qizk*yr - qiyk*zr
                        qikry = qixk*zr - qizk*xr
                        qikrz = qiyk*xr - qixk*yr
                        qkirx = qkzi*yr - qkyi*zr
                        qkiry = qkxi*zr - qkzi*xr
                        qkirz = qkyi*xr - qkxi*yr
                        diqkx = dix*qkxx + diy*qkxy + diz*qkxz
                        diqky = dix*qkxy + diy*qkyy + diz*qkyz
                        diqkz = dix*qkxz + diy*qkyz + diz*qkzz
                        dkqix = dkx*qixx + dky*qixy + dkz*qixz
                        dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz
                        dkqiz = dkx*qixz + dky*qiyz + dkz*qizz
                        diqkrx = diqkz*yr - diqky*zr
                        diqkry = diqkx*zr - diqkz*xr
                        diqkrz = diqky*xr - diqkx*yr
                        dkqirx = dkqiz*yr - dkqiy*zr
                        dkqiry = dkqix*zr - dkqiz*xr
                        dkqirz = dkqiy*xr - dkqix*yr
                        dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy
     &                          - 2.0d0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz
     &                                  -qixz*qkxy-qiyz*qkyy-qizz*qkyz)
                        dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz
     &                          - 2.0d0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz
     &                                  -qixx*qkxz-qixy*qkyz-qixz*qkzz)
                        dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix
     &                          - 2.0d0*(qixx*qkxy+qixy*qkyy+qixz*qkyz
     &                                  -qixy*qkxx-qiyy*qkxy-qiyz*qkxz)
c
c     get reciprocal distance terms for this interaction
c
                        rr1 = 1.0d0 / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
                        rr11 = 9.0d0 * rr9 / r2
c
c     get damping coefficients for the Pauli repulsion energy
c
                        call damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                                   11,dmpi,dmpk,dmpik)                  
c
c     compute the Pauli repulsion energy for this interaction
c
                        term1 = vali*valk
                        term2 = valk*dir - vali*dkr + dik
                        term3 = vali*qkr + valk*qir - dir*dkr
     &                             + 2.0d0*(dkqi-diqk+qiqk)
                        term4 = dir*qkr - dkr*qir - 4.0d0*qik
                        term5 = qir*qkr
                        eterm = term1*dmpik(1) + term2*dmpik(3)
     &                             + term3*dmpik(5) + term4*dmpik(7)
     &                             + term5*dmpik(9)
c
c     compute the Pauli repulsion energy for this interaction
c
                        sizik = sizi * sizk
                        e = sizik * rscale(k) * eterm * rr1
c
c     calculate intermediate terms for force and torque
c
                        de = term1*dmpik(3) + term2*dmpik(5)
     &                          + term3*dmpik(7) + term4*dmpik(9)
     &                          + term5*dmpik(11)
                        term1 = -valk*dmpik(3) + dkr*dmpik(5)
     &                             - qkr*dmpik(7)
                        term2 = vali*dmpik(3) + dir*dmpik(5)
     &                             + qir*dmpik(7)
                        term3 = 2.0d0 * dmpik(5)
                        term4 = 2.0d0 * (-valk*dmpik(5) + dkr*dmpik(7)
     &                                      - qkr*dmpik(9))
                        term5 = 2.0d0 * (-vali*dmpik(5) - dir*dmpik(7)
     &                                      - qir*dmpik(9))
                        term6 = 4.0d0 * dmpik(7)
c     
c     compute the force components for this interaction
c     
                        frcx = de*xr + term1*dix + term2*dkx
     &                            + term3*(diqkx-dkqix) + term4*qix
     &                            + term5*qkx + term6*(qixk+qkxi)
                        frcy = de*yr + term1*diy + term2*dky
     &                            + term3*(diqky-dkqiy) + term4*qiy
     &                            + term5*qky + term6*(qiyk+qkyi)
                        frcz = de*zr + term1*diz + term2*dkz
     &                            + term3*(diqkz-dkqiz) + term4*qiz
     &                            + term5*qkz + term6*(qizk+qkzi)
                        frcx = frcx*rr1 + eterm*rr3*xr
                        frcy = frcy*rr1 + eterm*rr3*yr
                        frcz = frcz*rr1 + eterm*rr3*zr
                        frcx = sizik * frcx
                        frcy = sizik * frcy
                        frcz = sizik * frcz
c
c     compute the torque components for this interaction
c
                        ttri(1) = -dmpik(3)*dikx + term1*dirx
     &                               + term3*(dqikx+dkqirx)
     &                               - term4*qirx - term6*(qikrx+qikx)
                        ttri(2) = -dmpik(3)*diky + term1*diry
     &                               + term3*(dqiky+dkqiry)
     &                               - term4*qiry - term6*(qikry+qiky)
                        ttri(3) = -dmpik(3)*dikz + term1*dirz
     &                               + term3*(dqikz+dkqirz)
     &                               - term4*qirz - term6*(qikrz+qikz)
                        ttrk(1) = dmpik(3)*dikx + term2*dkrx
     &                               - term3*(dqikx+diqkrx)
     &                               - term5*qkrx - term6*(qkirx-qikx)
                        ttrk(2) = dmpik(3)*diky + term2*dkry
     &                               - term3*(dqiky+diqkry)
     &                               - term5*qkry - term6*(qkiry-qiky)
                        ttrk(3) = dmpik(3)*dikz + term2*dkrz
     &                               - term3*(dqikz+diqkrz)
     &                               - term5*qkrz - term6*(qkirz-qikz)
                        ttri(1) = sizik * ttri(1) * rr1
                        ttri(2) = sizik * ttri(2) * rr1
                        ttri(3) = sizik * ttri(3) * rr1
                        ttrk(1) = sizik * ttrk(1) * rr1
                        ttrk(2) = sizik * ttrk(2) * rr1
                        ttrk(3) = sizik * ttrk(3) * rr1
c
c     scale the interaction based on its group membership
c
                        if (use_group) then
                           e = fgrp * e
                           frcx = fgrp * frcx
                           frcy = fgrp * frcy
                           frcz = fgrp * frcz
                           do j = 1, 3
                              ttri(j) = fgrp * ttri(j)
                              ttrk(j) = fgrp * ttrk(j)
                           end do
                        end if
c
c     use energy switching if near the cutoff distance
c
                        if (r2 .gt. cut2) then
                           r3 = r2 * r
                           r4 = r2 * r2
                           r5 = r2 * r3
                           taper = c5*r5 + c4*r4 + c3*r3
     &                                + c2*r2 + c1*r + c0
                           dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                                 + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                           dtaper = dtaper * e * rr1
                           frcx = frcx*taper - dtaper*xr
                           frcy = frcy*taper - dtaper*yr
                           frcz = frcz*taper - dtaper*zr
                           do j = 1, 3
                              ttri(j) = ttri(j) * taper
                              ttrk(j) = ttrk(j) * taper
                           end do
                        end if
c
c     increment the overall Pauli repulsion energy component
c
                        if (i .eq. k) then
                           frcx = 0.5d0 * frcx
                           frcy = 0.5d0 * frcy
                           frcz = 0.5d0 * frcz
                           do j = 1, 3
                              ttri(j) = 0.5d0 * ttri(j)
                              ttrk(j) = 0.5d0 * ttrk(j)
                           end do
                        end if
c
c     increment force-based gradient and torque on first site
c
                        der(1,i) = der(1,i) + frcx
                        der(2,i) = der(2,i) + frcy
                        der(3,i) = der(3,i) + frcz
                        ter(1,i) = ter(1,i) + ttri(1)
                        ter(2,i) = ter(2,i) + ttri(2)
                        ter(3,i) = ter(3,i) + ttri(3)
c
c     increment force-based gradient and torque on second site
c
                        der(1,k) = der(1,k) - frcx
                        der(2,k) = der(2,k) - frcy
                        der(3,k) = der(3,k) - frcz
                        ter(1,k) = ter(1,k) + ttrk(1)
                        ter(2,k) = ter(2,k) + ttrk(2)
                        ter(3,k) = ter(3,k) + ttrk(3)
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               rscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               rscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               rscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               rscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     resolve site torques then increment forces and virial
c
      do ii = 1, npole
         i = ipole(ii)
         call torque (ii,ter(1,i),fix,fiy,fiz,der)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rscale)
      deallocate (ter)
      return
      end
