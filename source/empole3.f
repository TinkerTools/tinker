c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
      subroutine empole3
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'mpole.i'
      include 'potent.i'
      integer i,ii
c
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole3d
         else
            call empole3c
         end if
      else
         if (use_mlist) then
            call empole3b
         else
            call empole3a
         end if
      end if
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = 0.0d0
         end do
      end if
      if (.not. use_polar) then
         ep = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aep(ii) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole and dipole
c     polarizability interaction energy using a double loop,
c     and partitions the energy among the atoms
c
c
      subroutine empole3a
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,ei,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch ('MPOLE')
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do k = i+1, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                     end if
                  end if
                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                  fm = f * mscale(kk)
                  fp = f * pscale(kk)
                  e = fm * e
                  ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     e = e * fgrp
c                    ei = ei * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  muse = (use_mpole .and. mscale(kk).ne.0.0d0)
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
                  if (muse)  nem = nem + 1
                  if (puse)  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e + ei
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (muse .or. puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                9x,'Distance',6x,'Energies',
     &                                ' (MPole, Polar)',/)
                        end if
                        write (iout,20)  ii,name(ii),kk,name(kk),r,e,ei
   20                   format (' M-Pole',5x,i5,'-',a3,1x,i5,
     &                             '-',a3,5x,f8.4,4x,2f12.4)
                     end if
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(k)
         pdi = pdamp(i)
         pti = thole(i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do k = i, npole
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
                     ukx = uind(1,k)
                     uky = uind(2,k)
                     ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                     sc(2) = dix*dkx + diy*dky + diz*dkz
                     sc(3) = dix*xr + diy*yr + diz*zr
                     sc(4) = dkx*xr + dky*yr + dkz*zr
                     sc(5) = qix*xr + qiy*yr + qiz*zr
                     sc(6) = qkx*xr + qky*yr + qkz*zr
                     sc(7) = qix*dkx + qiy*dky + qiz*dkz
                     sc(8) = qkx*dix + qky*diy + qkz*diz
                     sc(9) = qix*qkx + qiy*qky + qiz*qkz
                     sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                     sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                           + diy*uky + uiz*dkz + diz*ukz
                     sci(3) = uix*xr + uiy*yr + uiz*zr
                     sci(4) = ukx*xr + uky*yr + ukz*zr
                     sci(7) = qix*ukx + qiy*uky + qiz*ukz
                     sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                     gl(0) = ci*ck
                     gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                     gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                          + 2.0d0*(sc(7)-sc(8)+sc(10))
                     gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                     gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                     gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                     gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                           - sc(3)*sci(4)
                     gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                     rr1 = 1.0d0 / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                           scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                             *expdamp
                        end if
                     end if
                     e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                     ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                       + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                     fm = f
                     fp = f
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           fm = fm * mscale(kk)
                           fp = fp * pscale(kk)
                        end if
                     end if
                     e = fm * e
                     ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                     if (use_group) then
                        e = e * fgrp
c                       ei = ei * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        ei = 0.5d0 * ei
                     end if
                     nem = nem + 1
                     nep = nep + 1
                     em = em + e
                     ep = ep + ei
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     aep(ii) = aep(ii) + 0.5d0*ei
                     aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                     einter = einter + e + ei
c
c     print a message if the energy of this interaction is large
c
                     huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                9x,'Distance',6x,'Energies',
     &                                ' (MPole, Polar)',/)
                        end if
                        write (iout,40)  ii,name(ii),kk,name(kk),r,e,ei
   40                   format (' M-Pole',5x,i5,'-',a3,1x,i5,'-',a3,
     &                             1x,'(X)',1x,f8.4,4x,2f12.4)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using a neighbor list,
c     and partitions the energy among the atoms
c
c
      subroutine empole3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,ei,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch ('MPOLE')
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                     end if
                  end if
                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                  fm = f * mscale(kk)
                  fp = f * pscale(kk)
                  e = fm * e
                  ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     e = e * fgrp
c                    ei = ei * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  muse = (use_mpole .and. mscale(kk).ne.0.0d0)
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
                  if (muse)  nem = nem + 1
                  if (puse)  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e + ei
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (muse .or. puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',13x,'Atom Names',
     &                                9x,'Distance',6x,'Energies',
     &                                ' (MPole, Polar)',/)
                        end if
                        write (iout,20)  ii,name(ii),kk,name(kk),r,e,ei
   20                   format (' M-Pole',5x,i5,'-',a3,1x,i5,
     &                             '-',a3,5x,f8.4,4x,2f12.4)
                     end if
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3c  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3c" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and double loop, and partitions the energy
c     among the atoms
c
c
      subroutine empole3c
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      integer i,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the real space part of the Ewald summation
c
      call ereal3c (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
         nem = nem + 1
         nep = nep + 1
         em = em + e
         ep = ep + ei
         aem(i) = aem(i) + e
         aep(i) = aep(i) + ei
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
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
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         nem = nem + 1
         nep = nep + 1
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3c  --  real space mpole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3c (eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      integer i,j,k,m
      integer ii,kk
      real*8 e,ei,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical header,huge
      logical muse,puse
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
      header = .true.
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = pscale(kk)
               scale5 = pscale(kk)
               scale7 = pscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     scale7 = scale7 * (1.0d0-(1.0d0-damp
     &                                    +0.6d0*damp**2)*expdamp)
                  end if
               end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3)
     &                    + gli(2)*rr5*(1.0d0-scale5)
     &                    + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
               e = e - efix*(1.0d0-mscale(kk))
               ei = ei - eifix
               e = f * e
               ei = 0.5d0 * f * ei
c
c     increment the overall multipole and polarization energies
c
               muse = use_mpole
               puse = use_polar
               if (muse)  nem = nem + 1
               if (puse)  nep = nep + 1
               em = em + e
               ep = ep + ei
               aem(ii) = aem(ii) + 0.5d0*e
               aem(kk) = aem(kk) + 0.5d0*e
               aep(ii) = aep(ii) + 0.5d0*ei
               aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intramolecular energy
c
               efix = f * efix * mscale(kk)
               eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + efix + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             9x,'Distance',6x,'Energies',
     &                             ' (MPole, Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   20                format (' M-Pole',5x,i5,'-',a3,1x,i5,
     &                          '-',a3,5x,f8.4,4x,2f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do k = i, npole
            kk = ipole(k)
            do j = 1, ncell
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call imager (xr,yr,zr,j)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     calculate the error function damping terms
c
                  ralpha = aewald * r
                  bn(0) = erfc(ralpha) / r
                  alsq2 = 2.0d0 * aewald**2
                  alsq2n = 0.0d0
                  if (aewald .gt. 0.0d0)
     &               alsq2n = 1.0d0 / (sqrtpi*aewald)
                  exp2a = exp(-ralpha**2)
                  do m = 1, 4
                     bfac = dble(m+m-1)
                     alsq2n = alsq2 * alsq2n
                     bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
                  end do
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                   + gl(3)*bn(3) + gl(4)*bn(4)
                  ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                        if (use_polymer .and. r2.le.polycut2) then
                           scale3 = scale3 * pscale(kk)
                           scale5 = scale5 * pscale(kk)
                           scale7 = scale7 * pscale(kk)
                        end if
                     end if
                  end if
                  efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                  eifix = gli(1)*rr3*(1.0d0-scale3)
     &                       + gli(2)*rr5*(1.0d0-scale5)
     &                       + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
                  if (use_polymer .and. r2.le.polycut2)
     &               e = e - efix*(1.0d0-mscale(kk))
                  ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
                  e = f * e
                  ei = 0.5d0 * f * ei
                  if (ii .eq. kk) then
                     e = 0.5d0 * e
                     ei = 0.5d0 * ei
                  end if
                  nem = nem + 1
                  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     print a message if the energy of this interaction is large
c
                  efix = f * efix * mscale(kk)
                  eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                       + gli(3)*rr7*scale7
                  eifix = 0.5d0 * f * eifix
                  huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             9x,'Distance',6x,'Energies',
     &                             ' (MPole, Polar)',/)
                     end if
                     write (iout,40)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   40                format (' M-Pole',5x,i5,'-',a3,1x,i5,'-',a3,
     &                          1x,'(X)',1x,f8.4,4x,2f12.4)
                  end if
               end if
            end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3d  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3d" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and a neighbor list, and partitions the
c     energy among the atoms
c
c
      subroutine empole3d
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      integer i,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip
c
c     compute the real space part of the Ewald summation
c
      call ereal3d (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
         nem = nem + 1
         nep = nep + 1
         em = em + e
         ep = ep + ei
         aem(i) = aem(i) + e
         aep(i) = aep(i) + ei
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
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
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         nem = nem + 1
         nep = nep + 1
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3d  --  real space mpole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c     using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3d (eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 e,ei,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8 mscale(maxatm)
      real*8 pscale(maxatm)
      logical header,huge
      logical muse,puse
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
      header = .true.
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      call switch ('EWALD')
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = pscale(kk)
               scale5 = pscale(kk)
               scale7 = pscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     scale7 = scale7 * (1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &                                               *expdamp)
                  end if
               end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3)
     &                    + gli(2)*rr5*(1.0d0-scale5)
     &                    + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
               e = e - efix*(1.0d0-mscale(kk))
               ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
               e = f * e
               ei = 0.5d0 * f * ei
               muse = use_mpole
               puse = use_polar
               if (muse)  nem = nem + 1
               if (puse)  nep = nep + 1
               em = em + e
               ep = ep + ei
               aem(ii) = aem(ii) + 0.5d0*e
               aem(kk) = aem(kk) + 0.5d0*e
               aep(ii) = aep(ii) + 0.5d0*ei
               aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intramolecular energy
c
               efix = f * efix * mscale(kk)
               eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + efix + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',13x,'Atom Names',
     &                             9x,'Distance',6x,'Energies',
     &                             ' (MPole, Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   20                format (' M-Pole',5x,i5,'-',a3,1x,i5,
     &                          '-',a3,5x,f8.4,4x,2f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
      return
      end
