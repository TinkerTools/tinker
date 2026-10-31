c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1  --  charge-charge energy & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine echarge1
      use dlmda
      use energi
      use extfld
      use limits
      use virial
      use warp
      implicit none
      integer i,j
      real*8 exf
      logical edrv
      logical lmdadrv
      character*6 mode
c
c
c     zero out the partial charge virial, which is kept apart from
c     the main virial for the lambda derivatives
c
      do i = 1, 3
         do j = 1, 3
            ecvir(j,i) = 0.0d0
         end do
      end do
c
c     choose the method for summing over pairwise interactions
c
      if (use_smooth) then
         call echarge1g
      else if (use_ewald) then
         if (use_clist) then
            call echarge1f
         else if (use_lights) then
            call echarge1e
         else
            call echarge1d
         end if
      else if (use_clist) then
         call echarge1c
      else if (use_lights) then
         call echarge1b
      else
         call echarge1a
      end if
c
c     get contribution from external electric field if used, along
c     with its lambda derivatives as "lmdadrv" decides
c
      if (use_exfld) then
         mode = 'CHARGE'
         edrv = lmdadrv (use_edlmda,deldlmda,d2eldlmda2)
         if (edrv) then
            call exfield4 (mode)
         else
            call exfield1 (mode,exf)
            ec = ec + exf
         end if
      end if
c
c     add the partial charge virial to main virial
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + ecvir(j,i)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine echarge1a  --  charge derivs via double loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "echarge1a" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a pairwise double loop
c
c
      subroutine echarge1a
      use atoms
      use bound
      use cell
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mutant
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in,ic
      integer kk,kn,kc
      real*8 e,fgrp,de,dc
      real*8 f,fi,fik
      real*8 r,r2,rb,rb2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc,dldedxc,dldedyc,dldedzc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      logical proceed,usei
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      character*6 mode
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     compute the charge interaction energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(i)
         ic = kion(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(i)
         usei = (use(i) .or. use(ic))
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(k)
            kc = kion(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk * cscale(kn)
                  e = fik / rb
                  de = -fik / rb2
                  dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e*dtaper + dtrans) / rc
                     de = de * taper
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the internal virial tensor components
c
                  vxx = xr*dedx + xc*dedxc
                  vyx = yr*dedx + yc*dedxc
                  vzx = zr*dedx + zc*dedxc
                  vyy = yr*dedy + yc*dedyc
                  vzy = zr*dedy + zc*dedyc
                  vzz = zr*dedz + zc*dedzc
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     dldedxc = dldc * xc
                     dldedyc = dldc * yc
                     dldedzc = dldc * zc
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfmdl,demvirdl)
                     call dlcpair (ic,kc,xc,yc,zc,dldedxc,dldedyc,
     &                            dldedzc,dfmdl,demvirdl)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica) then
         deallocate (cscale)
         deallocate (dq)
         return
      end if
c
c     calculate interaction energy with other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         ic = kion(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(i)
         usei = (use(i) .or. use(ic))
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(kk)
            kc = kion(kk)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 2, ncell
                  xc = xic - x(kc)
                  yc = yic - y(kc)
                  zc = zic - z(kc)
                  call imager (xc,yc,zc,j)
                  rc2 = xc*xc + yc*yc + zc*zc
                  if (rc2 .le. off2) then
                     xr = xc + xi - x(k) + x(kc)
                     yr = yc + yi - y(k) + y(kc)
                     zr = zc + zi - z(k) + z(kc)
                     r2 = xr*xr + yr*yr + zr*zr
                     r = sqrt(r2)
                     rb = r + ebuffer
                     rb2 = rb * rb
                     mutik = .false.
                     if (dodl) then
                        mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                     end if
                     mutd2 = (mutik .and. dod2)
                     fqk = fi * pchg(k)
                     if (mutik)  fqk = f
                     fik = fqk
                     if (use_polymer) then
                        if (r2 .le. polycut2)  fik = fik * cscale(kn)
                     end if
                     e = fik / rb
                     de = -fik / rb2
                     dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                     shift = fik / (0.5d0*(off+cut))
                     e = e - shift
                     if (rc2 .gt. cut2) then
                        rc = sqrt(rc2)
                        rc3 = rc2 * rc
                        rc4 = rc2 * rc2
                        rc5 = rc2 * rc3
                        rc6 = rc3 * rc3
                        rc7 = rc3 * rc4
                        taper = c5*rc5 + c4*rc4 + c3*rc3
     &                             + c2*rc2 + c1*rc + c0
                        dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                              + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                        trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                                  + f3*rc3 + f2*rc2 + f1*rc + f0)
                        dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                                  + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                                + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                        dc = (e*dtaper + dtrans) / rc
                        de = de * taper
                        e = e*taper + trans
                     end if
c
c     scale the interaction based on its group membership
c
                     if (use_group) then
                        e = e * fgrp
                        de = de * fgrp
                        dc = dc * fgrp
                     end if
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                     if (mutik) then
                        call dlcscale (pchg(i),pchg(k),dq(i),dq(k),
     &                                 mutd2,e,de,dc,dlambda,dlambda2,
     &                                 dlde,dldc)
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
                     dedxc = dc * xc
                     dedyc = dc * yc
                     dedzc = dc * zc
c
c     increment the energy and gradient values
c
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        if (mutik)  dlambda = 0.5d0 * dlambda
                        if (mutd2)  dlambda2 = 0.5d0 * dlambda2
                     end if
                     ec = ec + e
                     if (mutik)  demdl = demdl + dlambda
                     if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     dec(1,ic) = dec(1,ic) + dedxc
                     dec(2,ic) = dec(2,ic) + dedyc
                     dec(3,ic) = dec(3,ic) + dedzc
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                        dec(1,kc) = dec(1,kc) - dedxc
                        dec(2,kc) = dec(2,kc) - dedyc
                        dec(3,kc) = dec(3,kc) - dedzc
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr*dedx + xc*dedxc
                     vyx = yr*dedx + yc*dedxc
                     vzx = zr*dedx + zc*dedxc
                     vyy = yr*dedy + yc*dedyc
                     vzy = zr*dedy + zc*dedyc
                     vzz = zr*dedz + zc*dedzc
                     ecvir(1,1) = ecvir(1,1) + vxx
                     ecvir(2,1) = ecvir(2,1) + vyx
                     ecvir(3,1) = ecvir(3,1) + vzx
                     ecvir(1,2) = ecvir(1,2) + vyx
                     ecvir(2,2) = ecvir(2,2) + vyy
                     ecvir(3,2) = ecvir(3,2) + vzy
                     ecvir(1,3) = ecvir(1,3) + vzx
                     ecvir(2,3) = ecvir(2,3) + vzy
                     ecvir(3,3) = ecvir(3,3) + vzz
                     if (mutd2) then
                        dlde = dlde / r
                        dldedx = dlde * xr
                        dldedy = dlde * yr
                        dldedz = dlde * zr
                        dldedxc = dldc * xc
                        dldedyc = dldc * yc
                        dldedzc = dldc * zc
                        dfmdl(1,i) = dfmdl(1,i) + dldedx
                        dfmdl(2,i) = dfmdl(2,i) + dldedy
                        dfmdl(3,i) = dfmdl(3,i) + dldedz
                        dfmdl(1,ic) = dfmdl(1,ic) + dldedxc
                        dfmdl(2,ic) = dfmdl(2,ic) + dldedyc
                        dfmdl(3,ic) = dfmdl(3,ic) + dldedzc
                        if (i .ne. k) then
                           dfmdl(1,k) = dfmdl(1,k) - dldedx
                           dfmdl(2,k) = dfmdl(2,k) - dldedy
                           dfmdl(3,k) = dfmdl(3,k) - dldedz
                           dfmdl(1,kc) = dfmdl(1,kc) - dldedxc
                           dfmdl(2,kc) = dfmdl(2,kc) - dldedyc
                           dfmdl(3,kc) = dfmdl(3,kc) - dldedzc
                        end if
                        call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                               demvirdl)
                        call dlcvir (xc,yc,zc,dldedxc,dldedyc,dldedzc,
     &                               demvirdl)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge1b  --  method of lights charge derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge1b" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using the method of lights
c
c
      subroutine echarge1b
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use light
      use mutant
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in,ic
      integer kk,kn,kc
      integer kgy,kgz,kmap
      integer start,stop
      real*8 e,fgrp,de,dc
      real*8 f,fi,fik
      real*8 r,r2,rb,rb2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc,dldedxc,dldedyc,dldedzc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      logical unique,repeat
      character*6 mode
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     transfer the interaction site coordinates to sorting arrays
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         xsort(ii) = x(ic)
         ysort(ii) = y(ic)
         zsort(ii) = z(ic)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         ic = kion(i)
         xic = xsort(rgx(ii))
         yic = ysort(rgy(ii))
         zic = zsort(rgz(ii))
         xi = x(i) - x(ic)
         yi = y(i) - y(ic)
         zi = z(i) - z(ic)
         fi = f * pchg(i)
         usei = (use(i) .or. use(ic))
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(k)
            kc = kion(k)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - xsort(j)
               yc = yic - ysort(kgy)
               zc = zic - zsort(kgz)
               if (use_bounds) then
                  if (abs(xc) .gt. xcell2)  xc = xc - sign(xcell,xc)
                  if (abs(yc) .gt. ycell2)  yc = yc - sign(ycell,yc)
                  if (abs(zc) .gt. zcell2)  zc = zc - sign(zcell,zc)
                  if (monoclinic) then
                     xc = xc + zc*beta_cos
                     zc = zc * beta_sin
                  else if (triclinic) then
                     xc = xc + yc*gamma_cos + zc*beta_cos
                     yc = yc*gamma_sin + zc*beta_term
                     zc = zc * gamma_term
                  end if
               end if
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk
                  if (prime)  fik = fik * cscale(kn)
                  if (use_polymer) then
                     if (r2 .gt. polycut2)  fik = fqk
                  end if
                  e = fik / rb
                  de = -fik / rb2
                  dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e*dtaper + dtrans) / rc
                     de = de * taper
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the internal virial tensor components
c
                  vxx = xr*dedx + xc*dedxc
                  vyx = yr*dedx + yc*dedxc
                  vzx = zr*dedx + zc*dedxc
                  vyy = yr*dedy + yc*dedyc
                  vzy = zr*dedy + zc*dedyc
                  vzz = zr*dedz + zc*dedzc
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     dldedxc = dldc * xc
                     dldedyc = dldc * yc
                     dldedzc = dldc * zc
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfmdl,demvirdl)
                     call dlcpair (ic,kc,xc,yc,zc,dldedxc,dldedyc,
     &                            dldedzc,dfmdl,demvirdl)
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge1c  --  charge derivs via neighbor list  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge1c" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a pairwise neighbor list
c
c
      subroutine echarge1c
      use atoms
      use bound
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use group
      use mutant
      use neigh
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in,ic
      integer kk,kn,kc
      real*8 e,fgrp,de,dc
      real*8 f,fi,fik
      real*8 r,r2,rb,rb2
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xc,yc,zc
      real*8 xic,yic,zic
      real*8 dedx,dedy,dedz
      real*8 dedxc,dedyc,dedzc
      real*8 shift,taper,trans
      real*8 dtaper,dtrans
      real*8 rc,rc2,rc3,rc4
      real*8 rc5,rc6,rc7
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc,dldedxc,dldedyc,dldedzc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      real*8, allocatable :: dfl(:,:)
      logical proceed,usei
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      character*6 mode
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'CHARGE'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     the lambda derivative of the force is reduced over threads in
c     a local array that only has full size when it is needed
c
      if (dod2) then
         allocate (dfl(3,n))
      else
         allocate (dfl(3,1))
      end if
      dfl = 0.0d0
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,kion,use,
!$OMP& x,y,z,f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,
!$OMP& c1scale,c2scale,c3scale,c4scale,c5scale,use_group,use_bounds,
!$OMP& off,off2,cut,cut2,c0,c1,c2,c3,c4,c5,f0,f1,f2,f3,f4,f5,
!$OMP& f6,f7,ebuffer,dq,dodl,dod2,mutg,live,doskip)
!$OMP& firstprivate(cscale)
!$OMP& shared (ec,dec,ecvir,demdl,d2emdl2,dfl,demvirdl)
!$OMP DO reduction(+:ec,dec,ecvir,demdl,d2emdl2,dfl,demvirdl)
c
c     compute the charge interaction energy and first derivatives
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         ic = kion(i)
         xic = x(ic)
         yic = y(ic)
         zic = z(ic)
         xi = x(i) - xic
         yi = y(i) - yic
         zi = z(i) - zic
         fi = f * pchg(i)
         usei = (use(i) .or. use(ic))
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nelst(i)
            k = elst(kk,i)
            kn = jion(k)
            kc = kion(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k) .or. use(kc))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xc = xic - x(kc)
               yc = yic - y(kc)
               zc = zic - z(kc)
               if (use_bounds)  call image (xc,yc,zc)
               rc2 = xc*xc + yc*yc + zc*zc
               if (rc2 .le. off2) then
                  xr = xc + xi - x(k) + x(kc)
                  yr = yc + yi - y(k) + y(kc)
                  zr = zc + zi - z(k) + z(kc)
                  r2 = xr*xr + yr*yr + zr*zr
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk * cscale(kn)
                  e = fik / rb
                  de = -fik / rb2
                  dc = 0.0d0
c
c     use shifted energy switching if near the cutoff distance
c
                  shift = fik / (0.5d0*(off+cut))
                  e = e - shift
                  if (rc2 .gt. cut2) then
                     rc = sqrt(rc2)
                     rc3 = rc2 * rc
                     rc4 = rc2 * rc2
                     rc5 = rc2 * rc3
                     rc6 = rc3 * rc3
                     rc7 = rc3 * rc4
                     taper = c5*rc5 + c4*rc4 + c3*rc3
     &                          + c2*rc2 + c1*rc + c0
                     dtaper = 5.0d0*c5*rc4 + 4.0d0*c4*rc3
     &                           + 3.0d0*c3*rc2 + 2.0d0*c2*rc + c1
                     trans = fik * (f7*rc7 + f6*rc6 + f5*rc5 + f4*rc4
     &                               + f3*rc3 + f2*rc2 + f1*rc + f0)
                     dtrans = fik * (7.0d0*f7*rc6 + 6.0d0*f6*rc5
     &                               + 5.0d0*f5*rc4 + 4.0d0*f4*rc3
     &                             + 3.0d0*f3*rc2 + 2.0d0*f2*rc + f1)
                     dc = (e*dtaper + dtrans) / rc
                     de = de * taper
                     e = e*taper + trans
                  end if
c
c     scale the interaction based on its group membership
c
                  if (use_group) then
                     e = e * fgrp
                     de = de * fgrp
                     dc = dc * fgrp
                  end if
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
                  dedxc = dc * xc
                  dedyc = dc * yc
                  dedzc = dc * zc
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,ic) = dec(1,ic) + dedxc
                  dec(2,ic) = dec(2,ic) + dedyc
                  dec(3,ic) = dec(3,ic) + dedzc
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
                  dec(1,kc) = dec(1,kc) - dedxc
                  dec(2,kc) = dec(2,kc) - dedyc
                  dec(3,kc) = dec(3,kc) - dedzc
c
c     increment the internal virial tensor components
c
                  vxx = xr*dedx + xc*dedxc
                  vyx = yr*dedx + yc*dedxc
                  vzx = zr*dedx + zc*dedxc
                  vyy = yr*dedy + yc*dedyc
                  vzy = zr*dedy + zc*dedyc
                  vzz = zr*dedz + zc*dedzc
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     dldedxc = dldc * xc
                     dldedyc = dldc * yc
                     dldedzc = dldc * zc
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfl,demvirdl)
                     call dlcpair (ic,kc,xc,yc,zc,dldedxc,dldedyc,
     &                            dldedzc,dfl,demvirdl)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add the lambda derivative of the force to its global array
c
      if (dod2) then
         do i = 1, n
            dfmdl(1,i) = dfmdl(1,i) + dfl(1,i)
            dfmdl(2,i) = dfmdl(2,i) + dfl(2,i)
            dfmdl(3,i) = dfmdl(3,i) + dfl(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      deallocate (dfl)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine echarge1d  --  double loop Ewald charge derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "echarge1d" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation
c
c
      subroutine echarge1d
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use ewald
      use group
      use math
      use mutant
      use pme
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in
      integer kk,kn
      real*8 e,de,dc,f
      real*8 fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2
      real*8 fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erfc,erfterm
      real*8 scale
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      logical proceed,usei
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     compute the Ewald self, background charge and cell dipole
c     terms, along with their lambda derivatives
c
      call ecself1 (dodl,dod2,dq)
c
c     compute reciprocal space Ewald energy and first derivatives
c
      call ecrecip1 (dodl,dod2,dq)
c
c     compute the real space Ewald energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(i)
         usei = use(i)
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  scale = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scale)
                  de = -fik * ((erfterm+scale)/rb2
     &                    + (2.0d0*aewald/rootpi)*exp(-rew**2)/rb)
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     dc = 0.0d0
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfmdl,demvirdl)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica) then
         deallocate (cscale)
         deallocate (dq)
         return
      end if
c
c     calculate real space portion involving other unit cells
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(i)
         usei = use(i)
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii, nion
            k = iion(kk)
            kn = jion(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 2, ncell
                  xr = xi - x(k)
                  yr = yi - y(k)
                  zr = zi - z(k)
c
c     find appropriate image and check the real space cutoff
c
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rb = r + ebuffer
                     rb2 = rb * rb
                     mutik = .false.
                     if (dodl) then
                        mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                     end if
                     mutd2 = (mutik .and. dod2)
                     fqk = fi * pchg(k)
                     if (mutik)  fqk = f
                     fik = fqk
                     rew = aewald * r
                     erfterm = erfc (rew)
                     scale = 1.0d0
                     if (use_group)  scale = scale * fgrp
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           scale = scale * cscale(kn)
                        end if
                     end if
                     scale = scale - 1.0d0
                     e = (fik/rb) * (erfterm+scale)
                     de = -fik * ((erfterm+scale)/rb2
     &                       + (2.0d0*aewald/rootpi)*exp(-rew**2)/rb)
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                     if (mutik) then
                        dc = 0.0d0
                        call dlcscale (pchg(i),pchg(k),dq(i),dq(k),
     &                                 mutd2,e,de,dc,dlambda,dlambda2,
     &                                 dlde,dldc)
                     end if
c
c     form the chain rule terms for derivative expressions
c
                     de = de / r
                     dedx = de * xr
                     dedy = de * yr
                     dedz = de * zr
c
c     increment the energy and gradient values
c
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        if (mutik)  dlambda = 0.5d0 * dlambda
                        if (mutd2)  dlambda2 = 0.5d0 * dlambda2
                     end if
                     ec = ec + e
                     if (mutik)  demdl = demdl + dlambda
                     if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                     dec(1,i) = dec(1,i) + dedx
                     dec(2,i) = dec(2,i) + dedy
                     dec(3,i) = dec(3,i) + dedz
                     if (i .ne. k) then
                        dec(1,k) = dec(1,k) - dedx
                        dec(2,k) = dec(2,k) - dedy
                        dec(3,k) = dec(3,k) - dedz
                     end if
c
c     increment the internal virial tensor components
c
                     vxx = xr * dedx
                     vyx = yr * dedx
                     vzx = zr * dedx
                     vyy = yr * dedy
                     vzy = zr * dedy
                     vzz = zr * dedz
                     ecvir(1,1) = ecvir(1,1) + vxx
                     ecvir(2,1) = ecvir(2,1) + vyx
                     ecvir(3,1) = ecvir(3,1) + vzx
                     ecvir(1,2) = ecvir(1,2) + vyx
                     ecvir(2,2) = ecvir(2,2) + vyy
                     ecvir(3,2) = ecvir(3,2) + vzy
                     ecvir(1,3) = ecvir(1,3) + vzx
                     ecvir(2,3) = ecvir(2,3) + vzy
                     ecvir(3,3) = ecvir(3,3) + vzz
                     if (mutd2) then
                        dlde = dlde / r
                        dldedx = dlde * xr
                        dldedy = dlde * yr
                        dldedz = dlde * zr
                        dfmdl(1,i) = dfmdl(1,i) + dldedx
                        dfmdl(2,i) = dfmdl(2,i) + dldedy
                        dfmdl(3,i) = dfmdl(3,i) + dldedz
                        if (i .ne. k) then
                           dfmdl(1,k) = dfmdl(1,k) - dldedx
                           dfmdl(2,k) = dfmdl(2,k) - dldedy
                           dfmdl(3,k) = dfmdl(3,k) - dldedz
                        end if
                        call dlcvir (xr,yr,zr,dldedx,dldedy,dldedz,
     &                               demvirdl)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine echarge1e  --  Ewald charge derivs via lights  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "echarge1e" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and the method of lights
c
c
      subroutine echarge1e
      use atoms
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use ewald
      use group
      use light
      use math
      use mutant
      use pme
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in,ic
      integer kk,kn
      integer kgy,kgz,kmap
      integer start,stop
      real*8 e,de,dc,f
      real*8 fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2
      real*8 fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erfc,erfterm
      real*8 scale
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      logical proceed,usei,prime
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      logical unique,repeat
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
      allocate (xsort(8*n))
      allocate (ysort(8*n))
      allocate (zsort(8*n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     compute the Ewald self, background charge and cell dipole
c     terms, along with their lambda derivatives
c
      call ecself1 (dodl,dod2,dq)
c
c     compute reciprocal space Ewald energy and first derivatives
c
      call ecrecip1 (dodl,dod2,dq)
c
c     compute the real space portion of the Ewald summation;
c     transfer the interaction site coordinates to sorting arrays
c
      do ii = 1, nion
         i = iion(ii)
         ic = kion(i)
         xsort(ii) = x(ic)
         ysort(ii) = y(ic)
         zsort(ii) = z(ic)
      end do
c
c     use the method of lights to generate neighbors
c
      unique = .true.
      call lights (off,nion,xsort,ysort,zsort,unique)
c
c     loop over all atoms computing the interactions
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         xi = xsort(rgx(ii))
         yi = ysort(rgy(ii))
         zi = zsort(rgz(ii))
         fi = f * pchg(i)
         usei = (use(i))
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     loop over method of lights neighbors of current atom
c
         if (kbx(ii) .le. kex(ii)) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = kex(ii)
         else
            repeat = .true.
            start = 1
            stop = kex(ii)
         end if
   10    continue
         do j = start, stop
            kk = locx(j)
            kgy = rgy(kk)
            if (kby(ii) .le. key(ii)) then
               if (kgy.lt.kby(ii) .or. kgy.gt.key(ii))  goto 20
            else
               if (kgy.lt.kby(ii) .and. kgy.gt.key(ii))  goto 20
            end if
            kgz = rgz(kk)
            if (kbz(ii) .le. kez(ii)) then
               if (kgz.lt.kbz(ii) .or. kgz.gt.kez(ii))  goto 20
            else
               if (kgz.lt.kbz(ii) .and. kgz.gt.kez(ii))  goto 20
            end if
            kmap = kk - ((kk-1)/nion)*nion
            k = iion(kmap)
            kn = jion(k)
            prime = (kk .le. nion)
c
c     decide whether to compute the current interaction
c
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               if (use_bounds) then
                  if (abs(xr) .gt. xcell2)  xr = xr - sign(xcell,xr)
                  if (abs(yr) .gt. ycell2)  yr = yr - sign(ycell,yr)
                  if (abs(zr) .gt. zcell2)  zr = zr - sign(zcell,zr)
                  if (monoclinic) then
                     xr = xr + zr*beta_cos
                     zr = zr * beta_sin
                  else if (triclinic) then
                     xr = xr + yr*gamma_cos + zr*beta_cos
                     yr = yr*gamma_sin + zr*beta_term
                     zr = zr * gamma_term
                  end if
               end if
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = 1.0d0
                  if (prime)  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk
                  scale = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scale)
                  de = -fik * ((erfterm+scale)/rb2
     &                    + (2.0d0*aewald/rootpi)*exp(-rew**2)/rb)
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     dc = 0.0d0
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfmdl,demvirdl)
                  end if
               end if
            end if
   20       continue
         end do
         if (repeat) then
            repeat = .false.
            start = kbx(ii) + 1
            stop = nlight
            goto 10
         end if
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine echarge1f  --  Ewald charge derivs via list  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "echarge1f" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     using a particle mesh Ewald summation and a neighbor list
c
c
      subroutine echarge1f
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use couple
      use deriv
      use dlmda
      use energi
      use ewald
      use group
      use math
      use mutant
      use neigh
      use pme
      use shunt
      use usage
      use virial
      implicit none
      integer i,j,k
      integer ii,in
      integer kk,kn
      real*8 e,de,dc,f
      real*8 fi,fik
      real*8 r,r2,rew
      real*8 rb,rb2
      real*8 fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 erfc,erfterm
      real*8 scale
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 fqk
      real*8 dlambda,dlambda2
      real*8 dlde,dldedx,dldedy,dldedz
      real*8 dldc
      real*8, allocatable :: cscale(:)
      real*8, allocatable :: dq(:)
      real*8, allocatable :: dfl(:,:)
      logical proceed,usei
      logical livei,doskip
      logical dodl,dod2,mutik,mutd2
      logical live(0:2)
      character*6 mode
      external erfc
c
c
c     zero out the Ewald summation energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     set the lambda derivative flags, the partial charge derivatives
c     in an array that only has full size when it is needed, and the
c     atom groups whose pairs carry energy or a derivative
c
      call chgdlflags (dodl,dod2)
      if (dodl) then
         allocate (dq(n))
         call chgdlmda (dq)
      else
         allocate (dq(1))
      end if
      call eleclive (live,doskip)
c
c     the lambda derivative of the force is reduced over threads in
c     a local array that only has full size when it is needed
c
      if (dod2) then
         allocate (dfl(3,n))
      else
         allocate (dfl(3,1))
      end if
      dfl = 0.0d0
c
c     compute the Ewald self, background charge and cell dipole
c     terms, along with their lambda derivatives
c
      call ecself1 (dodl,dod2,dq)
c
c     compute reciprocal space Ewald energy and first derivatives
c
      call ecrecip1 (dodl,dod2,dq)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nion,iion,jion,use,x,y,z,
!$OMP& f,pchg,nelst,elst,n12,n13,n14,n15,i12,i13,i14,i15,c1scale,
!$OMP& c2scale,c3scale,c4scale,c5scale,use_group,off2,aewald,ebuffer,
!$OMP& dq,dodl,dod2,mutg,live,doskip)
!$OMP& firstprivate(cscale)
!$OMP& shared (ec,dec,ecvir,demdl,d2emdl2,dfl,demvirdl)
!$OMP DO reduction(+:ec,dec,ecvir,demdl,d2emdl2,dfl,demvirdl)
c
c     compute the real space Ewald energy and first derivatives
c
      do ii = 1, nion
         i = iion(ii)
         in = jion(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(i)
         usei = use(i)
         livei = live(mutg(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = 1, nelst(i)
            k = elst(kk,i)
            kn = jion(k)
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. use(k))
            if (proceed .and. doskip)
     &         proceed = (livei .and. live(mutg(k)))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
c
c     find energy for interactions within real space cutoff
c
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  rb = r + ebuffer
                  rb2 = rb * rb
                  mutik = .false.
                  if (dodl) then
                     mutik = (dq(i).ne.0.0d0 .or. dq(k).ne.0.0d0)
                  end if
                  mutd2 = (mutik .and. dod2)
                  fqk = fi * pchg(k)
                  if (mutik)  fqk = f
                  fik = fqk
                  rew = aewald * r
                  erfterm = erfc (rew)
                  scale = cscale(kn)
                  if (use_group)  scale = scale * fgrp
                  scale = scale - 1.0d0
                  e = (fik/rb) * (erfterm+scale)
                  de = -fik * ((erfterm+scale)/rb2
     &                    + (2.0d0*aewald/rootpi)*exp(-rew**2)/rb)
c
c     find the lambda derivatives from the pair terms for a unit
c     charge product, then scale them by the charge product
c
                  if (mutik) then
                     dc = 0.0d0
                     call dlcscale (pchg(i),pchg(k),dq(i),dq(k),mutd2,
     &                              e,de,dc,dlambda,dlambda2,dlde,dldc)
                  end if
c
c     form the chain rule terms for derivative expressions
c
                  de = de / r
                  dedx = de * xr
                  dedy = de * yr
                  dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
                  ec = ec + e
                  if (mutik)  demdl = demdl + dlambda
                  if (mutd2)  d2emdl2 = d2emdl2 + dlambda2
                  dec(1,i) = dec(1,i) + dedx
                  dec(2,i) = dec(2,i) + dedy
                  dec(3,i) = dec(3,i) + dedz
                  dec(1,k) = dec(1,k) - dedx
                  dec(2,k) = dec(2,k) - dedy
                  dec(3,k) = dec(3,k) - dedz
c
c     increment the internal virial tensor components
c
                  vxx = xr * dedx
                  vyx = yr * dedx
                  vzx = zr * dedx
                  vyy = yr * dedy
                  vzy = zr * dedy
                  vzz = zr * dedz
                  ecvir(1,1) = ecvir(1,1) + vxx
                  ecvir(2,1) = ecvir(2,1) + vyx
                  ecvir(3,1) = ecvir(3,1) + vzx
                  ecvir(1,2) = ecvir(1,2) + vyx
                  ecvir(2,2) = ecvir(2,2) + vyy
                  ecvir(3,2) = ecvir(3,2) + vzy
                  ecvir(1,3) = ecvir(1,3) + vzx
                  ecvir(2,3) = ecvir(2,3) + vzy
                  ecvir(3,3) = ecvir(3,3) + vzz
                  if (mutd2) then
                     dlde = dlde / r
                     dldedx = dlde * xr
                     dldedy = dlde * yr
                     dldedz = dlde * zr
                     call dlcpair (i,k,xr,yr,zr,dldedx,dldedy,dldedz,
     &                            dfl,demvirdl)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add the lambda derivative of the force to its global array
c
      if (dod2) then
         do i = 1, n
            dfmdl(1,i) = dfmdl(1,i) + dfl(1,i)
            dfmdl(2,i) = dfmdl(2,i) + dfl(2,i)
            dfmdl(3,i) = dfmdl(3,i) + dfl(3,i)
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      deallocate (dq)
      deallocate (dfl)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine echarge1g  --  charge derivatives for smoothing  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "echarge1g" calculates the charge-charge interaction energy
c     and first derivatives with respect to Cartesian coordinates
c     for use with potential smoothing methods
c
c
      subroutine echarge1g
      use atoms
      use charge
      use chgpot
      use couple
      use deriv
      use energi
      use group
      use math
      use usage
      use warp
      implicit none
      integer i,j,k
      integer ii,in
      integer kk,kn
      real*8 e,de,fgrp
      real*8 r,r2,rb,rb2
      real*8 f,fi,fik
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 dedx,dedy,dedz
      real*8 erf,erfterm
      real*8 expcut,expterm
      real*8 wterm,width
      real*8 width2,width3
      real*8, allocatable :: cscale(:)
      logical proceed,usei
      external erf
c
c
c     zero out the charge interaction energy and derivatives
c
      ec = 0.0d0
      do i = 1, n
         dec(1,i) = 0.0d0
         dec(2,i) = 0.0d0
         dec(3,i) = 0.0d0
      end do
      if (nion .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         cscale(i) = 1.0d0
      end do
c
c     set the energy units conversion factor
c
      f = electric / dielec
      expcut = -50.0d0
c
c     set the extent of smoothing to be performed
c
      width = deform * diffc
      if (use_dem) then
         if (width .gt. 0.0d0)  width = 0.5d0 / sqrt(width)
      else if (use_gda) then
         wterm = sqrt(3.0d0/(2.0d0*diffc))
      end if
      width2 = width * width
      width3 = width * width2
c
c     compute the charge interaction energy and first derivatives
c
      do ii = 1, nion-1
         i = iion(ii)
         in = jion(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         fi = f * pchg(i)
         usei = (use(i))
c
c     set exclusion coefficients for connected atoms
c
         cscale(in) = c1scale
         do j = 1, n12(in)
            cscale(i12(j,in)) = c2scale
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = c3scale
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = c4scale
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = c5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kk = ii+1, nion
            k = iion(kk)
            kn = jion(k)
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed)  proceed = (usei .or. use(k))
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = xi - x(k)
               yr = yi - y(k)
               zr = zi - z(k)
               r2 = xr*xr + yr*yr + zr*zr
               r = sqrt(r2)
               rb = r + ebuffer
               rb2 = rb * rb
               fik = fi * pchg(k) * cscale(kn)
               e = fik / rb
               de = -fik / rb2
c
c     transform the potential function via smoothing
c
               if (use_dem) then
                  if (width .gt. 0.0d0) then
                     erfterm = erf(width*rb)
                     expterm = -rb2 * width2
                     if (expterm .gt. expcut) then
                        expterm = 2.0d0*fik*width*exp(expterm)
     &                               / (rootpi*rb)
                     else
                        expterm = 0.0d0
                     end if
                     e = e * erfterm
                     de = de*erfterm + expterm
                  end if
               else if (use_gda) then
                  width = m2(i) + m2(k)
                  if (width .gt. 0.0d0) then
                     width = wterm / sqrt(width)
                     erfterm = erf(width*rb)
                     expterm = -rb2 * width * width
                     if (expterm .gt. expcut) then
                        expterm = 2.0d0*fik*width*exp(expterm)
     &                               / (rootpi*rb)
                     else
                        expterm = 0.0d0
                     end if
                     e = e * erfterm
                     de = de*erfterm + expterm
                  end if
               else if (use_tophat) then
                  if (width .gt. rb) then
                     e = fik * (3.0d0*width2-rb2) / (2.0d0*width3)
                     de = -fik * rb / width3
                  end if
               else if (use_stophat) then
                  wterm = rb + width
                  e = fik / wterm
                  de = -fik / (wterm*wterm)
               end if
c
c     scale the interaction based on its group membership
c
               if (use_group) then
                  e = e * fgrp
                  de = de * fgrp
               end if
c
c     form the chain rule terms for derivative expressions
c
               de = de / r
               dedx = de * xr
               dedy = de * yr
               dedz = de * zr
c
c     increment the overall energy and derivative expressions
c
               ec = ec + e
               dec(1,i) = dec(1,i) + dedx
               dec(2,i) = dec(2,i) + dedy
               dec(3,i) = dec(3,i) + dedz
               dec(1,k) = dec(1,k) - dedx
               dec(2,k) = dec(2,k) - dedy
               dec(3,k) = dec(3,k) - dedz
c
c     increment the total intermolecular energy
c
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         cscale(in) = 1.0d0
         do j = 1, n12(in)
            cscale(i12(j,in)) = 1.0d0
         end do
         do j = 1, n13(in)
            cscale(i13(j,in)) = 1.0d0
         end do
         do j = 1, n14(in)
            cscale(i14(j,in)) = 1.0d0
         end do
         do j = 1, n15(in)
            cscale(i15(j,in)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ecself1  --  Ewald charge self & boundary terms  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ecself1" finds the Ewald self energy, the uniform background
c     charge correction and the cell dipole boundary correction for
c     partial charges, with their gradient and virial, along with
c     their lambda derivatives when "dodl" is set, where "dq" holds the
c     lambda derivative of each partial charge, and the second, force
c     and virial lambda derivatives are found when "dod2" is set
c
c     each term is a quadratic form in the partial charges, so its
c     first lambda derivative pairs the charges with their lambda
c     derivatives, and its second pairs the lambda derivatives alone
c
c
      subroutine ecself1 (dodl,dod2,dq)
      use atoms
      use boxes
      use charge
      use chgpot
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use virial
      implicit none
      integer i,j,ii
      real*8 e,de,f,fs,term
      real*8 sum,dsum
      real*8 dle,d2le
      real*8 xd,yd,zd
      real*8 dxd,dyd,dzd
      real*8 dedx,dedy,dedz
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 dq(*)
      logical dodl,dod2
c
c
c     set the conversion factor for the Ewald terms
c
      f = electric / dielec
      dle = 0.0d0
      d2le = 0.0d0
c
c     compute the Ewald self-energy term over all the atoms
c
      fs = -f * aewald / rootpi
      do ii = 1, nion
         i = iion(ii)
         e = fs * pchg(i)**2
         ec = ec + e
         if (dodl) then
            dle = dle + 2.0d0*fs*pchg(i)*dq(i)
            d2le = d2le + 2.0d0*fs*dq(i)**2
         end if
      end do
c
c     compute the uniform background charge correction term, where
c     the second lambda derivative is kept as the net charge vanishes
c
      fs = -0.5d0 * f * pi / (volbox*aewald**2)
      sum = 0.0d0
      dsum = 0.0d0
      do ii = 1, nion
         i = iion(ii)
         sum = sum + pchg(i)
         if (dodl)  dsum = dsum + dq(i)
      end do
      if (abs(sum) .gt. 1.0d-10) then
         e = fs * sum**2
         ec = ec + e
c
c     the background term scales as the inverse volume
c
         do j = 1, 3
            ecvir(j,j) = ecvir(j,j) - e
         end do
         if (dodl) then
            de = 2.0d0 * fs * sum * dsum
            dle = dle + de
            if (dod2) then
               do j = 1, 3
                  demvirdl(j,j) = demvirdl(j,j) - de
               end do
            end if
         end if
      end if
      if (dodl)  d2le = d2le + 2.0d0*fs*dsum**2
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         dxd = 0.0d0
         dyd = 0.0d0
         dzd = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            xd = xd + pchg(i)*x(i)
            yd = yd + pchg(i)*y(i)
            zd = zd + pchg(i)*z(i)
            if (dodl) then
               dxd = dxd + dq(i)*x(i)
               dyd = dyd + dq(i)*y(i)
               dzd = dzd + dq(i)*z(i)
            end if
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         ec = ec + e
         do ii = 1, nion
            i = iion(ii)
            de = 2.0d0 * term * pchg(i)
            dedx = de * xd
            dedy = de * yd
            dedz = de * zd
            dec(1,i) = dec(1,i) + dedx
            dec(2,i) = dec(2,i) + dedy
            dec(3,i) = dec(3,i) + dedz
         end do
c
c     boundary correction to virial from the charge positions, and
c     from the inverse volume dependence of the boundary energy
c
         vxx = 2.0d0*term*xd*xd - e
         vyx = 2.0d0*term*yd*xd
         vzx = 2.0d0*term*zd*xd
         vyy = 2.0d0*term*yd*yd - e
         vzy = 2.0d0*term*zd*yd
         vzz = 2.0d0*term*zd*zd - e
         ecvir(1,1) = ecvir(1,1) + vxx
         ecvir(2,1) = ecvir(2,1) + vyx
         ecvir(3,1) = ecvir(3,1) + vzx
         ecvir(1,2) = ecvir(1,2) + vyx
         ecvir(2,2) = ecvir(2,2) + vyy
         ecvir(3,2) = ecvir(3,2) + vzy
         ecvir(1,3) = ecvir(1,3) + vzx
         ecvir(2,3) = ecvir(2,3) + vzy
         ecvir(3,3) = ecvir(3,3) + vzz
c
c     the force of each charge follows the cell dipole, so its lambda
c     derivative takes the charge derivative and the dipole derivative
c
         if (dodl) then
            de = 2.0d0 * term * (xd*dxd+yd*dyd+zd*dzd)
            dle = dle + de
            d2le = d2le + 2.0d0*term*(dxd*dxd+dyd*dyd+dzd*dzd)
            if (dod2) then
               do ii = 1, nion
                  i = iion(ii)
                  dedx = 2.0d0 * term * (dq(i)*xd+pchg(i)*dxd)
                  dedy = 2.0d0 * term * (dq(i)*yd+pchg(i)*dyd)
                  dedz = 2.0d0 * term * (dq(i)*zd+pchg(i)*dzd)
                  dfmdl(1,i) = dfmdl(1,i) + dedx
                  dfmdl(2,i) = dfmdl(2,i) + dedy
                  dfmdl(3,i) = dfmdl(3,i) + dedz
               end do
               vxx = 4.0d0*term*xd*dxd - de
               vyx = 2.0d0*term*(yd*dxd+dyd*xd)
               vzx = 2.0d0*term*(zd*dxd+dzd*xd)
               vyy = 4.0d0*term*yd*dyd - de
               vzy = 2.0d0*term*(zd*dyd+dzd*yd)
               vzz = 4.0d0*term*zd*dzd - de
               demvirdl(1,1) = demvirdl(1,1) + vxx
               demvirdl(2,1) = demvirdl(2,1) + vyx
               demvirdl(3,1) = demvirdl(3,1) + vzx
               demvirdl(1,2) = demvirdl(1,2) + vyx
               demvirdl(2,2) = demvirdl(2,2) + vyy
               demvirdl(3,2) = demvirdl(3,2) + vzy
               demvirdl(1,3) = demvirdl(1,3) + vzx
               demvirdl(2,3) = demvirdl(2,3) + vzy
               demvirdl(3,3) = demvirdl(3,3) + vzz
            end if
         end if
      end if
c
c     increment the lambda derivatives of the charge energy
c
      if (dodl) then
         demdl = demdl + dle
         if (dod2)  d2emdl2 = d2emdl2 + d2le
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ecrecip1  --  PME recip charge energy & derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ecrecip1" evaluates the reciprocal space portion of the particle
c     mesh Ewald summation energy and gradient due to partial charges,
c     along with the lambda derivatives when "dodl" is set, where "dq"
c     holds the lambda derivative of each partial charge
c
c     the first lambda derivative is the potential of the charges at
c     each site times its charge derivative; when second derivatives
c     are requested by "dod2", the charge derivatives fill the
c     imaginary part of the grid, so that one transform carries both
c     densities, and each wavevector is visited with its opposite to
c     separate the two structure factors for the virial
c
c     literature reference:
c
c     U. Essmann, L. Perera, M. L Berkowitz, T. Darden, H. Lee and
c     L. G. Pedersen, "A Smooth Particle Mesh Ewald Method", Journal
c     of Chemical Physics, 103, 8577-8593 (1995)
c
c     W. Smith and D. Fincham, "The Ewald Sum in Truncated Octahedral
c     and Rhombic Dodecahedral Boundary Conditions", Molecular
c     Simulation, 10, 67-71 (1993)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine ecrecip1 (dodl,dod2,dq)
      use atoms
      use bound
      use boxes
      use charge
      use chgpot
      use deriv
      use dlmda
      use energi
      use ewald
      use math
      use pme
      use virial
      implicit none
      integer i,j,ii
      integer k1,k2,k3
      integer j1,j2,j3
      integer m1,m2,m3
      integer nf1,nf2,nf3
      integer nff,ntot
      integer ip,ipt,jpt,npt
      integer kp(3,2)
      real*8 e,term,eterm
      real*8 expterm
      real*8 vterm,pterm
      real*8 volterm
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 f,fi,denom
      real*8 hsq,struc2
      real*8 h1,h2,h3
      real*8 r1,r2,r3
      real*8 f1,f2,f3
      real*8 q1,q2,lq1,lq2
      real*8 dstruc2dl
      real*8 dldeterm,dldvterm
      real*8 dldvxx,dldvyy,dldvzz
      real*8 dldvxy,dldvxz,dldvyz
      real*8 dle,d2le,dfi
      real*8 qa(2,2)
      real*8 dq(*)
      real*8, allocatable :: fphi(:,:)
      real*8, allocatable :: dqphi(:,:)
      logical dodl,dod2
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid)) then
         if (size(qgrid) .ne. 2*ntot)  call fftclose
      end if
      if (.not. allocated(qgrid))  call fftsetup
c
c     setup spatial decomposition and B-spline coefficients
c
      call getchunk
      call moduli
      call bspline_fill
      call table_fill
c
c     assign PME grid and perform 3-D FFT forward transform, where
c     the lambda derivative charges fill the imaginary part of the
c     grid when second lambda derivatives are needed
c
      if (dod2) then
         call grid_pchg2 (dq)
      else
         call grid_pchg
      end if
      call fftfront (qgrid)
c
c     zero out the temporary virial accumulation variables
c
      vxx = 0.0d0
      vxy = 0.0d0
      vxz = 0.0d0
      vyy = 0.0d0
      vyz = 0.0d0
      vzz = 0.0d0
      dldvxx = 0.0d0
      dldvxy = 0.0d0
      dldvxz = 0.0d0
      dldvyy = 0.0d0
      dldvyz = 0.0d0
      dldvzz = 0.0d0
c
c     use scalar sum to get reciprocal space energy and virial
c
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      nff = nfft1 * nfft2
      ntot = nff * nfft3
      if (.not. dod2) then
         do i = 1, ntot-1
            k3 = i/nff + 1
            j = i - (k3-1)*nff
            k2 = j/nfft1 + 1
            k1 = j - (k2-1)*nfft1 + 1
            m1 = k1 - 1
            m2 = k2 - 1
            m3 = k3 - 1
            if (k1 .gt. nf1)  m1 = m1 - nfft1
            if (k2 .gt. nf2)  m2 = m2 - nfft2
            if (k3 .gt. nf3)  m3 = m3 - nfft3
            r1 = dble(m1)
            r2 = dble(m2)
            r3 = dble(m3)
            h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
            h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
            h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
            hsq = h1*h1 + h2*h2 + h3*h3
            term = -pterm * hsq
            expterm = 0.0d0
            if (term .gt. -50.0d0) then
               denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
               expterm = exp(term) / denom
               if (.not. use_bounds) then
                  expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
               else if (nonprism) then
                  if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
               end if
               struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
               eterm = 0.5d0 * f * expterm * struc2
               vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
               vxx = vxx + h1*h1*vterm - eterm
               vxy = vxy + h1*h2*vterm
               vxz = vxz + h1*h3*vterm
               vyy = vyy + h2*h2*vterm - eterm
               vyz = vyz + h2*h3*vterm
               vzz = vzz + h3*h3*vterm - eterm
            end if
            qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
            qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
         end do
c
c     visit each wavevector together with its opposite, keeping the
c     values at both before either is scaled in place, since the
c     packed charge and charge derivative transforms at one are
c     separated using the other
c
      else
         do i = 1, ntot-1
            k3 = i/nff + 1
            j = i - (k3-1)*nff
            k2 = j/nfft1 + 1
            k1 = j - (k2-1)*nfft1 + 1
            call pmeopp (k1,k2,k3,j1,j2,j3)
            ip = (j3-1)*nff + (j2-1)*nfft1 + j1 - 1
            if (ip .lt. i)  goto 10
            kp(1,1) = k1
            kp(2,1) = k2
            kp(3,1) = k3
            kp(1,2) = j1
            kp(2,2) = j2
            kp(3,2) = j3
            do j = 1, 2
               qa(j,1) = qgrid(j,k1,k2,k3)
               qa(j,2) = qgrid(j,j1,j2,j3)
            end do
            npt = 2
            if (ip .eq. i)  npt = 1
            do ipt = 1, npt
               jpt = 3 - ipt
               k1 = kp(1,ipt)
               k2 = kp(2,ipt)
               k3 = kp(3,ipt)
               m1 = k1 - 1
               m2 = k2 - 1
               m3 = k3 - 1
               if (k1 .gt. nf1)  m1 = m1 - nfft1
               if (k2 .gt. nf2)  m2 = m2 - nfft2
               if (k3 .gt. nf3)  m3 = m3 - nfft3
               r1 = dble(m1)
               r2 = dble(m2)
               r3 = dble(m3)
               h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
               h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
               h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
               hsq = h1*h1 + h2*h2 + h3*h3
               term = -pterm * hsq
               expterm = 0.0d0
               if (term .gt. -50.0d0) then
                  denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
                  expterm = exp(term) / denom
                  if (.not. use_bounds) then
                     expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
                  else if (nonprism) then
                     if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
                  end if
                  call pmesplit (qa(1,ipt),qa(1,jpt),q1,q2,lq1,lq2)
                  struc2 = q1**2 + q2**2
                  eterm = 0.5d0 * f * expterm * struc2
                  vterm = (2.0d0/hsq) * (1.0d0-term) * eterm
                  vxx = vxx + h1*h1*vterm - eterm
                  vxy = vxy + h1*h2*vterm
                  vxz = vxz + h1*h3*vterm
                  vyy = vyy + h2*h2*vterm - eterm
                  vyz = vyz + h2*h3*vterm
                  vzz = vzz + h3*h3*vterm - eterm
                  dstruc2dl = 2.0d0 * (q1*lq1 + q2*lq2)
                  dldeterm = 0.5d0 * f * expterm * dstruc2dl
                  dldvterm = (2.0d0/hsq) * (1.0d0-term) * dldeterm
                  dldvxx = dldvxx + h1*h1*dldvterm - dldeterm
                  dldvxy = dldvxy + h1*h2*dldvterm
                  dldvxz = dldvxz + h1*h3*dldvterm
                  dldvyy = dldvyy + h2*h2*dldvterm - dldeterm
                  dldvyz = dldvyz + h2*h3*dldvterm
                  dldvzz = dldvzz + h3*h3*dldvterm - dldeterm
               end if
               qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
               qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
            end do
   10       continue
         end do
      end if
c
c     increment the total internal virial tensor components
c
      ecvir(1,1) = ecvir(1,1) + vxx
      ecvir(2,1) = ecvir(2,1) + vxy
      ecvir(3,1) = ecvir(3,1) + vxz
      ecvir(1,2) = ecvir(1,2) + vxy
      ecvir(2,2) = ecvir(2,2) + vyy
      ecvir(3,2) = ecvir(3,2) + vyz
      ecvir(1,3) = ecvir(1,3) + vxz
      ecvir(2,3) = ecvir(2,3) + vyz
      ecvir(3,3) = ecvir(3,3) + vzz
      if (dod2) then
         demvirdl(1,1) = demvirdl(1,1) + dldvxx
         demvirdl(2,1) = demvirdl(2,1) + dldvxy
         demvirdl(3,1) = demvirdl(3,1) + dldvxz
         demvirdl(1,2) = demvirdl(1,2) + dldvxy
         demvirdl(2,2) = demvirdl(2,2) + dldvyy
         demvirdl(3,2) = demvirdl(3,2) + dldvyz
         demvirdl(1,3) = demvirdl(1,3) + dldvxz
         demvirdl(2,3) = demvirdl(2,3) + dldvyz
         demvirdl(3,3) = demvirdl(3,3) + dldvzz
      end if
c
c     account for zeroth grid point for nonperiodic system
c
      qgrid(1,1,1,1) = 0.0d0
      qgrid(2,1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qgrid(1,1,1,1) = expterm * qgrid(1,1,1,1)
         qgrid(2,1,1,1) = expterm * qgrid(2,1,1,1)
      end if
c
c     perform the 3-D FFT backward transformation
c
      call fftback (qgrid)
c
c     perform dynamic allocation of some local arrays
c
      allocate (fphi(4,n))
      if (dod2)  allocate (dqphi(4,n))
c
c     extract the partial charge electrostatic potential, along
c     with that of the charge derivatives as needed
c
      if (dod2) then
         call fphi_pchg2 (fphi,dqphi)
      else
         call fphi_pchg (fphi)
      end if
c
c     increment partial charge energy and gradient components
c
      e = 0.0d0
      do ii = 1, nion
         i = iion(ii)
         fi = f * pchg(i)
         e = e + fi*fphi(1,i)
         f1 = fi * dble(nfft1) * fphi(2,i)
         f2 = fi * dble(nfft2) * fphi(3,i)
         f3 = fi * dble(nfft3) * fphi(4,i)
         h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
         h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
         h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
         dec(1,i) = dec(1,i) + h1
         dec(2,i) = dec(2,i) + h2
         dec(3,i) = dec(3,i) + h3
      end do
      e = 0.5d0 * e
      ec = ec + e
c
c     the first lambda derivative is the charge derivative times
c     the potential, and the second and the force derivative follow
c     from the potential of the charge derivatives
c
      if (dodl) then
         dle = 0.0d0
         d2le = 0.0d0
         do ii = 1, nion
            i = iion(ii)
            dfi = f * dq(i)
            dle = dle + dfi*fphi(1,i)
            if (dod2) then
               fi = f * pchg(i)
               d2le = d2le + dfi*dqphi(1,i)
               f1 = dble(nfft1) * (dfi*fphi(2,i)+fi*dqphi(2,i))
               f2 = dble(nfft2) * (dfi*fphi(3,i)+fi*dqphi(3,i))
               f3 = dble(nfft3) * (dfi*fphi(4,i)+fi*dqphi(4,i))
               h1 = recip(1,1)*f1 + recip(1,2)*f2 + recip(1,3)*f3
               h2 = recip(2,1)*f1 + recip(2,2)*f2 + recip(2,3)*f3
               h3 = recip(3,1)*f1 + recip(3,2)*f2 + recip(3,3)*f3
               dfmdl(1,i) = dfmdl(1,i) + h1
               dfmdl(2,i) = dfmdl(2,i) + h2
               dfmdl(3,i) = dfmdl(3,i) + h3
            end if
         end do
         demdl = demdl + dle
         if (dod2)  d2emdl2 = d2emdl2 + d2le
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (fphi)
      if (dod2)  deallocate (dqphi)
      return
      end
