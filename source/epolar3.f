c
c
c     ##################################################
c     ##  COPYRIGHT (C) 2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved             ##
c     ##################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the induced dipole polarization energy,
c     and partitions the energy among atoms
c
c
      subroutine epolar3
      use limits
      implicit none
      logical pairwise
c
c
c     choose the method for summing over polarization interactions
c
      pairwise = .true.
      if (pairwise) then
         if (use_ewald) then
            if (use_mlist) then
               call epolar3d
            else
               call epolar3c
            end if
         else
            if (use_mlist) then
               call epolar3b
            else
               call epolar3a
            end if
         end if
      else
         call epolar3x
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3a  --  double loop polarization analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3a" calculates the induced dipole polarization energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine epolar3a
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpot
      use couple
      use energi
      use inform
      use iounit
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,f
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 sc3,sc5,sc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8, allocatable :: pscale(:)
      logical header,huge
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
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
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do k = i+1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kk)
               psr5 = rr5 * sc5 * pscale(kk)
               psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qri1*xr + qri2*yr + qri3*zr
               qrrk = qrk1*xr + qrk2*yr + qrk3*zr
               duik = uix*dkx + dix*ukx + uiy*dky
     &                   + diy*uky + uiz*dkz + diz*ukz
               quik = qri1*ukx + qri2*uky + qri3*ukz
     &                   - qrk1*uix - qrk2*uiy - qrk3*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contributions for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  nep = nep + 1
                  ep = ep + e
                  aep(ii) = aep(ii) + 0.5d0*e
                  aep(kk) = aep(kk) + 0.5d0*e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 100.0d0)
               if ((debug.and.e.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c
         do i = 1, npole
            ii = ipole(i)
            pdi = pdamp(i)
            pti = thole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
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
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
               do k = 1, np11(ii)
                   if (i14(j,ii) .eq. ip11(k,ii))
     &               pscale(i14(j,ii)) = p4scale * p41scale
               end do
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do k = i, npole
               kk = ipole(k)
               do jcell = 1, ncell
                  xr = x(kk) - xi
                  yr = y(kk) - yi
                  zr = z(kk) - zi
                  if (use_bounds)  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2)) then
                     pscale(kk) = 1.0d0
                  end if
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
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
c
c     apply Thole polarization damping to scale factors
c
                     sc3 = 1.0d0
                     sc5 = 1.0d0
                     sc7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           sc3 = 1.0d0 - expdamp
                           sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                           sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                        end if
                     end if
c
c     intermediates involving Thole damping and scale factors
c
                     psr3 = rr3 * sc3 * pscale(kk)
                     psr5 = rr5 * sc5 * pscale(kk)
                     psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
                     dri = dix*xr + diy*yr + diz*zr
                     drk = dkx*xr + dky*yr + dkz*zr
                     uri = uix*xr + uiy*yr + uiz*zr
                     urk = ukx*xr + uky*yr + ukz*zr
                     qri1 = qixx*xr + qixy*yr + qixz*zr
                     qri2 = qixy*xr + qiyy*yr + qiyz*zr
                     qri3 = qixz*xr + qiyz*yr + qizz*zr
                     qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
                     qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
                     qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
                     qrri = qri1*xr + qri2*yr + qri3*zr
                     qrrk = qrk1*xr + qrk2*yr + qrk3*zr
                     duik = uix*dkx + dix*ukx + uiy*dky
     &                         + diy*uky + uiz*dkz + diz*ukz
                     quik = qri1*ukx + qri2*uky + qri3*ukz
     &                         - qrk1*uix - qrk2*uiy - qrk3*uiz
c
c     calculate intermediate terms for polarization interaction
c
                     term1 = ck*uri - ci*urk + duik
                     term2 = 2.0d0*quik - uri*drk - dri*urk
                     term3 = uri*qrrk - urk*qrri
c
c     compute the energy contributions for this interaction
c
                     e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
                     if (e .ne. 0.0d0) then
                        nep = nep + 1
                        ep = ep + e
                        aep(ii) = aep(ii) + 0.5d0*e
                        aep(kk) = aep(kk) + 0.5d0*e
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 100.0d0)
                     if ((debug.and.e.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Polarization',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  ii,name(ii),kk,name(kk),r,e
   50                   format (' Polar',5x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epolar3b  --  polarization analysis via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epolar3b" calculates the induced dipole polarization energy
c     using a neighbor list, and partitions the energy among atoms
c
c
      subroutine epolar3b
      use sizes
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpot
      use couple
      use energi
      use inform
      use iounit
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,f
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 sc3,sc5,sc7
      real*8 psr3,psr5,psr7
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dri,drk,uri,urk
      real*8 qri1,qri2,qri3
      real*8 qrk1,qrk2,qrk3
      real*8 qrri,qrrk
      real*8 duik,quik
      real*8 term1,term2,term3
      real*8, allocatable :: pscale(:)
      logical header,huge
      character*6 mode
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     set arrays needed to scale interactions and store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,pdamp,thole,x,y,z,rpole,uind,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,ip11,p2scale,p3scale,p4scale,p5scale,
!$OMP& p41scale,nelst,elst,use_bounds,off2,f,name,verbose,debug,
!$OMP& header,iout)
!$OMP& firstprivate(pscale) shared (ep,nep,aep)
!$OMP DO reduction(+:ep,nep,aep) schedule(guided)
c
c     compute the dipole polarization energy component
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         xi = x(ii)
         yi = y(ii)
         zi = z(ii)
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
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
c
c     apply Thole polarization damping to scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     sc3 = 1.0d0 - expdamp
                     sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                     sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                    *expdamp
                  end if
               end if
c
c     intermediates involving Thole damping and scale factors
c
               psr3 = rr3 * sc3 * pscale(kk)
               psr5 = rr5 * sc5 * pscale(kk)
               psr7 = rr7 * sc7 * pscale(kk)
c
c     intermediates involving moments and distance separation
c
               dri = dix*xr + diy*yr + diz*zr
               drk = dkx*xr + dky*yr + dkz*zr
               uri = uix*xr + uiy*yr + uiz*zr
               urk = ukx*xr + uky*yr + ukz*zr
               qri1 = qixx*xr + qixy*yr + qixz*zr
               qri2 = qixy*xr + qiyy*yr + qiyz*zr
               qri3 = qixz*xr + qiyz*yr + qizz*zr
               qrk1 = qkxx*xr + qkxy*yr + qkxz*zr
               qrk2 = qkxy*xr + qkyy*yr + qkyz*zr
               qrk3 = qkxz*xr + qkyz*yr + qkzz*zr
               qrri = qri1*xr + qri2*yr + qri3*zr
               qrrk = qrk1*xr + qrk2*yr + qrk3*zr
               duik = uix*dkx + dix*ukx + uiy*dky
     &                   + diy*uky + uiz*dkz + diz*ukz
               quik = qri1*ukx + qri2*uky + qri3*ukz
     &                   - qrk1*uix - qrk2*uiy - qrk3*uiz
c
c     calculate intermediate terms for polarization interaction
c
               term1 = ck*uri - ci*urk + duik
               term2 = 2.0d0*quik - uri*drk - dri*urk
               term3 = uri*qrrk - urk*qrri
c
c     compute the energy contributions for this interaction
c
               e = term1*psr3 + term2*psr5 + term3*psr7
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  nep = nep + 1
                  ep = ep + e
                  aep(ii) = aep(ii) + 0.5d0*e
                  aep(kk) = aep(kk) + 0.5d0*e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 100.0d0)
               if ((debug.and.e.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  ii,name(ii),kk,name(kk),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3c  --  Ewald polarization analysis; loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3c" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a double loop, and
c     partitions the energy among atoms
c
c
      subroutine epolar3c
      implicit none
c
c
c     get the polarization energy from the field-based routine
c
      call epolar3x
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3d  --  Ewald polarization analysis; list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3d" calculates the induced dipole polarization energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among atoms
c
c
      subroutine epolar3d
      implicit none
c
c
c     get the polarization energy from the field-based routine
c
      call epolar3x
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3x  --  single-loop polarization analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3x" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar3x
      use sizes
      use action
      use analyz
      use atoms
      use atomid
      use boxes
      use chgpot
      use energi
      use ewald
      use inform
      use iounit
      use limits
      use math
      use mpole
      use polar
      use potent
      use units
      implicit none
      integer i,j,ii
      real*8 e,f,fi
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      logical header,huge
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
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
c     print header information if debug output was requested
c
      header = .true.
      if (debug) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',9x,'Atom Name',24x,'Alpha',8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,polarity,f,uind,udirp,name,verbose,
!$OMP& debug,header,iout,ep,nep,aep)
!$OMP DO reduction(+:ep,nep,aep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            ii = ipole(i)
            fi = f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*uind(j,i)*udirp(j,i)
            end do
            nep = nep + 1
            ep = ep + e
            aep(ii) = aep(ii) + e
c
c     print a message if the energy for this site is large
c
            huge = (abs(e) .gt. 100.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Dipole Polarization',
     &                       ' Interactions :',
     &                    //,' Type',9x,'Atom Name',24x,'Alpha',
     &                       8x,'Energy',/)
               end if
               write (iout,30)  ii,name(ii),polarity(i),e
   30          format (' Polar',5x,i7,'-',a3,16x,2f14.4)
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
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
            nep = nep + 1
            ep = ep + e
            do i = 1, npole
               ii = ipole(i)
               aep(ii) = aep(ii) + e/dble(npole)
            end do
         end if
      end if
      return
      end
