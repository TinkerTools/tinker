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
c     choose the method to sum over polarization interactions
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
         call epolar3e
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
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use extfld
      use inform
      use inter
      use iounit
      use mplpot
      use molcul
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,f,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7)
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
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole ('MPOLE')
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
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
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  call damptholed (i,k,7,r,dmpik)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3 = dmpik(3) * rr3
                  rr5 = dmpik(5) * rr5
                  rr7 = dmpik(7) * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  ep = ep + e
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*e
                  aep(k) = aep(k) + 0.5d0*e
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + e
                  end if
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 10.0d0)
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
                  write (iout,30)  i,name(i),k,name(k),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
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
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
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
            if (use_chgpen) then
               corei = pcore(i)
               vali = pval(i)
               alphai = palpha(i)
            end if
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               do jcell = 2, ncell
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2)) then
                     pscale(k) = 1.0d0
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
                     diu = dix*ukx + diy*uky + diz*ukz
                     qiu = qix*ukx + qiy*uky + qiz*ukz
                     uir = uix*xr + uiy*yr + uiz*zr
                     dku = dkx*uix + dky*uiy + dkz*uiz
                     qku = qkx*uix + qky*uiy + qkz*uiz
                     ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
                     if (use_thole) then
                        call damptholed (i,k,7,r,dmpik)
                        scalek = pscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3 = dmpik(3) * rr3
                        rr5 = dmpik(5) * rr5
                        rr7 = dmpik(7) * rr7
                        term1 = ck*uir - ci*ukr + diu + dku
                        term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                        term3 = uir*qkr - ukr*qir
                        e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(k)
                        valk = pval(k)
                        alphak = palpha(k)
                        call dampdir (r,alphai,alphak,dmpi,dmpk)
                        scalek = pscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = dmpi(3) * rr3
                        rr5i = dmpi(5) * rr5
                        rr7i = dmpi(7) * rr7
                        rr3k = dmpk(3) * rr3
                        rr5k = dmpk(5) * rr5
                        rr7k = dmpk(7) * rr7
                        e = uir*(corek*rr3+valk*rr3k)
     &                         - ukr*(corei*rr3+vali*rr3i)
     &                         + diu*rr3i + dku*rr3k
     &                         + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                         - dkr*uir*rr5k - dir*ukr*rr5i
     &                         + qkr*uir*rr7k - qir*ukr*rr7i
                     end if
c
c     increment the overall polarization energy components
c
                     if (i .eq. k)  e = 0.5d0 * e
                     if (e .ne. 0.0d0) then
                        ep = ep + e
                        nep = nep + 1
                        aep(i) = aep(i) + 0.5d0*e
                        aep(k) = aep(k) + 0.5d0*e
                        einter = einter + e
                     end if
c
c     print a message if the energy of this interaction is large
c
                     huge = (abs(e) .gt. 10.0d0)
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
                        write (iout,50)  i,name(i),k,name(k),r,e
   50                   format (' Polar',5x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     increment polarization energy due to external field
c
      if (use_exfld) then
         do i = 1, npole
            e = 0.0d0
            do j = 1, 3
               e = e - f*uind(j,i)*texfld(j)
            end do
            ep = ep + e
            nep = nep + 1
            aep(i) = aep(i) + e
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
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use extfld
      use inform
      use inter
      use iounit
      use molcul
      use mplpot
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
      real*8 e,f,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7)
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
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole ('MPOLE')
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
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
!$OMP& shared(npole,ipole,rpole,x,y,z,pcore,pval,palpha,uind,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP& p2scale,p3scale,p4scale,p5scale,p2iscale,p3iscale,p4iscale,
!$OMP& p5iscale,nelst,elst,use_thole,use_chgpen,use_bounds,off2,f,
!$OMP& texfld,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(pscale) shared (ep,nep,aep,einter)
!$OMP DO reduction(+:ep,nep,aep,einter) schedule(guided)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  call damptholed (i,k,7,r,dmpik)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3 = dmpik(3) * rr3
                  rr5 = dmpik(5) * rr5
                  rr7 = dmpik(7) * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  e = term1*rr3 + term2*rr5 + term3*rr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     increment the overall polarization energy components
c
               if (e .ne. 0.0d0) then
                  ep = ep + e
                  nep = nep + 1
                  aep(i) = aep(i) + 0.5d0*e
                  aep(k) = aep(k) + 0.5d0*e
                  if (molcule(i) .ne. molcule(k))
     &               einter = einter + e
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(e) .gt. 10.0d0)
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
                  write (iout,30)  i,name(i),k,name(k),r,e
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
c
c     increment polarization energy due to external field
c
      if (use_exfld) then
!$OMP    DO reduction(+:ep,nep,aep) schedule(guided)
         do i = 1, npole
            e = 0.0d0
            do j = 1, 3
               e = e - f*uind(j,i)*texfld(j)
            end do
            ep = ep + e
            nep = nep + 1
            aep(i) = aep(i) + e
         end do
!$OMP    END DO
      end if
c
c     OpenMP directives for the major loop structure
c
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
c     "epolar3c" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a double loop
c
c
      subroutine epolar3c
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole ('MPOLE')
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip3
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
         nep = nep + 1
         aep(i) = aep(i) + e
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,i) + rpole(1,i)*x(i)
            yd = yd + rpole(3,i) + rpole(1,i)*y(i)
            zd = zd + rpole(4,i) + rpole(1,i)*z(i)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         nep = nep + 1
         do ii = 1, npole
            i = ipole(ii)
            aep(i) = aep(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epreal3c  --  real space polar analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epreal3c" calculates the induced dipole polarization energy and
c     analysis using particle mesh Ewald summation and a double loop
c
c
      subroutine epreal3c
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use energi
      use extfld
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,jcell
      real*8 e,efull
      real*8 f,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 sr3,sr5,sr7
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7),dmpe(7)
      real*8, allocatable :: pscale(:)
      logical header,huge
      character*6 mode
c
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate real space Ewald error function damping
c
               call dampewald (7,r,r2,f,dmpe)
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  call damptholed (i,k,7,r,dmpik)
                  scalek = pscale(k)
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  sr3 = scalek * dmpik(3) * rr3
                  sr5 = scalek * dmpik(5) * rr5
                  sr7 = scalek * dmpik(7) * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  efull = term1*sr3 + term2*sr5 + term3*sr7
                  sr3 = dmpe(3) - rr3 + sr3
                  sr5 = dmpe(5) - rr5 + sr5
                  sr7 = dmpe(7) - rr7 + sr7
                  e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  efull = uir*(corek*rr3+valk*rr3k)
     &                       - ukr*(corei*rr3+vali*rr3i)
     &                       + diu*rr3i + dku*rr3k
     &                       + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                       - dkr*uir*rr5k - dir*ukr*rr5i
     &                       + qkr*uir*rr7k - qir*ukr*rr7i
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpe(3) - rr3 + rr3i
                  rr5i = dmpe(5) - rr5 + rr5i
                  rr7i = dmpe(7) - rr7 + rr7i
                  rr3k = dmpe(3) - rr3 + rr3k
                  rr5k = dmpe(5) - rr5 + rr5k
                  rr7k = dmpe(7) - rr7 + rr7k
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
               aep(i) = aep(i) + 0.5d0*e
               aep(k) = aep(k) + 0.5d0*e
               if (efull .ne. 0.0d0) then
                  nep = nep + 1
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 10.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
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
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
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
            if (use_chgpen) then
               corei = pcore(i)
               vali = pval(i)
               alphai = palpha(i)
            end if
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = i, npole
               k = ipole(kk)
               do jcell = 2, ncell
                  xr = x(k) - xi
                  yr = y(k) - yi
                  zr = z(k) - zi
                  call imager (xr,yr,zr,jcell)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2)) then
                     pscale(k) = 1.0d0
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
                     diu = dix*ukx + diy*uky + diz*ukz
                     qiu = qix*ukx + qiy*uky + qiz*ukz
                     uir = uix*xr + uiy*yr + uiz*zr
                     dku = dkx*uix + dky*uiy + dkz*uiz
                     qku = qkx*uix + qky*uiy + qkz*uiz
                     ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate real space Ewald error function damping
c
                     call dampewald (7,r,r2,f,dmpe)
c
c     find the energy value for Thole polarization damping
c
                     if (use_thole) then
                        call damptholed (i,k,7,r,dmpik)
                        scalek = pscale(k)
                        rr3 = f / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        sr3 = scalek * dmpik(3) * rr3
                        sr5 = scalek * dmpik(5) * rr5
                        sr7 = scalek * dmpik(7) * rr7
                        term1 = ck*uir - ci*ukr + diu + dku
                        term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                        term3 = uir*qkr - ukr*qir
                        efull = term1*sr3 + term2*sr5 + term3*sr7
                        sr3 = dmpe(3) - rr3 + sr3
                        sr5 = dmpe(5) - rr5 + sr5
                        sr7 = dmpe(7) - rr7 + sr7
                        e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(k)
                        valk = pval(k)
                        alphak = palpha(k)
                        call dampdir (r,alphai,alphak,dmpi,dmpk)
                        scalek = pscale(k)
                        rr3 = f * scalek / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = dmpi(3) * rr3
                        rr5i = dmpi(5) * rr5
                        rr7i = dmpi(7) * rr7
                        rr3k = dmpk(3) * rr3
                        rr5k = dmpk(5) * rr5
                        rr7k = dmpk(7) * rr7
                        efull = uir*(corek*rr3+valk*rr3k)
     &                             - ukr*(corei*rr3+vali*rr3i)
     &                             + diu*rr3i + dku*rr3k
     &                             + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                             - dkr*uir*rr5k - dir*ukr*rr5i
     &                             + qkr*uir*rr7k - qir*ukr*rr7i
                        rr3 = f / (r*r2)
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr3i = dmpe(3) - rr3 + rr3i
                        rr5i = dmpe(5) - rr5 + rr5i
                        rr7i = dmpe(7) - rr7 + rr7i
                        rr3k = dmpe(3) - rr3 + rr3k
                        rr5k = dmpe(5) - rr5 + rr5k
                        rr7k = dmpe(7) - rr7 + rr7k
                        rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                        e = uir*(corek*rr3+valk*rr3k)
     &                         - ukr*(corei*rr3+vali*rr3i)
     &                         + diu*rr3i + dku*rr3k
     &                         + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                         - dkr*uir*rr5k - dir*ukr*rr5i
     &                         + qkr*uir*rr7k - qir*ukr*rr7i
                     end if
c
c     compute the energy contribution for this interaction
c
                     if (i .eq. k) then
                        e = 0.5d0 * e
                        efull = 0.5d0 * efull
                     end if
                     ep = ep + e
                     aep(i) = aep(i) + 0.5d0*e
                     aep(k) = aep(k) + 0.5d0*e
                     if (efull .ne. 0.0d0) then
                        nep = nep + 1
                        if (molcule(i) .ne. molcule(k)) then
                           einter = einter + efull
                        end if
                     end if
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 10.0d0)
                     if ((debug.and.efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Dipole Polarization',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,efull
   50                   format (' Polar',5x,2(i7,'-',a3),9x,
     &                             f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     increment polarization energy due to external field
c
      if (use_exfld) then
         do i = 1, npole
            e = 0.0d0
            do j = 1, 3
               e = e - f*uind(j,i)*texfld(j)
            end do
            ep = ep + e
            nep = nep + 1
            aep(i) = aep(i) + e
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar3d  --  Ewald polarization analysis; list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar3d" calculates the polarization energy and analysis with
c     respect to Cartesian coordinates using particle mesh Ewald and
c     a neighbor list
c
c
      subroutine epolar3d
      use action
      use analyz
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
c
c
c     zero out the dipole polarization energy and components
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole ('MPOLE')
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     compute the real space part of the Ewald summation
c
      call epreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call eprecip3
c
c     compute the Ewald self-energy term over all the atoms
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * term * uii / 3.0d0
         ep = ep + e
         nep = nep + 1
         aep(i) = aep(i) + e
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
         do ii = 1, npole
            i = ipole(ii)
            xd = xd + rpole(2,i) + rpole(1,i)*x(i)
            yd = yd + rpole(3,i) + rpole(1,i)*y(i)
            zd = zd + rpole(4,i) + rpole(1,i)*z(i)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         nep = nep + 1
         do ii = 1, npole
            i = ipole(ii)
            aep(i) = aep(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epreal3d  --  real space polar analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epreal3d" calculates the induced dipole polarization energy
c     and analysis using particle mesh Ewald and a neighbor list
c
c
      subroutine epreal3d
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use extfld
      use inform
      use inter
      use iounit
      use math
      use mplpot
      use molcul
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
      real*8 e,efull
      real*8 f,scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 sr3,sr5,sr7
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 dir,diu,qiu,uir
      real*8 dkr,dku,qku,ukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(7),dmpe(7)
      real*8, allocatable :: pscale(:)
      logical header,huge
      character*6 mode
c
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
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,rpole,uind,x,y,z,pcore,pval,palpha,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,
!$OMP& p2scale,p3scale,p4scale,p5scale,p2iscale,p3iscale,p4iscale,
!$OMP& p5iscale,nelst,elst,use_thole,use_chgpen,use_bounds,off2,f,
!$OMP& texfld,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(pscale) shared (ep,nep,aep,einter)
!$OMP DO reduction(+:ep,nep,aep,einter) schedule(guided)
c
c     compute the dipole polarization energy component
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(i)
            vali = pval(i)
            alphai = palpha(i)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            do k = 1, np11(i)
               if (i12(j,i) .eq. ip11(k,i))
     &            pscale(i12(j,i)) = p2iscale
            end do
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            do k = 1, np11(i)
               if (i13(j,i) .eq. ip11(k,i))
     &            pscale(i13(j,i)) = p3iscale
            end do
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            do k = 1, np11(i)
               if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4iscale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
            do k = 1, np11(i)
               if (i15(j,i) .eq. ip11(k,i))
     &            pscale(i15(j,i)) = p5iscale
            end do
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
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
               diu = dix*ukx + diy*uky + diz*ukz
               qiu = qix*ukx + qiy*uky + qiz*ukz
               uir = uix*xr + uiy*yr + uiz*zr
               dku = dkx*uix + dky*uiy + dkz*uiz
               qku = qkx*uix + qky*uiy + qkz*uiz
               ukr = ukx*xr + uky*yr + ukz*zr
c
c     calculate real space Ewald error function damping
c
               call dampewald (7,r,r2,f,dmpe)
c
c     find the energy value for Thole polarization damping
c
               if (use_thole) then
                  call damptholed (i,k,7,r,dmpik)
                  scalek = pscale(k)
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  sr3 = scalek * dmpik(3) * rr3
                  sr5 = scalek * dmpik(5) * rr5
                  sr7 = scalek * dmpik(7) * rr7
                  term1 = ck*uir - ci*ukr + diu + dku
                  term2 = 2.0d0*(qiu-qku) - uir*dkr - dir*ukr
                  term3 = uir*qkr - ukr*qir
                  efull = term1*sr3 + term2*sr5 + term3*sr7
                  sr3 = dmpe(3) - rr3 + sr3
                  sr5 = dmpe(5) - rr5 + sr5
                  sr7 = dmpe(7) - rr7 + sr7
                  e = term1*sr3 + term2*sr5 + term3*sr7
c
c     find the energy value for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(k)
                  valk = pval(k)
                  alphak = palpha(k)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = pscale(k)
                  rr3 = f * scalek / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpi(3) * rr3
                  rr5i = dmpi(5) * rr5
                  rr7i = dmpi(7) * rr7
                  rr3k = dmpk(3) * rr3
                  rr5k = dmpk(5) * rr5
                  rr7k = dmpk(7) * rr7
                  efull = uir*(corek*rr3+valk*rr3k)
     &                       - ukr*(corei*rr3+vali*rr3i)
     &                       + diu*rr3i + dku*rr3k
     &                       + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                       - dkr*uir*rr5k - dir*ukr*rr5i
     &                       + qkr*uir*rr7k - qir*ukr*rr7i
                  rr3 = f / (r*r2)
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr3i = dmpe(3) - rr3 + rr3i
                  rr5i = dmpe(5) - rr5 + rr5i
                  rr7i = dmpe(7) - rr7 + rr7i
                  rr3k = dmpe(3) - rr3 + rr3k
                  rr5k = dmpe(5) - rr5 + rr5k
                  rr7k = dmpe(7) - rr7 + rr7k
                  rr3 = dmpe(3) - (1.0d0-scalek)*rr3
                  e = uir*(corek*rr3+valk*rr3k)
     &                   - ukr*(corei*rr3+vali*rr3i)
     &                   + diu*rr3i + dku*rr3k
     &                   + 2.0d0*(qiu*rr5i-qku*rr5k)
     &                   - dkr*uir*rr5k - dir*ukr*rr5i
     &                   + qkr*uir*rr7k - qir*ukr*rr7i
               end if
c
c     compute the energy contribution for this interaction
c
               ep = ep + e
               aep(i) = aep(i) + 0.5d0*e
               aep(k) = aep(k) + 0.5d0*e
               if (efull .ne. 0.0d0) then
                  nep = nep + 1
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     print message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 10.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Dipole Polarization',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Polar',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
c
c     increment polarization energy due to external field
c
      if (use_exfld) then
!$OMP    DO reduction(+:ep,nep,aep) schedule(guided)
         do i = 1, npole
            e = 0.0d0
            do j = 1, 3
               e = e - f*uind(j,i)*texfld(j)
            end do
            ep = ep + e
            nep = nep + 1
            aep(i) = aep(i) + e
         end do
!$OMP    END DO
      end if
c
c     OpenMP directives for the major loop structure
c
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
c     ##  subroutine epolar3e  --  single-loop polarization analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epreal3e" calculates the induced dipole polarization energy
c     and analysis from the induced dipoles times the electric field
c
c
      subroutine epolar3e
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
      use polpot
      use potent
      use units
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
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
      if (.not. use_mpole)  call rotpole ('MPOLE')
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
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
c     set the energy unit conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private (i,j,fi,e,ii,huge)
!$OMP DO reduction(+:ep,nep,aep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do ii = 1, npole
         i = ipole(ii)
         if (douind(i)) then
            fi = f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*uind(j,i)*udirp(j,i)
            end do
            nep = nep + 1
            ep = ep + e
            aep(i) = aep(i) + e
c
c     print a message if the energy for this site is large
c
            huge = (abs(e) .gt. 10.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Dipole Polarization',
     &                       ' Interactions :',
     &                    //,' Type',9x,'Atom Name',24x,'Alpha',
     &                       8x,'Energy',/)
               end if
               write (iout,30)  i,name(i),polarity(i),e
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
            do ii = 1, npole
               i = ipole(ii)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               uix = uind(1,i)
               uiy = uind(2,i)
               uiz = uind(3,i)
               xd = xd + dix + rpole(1,i)*x(i)
               yd = yd + diy + rpole(1,i)*y(i)
               zd = zd + diz + rpole(1,i)*z(i)
               xu = xu + uix
               yu = yu + uiy
               zu = zu + uiz
            end do
            term = (2.0d0/3.0d0) * f * (pi/volbox)
            e = term * (xd*xu+yd*yu+zd*zu)
            nep = nep + 1
            ep = ep + e
            do ii = 1, npole
               i = ipole(ii)
               aep(i) = aep(i) + e/dble(npole)
            end do
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine eprecip3  --  PME recip polarization analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "eprecip3" evaluates the reciprocal space portion of particle
c     mesh Ewald summation energy due to dipole polarization, and
c     partitions the energy among the atoms
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine eprecip3
      use analyz
      use atoms
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use polar
      use polpot
      use potent
      implicit none
      integer i,j,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 f,h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 a(3,3)
      real*8, allocatable :: fuind(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = 0.5d0 * electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (.not.use_mpole .or. aewald.ne.aeewald) then
         if (allocated(cmp)) then
            if (size(cmp) .lt. 10*n)  deallocate (cmp)
         end if
         if (allocated(fmp)) then
            if (size(fmp) .lt. 10*n)  deallocate (fmp)
         end if
         if (allocated(fphi)) then
            if (size(fphi) .lt. 20*n)  deallocate (fphi)
         end if
         if (.not. allocated(cmp))  allocate (cmp(10,n))
         if (.not. allocated(fmp))  allocate (fmp(10,n))
         if (.not. allocated(fphi))  allocate (fphi(20,n))
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
c     assign only the permanent multipoles to the PME grid
c     and perform the 3-D FFT forward transformation
c
         do ii = 1, npole
            i = ipole(ii)
            cmp(1,i) = rpole(1,i)
            cmp(2,i) = rpole(2,i)
            cmp(3,i) = rpole(3,i)
            cmp(4,i) = rpole(4,i)
            cmp(5,i) = rpole(5,i)
            cmp(6,i) = rpole(9,i)
            cmp(7,i) = rpole(13,i)
            cmp(8,i) = 2.0d0 * rpole(6,i)
            cmp(9,i) = 2.0d0 * rpole(7,i)
            cmp(10,i) = 2.0d0 * rpole(10,i)
         end do
         call cmp_to_fmp (cmp,fmp)
         call grid_mpole (fmp)
         call fftfront
c
c     make the scalar summation over reciprocal lattice
c
         pterm = (pi/aewald)**2
         volterm = pi * volbox
         nf1 = (nfft1+1) / 2
         nf2 = (nfft2+1) / 2
         nf3 = (nfft3+1) / 2
         nff = nfft1 * nfft2
         ntot = nff * nfft3
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
            end if
            qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
            qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
         end do
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
c     perform 3-D FFT backward transform and get potential
c
         call fftback
         call fphi_mpole (fphi)
      end if
c
c     set matrix for Cartesian to fractional induced dipoles
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,n))
c
c     increment the induced dipole polarization energy
c
      e = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            fuind(j,i) = a(j,1)*uind(1,i) + a(j,2)*uind(2,i)
     &                      + a(j,3)*uind(3,i)
            term = f * fuind(j,i) * fphi(j+1,i)
            e = e + term
            aep(i) = aep(i) + term
         end do
      end do
      ep = ep + e
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      return
      end
