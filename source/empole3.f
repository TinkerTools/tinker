c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3  --  atomic multipole energy & analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3" calculates the electrostatic energy due to atomic
c     multipole interactions, and partitions the energy among atoms
c
c
      subroutine empole3
      use limits
      implicit none
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
c     "empole3a" calculates the atomic multipole interaction energy
c     using a double loop, and partitions the energy among atoms
c
c
      subroutine empole3a
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
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 rr1i,rr3i,rr5i
      real*8 rr1k,rr3k,rr5k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
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
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9)
      real*8, allocatable :: mscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
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
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole-1
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = abs(yaxis(kk))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
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
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(k) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     find damped multipole intermediates and energy value
c
                  if (use_chgpen) then
                     corek = pcore(kk)
                     valk = pval(kk)
                     alphak = palpha(kk)
                     term1 = corei*corek
                     term1i = corek*vali
                     term2i = corek*dir
                     term3i = corek*qir
                     term1k = corei*valk
                     term2k = -corei*dkr
                     term3k = corei*qkr
                     term1ik = vali*valk
                     term2ik = valk*dir - vali*dkr + dik
                     term3ik = vali*qkr + valk*qir - dir*dkr
     &                            + 2.0d0*(dkqi-diqk+qiqk)
                     term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                     term5ik = qir*qkr
                     call damppole (r,9,alphai,alphak,
     &                               dmpi,dmpk,dmpik)
                     rr1i = dmpi(1)*rr1
                     rr3i = dmpi(3)*rr3
                     rr5i = dmpi(5)*rr5
                     rr1k = dmpk(1)*rr1
                     rr3k = dmpk(3)*rr3
                     rr5k = dmpk(5)*rr5
                     rr1ik = dmpik(1)*rr1
                     rr3ik = dmpik(3)*rr3
                     rr5ik = dmpik(5)*rr5
                     rr7ik = dmpik(7)*rr7
                     rr9ik = dmpik(9)*rr9
                     e = term1*rr1 + term1i*rr1i
     &                      + term1k*rr1k + term1ik*rr1ik
     &                      + term2i*rr3i + term2k*rr3k
     &                      + term2ik*rr3ik + term3i*rr5i
     &                      + term3k*rr5k + term3ik*rr5ik
     &                      + term4ik*rr7ik + term5ik*rr9ik
c
c     find standard multipole intermediates and energy value
c
                  else
                     term1 = ci*ck
                     term2 = ck*dir - ci*dkr + dik
                     term3 = ci*qkr + ck*qir - dir*dkr
     &                          + 2.0d0*(dkqi-diqk+qiqk)
                     term4 = dir*qkr - dkr*qir - 4.0d0*qik
                     term5 = qir*qkr
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                  end if
c
c     increment the overall multipole energy components
c
                  if (use_group)  e = e * fgrp
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(i) = aem(i) + 0.5d0*e
                     aem(k) = aem(k) + 0.5d0*e
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Mpole',5x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            iz = zaxis(ii)
            ix = xaxis(ii)
            iy = abs(yaxis(ii))
            xi = x(i)
            yi = y(i)
            zi = z(i)
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
            if (use_chgpen) then
               corei = pcore(ii)
               vali = pval(ii)
               alphai = palpha(ii)
            end if
            usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = m2scale
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = m3scale
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = m4scale
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = m5scale
            end do
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               kz = zaxis(kk)
               kx = xaxis(kk)
               ky = abs(yaxis(kk))
               usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
               if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
               proceed = .true.
               if (proceed)  proceed = (usei .or. usek)
               if (proceed) then
                  do j = 2, ncell
                     xr = x(k) - xi
                     yr = y(k) - yi
                     zr = z(k) - zi
                     call imager (xr,yr,zr,j)
                     r2 = xr*xr + yr* yr + zr*zr
                     if (.not. (use_polymer .and. r2.le.polycut2))
     &                  mscale(k) = 1.0d0
                     if (r2 .le. off2) then
                        r = sqrt(r2)
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
c     get reciprocal distance terms for this interaction
c
                        rr1 = f * mscale(k) / r
                        rr3 = rr1 / r2
                        rr5 = 3.0d0 * rr3 / r2
                        rr7 = 5.0d0 * rr5 / r2
                        rr9 = 7.0d0 * rr7 / r2
c
c     find damped multipole intermediates and energy value
c
                        if (use_chgpen) then
                           corek = pcore(kk)
                           valk = pval(kk)
                           alphak = palpha(kk)
                           term1 = corei*corek
                           term1i = corek*vali
                           term2i = corek*dir
                           term3i = corek*qir
                           term1k = corei*valk
                           term2k = -corei*dkr
                           term3k = corei*qkr
                           term1ik = vali*valk
                           term2ik = valk*dir - vali*dkr + dik
                           term3ik = vali*qkr + valk*qir - dir*dkr
     &                                  + 2.0d0*(dkqi-diqk+qiqk)
                           term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                           term5ik = qir*qkr
                           call damppole (r,9,alphai,alphak,
     &                                     dmpi,dmpk,dmpik)
                           rr1i = dmpi(1)*rr1
                           rr3i = dmpi(3)*rr3
                           rr5i = dmpi(5)*rr5
                           rr1k = dmpk(1)*rr1
                           rr3k = dmpk(3)*rr3
                           rr5k = dmpk(5)*rr5
                           rr1ik = dmpik(1)*rr1
                           rr3ik = dmpik(3)*rr3
                           rr5ik = dmpik(5)*rr5
                           rr7ik = dmpik(7)*rr7
                           rr9ik = dmpik(9)*rr9
                           e = term1*rr1 + term1i*rr1i
     &                            + term1k*rr1k + term1ik*rr1ik
     &                            + term2i*rr3i + term2k*rr3k
     &                            + term2ik*rr3ik + term3i*rr5i
     &                            + term3k*rr5k + term3ik*rr5ik
     &                            + term4ik*rr7ik + term5ik*rr9ik
c
c     find standard multipole intermediates and energy value
c
                        else
                           term1 = ci*ck
                           term2 = ck*dir - ci*dkr + dik
                           term3 = ci*qkr + ck*qir - dir*dkr
     &                                + 2.0d0*(dkqi-diqk+qiqk)
                           term4 = dir*qkr - dkr*qir - 4.0d0*qik
                           term5 = qir*qkr
                           e = term1*rr1 + term2*rr3 + term3*rr5
     &                            + term4*rr7 + term5*rr9
                        end if
c
c     increment the overall multipole energy components
c
                        if (use_group)  e = e * fgrp
                        if (i .eq. k)  e = 0.5d0 * e
                        if (e .ne. 0.0d0) then
                           nem = nem + 1
                           em = em + e
                           aem(i) = aem(i) + 0.5d0*e
                           aem(k) = aem(k) + 0.5d0*e
                           einter = einter + e
                        end if
c
c     print message if the energy of this interaction is large
c
                        huge = (abs(e) .gt. 100.0d0)
                        if ((debug.and.e.ne.0.0d0)
     &                        .or. (verbose.and.huge)) then
                           if (header) then
                              header = .false.
                              write (iout,40)
   40                         format (/,' Individual Atomic Multipole',
     &                                   ' Interactions :',
     &                                //,' Type',14x,'Atom Names',
     &                                   15x,'Distance',8x,'Energy',/)
                           end if
                           write (iout,50)  i,name(i),k,name(k),r,e
   50                      format (' Mpole',5x,2(i7,'-',a3),1x,
     &                                '(XTAL)',2x,f10.4,2x,f12.4)
                        end if
                     end if
                  end do
               end if
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
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
c     "empole3b" calculates the atomic multipole interaction energy
c     using a neighbor list, and partitions the energy among the atoms
c
c
      subroutine empole3b
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use group
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      use usage
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,f,fgrp
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 rr1i,rr3i,rr5i
      real*8 rr1k,rr3k,rr5k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
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
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9)
      real*8, allocatable :: mscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      character*6 mode
c
c
c     zero out total atomic multipole energy and partitioning
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
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
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,xaxis,yaxis,zaxis,rpole,pcore,pval,
!$OMP& palpha,use,n12,i12,n13,i13,n14,i14,n15,i15,m2scale,m3scale,
!$OMP& m4scale,m5scale,f,nelst,elst,use_chgpen,use_group,use_intra,
!$OMP& use_bounds,off2,molcule,name,verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,nem,aem,einter)
!$OMP DO reduction(+:em,nem,aem,einter) schedule(guided)
c
c     calculate the multipole interaction energy term
c
      do ii = 1, npole
         i = ipole(ii)
         iz = zaxis(ii)
         ix = xaxis(ii)
         iy = abs(yaxis(ii))
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
         usei = (use(i) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
         end do
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            kz = zaxis(kk)
            kx = xaxis(kk)
            ky = abs(yaxis(kk))
            usek = (use(k) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (proceed) then
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               if (use_bounds)  call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
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
c     get reciprocal distance terms for this interaction
c
                  rr1 = f * mscale(k) / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
c
c     find damped multipole intermediates and energy value
c
                  if (use_chgpen) then
                     corek = pcore(kk)
                     valk = pval(kk)
                     alphak = palpha(kk)
                     term1 = corei*corek
                     term1i = corek*vali
                     term2i = corek*dir
                     term3i = corek*qir
                     term1k = corei*valk
                     term2k = -corei*dkr
                     term3k = corei*qkr
                     term1ik = vali*valk
                     term2ik = valk*dir - vali*dkr + dik
                     term3ik = vali*qkr + valk*qir - dir*dkr
     &                            + 2.0d0*(dkqi-diqk+qiqk)
                     term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                     term5ik = qir*qkr
                     call damppole (r,9,alphai,alphak,
     &                               dmpi,dmpk,dmpik)
                     rr1i = dmpi(1)*rr1
                     rr3i = dmpi(3)*rr3
                     rr5i = dmpi(5)*rr5
                     rr1k = dmpk(1)*rr1
                     rr3k = dmpk(3)*rr3
                     rr5k = dmpk(5)*rr5
                     rr1ik = dmpik(1)*rr1
                     rr3ik = dmpik(3)*rr3
                     rr5ik = dmpik(5)*rr5
                     rr7ik = dmpik(7)*rr7
                     rr9ik = dmpik(9)*rr9
                     e = term1*rr1 + term1i*rr1i
     &                      + term1k*rr1k + term1ik*rr1ik
     &                      + term2i*rr3i + term2k*rr3k
     &                      + term2ik*rr3ik + term3i*rr5i
     &                      + term3k*rr5k + term3ik*rr5ik
     &                      + term4ik*rr7ik + term5ik*rr9ik
c
c     find standard multipole intermediates and energy value
c
                  else
                     term1 = ci*ck
                     term2 = ck*dir - ci*dkr + dik
                     term3 = ci*qkr + ck*qir - dir*dkr
     &                          + 2.0d0*(dkqi-diqk+qiqk)
                     term4 = dir*qkr - dkr*qir - 4.0d0*qik
                     term5 = qir*qkr
                     e = term1*rr1 + term2*rr3 + term3*rr5
     &                      + term4*rr7 + term5*rr9
                  end if
c
c     increment the overall multipole energy components
c
                  if (use_group)  e = e * fgrp
                  if (e .ne. 0.0d0) then
                     nem = nem + 1
                     em = em + e
                     aem(i) = aem(i) + 0.5d0*e
                     aem(k) = aem(k) + 0.5d0*e
                     if (molcule(i) .ne. molcule(k)) then
                        einter = einter + e
                     end if
                  end if
c
c     print message if the energy of this interaction is large
c
                  huge = (abs(e) .gt. 100.0d0)
                  if ((debug.and.e.ne.0.0d0)
     &                  .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,20)
   20                   format (/,' Individual Atomic Multipole',
     &                             ' Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',8x,'Energy',/)
                     end if
                     write (iout,30)  i,name(i),k,name(k),r,e
   30                format (' Mpole',5x,2(i7,'-',a3),9x,
     &                          f10.4,2x,f12.4)
                  end if
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
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
      deallocate (mscale)
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
c     "empole3c" calculates the atomic multipole interaction energy
c     using a particle mesh Ewald summation and double loop, and
c     partitions the energy among the atoms
c
c
      subroutine empole3c
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
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
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
c     compute the real space part of the Ewald summation
c
      call emreal3c
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip3
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
         nem = nem + 1
         aem(i) = aem(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, npole
            i = ipole(ii)
            dix = rpole(2,ii)
            diy = rpole(3,ii)
            diz = rpole(4,ii)
            xd = xd + dix + rpole(1,ii)*x(i)
            yd = yd + diy + rpole(1,ii)*y(i)
            zd = zd + diz + rpole(1,ii)*z(i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do ii = 1, npole
            i = ipole(ii)
            aem(i) = aem(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3c  --  real space mpole analysis via loop  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and partitions
c     the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  [newsletters are available at
c     https://www.ccp5.ac.uk/infoweb/newsletters]
c
c
      subroutine emreal3c
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
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk
      integer jcell
      real*8 e,efull,f
      real*8 scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 rr1i,rr3i,rr5i
      real*8 rr1k,rr3k,rr5k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
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
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole-1
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
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
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
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
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,9,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  rr1i = dmpi(1)*rr1
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr1k = dmpk(1)*rr1
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr1ik = dmpik(1)*rr1
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if
c
c     compute the full undamped energy for this interaction
c
               efull = mscale(k) * e
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     compute the energy contribution for this interaction
c
               if (use_chgpen) then
                  scalek = mscale(k)
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
               else
                  scalek = 1.0d0 - mscale(k)
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if
c
c     increment the overall multipole energy component
c
               em = em + e
               aem(i) = aem(i) + 0.5d0*e
               aem(k) = aem(k) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Mpole',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction energy with other unit cells
c
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
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
            if (use_chgpen) then
               corei = pcore(ii)
               vali = pval(ii)
               alphai = palpha(ii)
            end if
c
c     set exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = m2scale
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = m3scale
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = m4scale
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = m5scale
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
                  r2 = xr*xr + yr* yr + zr*zr
                  if (.not. (use_polymer .and. r2.le.polycut2))
     &               mscale(k) = 1.0d0
                  if (r2 .le. off2) then
                     r = sqrt(r2)
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
     &                         + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
                     rr1 = f / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
                     call dampewald (9,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
                     if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        term1 = corei*corek
                        term1i = corek*vali
                        term2i = corek*dir
                        term3i = corek*qir
                        term1k = corei*valk
                        term2k = -corei*dkr
                        term3k = corei*qkr
                        term1ik = vali*valk
                        term2ik = valk*dir - vali*dkr + dik
                        term3ik = vali*qkr + valk*qir - dir*dkr
     &                               + 2.0d0*(dkqi-diqk+qiqk)
                        term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                        term5ik = qir*qkr
                        call damppole (r,9,alphai,alphak,
     &                                  dmpi,dmpk,dmpik)
                        rr1i = dmpi(1)*rr1
                        rr3i = dmpi(3)*rr3
                        rr5i = dmpi(5)*rr5
                        rr1k = dmpk(1)*rr1
                        rr3k = dmpk(3)*rr3
                        rr5k = dmpk(5)*rr5
                        rr1ik = dmpik(1)*rr1
                        rr3ik = dmpik(3)*rr3
                        rr5ik = dmpik(5)*rr5
                        rr7ik = dmpik(7)*rr7
                        rr9ik = dmpik(9)*rr9
                        e = term1*rr1 + term1i*rr1i
     &                         + term1k*rr1k + term1ik*rr1ik
     &                         + term2i*rr3i + term2k*rr3k
     &                         + term2ik*rr3ik + term3i*rr5i
     &                         + term3k*rr5k + term3ik*rr5ik
     &                         + term4ik*rr7ik + term5ik*rr9ik
c
c     find standard multipole intermediates and energy value
c
                     else
                        term1 = ci*ck
                        term2 = ck*dir - ci*dkr + dik
                        term3 = ci*qkr + ck*qir - dir*dkr
     &                             + 2.0d0*(dkqi-diqk+qiqk)
                        term4 = dir*qkr - dkr*qir - 4.0d0*qik
                        term5 = qir*qkr
                        e = term1*rr1 + term2*rr3 + term3*rr5
     &                         + term4*rr7 + term5*rr9
                     end if
c
c     compute the full undamped energy for this interaction
c
                     efull = mscale(k) * e
                     if (efull .ne. 0.0d0) then
                        nem = nem + 1
                        einter = einter + efull
                     end if
c
c     compute the energy contribution for this interaction
c
                     if (use_chgpen) then
                        scalek = mscale(k)
                        rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                        rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                        rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                        rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                        rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                        rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                        rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                        rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                        rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                        rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                        rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                        rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                        e = term1*rr1 + term1i*rr1i
     &                         + term1k*rr1k + term1ik*rr1ik
     &                         + term2i*rr3i + term2k*rr3k
     &                         + term2ik*rr3ik + term3i*rr5i
     &                         + term3k*rr5k + term3ik*rr5ik
     &                         + term4ik*rr7ik + term5ik*rr9ik
                     else
                        scalek = 1.0d0 - mscale(k)
                        rr1 = dmpe(1) - scalek*rr1
                        rr3 = dmpe(3) - scalek*rr3
                        rr5 = dmpe(5) - scalek*rr5
                        rr7 = dmpe(7) - scalek*rr7
                        rr9 = dmpe(9) - scalek*rr9
                        e = term1*rr1 + term2*rr3 + term3*rr5
     &                               + term4*rr7 + term5*rr9
                     end if
c
c     increment the overall multipole energy component
c
                     if (i .eq. k)  e = 0.5d0 * e
                     em = em + e
                     aem(i) = aem(i) + 0.5d0*e
                     aem(k) = aem(k) + 0.5d0*e
c
c     print message if the energy of this interaction is large
c
                     huge = (abs(efull) .gt. 100.0d0)
                     if ((debug .and. efull.ne.0.0d0)
     &                     .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,40)
   40                      format (/,' Individual Atomic Multipole',
     &                                ' Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',8x,'Energy',/)
                        end if
                        write (iout,50)  i,name(i),k,name(k),r,efull
   50                   format (' Mpole',5x,2(i7,'-',a3),1x,
     &                             '(XTAL)',2x,f10.4,2x,f12.4)
                     end if
                  end if
               end do
            end do
c
c     reset exclusion coefficients for connected atoms
c
            do j = 1, n12(i)
               mscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               mscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               mscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               mscale(i15(j,i)) = 1.0d0
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
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
c     "empole3d" calculates the atomic multipole interaction energy
c     using particle mesh Ewald summation and a neighbor list, and
c     partitions the energy among the atoms
c
c
      subroutine empole3d
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
      implicit none
      integer i,ii
      real*8 e,f
      real*8 term,fterm
      real*8 cii,dii,qii
      real*8 xd,yd,zd
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      em = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bseorder
      aewald = aeewald
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
c     compute the real space part of the Ewald summation
c
      call emreal3d
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip3
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / rootpi
      do ii = 1, npole
         i = ipole(ii)
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
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
     &            + qixx*qixx + qiyy*qiyy + qizz*qizz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         em = em + e
         nem = nem + 1
         aem(i) = aem(i) + e
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         do ii = 1, npole
            i = ipole(ii)
            dix = rpole(2,ii)
            diy = rpole(3,ii)
            diz = rpole(4,ii)
            xd = xd + dix + rpole(1,ii)*x(i)
            yd = yd + diy + rpole(1,ii)*y(i)
            zd = zd + diz + rpole(1,ii)*z(i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         e = term * (xd*xd+yd*yd+zd*zd)
         em = em + e
         nem = nem + 1
         do ii = 1, npole
            i = ipole(ii)
            aem(i) = aem(i) + e/dble(npole)
         end do
      end if
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emreal3d  --  real space mpole analysis via list  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emreal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions, and partitions
c     the energy among the atoms using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  [newsletters are available at
c     https://www.ccp5.ac.uk/infoweb/newsletters]
c
c
      subroutine emreal3d
      use action
      use analyz
      use atomid
      use atoms
      use bound
      use chgpen
      use chgpot
      use couple
      use energi
      use inform
      use inter
      use iounit
      use math
      use molcul
      use mplpot
      use mpole
      use neigh
      use potent
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 e,efull,f
      real*8 scalek
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 rr1i,rr3i,rr5i
      real*8 rr1k,rr3k,rr5k
      real*8 rr1ik,rr3ik,rr5ik
      real*8 rr7ik,rr9ik
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
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 term1,term2,term3
      real*8 term4,term5
      real*8 term1i,term2i,term3i
      real*8 term1k,term2k,term3k
      real*8 term1ik,term2ik,term3ik
      real*8 term4ik,term5ik
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9),dmpe(9)
      real*8, allocatable :: mscale(:)
      logical header,huge
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         mscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug .and. npole.ne.0) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Atomic Multipole Interactions :',
     &           //,' Type',14x,'Atom Names',15x,'Distance',
     &              8x,'Energy',/)
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,x,y,z,rpole,pcore,pval,palpha,n12,i12,
!$OMP& n13,i13,n14,i14,n15,i15,m2scale,m3scale,m4scale,m5scale,
!$OMP& nelst,elst,use_chgpen,use_bounds,f,off2,molcule,name,
!$OMP& verbose,debug,header,iout)
!$OMP& firstprivate(mscale) shared (em,nem,aem,einter)
!$OMP DO reduction(+:em,nem,aem,einter) schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole
         i = ipole(ii)
         xi = x(i)
         yi = y(i)
         zi = z(i)
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
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = m2scale
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = m3scale
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = m4scale
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = m5scale
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
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
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
     &                   + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     calculate real space Ewald error function damping
c
               call dampewald (9,r,r2,f,dmpe)
c
c     find damped multipole intermediates and energy value
c
               if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  term1 = corei*corek
                  term1i = corek*vali
                  term2i = corek*dir
                  term3i = corek*qir
                  term1k = corei*valk
                  term2k = -corei*dkr
                  term3k = corei*qkr
                  term1ik = vali*valk
                  term2ik = valk*dir - vali*dkr + dik
                  term3ik = vali*qkr + valk*qir - dir*dkr
     &                         + 2.0d0*(dkqi-diqk+qiqk)
                  term4ik = dir*qkr - dkr*qir - 4.0d0*qik
                  term5ik = qir*qkr
                  call damppole (r,9,alphai,alphak,
     &                            dmpi,dmpk,dmpik)
                  rr1i = dmpi(1)*rr1
                  rr3i = dmpi(3)*rr3
                  rr5i = dmpi(5)*rr5
                  rr1k = dmpk(1)*rr1
                  rr3k = dmpk(3)*rr3
                  rr5k = dmpk(5)*rr5
                  rr1ik = dmpik(1)*rr1
                  rr3ik = dmpik(3)*rr3
                  rr5ik = dmpik(5)*rr5
                  rr7ik = dmpik(7)*rr7
                  rr9ik = dmpik(9)*rr9
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
c
c     find standard multipole intermediates and energy value
c
               else
                  term1 = ci*ck
                  term2 = ck*dir - ci*dkr + dik
                  term3 = ci*qkr + ck*qir - dir*dkr
     &                       + 2.0d0*(dkqi-diqk+qiqk)
                  term4 = dir*qkr - dkr*qir - 4.0d0*qik
                  term5 = qir*qkr
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if
c
c     compute the full undamped energy for this interaction
c
               efull = mscale(k) * e
               if (efull .ne. 0.0d0) then
                  nem = nem + 1
                  if (molcule(i) .ne. molcule(k)) then
                     einter = einter + efull
                  end if
               end if
c
c     compute the energy contribution for this interaction
c
               if (use_chgpen) then
                  scalek = mscale(k)
                  rr1i = dmpe(1) - (1.0d0-scalek*dmpi(1))*rr1
                  rr3i = dmpe(3) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = dmpe(5) - (1.0d0-scalek*dmpi(5))*rr5
                  rr1k = dmpe(1) - (1.0d0-scalek*dmpk(1))*rr1
                  rr3k = dmpe(3) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = dmpe(5) - (1.0d0-scalek*dmpk(5))*rr5
                  rr1ik = dmpe(1) - (1.0d0-scalek*dmpik(1))*rr1
                  rr3ik = dmpe(3) - (1.0d0-scalek*dmpik(3))*rr3
                  rr5ik = dmpe(5) - (1.0d0-scalek*dmpik(5))*rr5
                  rr7ik = dmpe(7) - (1.0d0-scalek*dmpik(7))*rr7
                  rr9ik = dmpe(9) - (1.0d0-scalek*dmpik(9))*rr9
                  rr1 = dmpe(1) - (1.0d0-scalek)*rr1
                  e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik
     &                   + term1i*rr1i + term1k*rr1k + term1ik*rr1ik
     &                   + term2i*rr3i + term2k*rr3k + term2ik*rr3ik
     &                   + term3i*rr5i + term3k*rr5k + term3ik*rr5ik
               else
                  scalek = 1.0d0 - mscale(k)
                  rr1 = dmpe(1) - scalek*rr1
                  rr3 = dmpe(3) - scalek*rr3
                  rr5 = dmpe(5) - scalek*rr5
                  rr7 = dmpe(7) - scalek*rr7
                  rr9 = dmpe(9) - scalek*rr9
                  e = term1*rr1 + term2*rr3 + term3*rr5
     &                   + term4*rr7 + term5*rr9
               end if
c
c     compute the energy contribution for this interaction
c
               em = em + e
               aem(i) = aem(i) + 0.5d0*e
               aem(k) = aem(k) + 0.5d0*e
c
c     print a message if the energy of this interaction is large
c
               huge = (abs(efull) .gt. 100.0d0)
               if ((debug.and.efull.ne.0.0d0)
     &               .or. (verbose.and.huge)) then
                  if (header) then
                     header = .false.
                     write (iout,20)
   20                format (/,' Individual Atomic Multipole',
     &                          ' Interactions :',
     &                       //,' Type',14x,'Atom Names',
     &                          15x,'Distance',8x,'Energy',/)
                  end if
                  write (iout,30)  i,name(i),k,name(k),r,efull
   30             format (' Mpole',5x,2(i7,'-',a3),9x,
     &                       f10.4,2x,f12.4)
               end if
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            mscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            mscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            mscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            mscale(i15(j,i)) = 1.0d0
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
      deallocate (mscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine emrecip3  --  PME reciprocal multipole analysis  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "emrecip3" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions, and
c     partitions the energy among the atoms
c
c     literature references:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     W. Smith and D. Fincham, "The Ewald Sum in Truncated Octahedral
c     and Rhombic Dodecahedral Boundary Conditions", Molecular
c     Simulation, 10, 67-71 (1993)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip3
      use analyz
      use bound
      use boxes
      use chgpot
      use energi
      use ewald
      use math
      use mpole
      use mrecip
      use pme
      use potent
      implicit none
      integer i,j,k
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 f,h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
      f = 0.5d0 * electric / dielec
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(cmp)) then
         if (size(cmp) .lt. 10*npole)  deallocate (cmp)
      end if
      if (allocated(fmp)) then
         if (size(fmp) .lt. 10*npole)  deallocate (fmp)
      end if
      if (allocated(fphi)) then
         if (size(fphi) .lt. 20*npole)  deallocate (fphi)
      end if
      if (.not. allocated(cmp))  allocate (cmp(10,npole))
      if (.not. allocated(fmp))  allocate (fmp(10,npole))
      if (.not. allocated(fphi))  allocate (fphi(20,npole))
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
c     copy the multipole moments into local storage areas
c
      do i = 1, npole
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
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
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
c
c     sum over multipoles and increment total multipole energy
c
      e = 0.0d0
      do i = 1, npole
         j = ipole(i)
         do k = 1, 10
            term = f * fmp(k,i) * fphi(k,i)
            e = e + term
            aem(j) = aem(j) + term
         end do
      end do
      em = em + e
      return
      end
