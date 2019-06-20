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
      use atoms
      use inform
      use iounit
      use molcul
      use mpole
      use polar
      use polpot
      use potent
      implicit none
      integer i
      real*8 addu,malpha
      real*8 external
      real*8 exfield(3)
      real*8 umol(3)
      real*8 umol0(3)
      real*8 dalpha(3)
      real*8 alpha(3,3)
      real*8 valpha(3,3)
      character*40 fstr
c
c
c     get the coordinates and required force field parameters
c
      call initial
      call getxyz
      call attach
      call field
      call molecule
      call katom
      call kmpole
      call kpolar
      call mutate
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
      fstr = ' Additive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Additive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,20)  fstr(1:36),addu
   20    format (/,a36,f19.8)
      else if (digits .ge. 6) then
         write (iout,30)  fstr(1:36),addu
   30    format (/,a36,f17.6)
      else
         write (iout,40)  fstr(1:36),addu
   40    format (/,a36,f15.4)
      end if
c
c     find induced dipoles in absence of an external field
c
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      call moluind (exfield,umol0)
c
c     compute each column of the polarizability tensor
c
      external = 0.01d0
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(1) = external
      call moluind (exfield,umol)
      alpha(1,1) = (umol(1)-umol0(1)) / exfield(1)
      alpha(2,1) = (umol(2)-umol0(2)) / exfield(1)
      alpha(3,1) = (umol(3)-umol0(3)) / exfield(1)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(2) = external
      call moluind (exfield,umol)
      alpha(1,2) = (umol(1)-umol0(1)) / exfield(2)
      alpha(2,2) = (umol(2)-umol0(2)) / exfield(2)
      alpha(3,2) = (umol(3)-umol0(3)) / exfield(2)
      do i = 1, 3
         exfield(i) = 0.0d0
      end do
      exfield(3) = external
      call moluind (exfield,umol)
      alpha(1,3) = (umol(1)-umol0(1)) / exfield(3)
      alpha(2,3) = (umol(2)-umol0(2)) / exfield(3)
      alpha(3,3) = (umol(3)-umol0(3)) / exfield(3)
c
c     print out the full polarizability tensor
c
      fstr = ' Total Polarizability Tensor :    '
      if (nmol .eq. 1)  fstr = ' Molecular Polarizability Tensor :'
      write (iout,50)  fstr(1:34)
   50 format (/,a34,/)
      if (digits .ge. 8) then
         write (iout,60)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   60    format (15x,3f16.8,/,15x,3f16.8,/,15x,3f16.8)
      else if (digits .ge. 6) then
         write (iout,70)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   70    format (15x,3f14.6,/,15x,3f14.6,/,15x,3f14.6)
      else
         write (iout,80)  alpha(1,1),alpha(1,2),alpha(1,3),
     &                    alpha(2,1),alpha(2,2),alpha(2,3),
     &                    alpha(3,1),alpha(3,2),alpha(3,3)
   80    format (15x,3f12.4,/,15x,3f12.4,/,15x,3f12.4)
      end if
c
c     diagonalize the tensor and get molecular polarizability
c
      call jacobi (3,alpha,dalpha,valpha)
      if (nmol .eq. 1)  fstr = ' Polarizability Tensor Eigenvalues :'
      write (iout,90)  fstr(1:36)
   90 format (/,a36,/)
      if (digits .ge. 8) then
         write (iout,100)  dalpha(1),dalpha(2),dalpha(3)
  100    format (15x,3f16.8)
      else if (digits .ge. 6) then
         write (iout,110)  dalpha(1),dalpha(2),dalpha(3)
  110    format (15x,3f14.6)
      else
         write (iout,120)  dalpha(1),dalpha(2),dalpha(3)
  120    format (15x,3f12.4)
      end if
      malpha = (dalpha(1)+dalpha(2)+dalpha(3)) / 3.0d0
      fstr = ' Interactive Total Polarizability :    '
      if (nmol .eq. 1)  fstr = ' Interactive Molecular Polarizability :'
      if (digits .ge. 8) then
         write (iout,130)  fstr(1:39),malpha
  130    format (/,a39,f16.8)
      else if (digits .ge. 6) then
         write (iout,140)  fstr(1:39),malpha
  140    format (/,a39,f14.6)
      else
         write (iout,150)  fstr(1:39),malpha
  150    format (/,a39,f12.4)
      end if
c
c     perform any final tasks before program exit
c
      call final
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
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polpot
      use units
      implicit none
      integer i,j,iter
      integer maxiter
      real*8 eps,epsold
      real*8 polmin
      real*8 a,b,sum
      real*8 umol(3)
      real*8 exfield(3)
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: vec(:,:)
      logical done
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (poli(npole))
      allocate (field(3,npole))
      allocate (rsd(3,npole))
      allocate (zrsd(3,npole))
      allocate (conj(3,npole))
      allocate (vec(3,npole))
c
c     set induced dipoles to polarizability times direct field
c
      call dfield (field)
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarity(i) * field(j,i)
         end do
      end do
c
c     increment induced dipoles to account for external field
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = uind(j,i) + polarity(i)*exfield(j)
         end do
      end do
c
c     compute mutual induced dipole moments via CG algorithm
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         eps = 100.0d0
         polmin = 0.00000001d0
         call ufield (field)
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = field(j,i)
               zrsd(j,i) = rsd(j,i) * poli(i)
               conj(j,i) = zrsd(j,i)
            end do
         end do
c
c     iterate the mutual induced dipoles and check convergence
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  uind(j,i) = conj(j,i)
               end do
            end do
            call ufield (field)
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
               end do
            end do
            a = 0.0d0
            sum = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
               end do
            end do
            b = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  b = b + rsd(j,i)*zrsd(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            eps = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  eps = eps + rsd(j,i)*rsd(j,i)
               end do
            end do
            eps = debye * sqrt(eps/dble(npolar))
            epsold = eps
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debye)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. maxiter)  done = .true.
         end do
c
c     print a warning if induced dipoles failed to converge
c
         if (eps .gt. poleps) then
            write (iout,30)
   30       format (/,' MOLUIND  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
         end if
      end if
c
c     sum up the total molecular induced dipole components
c
      do j = 1, 3
         umol(j) = 0.0d0
      end do
      do i = 1, npole
         umol(1) = umol(1) + uind(1,i)
         umol(2) = umol(2) + uind(2,i)
         umol(3) = umol(3) + uind(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (poli)
      deallocate (field)
      deallocate (rsd)
      deallocate (zrsd)
      deallocate (conj)
      deallocate (vec)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dfield  --  electric field from induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dfield" finds the field at each polarizable site due to the
c     induced dipoles at the other polarizable sites
c
c
      subroutine dfield (field)
      use atoms
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      character*6 mode
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     find the electrostatic field due to permanent multipoles
c
      do ii = 1, npole-1
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
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               dscale(i12(j,i)) = pscale(i12(j,i))
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
            end do
         else
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
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
            end do
         end if
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            r2 = xr*xr + yr* yr + zr*zr
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
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               damp = pdi * pdamp(kk)
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(kk))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
               rr3 = scale3 / (r*r2)
               rr5 = 3.0d0 * scale5 / (r*r2*r2)
               rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
               fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dkx + 2.0d0*rr5*qkx
               fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dky + 2.0d0*rr5*qky
               fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dkz + 2.0d0*rr5*qkz
               fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*dix - 2.0d0*rr5*qix
               fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*diy - 2.0d0*rr5*qiy
               fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*diz - 2.0d0*rr5*qiz
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kk)
               valk = pval(kk)
               alphak = palpha(kk)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               rr3 = 1.0d0 / (r*r2)
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr3i = dmpi(3) * rr3
               rr5i = dmpi(5) * rr5
               rr7i = dmpi(7) * rr7
               rr3k = dmpk(3) * rr3
               rr5k = dmpk(5) * rr5
               rr7k = dmpk(7) * rr7
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx        
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
            end if
            do j = 1, 3
               field(j,ii) = field(j,ii) + fid(j)*dscale(k)
               field(j,kk) = field(j,kk) + fkd(j)*dscale(k)
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
            end do
         else
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
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
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
c     induced dipoles at the other polarizable sites
c
c
      subroutine ufield (field)
      use atoms
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 r,r2,rr3,rr5
      real*8 xr,yr,zr
      real*8 uix,uiy,uiz
      real*8 ukx,uky,ukz
      real*8 uir,ukr
      real*8 pdi,pti
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 pgamma,damp
      real*8 expdamp
      real*8 scale3,scale5
      real*8 fi(3),fk(3)
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 field(3,*)
c
c
c     zero out the value of the electric field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
c
c     initialize connected atom exclusion coefficients
c
      do i = 1, n
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
      end do
c
c     loop over pairs of sites incrementing the electric field
c
      do ii = 1, npole-1
         i = ipole(ii)
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
         else if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
         do j = 1, np11(i)
            uscale(ip11(j,i)) = u1scale
         end do
         do j = 1, np12(i)
            uscale(ip12(j,i)) = u2scale
         end do
         do j = 1, np13(i)
            uscale(ip13(j,i)) = u3scale
         end do
         do j = 1, np14(i)
            uscale(ip14(j,i)) = u4scale
         end do
         do j = 1, n12(i)
            wscale(i12(j,i)) = w2scale
         end do
         do j = 1, n13(i)
            wscale(i13(j,i)) = w3scale
         end do
         do j = 1, n14(i)
            wscale(i14(j,i)) = w4scale
         end do
         do j = 1, n15(i)
            wscale(i15(j,i)) = w5scale
         end do
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            r2 = xr*xr + yr* yr + zr*zr
            ukx = uind(1,kk)
            uky = uind(2,kk)
            ukz = uind(3,kk)
            r = sqrt(r2)
            uir = uix*xr + uiy*yr + uiz*zr
            ukr = ukx*xr + uky*yr + ukz*zr
c
c     adjust the field to account for polarization damping
c
            if (use_thole) then
               scale3 = uscale(k)
               scale5 = uscale(k)
               damp = pdi * pdamp(kk)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(kk))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                  end if
               end if
            else if (use_chgpen) then
               corek = pcore(kk)
               valk = pval(kk)
               alphak = palpha(kk)
               call dampmut (r,alphai,alphak,dmpik)
               scale3 = wscale(k) * dmpik(3)
               scale5 = wscale(k) * dmpik(5)
            end if
            rr3 = scale3 / (r*r2)
            rr5 = 3.0d0 * scale5 / (r*r2*r2)
            fi(1) = -rr3*ukx + rr5*ukr*xr
            fi(2) = -rr3*uky + rr5*ukr*yr
            fi(3) = -rr3*ukz + rr5*ukr*zr
            fk(1) = -rr3*uix + rr5*uir*xr
            fk(2) = -rr3*uiy + rr5*uir*yr
            fk(3) = -rr3*uiz + rr5*uir*zr
c
c     increment the field at each site due to this interaction
c
            do j = 1, 3
               field(j,ii) = field(j,ii) + fi(j)
               field(j,kk) = field(j,kk) + fk(j)
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         do j = 1, np11(i)
            uscale(ip11(j,i)) = 1.0d0
         end do
         do j = 1, np12(i)
            uscale(ip12(j,i)) = 1.0d0
         end do
         do j = 1, np13(i)
            uscale(ip13(j,i)) = 1.0d0
         end do
         do j = 1, np14(i)
            uscale(ip14(j,i)) = 1.0d0
         end do
         do j = 1, n12(i)
            wscale(i12(j,i)) = 1.0d0
         end do
         do j = 1, n13(i)
            wscale(i13(j,i)) = 1.0d0
         end do
         do j = 1, n14(i)
            wscale(i14(j,i)) = 1.0d0
         end do
         do j = 1, n15(i)
            wscale(i15(j,i)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
      return
      end
