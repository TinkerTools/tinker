c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce  --  evaluate induced dipole moments  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce" computes the induced dipole moments at polarizable
c     sites due to direct or mutual polarization
c
c     assumes multipole components have already been rotated into
c     the global coordinate frame; computes induced dipoles based
c     on full system, use of active or inactive atoms is ignored
c
c
      subroutine induce
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use polpot
      use potent
      use solpot
      use units
      use uprior
      implicit none
      integer i,j,k,ii
      real*8 norm
      logical header
c
c
c     choose the method for computation of induced dipoles
c
      if (solvtyp(1:2) .eq. 'PB') then
         call induce0d
      else if (solvtyp(1:2) .eq. 'GK') then
         call induce0c
      else if (poltyp .eq. 'TCG') then
         call induce0b
      else
         call induce0a
      end if
c
c     update the lists of previous induced dipole values
c
      if (use_pred) then
         nualt = min(nualt+1,maxualt)
         do ii = 1, npole
            do j = 1, 3
               do k = nualt, 2, -1
                  udalt(k,j,ii) = udalt(k-1,j,ii)
                  upalt(k,j,ii) = upalt(k-1,j,ii)
               end do
               udalt(1,j,ii) = uind(j,ii)
               upalt(1,j,ii) = uinp(j,ii)
               if (use_solv) then
                  do k = nualt, 2, -1
                     usalt(k,j,ii) = usalt(k-1,j,ii)
                     upsalt(k,j,ii) = upsalt(k-1,j,ii)
                  end do
                  usalt(1,j,ii) = uinds(j,ii)
                  upsalt(1,j,ii) = uinps(j,ii)
               end if
            end do
         end do
      end if
c
c     print out a list of the final induced dipole moments
c
      if (debug .and. use_polar) then
         header = .true.
         do ii = 1, npole
            i = ipole(ii)
            if (polarity(ii) .ne. 0.0d0) then
               if (header) then
                  header = .false.
                  if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
                     write (iout,10)
   10                format (/,' Vacuum Induced Dipole Moments',
     &                          ' (Debye) :')
                  else
                     write (iout,20)
   20                format (/,' Induced Dipole Moments (Debye) :')
                  end if
                  write (iout,30)
   30             format (/,4x,'Atom',15x,'X',12x,'Y',12x,'Z',
     &                       11x,'Total',/)
               end if
               norm = sqrt(uind(1,ii)**2+uind(2,ii)**2+uind(3,ii)**2)
               write (iout,40)  i,(debye*uind(j,ii),j=1,3),debye*norm
   40          format (i8,5x,3f13.4,1x,f13.4)
            end if
         end do
         header = .true.
         if (solvtyp.eq.'GK' .or. solvtyp.eq.'PB') then
            do ii = 1, npole
               i = ipole(ii)
               if (polarity(ii) .ne. 0.0d0) then
                  if (header) then
                     header = .false.
                     write (iout,50)
   50                format (/,' SCRF Induced Dipole Moments',
     &                          ' (Debye) :')
                     write (iout,60)
   60                format (/,4x,'Atom',15x,'X',12x,'Y',12x,'Z',
     &                          11x,'Total',/)
                  end if
                  norm = sqrt(uinds(1,ii)**2+uinds(2,ii)**2
     &                           +uinds(3,ii)**2)
                  write (iout,70)  i,(debye*uinds(j,ii),j=1,3),
     &                             debye*norm
   70             format (i8,5x,3f13.4,1x,f13.4)
               end if
            end do
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0a  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0a" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0a
      use atoms
      use ielscf
      use inform
      use iounit
      use limits
      use mpole
      use polar
      use polopt
      use polpcg
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 udsum,upsum
      real*8 a,ap,b,bp
      real*8 sum,sump,term
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: usum(:,:)
      real*8, allocatable :: usump(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
c
c     get the electrostatic field due to permanent multipoles
c
      if (use_ewald) then
         call dfield0c (field,fieldp)
      else if (use_mlist) then
         call dfield0b (field,fieldp)
      else
         call dfield0a (field,fieldp)
      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         if (douind(ipole(i))) then
            do j = 1, 3
               udir(j,i) = polarity(i) * field(j,i)
               udirp(j,i) = polarity(i) * fieldp(j,i)
               if (pcgguess) then
                  uind(j,i) = udir(j,i)
                  uinp(j,i) = udirp(j,i)
               end if
            end do
         end if
      end do
c
c     get induced dipoles via the OPT extrapolation method
c
      if (poltyp .eq. 'OPT') then
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uopt(0,j,i) = udir(j,i)
                  uoptp(0,j,i) = udirp(j,i)
               end do
            end if
         end do
         do k = 1, optorder
            optlevel = k - 1
            if (use_ewald) then
               call ufield0c (field,fieldp)
            else if (use_mlist) then
               call ufield0b (field,fieldp)
            else
               call ufield0a (field,fieldp)
            end if
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uopt(k,j,i) = polarity(i) * field(j,i)
                     uoptp(k,j,i) = polarity(i) * fieldp(j,i)
                     uind(j,i) = uopt(k,j,i)
                     uinp(j,i) = uoptp(k,j,i)
                  end do
               end if
            end do
         end do
         allocate (usum(3,n))
         allocate (usump(3,n))
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uind(j,i) = 0.0d0
                  uinp(j,i) = 0.0d0
                  usum(j,i) = 0.0d0
                  usump(j,i) = 0.0d0
                  do k = 0, optorder
                     usum(j,i) = usum(j,i) + uopt(k,j,i)
                     usump(j,i) = usump(j,i) + uoptp(k,j,i)
                     uind(j,i) = uind(j,i) + copt(k)*usum(j,i)
                     uinp(j,i) = uinp(j,i) + copt(k)*usump(j,i)
                  end do
               end do
            end if
         end do
         deallocate (usum)
         deallocate (usump)
      end if
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 100
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimate induced dipoles using a polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            call ulspred
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
               end do
            end do
         end if
c
c     estimate induced dipoles via inertial extended Lagrangian
c
         if (use_ielscf) then
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uaux(j,i)
                  uinp(j,i) = upaux(j,i)
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
c
c     get the electrostatic field due to induced dipoles
c
         if (use_ewald) then
            call ufield0c (field,fieldp)
         else if (use_mlist) then
            call ufield0b (field,fieldp)
         else
            call ufield0a (field,fieldp)
         end if

c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            if (douind(ipole(i))) then
               poli(i) = max(polmin,polarity(i))
               do j = 1, 3
                  if (pcgguess) then
                     rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                             + field(j,i)
                     rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                             + fieldp(j,i)
                  else
                     rsd(j,i) = udir(j,i) / poli(i)
                     rsdp(j,i) = udirp(j,i) / poli(i)
                  end if
                  zrsd(j,i) = rsd(j,i)
                  zrsdp(j,i) = rsdp(j,i)
               end do
            else
               do j = 1, 3
                  rsd(j,i) = 0.0d0
                  rsdp(j,i) = 0.0d0
                  zrsd(j,i) = 0.0d0
                  zrsdp(j,i) = 0.0d0
               end do
            end if
         end do
         if (pcgprec) then
            mode = 'BUILD'
            if (use_mlist) then
               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
               mode = 'APPLY'
               call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
            else
               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
               mode = 'APPLY'
               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
            end if
         end if
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  conj(j,i) = zrsd(j,i)
                  conjp(j,i) = zrsdp(j,i)
               end do
            end if
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     vec(j,i) = uind(j,i)
                     vecp(j,i) = uinp(j,i)
                     uind(j,i) = conj(j,i)
                     uinp(j,i) = conjp(j,i)
                  end do
               end if
            end do
            if (use_ewald) then
               call ufield0c (field,fieldp)
            else if (use_mlist) then
               call ufield0b (field,fieldp)
            else
               call ufield0a (field,fieldp)
            end if
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = vec(j,i)
                     uinp(j,i) = vecp(j,i)
                     vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                     vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                  end do
               end if
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     a = a + conj(j,i)*vec(j,i)
                     ap = ap + conjp(j,i)*vecp(j,i)
                     sum = sum + rsd(j,i)*zrsd(j,i)
                     sump = sump + rsdp(j,i)*zrsdp(j,i)
                  end do
               end if
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = uind(j,i) + a*conj(j,i)
                     uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                     rsd(j,i) = rsd(j,i) - a*vec(j,i)
                     rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                     zrsd(j,i) = rsd(j,i)
                     zrsdp(j,i) = rsdp(j,i)
                  end do
               end if
            end do
            if (pcgprec) then
               if (use_mlist) then
                  call uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
               else
                  call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
               end if
            end if
            b = 0.0d0
            bp = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     b = b + rsd(j,i)*zrsd(j,i)
                     bp = bp + rsdp(j,i)*zrsdp(j,i)
                  end do
               end if
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     conj(j,i) = zrsd(j,i) + b*conj(j,i)
                     conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                     epsd = epsd + rsd(j,i)*rsd(j,i)
                     epsp = epsp + rsdp(j,i)*rsdp(j,i)
                  end do
               end if
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp)
            eps = debye * sqrt(eps/dble(npolar))
            if (debug) then
               if (iter .eq. 1) then
                  write (iout,10)
   10             format (/,' Determination of SCF Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',7x,'RMS Residual (Debye)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
c
c     apply a "peek" iteration to the mutual induced dipoles
c
            if (done) then
               do i = 1, npole
                  if (douind(ipole(i))) then
                     term = pcgpeek * poli(i)
                     do j = 1, 3
                        uind(j,i) = uind(j,i) + term*rsd(j,i)
                        uinp(j,i) = uinp(j,i) + term*rsdp(j,i)
                     end do
                  end if
               end do
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (conj)
         deallocate (conjp)
         deallocate (vec)
         deallocate (vecp)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',4x,'Iterations',i5,
     &                 7x,'RMS Residual',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dfield0a  --  direct induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dfield0a" computes the direct electrostatic field due to
c     permanent multipole moments via a double loop
c
c
      subroutine dfield0a (field,fieldp)
      use atoms
      use bound
      use cell
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
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
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
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
            ddi = dirdamp(ii)
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
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  if (damp .ne. 0.0d0) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .ne. 0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
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
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx        
                  fid(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fid(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkd(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkd(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkd(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
               end if
               do j = 1, 3
                  field(j,ii) = field(j,ii) + fid(j)*dscale(k)
                  field(j,kk) = field(j,kk) + fkd(j)*dscale(k)
                  fieldp(j,ii) = fieldp(j,ii) + fid(j)*pscale(k)
                  fieldp(j,kk) = fieldp(j,kk) + fkd(j)*pscale(k)
               end do
            end if
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
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
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
            if (use_thole) then
               pdi = pdamp(ii)
               pti = thole(ii)
               ddi = dirdamp(ii)
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
     &                  pscale(i12(j,i)) = p2iscale
                  end do
                  dscale(i12(j,i)) = pscale(i12(j,i))
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
                  dscale(i13(j,i)) = pscale(i13(j,i))
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
                  dscale(i14(j,i)) = pscale(i14(j,i))
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
                  dscale(i15(j,i)) = pscale(i15(j,i))
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  pscale(i12(j,i)) = p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
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
            do kk = ii, npole
               k = ipole(kk)
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
               do m = 2, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
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
                           pgamma = min(ddi,dirdamp(kk))
                           if (pgamma .ne. 0.0d0) then
                              damp = pgamma * (r/damp)**(1.5d0)
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp) 
                                 scale3 = 1.0d0 - expdamp 
                                 scale5 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.5d0*damp)
                                 scale7 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.65d0*damp
     &                                                  +0.15d0*damp**2)
                              end if
                           else
                              pgamma = min(pti,thole(kk))
                              damp = pgamma * (r/damp)**3
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp)
                                 scale3 = 1.0d0 - expdamp
                                 scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                                 scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                               +0.6d0*damp**2)
                              end if
                           end if
                        end if
                        rr3 = scale3 / (r*r2)
                        rr5 = 3.0d0 * scale5 / (r*r2*r2)
                        rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                        fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                              - rr3*dkx + 2.0d0*rr5*qkx
                        fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                              - rr3*dky + 2.0d0*rr5*qky
                        fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                              - rr3*dkz + 2.0d0*rr5*qkz
                        fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                              - rr3*dix - 2.0d0*rr5*qix
                        fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                              - rr3*diy - 2.0d0*rr5*qiy
                        fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                              - rr3*diz - 2.0d0*rr5*qiz
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
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dkx + 2.0d0*rr5k*qkx        
                        fid(2) = -yr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr+rr7k*qkr)
     &                              - rr3k*dky + 2.0d0*rr5k*qky
                        fid(3) = -zr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr+rr7k*qkr)
     &                              - rr3k*dkz + 2.0d0*rr5k*qkz
                        fkd(1) = xr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*dix - 2.0d0*rr5i*qix
                        fkd(2) = yr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diy - 2.0d0*rr5i*qiy
                        fkd(3) = zr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diz - 2.0d0*rr5i*qiz
                     end if
                     do j = 1, 3
                        fip(j) = fid(j)
                        fkp(j) = fkd(j)
                     end do
                     if (use_polymer .and. r2.le.polycut2) then
                        do j = 1, 3
                           fid(j) = fid(j) * dscale(k)
                           fip(j) = fip(j) * pscale(k)
                           fkd(j) = fkd(j) * dscale(k)
                           fkp(j) = fkp(j) * pscale(k)
                        end do
                     end if
                     do j = 1, 3
                        field(j,ii) = field(j,ii) + fid(j)
                        fieldp(j,ii) = fieldp(j,ii) + fip(j)
                        if (i .ne. k) then
                           field(j,kk) = field(j,kk) + fkd(j)
                           fieldp(j,kk) = fieldp(j,kk) + fkp(j)
                        end if
                     end do
                  end if
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ufield0a  --  mutual induction via double loop  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ufield0a" computes the mutual electrostatic field due to
c     induced dipole moments via a double loop
c
c
      subroutine ufield0a (field,fieldp)
      use atoms
      use bound
      use cell
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 dix,diy,diz
      real*8 pix,piy,piz
      real*8 dkx,dky,dkz
      real*8 pkx,pky,pkz
      real*8 dir,pir
      real*8 dkr,pkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
         end do
      end do
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
      end do
c
c     find the electrostatic field due to mutual induced dipoles
c
      do ii = 1, npole-1
         i = ipole(ii)
         dix = uind(1,ii)
         diy = uind(2,ii)
         diz = uind(3,ii)
         pix = uinp(1,ii)
         piy = uinp(2,ii)
         piz = uinp(3,ii)
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               dkx = uind(1,kk)
               dky = uind(2,kk)
               dkz = uind(3,kk)
               pkx = uinp(1,kk)
               pky = uinp(2,kk)
               pkz = uinp(3,kk)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               pir = pix*xr + piy*yr + piz*zr
               pkr = pkx*xr + pky*yr + pkz*zr
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  scale3 = uscale(k)
                  scale5 = uscale(k)
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0+damp))
                     end if
                  end if
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampmut (r,alphai,alphak,dmpik)
                  scale3 = wscale(k) * dmpik(3)
                  scale5 = wscale(k) * dmpik(5)
               end if
               rr3 = -scale3 / (r*r2)
               rr5 = 3.0d0 * scale5 / (r*r2*r2)
               fid(1) = rr3*dkx + rr5*dkr*xr
               fid(2) = rr3*dky + rr5*dkr*yr
               fid(3) = rr3*dkz + rr5*dkr*zr
               fkd(1) = rr3*dix + rr5*dir*xr
               fkd(2) = rr3*diy + rr5*dir*yr
               fkd(3) = rr3*diz + rr5*dir*zr
               fip(1) = rr3*pkx + rr5*pkr*xr
               fip(2) = rr3*pky + rr5*pkr*yr
               fip(3) = rr3*pkz + rr5*pkr*zr
               fkp(1) = rr3*pix + rr5*pir*xr
               fkp(2) = rr3*piy + rr5*pir*yr
               fkp(3) = rr3*piz + rr5*pir*zr
               do j = 1, 3
                  field(j,ii) = field(j,ii) + fid(j)
                  field(j,kk) = field(j,kk) + fkd(j)
                  fieldp(j,ii) = fieldp(j,ii) + fip(j)
                  fieldp(j,kk) = fieldp(j,kk) + fkp(j)
               end do
            end if
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
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do ii = 1, npole
            i = ipole(ii)
            dix = uind(1,ii)
            diy = uind(2,ii)
            diz = uind(3,ii)
            pix = uinp(1,ii)
            piy = uinp(2,ii)
            piz = uinp(3,ii)
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
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               dkx = uind(1,kk)
               dky = uind(2,kk)
               dkz = uind(3,kk)
               pkx = uinp(1,kk)
               pky = uinp(2,kk)
               pkz = uinp(3,kk)
               do m = 2, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
c
c     intermediates involving moments and separation distance
c
                     dir = dix*xr + diy*yr + diz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     pir = pix*xr + piy*yr + piz*zr
                     pkr = pkx*xr + pky*yr + pkz*zr
c
c     find the field components for Thole polarization damping
c
                     if (use_thole) then
                        scale3 = uscale(k)
                        scale5 = uscale(k)
                        damp = pdi * pdamp(kk)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(pti,thole(kk))
                           damp = pgamma * (r/damp)**3
                           if (damp .lt. 50.0d0) then
                              expdamp = exp(damp)
                              scale3 = scale3 * (1.0d0-expdamp)
                              scale5 = scale5 * (1.0d0-expdamp
     &                                              *(1.0d0+damp))
                           end if
                        end if
c
c     find the field components for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        call dampmut (r,alphai,alphak,dmpik)
                        scale3 = wscale(k) * dmpik(3)
                        scale5 = wscale(k) * dmpik(5)
                     end if
                     rr3 = -scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     fid(1) = rr3*dkx + rr5*dkr*xr
                     fid(2) = rr3*dky + rr5*dkr*yr
                     fid(3) = rr3*dkz + rr5*dkr*zr
                     fkd(1) = rr3*dix + rr5*dir*xr
                     fkd(2) = rr3*diy + rr5*dir*yr
                     fkd(3) = rr3*diz + rr5*dir*zr
                     fip(1) = rr3*pkx + rr5*pkr*xr
                     fip(2) = rr3*pky + rr5*pkr*yr
                     fip(3) = rr3*pkz + rr5*pkr*zr
                     fkp(1) = rr3*pix + rr5*pir*xr
                     fkp(2) = rr3*piy + rr5*pir*yr
                     fkp(3) = rr3*piz + rr5*pir*zr
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           do j = 1, 3
                              fid(j) = fid(j) * uscale(k)
                              fkd(j) = fkd(j) * uscale(k)
                              fip(j) = fip(j) * uscale(k)
                              fkp(j) = fkp(j) * uscale(k)
                           end do
                        end if
                     end if
                     do j = 1, 3
                        field(j,ii) = field(j,ii) + fid(j)
                        fieldp(j,ii) = fieldp(j,ii) + fip(j)
                        if (ii .ne. kk) then
                           field(j,kk) = field(j,kk) + fkd(j)
                           fieldp(j,kk) = fieldp(j,kk) + fkp(j)
                        end if
                     end do
                  end if
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0b  --  direct induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0b" computes the direct electrostatic field due to
c     permanent multipole moments via a pair list
c
c
      subroutine dfield0b (field,fieldp)
      use atoms
      use bound
      use chgpen
      use couple
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
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
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 fid(3),fkd(3)
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      character*6 mode
c
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,rpole,x,y,z,pdamp,thole,dirdamp,pcore,pval,
!$OMP& palpha,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,np13,
!$OMP& ip13,np14,ip14,p2scale,p3scale,p4scale,p5scale,p2iscale,p3iscale,
!$OMP& p4iscale,p5iscale,d1scale,d2scale,d3scale,d4scale,nelst,elst,
!$OMP& dpequal,use_thole,use_chgpen,use_bounds,off2,field,fieldp)
!$OMP& firstprivate(dscale,pscale) shared (fieldt,fieldtp)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
c
c     find the electrostatic field due to permanent multipoles
c
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
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
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
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
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
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  if (damp .ne. 0.0d0) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .ne. 0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
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
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx        
                  fid(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fid(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkd(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkd(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkd(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
               end if
               do j = 1, 3
                  fieldt(j,ii) = fieldt(j,ii) + fid(j)*dscale(k)
                  fieldt(j,kk) = fieldt(j,kk) + fkd(j)*dscale(k)
                  fieldtp(j,ii) = fieldtp(j,ii) + fid(j)*pscale(k)
                  fieldtp(j,kk) = fieldtp(j,kk) + fkd(j)*pscale(k)
               end do
            end if
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
!$OMP END DO
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = fieldt(j,ii)
            fieldp(j,ii) = fieldtp(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0b  --  mutual induction via pair list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0b" computes the mutual electrostatic field due to
c     induced dipole moments via a pair list
c
c
      subroutine ufield0b (field,fieldp)
      use atoms
      use bound
      use chgpen
      use couple
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 dix,diy,diz
      real*8 pix,piy,piz
      real*8 dkx,dky,dkz
      real*8 pkx,pky,pkz
      real*8 dir,pir
      real*8 dkr,pkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      character*6 mode
c
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,uind,uinp,x,y,z,pdamp,thole,pcore,pval,
!$OMP& palpha,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,u1scale,u2scale,u3scale,u4scale,w2scale,
!$OMP& w3scale,w4scale,w5scale,nelst,elst,use_thole,use_chgpen,
!$OMP& use_bounds,off2,field,fieldp)
!$OMP& firstprivate(uscale,wscale) shared (fieldt,fieldtp)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
c
c     find the electrostatic field due to mutual induced dipoles
c
      do ii = 1, npole
         i = ipole(ii)
         dix = uind(1,ii)
         diy = uind(2,ii)
         diz = uind(3,ii)
         pix = uinp(1,ii)
         piy = uinp(2,ii)
         piz = uinp(3,ii)
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
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               dkx = uind(1,kk)
               dky = uind(2,kk)
               dkz = uind(3,kk)
               pkx = uinp(1,kk)
               pky = uinp(2,kk)
               pkz = uinp(3,kk)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               pir = pix*xr + piy*yr + piz*zr
               pkr = pkx*xr + pky*yr + pkz*zr
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  scale3 = uscale(k)
                  scale5 = uscale(k)
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0+damp))
                     end if
                  end if
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampmut (r,alphai,alphak,dmpik)
                  scale3 = wscale(k) * dmpik(3)
                  scale5 = wscale(k) * dmpik(5)
               end if
               rr3 = -scale3 / (r*r2)
               rr5 = 3.0d0 * scale5 / (r*r2*r2)
               fid(1) = rr3*dkx + rr5*dkr*xr
               fid(2) = rr3*dky + rr5*dkr*yr
               fid(3) = rr3*dkz + rr5*dkr*zr
               fkd(1) = rr3*dix + rr5*dir*xr
               fkd(2) = rr3*diy + rr5*dir*yr
               fkd(3) = rr3*diz + rr5*dir*zr
               fip(1) = rr3*pkx + rr5*pkr*xr
               fip(2) = rr3*pky + rr5*pkr*yr
               fip(3) = rr3*pkz + rr5*pkr*zr
               fkp(1) = rr3*pix + rr5*pir*xr
               fkp(2) = rr3*piy + rr5*pir*yr
               fkp(3) = rr3*piz + rr5*pir*zr
               do j = 1, 3
                  fieldt(j,ii) = fieldt(j,ii) + fid(j)
                  fieldt(j,kk) = fieldt(j,kk) + fkd(j)
                  fieldtp(j,ii) = fieldtp(j,ii) + fip(j)
                  fieldtp(j,kk) = fieldtp(j,kk) + fkp(j)
               end do
            end if
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
!$OMP END DO
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = fieldt(j,ii)
            fieldp(j,ii) = fieldtp(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0c  --  direct induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0c" computes the mutual electrostatic field due to
c     permanent multipole moments via Ewald summation
c
c
      subroutine dfield0c (field,fieldp)
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use pme
      use polar
      implicit none
      integer i,j,ii
      real*8 term
      real*8 ucell(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
         end do
      end do
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
c
c     get the reciprocal space part of the permanent field
c
      call udirect1 (field)
      do ii = 1, npole
         do j = 1, 3
            fieldp(j,ii) = field(j,ii)
         end do
      end do
c
c     get the real space portion of the permanent field
c
      if (use_mlist) then
         call udirect2b (field,fieldp)
      else
         call udirect2a (field,fieldp)
      end if
c
c     get the self-energy portion of the permanent field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + term*rpole(j+1,ii)
            fieldp(j,ii) = fieldp(j,ii) + term*rpole(j+1,ii)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do ii = 1, npole
            i = ipole(ii)
            ucell(1) = ucell(1) + rpole(2,ii) + rpole(1,ii)*x(i)
            ucell(2) = ucell(2) + rpole(3,ii) + rpole(1,ii)*y(i)
            ucell(3) = ucell(3) + rpole(4,ii) + rpole(1,ii)*z(i)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do ii = 1, npole
            do j = 1, 3
               field(j,ii) = field(j,ii) - term*ucell(j)
               fieldp(j,ii) = fieldp(j,ii) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine udirect1  --  Ewald recip direct induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1" computes the reciprocal space contribution of the
c     permanent atomic multipole moments to the field
c
c     note that cmp, fmp, cphi and fphi should not be made global
c     since corresponding values in empole and epolar are different
c
c
      subroutine udirect1 (field)
      use bound
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use polpot
      implicit none
      integer i,j,k,ii
      integer k1,k2,k3
      integer m1,m2,m3
      integer ntot,nff
      integer nf1,nf2,nf3
      real*8 r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 field(3,*)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: cphi(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (cphi(10,npole))
      allocate (fphi(20,npole))
c
c     perform dynamic allocation of some global arrays
c
      ntot = nfft1 * nfft2 * nfft3
      if (allocated(qgrid)) then
         if (size(qgrid) .ne. 2*ntot)  call fftclose
      end if
      if (allocated(qfac)) then
         if (size(qfac) .ne. ntot)  deallocate (qfac)
      end if
      if (.not. allocated(qgrid))  call fftsetup
      if (.not. allocated(qfac))  allocate (qfac(nfft1,nfft2,nfft3))
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
      do ii = 1, npole
         cmp(1,ii) = rpole(1,ii)
         cmp(2,ii) = rpole(2,ii)
         cmp(3,ii) = rpole(3,ii)
         cmp(4,ii) = rpole(4,ii)
         cmp(5,ii) = rpole(5,ii)
         cmp(6,ii) = rpole(9,ii)
         cmp(7,ii) = rpole(13,ii)
         cmp(8,ii) = 2.0d0 * rpole(6,ii)
         cmp(9,ii) = 2.0d0 * rpole(7,ii)
         cmp(10,ii) = 2.0d0 * rpole(10,ii)
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
      qfac(1,1,1) = 0.0d0
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
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for zeroth grid point for nonperiodic system
c
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_mpole (fphi)
c
c     convert the field from fractional to Cartesian
c
      call fphi_to_cphi (fphi,cphi)
c
c     increment the field at each multipole site
c
      do ii = 1, npole
         field(1,ii) = field(1,ii) - cphi(2,ii)
         field(2,ii) = field(2,ii) - cphi(3,ii)
         field(3,ii) = field(3,ii) - cphi(4,ii)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (cmp)
      deallocate (fmp)
      deallocate (cphi)
      deallocate (fphi)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2a  --  Ewald real direct field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2a" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a double loop
c
c
      subroutine udirect2a (field,fieldp)
      use atoms
      use boxes
      use bound
      use cell
      use chgpen
      use couple
      use ewald
      use math
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3),bcn(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute real space Ewald field due to permanent multipoles
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
            ddi = dirdamp(ii)
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
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = 3.0d0 * rr2 * rr3
               rr7 = 5.0d0 * rr2 * rr5
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
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .ne. 0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  scalek = dscale(k)
                  bcn(1) = bn(1) - (1.0d0-scalek*scale3)*rr3
                  bcn(2) = bn(2) - (1.0d0-scalek*scale5)*rr5
                  bcn(3) = bn(3) - (1.0d0-scalek*scale7)*rr7
                  fid(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                  fid(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dky + 2.0d0*bcn(2)*qky
                  fid(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                  fkd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*dix - 2.0d0*bcn(2)*qix
                  fkd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                  fkd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diz - 2.0d0*bcn(2)*qiz
                  scalek = pscale(k)
                  bcn(1) = bn(1) - (1.0d0-scalek*scale3)*rr3
                  bcn(2) = bn(2) - (1.0d0-scalek*scale5)*rr5
                  bcn(3) = bn(3) - (1.0d0-scalek*scale7)*rr7
                  fip(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                  fip(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dky + 2.0d0*bcn(2)*qky
                  fip(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                  fkp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*dix - 2.0d0*bcn(2)*qix
                  fkp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                  fkp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  fid(1) = -xr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx
                  fid(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fid(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkd(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkd(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkd(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
                  scalek = pscale(k)
                  rr3 = rr2 * rr1
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  fip(1) = -xr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx
                  fip(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fip(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkp(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkp(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkp(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
               end if
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,ii) = field(j,ii) + fid(j)
                  field(j,kk) = field(j,kk) + fkd(j)
                  fieldp(j,ii) = fieldp(j,ii) + fip(j)
                  fieldp(j,kk) = fieldp(j,kk) + fkp(j)
               end do
            end if
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
c
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
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
            if (use_thole) then
               pdi = pdamp(ii)
               pti = thole(ii)
               ddi = dirdamp(ii)
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
     &                  pscale(i12(j,i)) = p2iscale
                  end do
                  dscale(i12(j,i)) = pscale(i12(j,i))
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
                  dscale(i13(j,i)) = pscale(i13(j,i))
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
                  dscale(i14(j,i)) = pscale(i14(j,i))
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
                  dscale(i15(j,i)) = pscale(i15(j,i))
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  pscale(i12(j,i)) = p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
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
            do kk = ii, npole
               k = ipole(kk)
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
               do m = 2, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c
c     calculate the error function damping factors
c
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = 3.0d0 * rr2 * rr3
                     rr7 = 5.0d0 * rr2 * rr5
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
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     find the field components for Thole polarization damping
c
                     if (use_thole) then
                        scale3 = 1.0d0
                        scale5 = 1.0d0
                        scale7 = 1.0d0
                        damp = pdi * pdamp(kk)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(ddi,dirdamp(kk))
                           if (pgamma .ne. 0.0d0) then
                              damp = pgamma * (r/damp)**(1.5d0)
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp) 
                                 scale3 = 1.0d0 - expdamp 
                                 scale5 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.5d0*damp)
                                 scale7 = 1.0d0 - expdamp
     &                                               *(1.0d0+0.65d0*damp
     &                                                  +0.15d0*damp**2)
                              end if
                           else
                              pgamma = min(pti,thole(kk))
                              damp = pgamma * (r/damp)**3
                              if (damp .lt. 50.0d0) then
                                 expdamp = exp(-damp)
                                 scale3 = 1.0d0 - expdamp
                                 scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                                 scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                               +0.6d0*damp**2)
                              end if
                           end if
                        end if
                        dsc3 = scale3
                        dsc5 = scale5
                        dsc7 = scale7
                        psc3 = scale3
                        psc5 = scale5
                        psc7 = scale7
                        if (use_polymer) then
                           if (r2 .le. polycut2) then
                              dsc3 = scale3 * dscale(k)
                              dsc5 = scale5 * dscale(k)
                              dsc7 = scale7 * dscale(k)
                              psc3 = scale3 * pscale(k)
                              psc5 = scale5 * pscale(k)
                              psc7 = scale7 * pscale(k)
                           end if
                        end if
                        bcn(1) = bn(1) - (1.0d0-dsc3)*rr3
                        bcn(2) = bn(2) - (1.0d0-dsc5)*rr5
                        bcn(3) = bn(3) - (1.0d0-dsc7)*rr7
                        fid(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                        fid(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dky + 2.0d0*bcn(2)*qky
                        fid(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                        fkd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*dix - 2.0d0*bcn(2)*qix
                        fkd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                        fkd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*diz - 2.0d0*bcn(2)*qiz
                        bcn(1) = bn(1) - (1.0d0-psc3)*rr3
                        bcn(2) = bn(2) - (1.0d0-psc5)*rr5
                        bcn(3) = bn(3) - (1.0d0-psc7)*rr7
                        fip(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                        fip(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dky + 2.0d0*bcn(2)*qky
                        fip(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                              - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                        fkp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*dix - 2.0d0*bcn(2)*qix
                        fkp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                        fkp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                              - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find the field components for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        call dampdir (r,alphai,alphak,dmpi,dmpk)
                        scalek = 1.0d0
                        if (use_polymer) then
                           if (r2 .le. polycut2) then
                              scalek = dscale(k)
                           end if
                        end if
                        rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                        rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                        rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                        rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                        rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                        rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                        rr3 = bn(1) - (1.0d0-scalek)*rr3
                        fid(1) = -xr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dkx + 2.0d0*rr5k*qkx
                        fid(2) = -yr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dky + 2.0d0*rr5k*qky
                        fid(3) = -zr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dkz + 2.0d0*rr5k*qkz
                        fkd(1) = xr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*dix - 2.0d0*rr5i*qix
                        fkd(2) = yr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diy - 2.0d0*rr5i*qiy
                        fkd(3) = zr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diz - 2.0d0*rr5i*qiz
                        scalek = 1.0d0
                        if (use_polymer) then
                           if (r2 .le. polycut2) then
                              scalek = pscale(k)
                           end if
                        end if
                        rr3 = rr2 * rr1
                        rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                        rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                        rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                        rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                        rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                        rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                        rr3 = bn(1) - (1.0d0-scalek)*rr3
                        fip(1) = -xr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dkx + 2.0d0*rr5k*qkx
                        fip(2) = -yr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dky + 2.0d0*rr5k*qky
                        fip(3) = -zr*(rr3*corek + rr3k*valk
     &                              - rr5k*dkr + rr7k*qkr)
     &                              - rr3k*dkz + 2.0d0*rr5k*qkz
                        fkp(1) = xr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*dix - 2.0d0*rr5i*qix
                        fkp(2) = yr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diy - 2.0d0*rr5i*qiy
                        fkp(3) = zr*(rr3*corei + rr3i*vali
     &                              + rr5i*dir + rr7i*qir)
     &                              - rr3i*diz - 2.0d0*rr5i*qiz
                     end if
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,ii) = field(j,ii) + fid(j)
                        fieldp(j,ii) = fieldp(j,ii) + fid(j)
                        if (i .ne. k) then
                           field(j,kk) = field(j,kk) + fkp(j)
                           fieldp(j,kk) = fieldp(j,kk) + fkp(j)
                        end if
                     end do
                  end if
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2b  --  Ewald real direct field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2b" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a neighbor list
c
c
      subroutine udirect2b (field,fieldp)
      use atoms
      use boxes
      use bound
      use chgpen
      use couple
      use ewald
      use math
      use mplpot
      use mpole
      use neigh
      use openmp
      use polar
      use polgrp
      use polpot
      use shunt
      use tarray
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      integer nlocal,nchunk
      integer tid,maxlocal
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 rr3ik,rr5ik
      real*8 erfc,bfac,exp2a
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7,scalek
      real*8 bn(0:3),bcn(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpi(7),dmpk(7)
      real*8 dmpik(5)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: dlocal(:,:)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     values for storage of mutual polarization intermediates
c
      nchunk = int(0.5d0*dble(npole)/dble(nthread)) + 1
      maxlocal = int(dble(npole)*dble(maxelst)/dble(nthread))
      nlocal = 0
      ntpair = 0
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (toffset(0:nthread-1))
      if (poltyp .ne. 'DIRECT') then
         allocate (ilocal(2,maxlocal))
         allocate (dlocal(6,maxlocal))
      end if
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         wscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,rpole,x,y,z,pdamp,
!$OMP& thole,dirdamp,pcore,pval,palpha,p2scale,p3scale,p4scale,p5scale,
!$OMP& p2iscale,p3iscale,p4iscale,p5iscale,w2scale,w3scale,w4scale,
!$OMP& w5scale,d1scale,d2scale,d3scale,d4scale,u1scale,u2scale,u3scale,
!$OMP& u4scale,n12,i12,n13,i13,n14,i14,n15,i15,np11,ip11,np12,ip12,
!$OMP& np13,ip13,np14,ip14,nelst,elst,dpequal,use_thole,use_chgpen,
!$OMP& use_bounds,off2,aewald,aesq2,aesq2n,poltyp,nchunk,ntpair,tindex,
!$OMP& tdipdip,toffset,field,fieldp,fieldt,fieldtp)
!$OMP& firstprivate(pscale,dscale,uscale,wscale,nlocal)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(static,nchunk)
c
c     compute the real space portion of the Ewald summation
c
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
         if (use_thole) then
            pdi = pdamp(ii)
            pti = thole(ii)
            ddi = dirdamp(ii)
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
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               dscale(i13(j,i)) = pscale(i13(j,i))
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               dscale(i14(j,i)) = pscale(i14(j,i))
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               dscale(i15(j,i)) = pscale(i15(j,i))
               wscale(i15(j,i)) = w5scale
            end do
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
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               pscale(i12(j,i)) = p2iscale
               end do
               wscale(i12(j,i)) = w2scale
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               pscale(i13(j,i)) = p3iscale
               end do
               wscale(i13(j,i)) = w3scale
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                   if (i14(j,i) .eq. ip11(k,i))
     &               pscale(i14(j,i)) = p4iscale
               end do
               wscale(i14(j,i)) = w4scale
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               pscale(i15(j,i)) = p5iscale
               end do
               wscale(i15(j,i)) = w5scale
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
               uscale(ip11(j,i)) = u1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
               uscale(ip12(j,i)) = u2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
               uscale(ip13(j,i)) = u3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
               uscale(ip14(j,i)) = u4scale
            end do
         end if
c
c     evaluate all sites within the cutoff distance
c
         do kkk = 1, nelst(ii)
            kk = elst(kkk,ii)
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = 3.0d0 * rr2 * rr3
               rr7 = 5.0d0 * rr2 * rr5
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
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 3
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .ne. 0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  scalek = dscale(k)
                  bcn(1) = bn(1) - (1.0d0-scalek*scale3)*rr3
                  bcn(2) = bn(2) - (1.0d0-scalek*scale5)*rr5
                  bcn(3) = bn(3) - (1.0d0-scalek*scale7)*rr7
                  fid(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                  fid(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dky + 2.0d0*bcn(2)*qky
                  fid(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                  fkd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*dix - 2.0d0*bcn(2)*qix
                  fkd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                  fkd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diz - 2.0d0*bcn(2)*qiz
                  scalek = pscale(k)
                  bcn(1) = bn(1) - (1.0d0-scalek*scale3)*rr3
                  bcn(2) = bn(2) - (1.0d0-scalek*scale5)*rr5
                  bcn(3) = bn(3) - (1.0d0-scalek*scale7)*rr7
                  fip(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
                  fip(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dky + 2.0d0*bcn(2)*qky
                  fip(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                        - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
                  fkp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*dix - 2.0d0*bcn(2)*qix
                  fkp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diy - 2.0d0*bcn(2)*qiy
                  fkp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                        - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find terms needed later to compute mutual polarization
c
                  if (poltyp .ne. 'DIRECT') then
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     damp = pdi * pdamp(kk)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                        end if
                     end if
                     scalek = uscale(k)
                     bcn(1) = bn(1) - (1.0d0-scalek*scale3)*rr3
                     bcn(2) = bn(2) - (1.0d0-scalek*scale5)*rr5
                     nlocal = nlocal + 1
                     ilocal(1,nlocal) = ii
                     ilocal(2,nlocal) = kk
                     dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                     dlocal(2,nlocal) = bcn(2)*xr*yr
                     dlocal(3,nlocal) = bcn(2)*xr*zr
                     dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                     dlocal(5,nlocal) = bcn(2)*yr*zr
                     dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr
                  end if
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampdir (r,alphai,alphak,dmpi,dmpk)
                  scalek = dscale(k)
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  fid(1) = -xr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx
                  fid(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fid(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkd(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkd(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkd(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
                  scalek = pscale(k)
                  rr3 = rr2 * rr1
                  rr3i = bn(1) - (1.0d0-scalek*dmpi(3))*rr3
                  rr5i = bn(2) - (1.0d0-scalek*dmpi(5))*rr5
                  rr7i = bn(3) - (1.0d0-scalek*dmpi(7))*rr7
                  rr3k = bn(1) - (1.0d0-scalek*dmpk(3))*rr3
                  rr5k = bn(2) - (1.0d0-scalek*dmpk(5))*rr5
                  rr7k = bn(3) - (1.0d0-scalek*dmpk(7))*rr7
                  rr3 = bn(1) - (1.0d0-scalek)*rr3
                  fip(1) = -xr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkx + 2.0d0*rr5k*qkx
                  fip(2) = -yr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dky + 2.0d0*rr5k*qky
                  fip(3) = -zr*(rr3*corek + rr3k*valk
     &                        - rr5k*dkr + rr7k*qkr)
     &                        - rr3k*dkz + 2.0d0*rr5k*qkz
                  fkp(1) = xr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*dix - 2.0d0*rr5i*qix
                  fkp(2) = yr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diy - 2.0d0*rr5i*qiy
                  fkp(3) = zr*(rr3*corei + rr3i*vali
     &                        + rr5i*dir + rr7i*qir)
     &                        - rr3i*diz - 2.0d0*rr5i*qiz
c
c     find terms needed later to compute mutual polarization
c
                  if (poltyp .ne. 'DIRECT') then
                     call dampmut (r,alphai,alphak,dmpik)
                     scalek = wscale(k)
                     rr3 = rr2 * rr1
                     rr3ik = bn(1) - (1.0d0-scalek*dmpik(3))*rr3
                     rr5ik = bn(2) - (1.0d0-scalek*dmpik(5))*rr5
                     nlocal = nlocal + 1
                     ilocal(1,nlocal) = ii
                     ilocal(2,nlocal) = kk
                     dlocal(1,nlocal) = -rr3ik + rr5ik*xr*xr
                     dlocal(2,nlocal) = rr5ik*xr*yr
                     dlocal(3,nlocal) = rr5ik*xr*zr
                     dlocal(4,nlocal) = -rr3ik + rr5ik*yr*yr
                     dlocal(5,nlocal) = rr5ik*yr*zr
                     dlocal(6,nlocal) = -rr3ik + rr5ik*zr*zr
                  end if
               end if
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  fieldt(j,ii) = fieldt(j,ii) + fid(j)
                  fieldt(j,kk) = fieldt(j,kk) + fkd(j)
                  fieldtp(j,ii) = fieldtp(j,ii) + fip(j)
                  fieldtp(j,kk) = fieldtp(j,kk) + fkp(j)
               end do
            end if
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               dscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               dscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               dscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               dscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
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
         else
            do j = 1, n12(i)
               pscale(i12(j,i)) = 1.0d0
               wscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               pscale(i13(j,i)) = 1.0d0
               wscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               pscale(i14(j,i)) = 1.0d0
               wscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               pscale(i15(j,i)) = 1.0d0
               wscale(i15(j,i)) = 1.0d0
            end do
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
               uscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
               uscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
               uscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
               uscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
!$OMP END DO
c
c     find offset into global arrays for the current thread
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = ntpair
      ntpair = ntpair + nlocal
!$OMP END CRITICAL
c
c     store terms used later to compute mutual polarization
c
      if (poltyp .ne. 'DIRECT') then
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
            tindex(1,m) = ilocal(1,i)
            tindex(2,m) = ilocal(2,i)
            do j = 1, 6
               tdipdip(j,m) = dlocal(j,i)
            end do
         end do
      end if
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + fieldt(j,ii)
            fieldp(j,ii) = fieldp(j,ii) + fieldtp(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (wscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (toffset)
      if (allocated(ilocal))  deallocate (ilocal)
      if (allocated(dlocal))  deallocate (dlocal)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0c  --  mutual induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0c" computes the mutual electrostatic field due to
c     induced dipole moments via Ewald summation
c
c
      subroutine ufield0c (field,fieldp)
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use pme
      use polar
      implicit none
      integer ii,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c
c
c     zero out the electrostatic field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
         end do
      end do
c
c     set grid size, spline order and Ewald coefficient
c
      nfft1 = nefft1
      nfft2 = nefft2
      nfft3 = nefft3
      bsorder = bsporder
      aewald = apewald
c
c     get the reciprocal space part of the mutual field
c
      call umutual1 (field,fieldp)
c
c     get the real space portion of the mutual field
c
      if (use_mlist) then
         call umutual2b (field,fieldp)
      else
         call umutual2a (field,fieldp)
      end if
c
c     get the self-energy portion of the mutual field
c
      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + term*uind(j,ii)
            fieldp(j,ii) = fieldp(j,ii) + term*uinp(j,ii)
         end do
      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do j = 1, 3
            ucell(j) = 0.0d0
            ucellp(j) = 0.0d0
         end do
         do ii = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind(j,ii)
               ucellp(j) = ucellp(j) + uinp(j,ii)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do ii = 1, npole
            do j = 1, 3
               field(j,ii) = field(j,ii) - term*ucell(j)
               fieldp(j,ii) = fieldp(j,ii) - term*ucellp(j)
            end do
         end do
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine umutual1  --  Ewald recip mutual induced field  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "umutual1" computes the reciprocal space contribution of the
c     induced atomic dipole moments to the field
c
c
      subroutine umutual1 (field,fieldp)
      use boxes
      use ewald
      use math
      use mpole
      use pme
      use polar
      use polopt
      use polpot
      implicit none
      integer i,j,k,ii
      real*8 term
      real*8 a(3,3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: fuinp(:,:)
      real*8, allocatable :: fdip_phi1(:,:)
      real*8, allocatable :: fdip_phi2(:,:)
      real*8, allocatable :: fdip_sum_phi(:,:)
      real*8, allocatable :: dipfield1(:,:)
      real*8, allocatable :: dipfield2(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (fuinp(3,npole))
      allocate (fdip_phi1(10,npole))
      allocate (fdip_phi2(10,npole))
      allocate (fdip_sum_phi(20,npole))
      allocate (dipfield1(3,npole))
      allocate (dipfield2(3,npole))
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
      do ii = 1, npole
         do j = 1, 3
            fuind(j,ii) = a(j,1)*uind(1,ii) + a(j,2)*uind(2,ii)
     &                       + a(j,3)*uind(3,ii)
            fuinp(j,ii) = a(j,1)*uinp(1,ii) + a(j,2)*uinp(2,ii)
     &                       + a(j,3)*uinp(3,ii)
         end do
      end do
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_uind (fuind,fuinp)
      call fftfront
c
c     complete the transformation of the PME grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get field
c
      call fftback
      call fphi_uind (fdip_phi1,fdip_phi2,fdip_sum_phi)
c
c     store fractional reciprocal potentials for OPT method
c
      if (poltyp .eq. 'OPT') then
         do ii = 1, npole
            do j = 1, 10
               fopt(optlevel,j,ii) = fdip_phi1(j,ii)
               foptp(optlevel,j,ii) = fdip_phi2(j,ii)
            end do
         end do
      end if
c
c     convert the dipole fields from fractional to Cartesian
c
      do i = 1, 3
         a(i,1) = dble(nfft1) * recip(i,1)
         a(i,2) = dble(nfft2) * recip(i,2)
         a(i,3) = dble(nfft3) * recip(i,3)
      end do
      do ii = 1, npole
         do j = 1, 3
            dipfield1(j,ii) = a(j,1)*fdip_phi1(2,ii)
     &                           + a(j,2)*fdip_phi1(3,ii)
     &                           + a(j,3)*fdip_phi1(4,ii)
            dipfield2(j,ii) = a(j,1)*fdip_phi2(2,ii)
     &                           + a(j,2)*fdip_phi2(3,ii)
     &                           + a(j,3)*fdip_phi2(4,ii)
         end do
      end do
c
c     increment the field at each multipole site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) - dipfield1(j,ii)
            fieldp(j,ii) = fieldp(j,ii) - dipfield2(j,ii)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (fuinp)
      deallocate (fdip_phi1)
      deallocate (fdip_phi2)
      deallocate (fdip_sum_phi)
      deallocate (dipfield1)
      deallocate (dipfield2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2a  --  Ewald real mutual field via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2a" computes the real space contribution of the induced
c     atomic dipole moments to the field via a double loop
c
c
      subroutine umutual2a (field,fieldp)
      use atoms
      use boxes
      use bound
      use cell
      use chgpen
      use couple
      use ewald
      use math
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use units
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 dix,diy,diz
      real*8 pix,piy,piz
      real*8 dkx,dky,dkz
      real*8 pkx,pky,pkz
      real*8 dir,dkr
      real*8 pir,pkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 bn(0:2)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      character*6 mode
      external erfc
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
      mode = 'EWALD'
      call switch (mode)
      aesq2 = 2.0 * aewald * aewald
      aesq2n = 0.0d0
      if (aewald .gt. 0.0d0)  aesq2n = 1.0d0 / (sqrtpi*aewald)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
      allocate (wscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do ii = 1, npole-1
         i = ipole(ii)
         dix = uind(1,ii)
         diy = uind(2,ii)
         diz = uind(3,ii)
         pix = uinp(1,ii)
         piy = uinp(2,ii)
         piz = uinp(3,ii)
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               dkx = uind(1,k)
               dky = uind(2,k)
               dkz = uind(3,k)
               pkx = uinp(1,k)
               pky = uinp(2,k)
               pkz = uinp(3,k)
c
c     intermediates involving moments and separation distance
c
               dir = dix*xr + diy*yr + diz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               pir = pix*xr + piy*yr + piz*zr
               pkr = pkx*xr + pky*yr + pkz*zr
c
c     calculate the real space Ewald error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) * rr1
               exp2a = exp(-ralpha**2)
               aefac = aesq2n
               do j = 1, 2
                  bfac = dble(j+j-1)
                  aefac = aesq2 * aefac
                  bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
               end do
c
c     find the field components for Thole polarization damping
c
               if (use_thole) then
                  scale3 = uscale(k)
                  scale5 = uscale(k)
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0+damp))
                     end if
                  end if
c
c     find the field components for charge penetration damping
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call dampmut (r,alphai,alphak,dmpik)
                  scale3 = wscale(k) * dmpik(3)
                  scale5 = wscale(k) * dmpik(5)
               end if
c
c     find the field terms for the current interaction
c
               rr3 = -bn(1) + (1.0d0-scale3)*rr3
               rr5 = bn(2) - 3.0d0*(1.0d0-scale5)*rr5
               fid(1) = rr3*dkx + rr5*dkr*xr
               fid(2) = rr3*dky + rr5*dkr*yr
               fid(3) = rr3*dkz + rr5*dkr*zr
               fkd(1) = rr3*dix + rr5*dir*xr
               fkd(2) = rr3*diy + rr5*dir*yr
               fkd(3) = rr3*diz + rr5*dir*zr
               fip(1) = rr3*pkx + rr5*pkr*xr
               fip(2) = rr3*pky + rr5*pkr*yr
               fip(3) = rr3*pkz + rr5*pkr*zr
               fkp(1) = rr3*pix + rr5*pir*xr
               fkp(2) = rr3*piy + rr5*pir*yr
               fkp(3) = rr3*piz + rr5*pir*zr
c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,ii) = field(j,ii) + fid(j)
                  field(j,kk) = field(j,kk) + fkd(j)
                  fieldp(j,ii) = fieldp(j,ii) + fip(j)
                  fieldp(j,kk) = fieldp(j,kk) + fkp(j)
               end do
            end if
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
c     periodic boundary for large cutoffs via replicates method
c
      if (use_replica) then
         do ii = 1, npole
            i = ipole(ii)
            dix = uind(1,ii)
            diy = uind(2,ii)
            diz = uind(3,ii)
            pix = uinp(1,ii)
            piy = uinp(2,ii)
            piz = uinp(3,ii)
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
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii, npole
               k = ipole(kk)
               dkx = uind(1,kk)
               dky = uind(2,kk)
               dkz = uind(3,kk)
               pkx = uinp(1,kk)
               pky = uinp(2,kk)
               pkz = uinp(3,kk)
               do m = 2, ncell
                  xr = x(k) - x(i)
                  yr = y(k) - y(i)
                  zr = z(k) - z(i)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     rr1 = 1.0d0 / r
                     rr2 = rr1 * rr1
                     rr3 = rr2 * rr1
                     rr5 = rr2 * rr3
c
c     intermediates involving moments and separation distance
c
                     dir = dix*xr + diy*yr + diz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     pir = pix*xr + piy*yr + piz*zr
                     pkr = pkx*xr + pky*yr + pkz*zr
c
c     calculate the real space Ewald error function terms
c
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) * rr1
                     exp2a = exp(-ralpha**2)
                     aefac = aesq2n
                     do j = 1, 2
                        bfac = dble(j+j-1)
                        aefac = aesq2 * aefac
                        bn(j) = (bfac*bn(j-1)+aefac*exp2a) * rr2
                     end do
c
c     find the field components for Thole polarization damping
c
                     if (use_thole) then
                        scale3 = uscale(k)
                        scale5 = uscale(k)
                        damp = pdi * pdamp(kk)
                        if (damp .ne. 0.0d0) then
                           pgamma = min(pti,thole(kk))
                           damp = pgamma * (r/damp)**3
                           if (damp .lt. 50.0d0) then
                              expdamp = exp(-damp)
                              scale3 = scale3 * (1.0d0-expdamp)
                              scale5 = scale5 * (1.0d0-expdamp
     &                                              *(1.0d0+damp))
                           end if
                        end if
c
c     find the field components for charge penetration damping
c
                     else if (use_chgpen) then
                        corek = pcore(kk)
                        valk = pval(kk)
                        alphak = palpha(kk)
                        call dampmut (r,alphai,alphak,dmpik)
                        scale3 = wscale(k) * dmpik(3)
                        scale5 = wscale(k) * dmpik(5)
                     end if
c
c     find the field terms for the current interaction
c
                     rr3 = -bn(1) + (1.0d0-scale3)*rr3
                     rr5 = bn(2) - 3.0d0*(1.0d0-scale5)*rr5
                     fid(1) = rr3*dkx + rr5*dkr*xr
                     fid(2) = rr3*dky + rr5*dkr*yr
                     fid(3) = rr3*dkz + rr5*dkr*zr
                     fkd(1) = rr3*dix + rr5*dir*xr
                     fkd(2) = rr3*diy + rr5*dir*yr
                     fkd(3) = rr3*diz + rr5*dir*zr
                     fip(1) = rr3*pkx + rr5*pkr*xr
                     fip(2) = rr3*pky + rr5*pkr*yr
                     fip(3) = rr3*pkz + rr5*pkr*zr
                     fkp(1) = rr3*pix + rr5*pir*xr
                     fkp(2) = rr3*piy + rr5*pir*yr
                     fkp(3) = rr3*piz + rr5*pir*zr
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,ii) = field(j,ii) + fid(j)
                        fieldp(j,ii) = fieldp(j,ii) + fip(j)
                        if (ii .ne. kk) then
                           field(j,kk) = field(j,kk) + fkd(j)
                           fieldp(j,kk) = fieldp(j,kk) + fkp(j)
                        end if
                     end do
                  end if
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
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (wscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2b  --  Ewald real mutual field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2b" computes the real space contribution of the induced
c     atomic dipole moments to the field via a neighbor list
c
c
      subroutine umutual2b (field,fieldp)
      use mpole
      use polar
      use tarray
      implicit none
      integer i,j,k,m,ii
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
c
c
c     check for multipoles and set cutoff coefficients
c
      if (npole .eq. 0)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,uind,uinp,ntpair,
!$OMP& tindex,tdipdip,field,fieldp,fieldt,fieldtp)
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
c
c     find the field terms for each pairwise interaction
c
      do m = 1, ntpair
         i = tindex(1,m)
         k = tindex(2,m)
         fid(1) = tdipdip(1,m)*uind(1,k) + tdipdip(2,m)*uind(2,k)
     &               + tdipdip(3,m)*uind(3,k)
         fid(2) = tdipdip(2,m)*uind(1,k) + tdipdip(4,m)*uind(2,k)
     &               + tdipdip(5,m)*uind(3,k)
         fid(3) = tdipdip(3,m)*uind(1,k) + tdipdip(5,m)*uind(2,k)
     &               + tdipdip(6,m)*uind(3,k)
         fkd(1) = tdipdip(1,m)*uind(1,i) + tdipdip(2,m)*uind(2,i)
     &               + tdipdip(3,m)*uind(3,i)
         fkd(2) = tdipdip(2,m)*uind(1,i) + tdipdip(4,m)*uind(2,i)
     &               + tdipdip(5,m)*uind(3,i)
         fkd(3) = tdipdip(3,m)*uind(1,i) + tdipdip(5,m)*uind(2,i)
     &               + tdipdip(6,m)*uind(3,i)
         fip(1) = tdipdip(1,m)*uinp(1,k) + tdipdip(2,m)*uinp(2,k)
     &               + tdipdip(3,m)*uinp(3,k)
         fip(2) = tdipdip(2,m)*uinp(1,k) + tdipdip(4,m)*uinp(2,k)
     &               + tdipdip(5,m)*uinp(3,k)
         fip(3) = tdipdip(3,m)*uinp(1,k) + tdipdip(5,m)*uinp(2,k)
     &               + tdipdip(6,m)*uinp(3,k)
         fkp(1) = tdipdip(1,m)*uinp(1,i) + tdipdip(2,m)*uinp(2,i)
     &               + tdipdip(3,m)*uinp(3,i)
         fkp(2) = tdipdip(2,m)*uinp(1,i) + tdipdip(4,m)*uinp(2,i)
     &               + tdipdip(5,m)*uinp(3,i)
         fkp(3) = tdipdip(3,m)*uinp(1,i) + tdipdip(5,m)*uinp(2,i)
     &               + tdipdip(6,m)*uinp(3,i)
c
c     increment the field at each site due to this interaction
c
         do j = 1, 3
            fieldt(j,i) = fieldt(j,i) + fid(j)
            fieldt(j,k) = fieldt(j,k) + fkd(j)
            fieldtp(j,i) = fieldtp(j,i) + fip(j)
            fieldtp(j,k) = fieldtp(j,k) + fkp(j)
         end do
      end do
!$OMP END DO
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + fieldt(j,ii)
            fieldp(j,ii) = fieldp(j,ii) + fieldtp(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine induce0c  --  Kirkwood SCRF induced dipoles  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "induce0c" computes the induced dipole moments at polarizable
c     sites for generalized Kirkwood SCRF and vacuum environments
c
c
      subroutine induce0c
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polopt
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      real*8, allocatable :: usum(:,:)
      real*8, allocatable :: usump(:,:)
      real*8, allocatable :: usums(:,:)
      real*8, allocatable :: usumps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles at each site; uind and uinp are
c     vacuum dipoles, uinds and uinps are SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .and. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0d (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     set SCRF induced dipoles to polarizability times direct field
c     plus the GK reaction field due to permanent multipoles
c
      do i = 1, npole
         if (douind(ipole(i))) then
            do j = 1, 3
               udir(j,i) = polarity(i) * field(j,i)
               udirp(j,i) = polarity(i) * fieldp(j,i)
               udirs(j,i) = polarity(i) * fields(j,i)
               udirps(j,i) = polarity(i) * fieldps(j,i)
               uind(j,i) = udir(j,i)
               uinp(j,i) = udirp(j,i)
               uinds(j,i) = udirs(j,i)
               uinps(j,i) = udirps(j,i)
            end do
         end if
      end do
c
c     get induced dipoles via the OPT extrapolation method
c
      if (poltyp .eq. 'OPT') then
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uopt(0,j,i) = udir(j,i)
                  uoptp(0,j,i) = udirp(j,i)
                  uopts(0,j,i) = udirs(j,i)
                  uoptps(0,j,i) = udirps(j,i)
               end do
            end if
         end do
         do k = 1, optorder
            call ufield0d (field,fieldp,fields,fieldps)
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uopt(k,j,i) = polarity(i) * field(j,i)
                     uoptp(k,j,i) = polarity(i) * fieldp(j,i)
                     uopts(k,j,i) = polarity(i) * fields(j,i)
                     uoptps(k,j,i) = polarity(i) * fieldps(j,i)
                     uind(j,i) = uopt(k,j,i)
                     uinp(j,i) = uoptp(k,j,i)
                     uinds(j,i) = uopts(k,j,i)
                     uinps(j,i) = uoptps(k,j,i)
                  end do
               end if
            end do
         end do
         allocate (usum(3,n))
         allocate (usump(3,n))
         allocate (usums(3,n))
         allocate (usumps(3,n))
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uind(j,i) = 0.0d0
                  uinp(j,i) = 0.0d0
                  uinds(j,i) = 0.0d0
                  uinps(j,i) = 0.0d0
                  usum(j,i) = 0.0d0
                  usump(j,i) = 0.0d0
                  usums(j,i) = 0.0d0
                  usumps(j,i) = 0.0d0
                  do k = 0, optorder
                     usum(j,i) = usum(j,i) + uopt(k,j,i)
                     usump(j,i) = usump(j,i) + uoptp(k,j,i)
                     usums(j,i) = usums(j,i) + uopts(k,j,i)
                     usumps(j,i) = usumps(j,i) + uoptps(k,j,i)
                     uind(j,i) = uind(j,i) + copt(k)*usum(j,i)
                     uinp(j,i) = uinp(j,i) + copt(k)*usump(j,i)
                     uinds(j,i) = uinds(j,i) + copt(k)*usums(j,i)
                     uinps(j,i) = uinps(j,i) + copt(k)*usumps(j,i)
                  end do
               end do
            end if
         end do
         deallocate (usum)
         deallocate (usump)
         deallocate (usums)
         deallocate (usumps)
      end if
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 100
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0d (field,fieldp,fields,fieldps)
         do i = 1, npole
            if (douind(ipole(i))) then
               poli(i) = max(polmin,polarity(i))
               do j = 1, 3
                  rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                          + field(j,i)
                  rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                           + fieldp(j,i)
                  rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                           + fields(j,i)
                  rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                            + fieldps(j,i)
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  conj(j,i) = zrsd(j,i)
                  conjp(j,i) = zrsdp(j,i)
                  conjs(j,i) = zrsds(j,i)
                  conjps(j,i) = zrsdps(j,i)
               end do
            end if
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     vec(j,i) = uind(j,i)
                     vecp(j,i) = uinp(j,i)
                     vecs(j,i) = uinds(j,i)
                     vecps(j,i) = uinps(j,i)
                     uind(j,i) = conj(j,i)
                     uinp(j,i) = conjp(j,i)
                     uinds(j,i) = conjs(j,i)
                     uinps(j,i) = conjps(j,i)
                  end do
               end if
            end do
            call ufield0d (field,fieldp,fields,fieldps)
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = vec(j,i)
                     uinp(j,i) = vecp(j,i)
                     uinds(j,i) = vecs(j,i)
                     uinps(j,i) = vecps(j,i)
                     vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                     vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                     vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                     vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
                  end do
               end if
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     a = a + conj(j,i)*vec(j,i)
                     ap = ap + conjp(j,i)*vecp(j,i)
                     as = as + conjs(j,i)*vecs(j,i)
                     aps = aps + conjps(j,i)*vecps(j,i)
                     sum = sum + rsd(j,i)*zrsd(j,i)
                     sump = sump + rsdp(j,i)*zrsdp(j,i)
                     sums = sums + rsds(j,i)*zrsds(j,i)
                     sumps = sumps + rsdps(j,i)*zrsdps(j,i)
                  end do
               end if
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = uind(j,i) + a*conj(j,i)
                     uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                     uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                     uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                     rsd(j,i) = rsd(j,i) - a*vec(j,i)
                     rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                     rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                     rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
                  end do
               end if
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     zrsd(j,i) = rsd(j,i) * poli(i)
                     zrsdp(j,i) = rsdp(j,i) * poli(i)
                     zrsds(j,i) = rsds(j,i) * poli(i)
                     zrsdps(j,i) = rsdps(j,i) * poli(i)
                     b = b + rsd(j,i)*zrsd(j,i)
                     bp = bp + rsdp(j,i)*zrsdp(j,i)
                     bs = bs + rsds(j,i)*zrsds(j,i)
                     bps = bps + rsdps(j,i)*zrsdps(j,i)
                  end do
               end if
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     conj(j,i) = zrsd(j,i) + b*conj(j,i)
                     conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                     conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                     conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                     epsd = epsd + rsd(j,i)*rsd(j,i)
                     epsp = epsp + rsdp(j,i)*rsdp(j,i)
                     epsds = epsds + rsds(j,i)*rsds(j,i)
                     epsps = epsps + rsdps(j,i)*rsdps(j,i)
                  end do
               end if
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
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
            if (iter .ge. politer)  done = .true.
c
c     apply a "peek" iteration to the mutual induced dipoles
c
            if (done) then
               do i = 1, npole
                  if (douind(ipole(i))) then
                     do j = 1, 3
                        uind(j,i) = uind(j,i) + poli(i)*rsd(j,i)
                        uinp(j,i) = uinp(j,i) + poli(i)*rsdp(j,i)
                        uinds(j,i) = uinds(j,i) + poli(i)*rsds(j,i)
                        uinps(j,i) = uinps(j,i) + poli(i)*rsdps(j,i)
                     end do
                  end if
               end do
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine dfield0d  --  generalized Kirkwood direct field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dfield0d" computes the direct electrostatic field due to
c     permanent multipole moments for use with with generalized
c     Kirkwood implicit solvation
c
c
      subroutine dfield0d (field,fieldp,fields,fieldps)
      use atoms
      use couple
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,uxi,uyi,uzi
      real*8 qxxi,qxyi,qxzi
      real*8 qyyi,qyzi,qzzi
      real*8 ck,uxk,uyk,uzk
      real*8 qxxk,qxyk,qxzk
      real*8 qyyk,qyzk,qzzk
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fc,fd,fq
      real*8 gf,gf2,gf3,gf5,gf7
      real*8 expterm,expc,expc1
      real*8 dexpc,expcdexpc
      real*8 a(0:3,0:2)
      real*8 gc(4),gux(10)
      real*8 guy(10),guz(10)
      real*8 gqxx(4),gqxy(4)
      real*8 gqxz(4),gqyy(4)
      real*8 gqyz(4),gqzz(4)
      real*8 fid(3),fkd(3)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
            fields(j,ii) = 0.0d0
            fieldps(j,ii) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factors for water
c
      dwater = 78.3d0
      fc = 1.0d0 * (1.0d0-dwater) / (1.0d0*dwater)
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
      fq = 3.0d0 * (1.0d0-dwater) / (2.0d0+3.0d0*dwater)
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
            fieldts(j,ii) = 0.0d0
            fieldtps(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,rpole,pdamp,thole,
!$OMP& dirdamp,rborn,n12,n13,n14,n15,np11,np12,np13,np14,i12,i13,i14,
!$OMP& i15,ip11,ip12,ip13,ip14,p2scale,p3scale,p4scale,p5scale,p2iscale,
!$OMP& p3iscale,p4iscale,p5iscale,d1scale,d2scale,d3scale,d4scale,
!$OMP& dpequal,use_intra,x,y,z,off2,fc,fd,fq,gkc,field,fieldp,fields,
!$OMP& fieldps)
!$OMP& firstprivate(dscale,pscale)
!$OMP& shared(fieldt,fieldtp,fieldts,fieldtps)
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps) schedule(guided)
c
c     find the field terms for each pairwise interaction
c
      do ii = 1, npole
         i = ipole(ii)
         ci = rpole(1,ii)
         uxi = rpole(2,ii)
         uyi = rpole(3,ii)
         uzi = rpole(4,ii)
         qxxi = rpole(5,ii)
         qxyi = rpole(6,ii)
         qxzi = rpole(7,ii)
         qyyi = rpole(9,ii)
         qyzi = rpole(10,ii)
         qzzi = rpole(13,ii)
         pdi = pdamp(ii)
         pti = thole(ii)
         ddi = dirdamp(ii)
         rbi = rborn(i)
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
         do kk = ii, npole
            k = ipole(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,kk)
                  uxk = rpole(2,kk)
                  uyk = rpole(3,kk)
                  uzk = rpole(4,kk)
                  qxxk = rpole(5,kk)
                  qxyk = rpole(6,kk)
                  qxzk = rpole(7,kk)
                  qyyk = rpole(9,kk)
                  qyzk = rpole(10,kk)
                  qzzk = rpole(13,kk)
                  rbk = rborn(k)
c
c     self-interactions for the solute field are skipped
c
                  if (i .ne. k) then
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(kk)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(ddi,dirdamp(kk))
                        if (pgamma .ne. 0.0d0) then
                           damp = pgamma * (r/damp)**(1.5d0)
                           if (damp .lt. 50.0d0) then
                              expdamp = exp(-damp) 
                              scale3 = 1.0d0 - expdamp 
                              scale5 = 1.0d0 - expdamp
     &                                            *(1.0d0+0.5d0*damp)
                              scale7 = 1.0d0 - expdamp
     &                                            *(1.0d0+0.65d0*damp
     &                                               +0.15d0*damp**2)
                           end if
                        else
                           pgamma = min(pti,thole(kk))
                           damp = pgamma * (r/damp)**3
                           if (damp .lt. 50.0d0) then
                              expdamp = exp(-damp)
                              scale3 = 1.0d0 - expdamp
                              scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                              scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                            +0.6d0*damp**2)
                           end if
                        end if
                     end if
                     rr3 = scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
                     dir = uxi*xr + uyi*yr + uzi*zr
                     qix = qxxi*xr + qxyi*yr + qxzi*zr
                     qiy = qxyi*xr + qyyi*yr + qyzi*zr
                     qiz = qxzi*xr + qyzi*yr + qzzi*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = uxk*xr + uyk*yr + uzk*zr
                     qkx = qxxk*xr + qxyk*yr + qxzk*zr
                     qky = qxyk*xr + qyyk*yr + qyzk*zr
                     qkz = qxzk*xr + qyzk*yr + qzzk*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uxk + 2.0d0*rr5*qkx
                     fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uyk + 2.0d0*rr5*qky
                     fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                           - rr3*uzk + 2.0d0*rr5*qkz
                     fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uxi - 2.0d0*rr5*qix
                     fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uyi - 2.0d0*rr5*qiy
                     fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                           - rr3*uzi - 2.0d0*rr5*qiz
                     do j = 1, 3
                        fieldt(j,ii) = fieldt(j,ii) + fid(j)*dscale(k)
                        fieldt(j,kk) = fieldt(j,kk) + fkd(j)*dscale(k)
                        fieldtp(j,ii) = fieldtp(j,ii) + fid(j)*pscale(k)
                        fieldtp(j,kk) = fieldtp(j,kk) + fkd(j)*pscale(k)
                     end do
                  end if
c
c     set the reaction potential auxiliary terms
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rb2)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  gf7 = gf5 * gf2
                  a(0,0) = gf
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  a(3,0) = -15.0d0 * gf7
c
c     set the reaction potential gradient auxiliary terms
c
                  expc1 = 1.0d0 - expc
                  a(0,1) = expc1 * a(1,0)
                  a(1,1) = expc1 * a(2,0)
                  a(2,1) = expc1 * a(3,0)
c
c     dipole second reaction potential gradient auxiliary term
c
                  expcdexpc = -expc * dexpc
                  a(1,2) = expc1*a(2,1) + expcdexpc*a(2,0)
c
c     multiply the auxiliary terms by dielectric functions
c
                  a(0,1) = fc * a(0,1)
                  a(1,0) = fd * a(1,0)
                  a(1,1) = fd * a(1,1)
                  a(1,2) = fd * a(1,2)
                  a(2,0) = fq * a(2,0)
                  a(2,1) = fq * a(2,1)
c
c     unweighted dipole reaction potential tensor
c
                  gux(1) = xr * a(1,0)
                  guy(1) = yr * a(1,0)
                  guz(1) = zr * a(1,0)
c
c     unweighted reaction potential gradient tensor
c
                  gc(2) = xr * a(0,1)
                  gc(3) = yr * a(0,1)
                  gc(4) = zr * a(0,1)
                  gux(2) = a(1,0) + xr2*a(1,1)
                  gux(3) = xr * yr * a(1,1)
                  gux(4) = xr * zr * a(1,1)
                  guy(2) = gux(3)
                  guy(3) = a(1,0) + yr2*a(1,1)
                  guy(4) = yr * zr * a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = a(1,0) + zr2*a(1,1)
                  gqxx(2) = xr * (2.0d0*a(2,0)+xr2*a(2,1))
                  gqxx(3) = yr * xr2*a(2,1)
                  gqxx(4) = zr * xr2*a(2,1)
                  gqyy(2) = xr * yr2*a(2,1)
                  gqyy(3) = yr * (2.0d0*a(2,0)+yr2*a(2,1))
                  gqyy(4) = zr * yr2 * a(2,1)
                  gqzz(2) = xr * zr2 * a(2,1)
                  gqzz(3) = yr * zr2 * a(2,1)
                  gqzz(4) = zr * (2.0d0*a(2,0)+zr2*a(2,1))
                  gqxy(2) = yr * (a(2,0)+xr2*a(2,1))
                  gqxy(3) = xr * (a(2,0)+yr2*a(2,1))
                  gqxy(4) = zr * xr * yr * a(2,1)
                  gqxz(2) = zr * (a(2,0)+xr2*a(2,1))
                  gqxz(3) = gqxy(4)
                  gqxz(4) = xr * (a(2,0)+zr2*a(2,1))
                  gqyz(2) = gqxy(4)
                  gqyz(3) = zr * (a(2,0)+yr2*a(2,1))
                  gqyz(4) = yr * (a(2,0)+zr2*a(2,1))
c
c     unweighted dipole second reaction potential gradient tensor
c
                  gux(5) = xr * (3.0d0*a(1,1)+xr2*a(1,2))
                  gux(6) = yr * (a(1,1)+xr2*a(1,2))
                  gux(7) = zr * (a(1,1)+xr2*a(1,2))
                  gux(8) = xr * (a(1,1)+yr2*a(1,2))
                  gux(9) = zr * xr * yr * a(1,2)
                  gux(10) = xr * (a(1,1)+zr2*a(1,2))
                  guy(5) = yr * (a(1,1)+xr2*a(1,2))
                  guy(6) = xr * (a(1,1)+yr2*a(1,2))
                  guy(7) = gux(9)
                  guy(8) = yr * (3.0d0*a(1,1)+yr2*a(1,2))
                  guy(9) = zr * (a(1,1)+yr2*a(1,2))
                  guy(10) = yr * (a(1,1)+zr2*a(1,2))
                  guz(5) = zr * (a(1,1)+xr2*a(1,2))
                  guz(6) = gux(9)
                  guz(7) = xr * (a(1,1)+zr2*a(1,2))
                  guz(8) = zr * (a(1,1)+yr2*a(1,2))
                  guz(9) = yr * (a(1,1)+zr2*a(1,2))
                  guz(10) = zr * (3.0d0*a(1,1)+zr2*a(1,2))
c
c     generalized Kirkwood permanent reaction field
c
                  fid(1) = uxk*gux(2) + uyk*gux(3) + uzk*gux(4)
     &                        + 0.5d0 * (ck*gux(1) + qxxk*gux(5)
     &                            + qyyk*gux(8) + qzzk*gux(10)
     &                            + 2.0d0*(qxyk*gux(6)+qxzk*gux(7)
     &                                         +qyzk*gux(9)))
     &                        + 0.5d0 * (ck*gc(2) + qxxk*gqxx(2)
     &                            + qyyk*gqyy(2) + qzzk*gqzz(2)
     &                            + 2.0d0*(qxyk*gqxy(2)+qxzk*gqxz(2)
     &                                         +qyzk*gqyz(2)))
                  fid(2) = uxk*guy(2) + uyk*guy(3) + uzk*guy(4)
     &                        + 0.5d0 * (ck*guy(1) + qxxk*guy(5)
     &                            + qyyk*guy(8) + qzzk*guy(10)
     &                            + 2.0d0*(qxyk*guy(6)+qxzk*guy(7)
     &                                         +qyzk*guy(9)))
     &                        + 0.5d0 * (ck*gc(3) + qxxk*gqxx(3)
     &                            + qyyk*gqyy(3) + qzzk*gqzz(3)
     &                            + 2.0d0*(qxyk*gqxy(3)+qxzk*gqxz(3)
     &                                         +qyzk*gqyz(3)))
                  fid(3) = uxk*guz(2) + uyk*guz(3) + uzk*guz(4)
     &                        + 0.5d0 * (ck*guz(1) + qxxk*guz(5)
     &                            + qyyk*guz(8) + qzzk*guz(10)
     &                            + 2.0d0*(qxyk*guz(6)+qxzk*guz(7)
     &                                         +qyzk*guz(9)))
     &                        + 0.5d0 * (ck*gc(4) + qxxk*gqxx(4)
     &                            + qyyk*gqyy(4) + qzzk*gqzz(4)
     &                            + 2.0d0*(qxyk*gqxy(4)+qxzk*gqxz(4)
     &                                         +qyzk*gqyz(4)))
                  fkd(1) = uxi*gux(2) + uyi*gux(3) + uzi*gux(4)
     &                        - 0.5d0 * (ci*gux(1) + qxxi*gux(5)
     &                            + qyyi*gux(8) + qzzi*gux(10)
     &                            + 2.0d0*(qxyi*gux(6)+qxzi*gux(7)
     &                                         +qyzi*gux(9)))
     &                        - 0.5d0 * (ci*gc(2) + qxxi*gqxx(2)
     &                            + qyyi*gqyy(2) + qzzi*gqzz(2)
     &                            + 2.0d0*(qxyi*gqxy(2)+qxzi*gqxz(2)
     &                                         +qyzi*gqyz(2)))
                  fkd(2) = uxi*guy(2) + uyi*guy(3) + uzi*guy(4)
     &                        - 0.5d0 * (ci*guy(1) + qxxi*guy(5)
     &                            + qyyi*guy(8) + qzzi*guy(10)
     &                            + 2.0d0*(qxyi*guy(6)+qxzi*guy(7)
     &                                         +qyzi*guy(9)))
     &                        - 0.5d0 * (ci*gc(3) + qxxi*gqxx(3)
     &                            + qyyi*gqyy(3) + qzzi*gqzz(3)
     &                            + 2.0d0*(qxyi*gqxy(3)+qxzi*gqxz(3)
     &                                         +qyzi*gqyz(3)))
                  fkd(3) = uxi*guz(2) + uyi*guz(3) + uzi*guz(4)
     &                        - 0.5d0 * (ci*guz(1) + qxxi*guz(5)
     &                            + qyyi*guz(8) + qzzi*guz(10)
     &                            + 2.0d0*(qxyi*guz(6)+qxzi*guz(7)
     &                                         +qyzi*guz(9)))
     &                        - 0.5d0 * (ci*gc(4) + qxxi*gqxx(4)
     &                            + qyyi*gqyy(4) + qzzi*gqzz(4)
     &                            + 2.0d0*(qxyi*gqxy(4)+qxzi*gqxz(4)
     &                                         +qyzi*gqyz(4)))
c
c     scale the self-field by half, such that it sums to one below
c
                  if (i .eq. k) then
                     do j = 1, 3
                        fid(j) = 0.5d0 * fid(j)
                        fkd(j) = 0.5d0 * fkd(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,ii) = fieldts(j,ii) + fid(j)
                     fieldts(j,kk) = fieldts(j,kk) + fkd(j)
                     fieldtps(j,ii) = fieldtps(j,ii) + fid(j)
                     fieldtps(j,kk) = fieldtps(j,kk) + fkd(j)
                  end do
               end if
            end if
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
!$OMP END DO
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + fieldt(j,ii)
            fieldp(j,ii) = fieldp(j,ii) + fieldtp(j,ii)
            fields(j,ii) = fields(j,ii) + fieldts(j,ii)
            fieldps(j,ii) = fieldps(j,ii) + fieldtps(j,ii)
         end do
      end do
!$OMP END DO
c
c     combine permanent multipole field and GK reaction field
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            fields(j,ii) = field(j,ii) + fields(j,ii)
            fieldps(j,ii) = fieldp(j,ii) + fieldps(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ufield0d  --  generalized Kirkwood mutual field  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ufield0d" computes the mutual electrostatic field due to
c     induced dipole moments for use with with generalized Kirkwood
c     implicit solvation
c
c
      subroutine ufield0d (field,fieldp,fields,fieldps)
      use atoms
      use gkstuf
      use group
      use mpole
      use polar
      use polgrp
      use polpot
      use shunt
      use solute
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,dukr
      real*8 puir,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 rb2,rbi,rbk
      real*8 dwater,fd
      real*8 gf,gf2,gf3,gf5
      real*8 expterm,expc
      real*8 expc1,dexpc
      real*8 a(0:3,0:2)
      real*8 gux(10),guy(10)
      real*8 guz(10)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8, allocatable :: uscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: fieldts(:,:)
      real*8, allocatable :: fieldtps(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
            fields(j,ii) = 0.0d0
            fieldps(j,ii) = 0.0d0
         end do
      end do
c
c     set dielectric constant and scaling factor for water
c
      dwater = 78.3d0
      fd = 2.0d0 * (1.0d0-dwater) / (1.0d0+2.0d0*dwater)
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
      allocate (fieldts(3,npole))
      allocate (fieldtps(3,npole))
c
c     initialize local variables for OpenMP calculation
c
      do ii = 1, npole
         do j = 1, 3
            fieldt(j,ii) = 0.0d0
            fieldtp(j,ii) = 0.0d0
            fieldts(j,ii) = 0.0d0
            fieldtps(j,ii) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,ipole,pdamp,thole,rborn,
!$OMP& uind,uinp,uinds,uinps,np11,np12,np13,np14,ip11,ip12,ip13,ip14,
!$OMP& u1scale,u2scale,u3scale,u4scale,use_intra,x,y,z,off2,fd,gkc,
!$OMP& field,fieldp,fields,fieldps)
!$OMP& firstprivate(uscale) shared(fieldt,fieldtp,fieldts,fieldtps)
!$OMP DO reduction(+:fieldt,fieldtp,fieldts,fieldtps) schedule(guided)
c
c     find the field terms for each pairwise interaction
c
      do ii = 1, npole
         i = ipole(ii)
         duix = uind(1,ii)
         duiy = uind(2,ii)
         duiz = uind(3,ii)
         puix = uinp(1,ii)
         puiy = uinp(2,ii)
         puiz = uinp(3,ii)
         duixs = uinds(1,ii)
         duiys = uinds(2,ii)
         duizs = uinds(3,ii)
         puixs = uinps(1,ii)
         puiys = uinps(2,ii)
         puizs = uinps(3,ii)
         pdi = pdamp(ii)
         pti = thole(ii)
         rbi = rborn(i)
c
c     set exclusion coefficients for connected atoms
c
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii, npole
            k = ipole(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,kk)
                  duky = uind(2,kk)
                  dukz = uind(3,kk)
                  pukx = uinp(1,kk)
                  puky = uinp(2,kk)
                  pukz = uinp(3,kk)
                  dukxs = uinds(1,kk)
                  dukys = uinds(2,kk)
                  dukzs = uinds(3,kk)
                  pukxs = uinps(1,kk)
                  pukys = uinps(2,kk)
                  pukzs = uinps(3,kk)
                  rbk = rborn(k)
                  if (i .ne. k) then
                     scale3 = uscale(k)
                     scale5 = uscale(k)
                     damp = pdi * pdamp(kk)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = scale3 * (1.0d0-expdamp)
                           scale5 = scale5 * (1.0d0-expdamp
     &                                           *(1.0d0+damp))
                        end if
                     end if
                     rr3 = -scale3 / (r*r2)
                     rr5 = 3.0d0 * scale5 / (r*r2*r2)
                     duir = xr*duix + yr*duiy + zr*duiz
                     dukr = xr*dukx + yr*duky + zr*dukz
                     puir = xr*puix + yr*puiy + zr*puiz
                     pukr = xr*pukx + yr*puky + zr*pukz
                     duirs = xr*duixs + yr*duiys + zr*duizs
                     dukrs = xr*dukxs + yr*dukys + zr*dukzs
                     puirs = xr*puixs + yr*puiys + zr*puizs
                     pukrs = xr*pukxs + yr*pukys + zr*pukzs
                     fid(1) = rr3*dukx + rr5*dukr*xr
                     fid(2) = rr3*duky + rr5*dukr*yr
                     fid(3) = rr3*dukz + rr5*dukr*zr
                     fkd(1) = rr3*duix + rr5*duir*xr
                     fkd(2) = rr3*duiy + rr5*duir*yr
                     fkd(3) = rr3*duiz + rr5*duir*zr
                     fip(1) = rr3*pukx + rr5*pukr*xr
                     fip(2) = rr3*puky + rr5*pukr*yr
                     fip(3) = rr3*pukz + rr5*pukr*zr
                     fkp(1) = rr3*puix + rr5*puir*xr
                     fkp(2) = rr3*puiy + rr5*puir*yr
                     fkp(3) = rr3*puiz + rr5*puir*zr
                     fids(1) = rr3*dukxs + rr5*dukrs*xr
                     fids(2) = rr3*dukys + rr5*dukrs*yr
                     fids(3) = rr3*dukzs + rr5*dukrs*zr
                     fkds(1) = rr3*duixs + rr5*duirs*xr
                     fkds(2) = rr3*duiys + rr5*duirs*yr
                     fkds(3) = rr3*duizs + rr5*duirs*zr
                     fips(1) = rr3*pukxs + rr5*pukrs*xr
                     fips(2) = rr3*pukys + rr5*pukrs*yr
                     fips(3) = rr3*pukzs + rr5*pukrs*zr
                     fkps(1) = rr3*puixs + rr5*puirs*xr
                     fkps(2) = rr3*puiys + rr5*puirs*yr
                     fkps(3) = rr3*puizs + rr5*puirs*zr
                     do j = 1, 3
                        fieldt(j,ii) = fieldt(j,ii) + fid(j)
                        fieldt(j,kk) = fieldt(j,kk) + fkd(j)
                        fieldtp(j,ii) = fieldtp(j,ii) + fip(j)
                        fieldtp(j,kk) = fieldtp(j,kk) + fkp(j)
                        fieldts(j,ii) = fieldts(j,ii) + fids(j)
                        fieldts(j,kk) = fieldts(j,kk) + fkds(j)
                        fieldtps(j,ii) = fieldtps(j,ii) + fips(j)
                        fieldtps(j,kk) = fieldtps(j,kk) + fkps(j)
                     end do
                  end if
c
c     unweighted dipole reaction potential gradient tensor
c
                  rb2 = rbi * rbk
                  expterm = exp(-r2/(gkc*rb2))
                  expc = expterm / gkc
                  dexpc = -2.0d0 / (gkc*rbi*rbk)
                  gf2 = 1.0d0 / (r2+rb2*expterm)
                  gf = sqrt(gf2)
                  gf3 = gf2 * gf
                  gf5 = gf3 * gf2
                  a(1,0) = -gf3
                  a(2,0) = 3.0d0 * gf5
                  expc1 = 1.0d0 - expc
                  a(1,1) = expc1 * a(2,0)
                  gux(2) = fd * (a(1,0) + xr2*a(1,1))
                  gux(3) = fd * xr*yr*a(1,1)
                  gux(4) = fd * xr*zr*a(1,1)
                  guy(2) = gux(3)
                  guy(3) = fd * (a(1,0) + yr2*a(1,1))
                  guy(4) = fd * yr*zr*a(1,1)
                  guz(2) = gux(4)
                  guz(3) = guy(4)
                  guz(4) = fd * (a(1,0) + zr2*a(1,1))
                  fids(1) = dukxs*gux(2) + dukys*guy(2) + dukzs*guz(2)
                  fids(2) = dukxs*gux(3) + dukys*guy(3) + dukzs*guz(3)
                  fids(3) = dukxs*gux(4) + dukys*guy(4) + dukzs*guz(4)
                  fkds(1) = duixs*gux(2) + duiys*guy(2) + duizs*guz(2)
                  fkds(2) = duixs*gux(3) + duiys*guy(3) + duizs*guz(3)
                  fkds(3) = duixs*gux(4) + duiys*guy(4) + duizs*guz(4)
                  fips(1) = pukxs*gux(2) + pukys*guy(2) + pukzs*guz(2)
                  fips(2) = pukxs*gux(3) + pukys*guy(3) + pukzs*guz(3)
                  fips(3) = pukxs*gux(4) + pukys*guy(4) + pukzs*guz(4)
                  fkps(1) = puixs*gux(2) + puiys*guy(2) + puizs*guz(2)
                  fkps(2) = puixs*gux(3) + puiys*guy(3) + puizs*guz(3)
                  fkps(3) = puixs*gux(4) + puiys*guy(4) + puizs*guz(4)
                  if (i .eq. k) then
                     do j = 1, 3
                        fids(j) = 0.5d0 * fids(j)
                        fkds(j) = 0.5d0 * fkds(j)
                        fips(j) = 0.5d0 * fips(j)
                        fkps(j) = 0.5d0 * fkps(j)
                     end do
                  end if
                  do j = 1, 3
                     fieldts(j,ii) = fieldts(j,ii) + fids(j)
                     fieldts(j,kk) = fieldts(j,kk) + fkds(j)
                     fieldtps(j,ii) = fieldtps(j,ii) + fips(j)
                     fieldtps(j,kk) = fieldtps(j,kk) + fkps(j)
                  end do
               end if
            end if
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
      end do
!$OMP END DO
c
c     add local to global variables for OpenMP calculation
c
!$OMP DO
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = field(j,ii) + fieldt(j,ii)
            fieldp(j,ii) = fieldp(j,ii) + fieldtp(j,ii)
            fields(j,ii) = fields(j,ii) + fieldts(j,ii)
            fieldps(j,ii) = fieldps(j,ii) + fieldtps(j,ii)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      deallocate (fieldts)
      deallocate (fieldtps)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine induce0d  --  Poisson-Boltzmann induced dipoles  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "induce0d" computes the induced dipole moments at polarizable
c     sites for Poisson-Boltzmann SCRF and vacuum environments
c
c
      subroutine induce0d
      use atoms
      use inform
      use iounit
      use mpole
      use polar
      use polopt
      use polpot
      use potent
      use units
      use uprior
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 polmin
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 epsds,epsps
      real*8 udsum,upsum
      real*8 ussum,upssum
      real*8 a,ap,as,aps
      real*8 b,bp,bs,bps
      real*8 sum,sump
      real*8 sums,sumps
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: fields(:,:)
      real*8, allocatable :: fieldps(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: rsds(:,:)
      real*8, allocatable :: rsdps(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: zrsds(:,:)
      real*8, allocatable :: zrsdps(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: conjs(:,:)
      real*8, allocatable :: conjps(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
      real*8, allocatable :: vecs(:,:)
      real*8, allocatable :: vecps(:,:)
      real*8, allocatable :: usum(:,:)
      real*8, allocatable :: usump(:,:)
      real*8, allocatable :: usums(:,:)
      real*8, allocatable :: usumps(:,:)
      logical done
      character*6 mode
c
c
c     zero out the induced dipoles; uind and uinp are vacuum dipoles,
c     uinds and uinps are Poisson-Boltzmann SCRF dipoles
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
            uinds(j,i) = 0.0d0
            uinps(j,i) = 0.0d0
         end do
      end do
      if (.not.use_polar .or. .not.use_solv)  return
c
c     set the switching function coefficients
c
      mode = 'MPOLE'
      call switch (mode)
c
c     perform dynamic allocation of some local arrays
c
      allocate (field(3,npole))
      allocate (fieldp(3,npole))
      allocate (fields(3,npole))
      allocate (fieldps(3,npole))
c
c     compute the direct induced dipole moment at each atom, and
c     another set that also includes RF due to permanent multipoles
c
      call dfield0e (field,fieldp,fields,fieldps)
c
c     set vacuum induced dipoles to polarizability times direct field;
c     SCRF induced dipoles are polarizability times direct field
c     plus the reaction field due to permanent multipoles
c
      do i = 1, npole
         if (douind(ipole(i))) then
            do j = 1, 3
               udir(j,i) = polarity(i) * field(j,i)
               udirp(j,i) = polarity(i) * fieldp(j,i)
               udirs(j,i) = polarity(i) * fields(j,i)
               udirps(j,i) = polarity(i) * fieldps(j,i)
               uind(j,i) = udir(j,i)
               uinp(j,i) = udirp(j,i)
               uinds(j,i) = udirs(j,i)
               uinps(j,i) = udirps(j,i)
            end do
         end if
      end do
c
c     get induced dipoles via the OPT extrapolation method
c
      if (poltyp .eq. 'OPT') then
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uopt(0,j,i) = udir(j,i)
                  uoptp(0,j,i) = udirp(j,i)
                  uopts(0,j,i) = udirs(j,i)
                  uoptps(0,j,i) = udirps(j,i)
               end do
            end if
         end do
         do k = 1, optorder
            call ufield0e (field,fieldp,fields,fieldps)
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uopt(k,j,i) = polarity(i) * field(j,i)
                     uoptp(k,j,i) = polarity(i) * fieldp(j,i)
                     uopts(k,j,i) = polarity(i) * fields(j,i)
                     uoptps(k,j,i) = polarity(i) * fieldps(j,i)
                     uind(j,i) = uopt(k,j,i)
                     uinp(j,i) = uoptp(k,j,i)
                     uinds(j,i) = uopts(k,j,i)
                     uinps(j,i) = uoptps(k,j,i)
                  end do
               end if
            end do
         end do
         allocate (usum(3,n))
         allocate (usump(3,n))
         allocate (usums(3,n))
         allocate (usumps(3,n))
         do i = 1, npole
            if (douind(ipole(i))) then
               do j = 1, 3
                  uind(j,i) = 0.0d0
                  uinp(j,i) = 0.0d0
                  uinds(j,i) = 0.0d0
                  uinps(j,i) = 0.0d0
                  usum(j,i) = 0.0d0
                  usump(j,i) = 0.0d0
                  usums(j,i) = 0.0d0
                  usumps(j,i) = 0.0d0
                  do k = 0, optorder
                     usum(j,i) = usum(j,i) + uopt(k,j,i)
                     usump(j,i) = usump(j,i) + uoptp(k,j,i)
                     usums(j,i) = usums(j,i) + uopts(k,j,i)
                     usumps(j,i) = usumps(j,i) + uoptps(k,j,i)
                     uind(j,i) = uind(j,i) + copt(k)*usum(j,i)
                     uinp(j,i) = uinp(j,i) + copt(k)*usump(j,i)
                     uinds(j,i) = uinds(j,i) + copt(k)*usums(j,i)
                     uinps(j,i) = uinps(j,i) + copt(k)*usumps(j,i)
                  end do
               end do
            end if
         end do
         deallocate (usum)
         deallocate (usump)
         deallocate (usums)
         deallocate (usumps)
      end if
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 100
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
c
         if (use_pred .and. nualt.eq.maxualt) then
            do i = 1, npole
               do j = 1, 3
                  udsum = 0.0d0
                  upsum = 0.0d0
                  ussum = 0.0d0
                  upssum = 0.0d0
                  do k = 1, nualt-1
                     udsum = udsum + bpred(k)*udalt(k,j,i)
                     upsum = upsum + bpredp(k)*upalt(k,j,i)
                     ussum = ussum + bpreds(k)*usalt(k,j,i)
                     upssum = upssum + bpredps(k)*upsalt(k,j,i)
                  end do
                  uind(j,i) = udsum
                  uinp(j,i) = upsum
                  uinds(j,i) = ussum
                  uinps(j,i) = upssum
               end do
            end do
         end if
c
c     perform dynamic allocation of some local arrays
c
         allocate (poli(npole))
         allocate (rsd(3,npole))
         allocate (rsdp(3,npole))
         allocate (rsds(3,npole))
         allocate (rsdps(3,npole))
         allocate (zrsd(3,npole))
         allocate (zrsdp(3,npole))
         allocate (zrsds(3,npole))
         allocate (zrsdps(3,npole))
         allocate (conj(3,npole))
         allocate (conjp(3,npole))
         allocate (conjs(3,npole))
         allocate (conjps(3,npole))
         allocate (vec(3,npole))
         allocate (vecp(3,npole))
         allocate (vecs(3,npole))
         allocate (vecps(3,npole))
c
c     set initial conjugate gradient residual and conjugate vector
c
         call ufield0e (field,fieldp,fields,fieldps)
         do i = 1, npole
            if (douind(ipole(i))) then
               poli(i) = max(polmin,polarity(i))
               do j = 1, 3
                  rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                          + field(j,i)
                  rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                           + fieldp(j,i)
                  rsds(j,i) = (udirs(j,i)-uinds(j,i))/poli(i)
     &                           + fields(j,i)
                  rsdps(j,i) = (udirps(j,i)-uinps(j,i))/poli(i)
     &                            + fieldps(j,i)
                  zrsd(j,i) = rsd(j,i) * poli(i)
                  zrsdp(j,i) = rsdp(j,i) * poli(i)
                  zrsds(j,i) = rsds(j,i) * poli(i)
                  zrsdps(j,i) = rsdps(j,i) * poli(i)
                  conj(j,i) = zrsd(j,i)
                  conjp(j,i) = zrsdp(j,i)
                  conjs(j,i) = zrsds(j,i)
                  conjps(j,i) = zrsdps(j,i)
               end do
            end if
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     vec(j,i) = uind(j,i)
                     vecp(j,i) = uinp(j,i)
                     vecs(j,i) = uinds(j,i)
                     vecps(j,i) = uinps(j,i)
                     uind(j,i) = conj(j,i)
                     uinp(j,i) = conjp(j,i)
                     uinds(j,i) = conjs(j,i)
                     uinps(j,i) = conjps(j,i)
                  end do
               end if
            end do
            call ufield0e (field,fieldp,fields,fieldps)
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = vec(j,i)
                     uinp(j,i) = vecp(j,i)
                     uinds(j,i) = vecs(j,i)
                     uinps(j,i) = vecps(j,i)
                     vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                     vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
                     vecs(j,i) = conjs(j,i)/poli(i) - fields(j,i)
                     vecps(j,i) = conjps(j,i)/poli(i) - fieldps(j,i)
                  end do
               end if
            end do
            a = 0.0d0
            ap = 0.0d0
            as = 0.0d0
            aps = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            sums = 0.0d0
            sumps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     a = a + conj(j,i)*vec(j,i)
                     ap = ap + conjp(j,i)*vecp(j,i)
                     as = as + conjs(j,i)*vecs(j,i)
                     aps = aps + conjps(j,i)*vecps(j,i)
                     sum = sum + rsd(j,i)*zrsd(j,i)
                     sump = sump + rsdp(j,i)*zrsdp(j,i)
                     sums = sums + rsds(j,i)*zrsds(j,i)
                     sumps = sumps + rsdps(j,i)*zrsdps(j,i)
                  end do
               end if
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            if (as .ne. 0.0d0)  as = sums / as
            if (aps .ne. 0.0d0)  aps = sumps / aps
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     uind(j,i) = uind(j,i) + a*conj(j,i)
                     uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                     uinds(j,i) = uinds(j,i) + as*conjs(j,i)
                     uinps(j,i) = uinps(j,i) + aps*conjps(j,i)
                     rsd(j,i) = rsd(j,i) - a*vec(j,i)
                     rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
                     rsds(j,i) = rsds(j,i) - as*vecs(j,i)
                     rsdps(j,i) = rsdps(j,i) - aps*vecps(j,i)
                  end do
               end if
            end do
            b = 0.0d0
            bp = 0.0d0
            bs = 0.0d0
            bps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     zrsd(j,i) = rsd(j,i) * poli(i)
                     zrsdp(j,i) = rsdp(j,i) * poli(i)
                     zrsds(j,i) = rsds(j,i) * poli(i)
                     zrsdps(j,i) = rsdps(j,i) * poli(i)
                     b = b + rsd(j,i)*zrsd(j,i)
                     bp = bp + rsdp(j,i)*zrsdp(j,i)
                     bs = bs + rsds(j,i)*zrsds(j,i)
                     bps = bps + rsdps(j,i)*zrsdps(j,i)
                  end do
               end if
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            if (sums .ne. 0.0d0)  bs = bs / sums
            if (sumps .ne. 0.0d0)  bps = bps / sumps
            epsd = 0.0d0
            epsp = 0.0d0
            epsds = 0.0d0
            epsps = 0.0d0
            do i = 1, npole
               if (douind(ipole(i))) then
                  do j = 1, 3
                     conj(j,i) = zrsd(j,i) + b*conj(j,i)
                     conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                     conjs(j,i) = zrsds(j,i) + bs*conjs(j,i)
                     conjps(j,i) = zrsdps(j,i) + bps*conjps(j,i)
                     epsd = epsd + rsd(j,i)*rsd(j,i)
                     epsp = epsp + rsdp(j,i)*rsdp(j,i)
                     epsds = epsds + rsds(j,i)*rsds(j,i)
                     epsps = epsps + rsdps(j,i)*rsdps(j,i)
                  end do
               end if
            end do
c
c     check the convergence of the mutual induced dipoles
c
            epsold = eps
            eps = max(epsd,epsp,epsds,epsps)
            eps = debye * sqrt(eps/dble(npolar))
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
            if (iter .ge. politer)  done = .true.
c
c     apply a "peek" iteration to the mutual induced dipoles
c
            if (done) then
               do i = 1, npole
                  if (douind(ipole(i))) then
                     do j = 1, 3
                        uind(j,i) = uind(j,i) + poli(i)*rsd(j,i)
                        uinp(j,i) = uinp(j,i) + poli(i)*rsdp(j,i)
                        uinds(j,i) = uinds(j,i) + poli(i)*rsds(j,i)
                        uinps(j,i) = uinps(j,i) + poli(i)*rsdps(j,i)
                     end do
                  end if
               end do
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (poli)
         deallocate (rsd)
         deallocate (rsdp)
         deallocate (rsds)
         deallocate (rsdps)
         deallocate (zrsd)
         deallocate (zrsdp)
         deallocate (zrsds)
         deallocate (zrsdps)
         deallocate (conj)
         deallocate (conjp)
         deallocate (conjs)
         deallocate (conjps)
         deallocate (vec)
         deallocate (vecp)
         deallocate (vecs)
         deallocate (vecps)
c
c     print the results from the conjugate gradient iteration
c
         if (debug) then
            write (iout,30)  iter,eps
   30       format (/,' Induced Dipoles :',6x,'Iterations',i5,
     &                 6x,'RMS Change',f15.10)
         end if
c
c     terminate the calculation if dipoles failed to converge
c
         if (iter.ge.maxiter .or. eps.gt.epsold) then
            write (iout,40)
   40       format (/,' INDUCE  --  Warning, Induced Dipoles',
     &                 ' are not Converged')
            call prterr
            call fatal
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (field)
      deallocate (fieldp)
      deallocate (fields)
      deallocate (fieldps)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0e  --  Poisson-Boltzmann direct field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0e" computes the direct electrostatic field due to
c     permanent multipole moments for use with in Poisson-Boltzmann
c
c
      subroutine dfield0e (field,fieldp,fields,fieldps)
      use atoms
      use couple
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,ddi
      real*8 pgamma
      real*8 fid(3),fkd(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
      end do
c
c     compute the direct electrostatic field at each atom, and
c     another field including RF due to permanent multipoles;
c     note self-interactions for the solute field are skipped
c
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
         pdi = pdamp(ii)
         pti = thole(ii)
         ddi = dirdamp(ii)
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
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
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
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(ddi,dirdamp(kk))
                     if (pgamma .ne. 0.0d0) then
                        damp = pgamma * (r/damp)**(1.5d0)
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp) 
                           scale3 = 1.0d0 - expdamp 
                           scale5 = 1.0d0 - expdamp*(1.0d0+0.5d0*damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+0.65d0*damp
     &                                         +0.15d0*damp**2)
                        end if
                     else
                        pgamma = min(pti,thole(kk))
                        damp = pgamma * (r/damp)**3
                        if (damp .lt. 50.0d0) then
                           expdamp = exp(-damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - expdamp*(1.0d0+damp)
                           scale7 = 1.0d0 - expdamp*(1.0d0+damp
     &                                         +0.6d0*damp**2)
                        end if
                     end if
                  end if
                  rr3 = scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
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
                  fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkx + 2.0d0*rr5*qkx
                  fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dky + 2.0d0*rr5*qky
                  fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                        - rr3*dkz + 2.0d0*rr5*qkz
                  fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*dix - 2.0d0*rr5*qix
                  fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diy - 2.0d0*rr5*qiy
                  fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                        - rr3*diz - 2.0d0*rr5*qiz
                  do j = 1, 3
                     field(j,ii) = field(j,ii) + fid(j)*dscale(k)
                     field(j,kk) = field(j,kk) + fkd(j)*dscale(k)
                     fieldp(j,ii) = fieldp(j,ii) + fid(j)*pscale(k)
                     fieldp(j,kk) = fieldp(j,kk) + fkd(j)*pscale(k)
                  end do
               end if
            end if
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
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
c
c     find the Poisson-Boltzmann reaction field at each site
c
      call pbempole
c
c     combine permanent multipole field and PB reaction field
c
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            fields(j,ii) = field(j,ii) + pbep(j,i)
            fieldps(j,ii) = fieldp(j,ii) + pbep(j,i)
         end do
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0e  --  Poisson-Boltzmann mutual field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0e" computes the mutual electrostatic field due to
c     induced dipole moments via a Poisson-Boltzmann solver
c
c
      subroutine ufield0e (field,fieldp,fields,fieldps)
      use atoms
      use group
      use mpole
      use pbstuf
      use polar
      use polgrp
      use polpot
      use shunt
      use solpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 xr2,yr2,zr2
      real*8 fgrp,r,r2
      real*8 rr3,rr5
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 duir,puir
      real*8 dukr,pukr
      real*8 duixs,duiys,duizs
      real*8 puixs,puiys,puizs
      real*8 dukxs,dukys,dukzs
      real*8 pukxs,pukys,pukzs
      real*8 duirs,puirs
      real*8 dukrs,pukrs
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8 fids(3),fkds(3)
      real*8 fips(3),fkps(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 fields(3,*)
      real*8 fieldps(3,*)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: indpole(:,:)
      real*8, allocatable :: inppole(:,:)
      logical proceed
c
c
c     zero out the value of the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            field(j,ii) = 0.0d0
            fieldp(j,ii) = 0.0d0
            fields(j,ii) = 0.0d0
            fieldps(j,ii) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (uscale(n))
c
c     set array needed to scale connected atom interactions
c
      do i = 1, n
         uscale(i) = 1.0d0
      end do
c
c     compute the mutual electrostatic field at each atom,
c     and another field including RF due to induced dipoles
c
      do ii = 1, npole
         i = ipole(ii)
         duix = uind(1,ii)
         duiy = uind(2,ii)
         duiz = uind(3,ii)
         puix = uinp(1,ii)
         puiy = uinp(2,ii)
         puiz = uinp(3,ii)
         duixs = uinds(1,ii)
         duiys = uinds(2,ii)
         duizs = uinds(3,ii)
         puixs = uinps(1,ii)
         puiys = uinps(2,ii)
         puizs = uinps(3,ii)
         pdi = pdamp(ii)
         pti = thole(ii)
c
c     set exclusion coefficients for connected atoms
c
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = ii+1, npole
            k = ipole(kk)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,i,k,0,0,0,0)
            if (proceed) then
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               xr2 = xr * xr
               yr2 = yr * yr
               zr2 = zr * zr
               r2 = xr2 + yr2 + zr2
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  dukx = uind(1,kk)
                  duky = uind(2,kk)
                  dukz = uind(3,kk)
                  pukx = uinp(1,kk)
                  puky = uinp(2,kk)
                  pukz = uinp(3,kk)
                  dukxs = uinds(1,kk)
                  dukys = uinds(2,kk)
                  dukzs = uinds(3,kk)
                  pukxs = uinps(1,kk)
                  pukys = uinps(2,kk)
                  pukzs = uinps(3,kk)
                  scale3 = uscale(k)
                  scale5 = uscale(k)
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = pgamma * (r/damp)**3
                     if (damp .lt. 50.0d0) then
                        expdamp = exp(-damp)
                        scale3 = scale3 * (1.0d0-expdamp)
                        scale5 = scale5 * (1.0d0-expdamp*(1.0d0+damp))
                     end if
                  end if
                  rr3 = -scale3 / (r*r2)
                  rr5 = 3.0d0 * scale5 / (r*r2*r2)
                  duir = xr*duix + yr*duiy + zr*duiz
                  dukr = xr*dukx + yr*duky + zr*dukz
                  puir = xr*puix + yr*puiy + zr*puiz
                  pukr = xr*pukx + yr*puky + zr*pukz
                  duirs = xr*duixs + yr*duiys + zr*duizs
                  dukrs = xr*dukxs + yr*dukys + zr*dukzs
                  puirs = xr*puixs + yr*puiys + zr*puizs
                  pukrs = xr*pukxs + yr*pukys + zr*pukzs
                  fid(1) = rr3*dukx + rr5*dukr*xr
                  fid(2) = rr3*duky + rr5*dukr*yr
                  fid(3) = rr3*dukz + rr5*dukr*zr
                  fkd(1) = rr3*duix + rr5*duir*xr
                  fkd(2) = rr3*duiy + rr5*duir*yr
                  fkd(3) = rr3*duiz + rr5*duir*zr
                  fip(1) = rr3*pukx + rr5*pukr*xr
                  fip(2) = rr3*puky + rr5*pukr*yr
                  fip(3) = rr3*pukz + rr5*pukr*zr
                  fkp(1) = rr3*puix + rr5*puir*xr
                  fkp(2) = rr3*puiy + rr5*puir*yr
                  fkp(3) = rr3*puiz + rr5*puir*zr
                  fids(1) = rr3*dukxs + rr5*dukrs*xr
                  fids(2) = rr3*dukys + rr5*dukrs*yr
                  fids(3) = rr3*dukzs + rr5*dukrs*zr
                  fkds(1) = rr3*duixs + rr5*duirs*xr
                  fkds(2) = rr3*duiys + rr5*duirs*yr
                  fkds(3) = rr3*duizs + rr5*duirs*zr
                  fips(1) = rr3*pukxs + rr5*pukrs*xr
                  fips(2) = rr3*pukys + rr5*pukrs*yr
                  fips(3) = rr3*pukzs + rr5*pukrs*zr
                  fkps(1) = rr3*puixs + rr5*puirs*xr
                  fkps(2) = rr3*puiys + rr5*puirs*yr
                  fkps(3) = rr3*puizs + rr5*puirs*zr
                  do j = 1, 3
                     field(j,ii) = field(j,ii) + fid(j)
                     field(j,kk) = field(j,kk) + fkd(j)
                     fieldp(j,ii) = fieldp(j,ii) + fip(j)
                     fieldp(j,kk) = fieldp(j,kk) + fkp(j)
                     fields(j,ii) = fields(j,ii) + fids(j)
                     fields(j,kk) = fields(j,kk) + fkds(j)
                     fieldps(j,ii) = fieldps(j,ii) + fips(j)
                     fieldps(j,kk) = fieldps(j,kk) + fkps(j)
                  end do
               end if
            end if
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
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (uscale)
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pbeuind))  allocate (pbeuind(3,n))
      if (.not. allocated(pbeuinp))  allocate (pbeuinp(3,n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (indpole(3,n))
      allocate (inppole(3,n))
c
c     zero out the PB reaction field at each atomic site
c
      do i = 1, n
         do j = 1, 3
            indpole(j,i) = 0.0d0
            inppole(j,i) = 0.0d0
            pbeuind(j,i) = 0.0d0
            pbeuinp(j,i) = 0.0d0
         end do
      end do
c
c     find the Poisson-Boltzmann reaction field at each site
c
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            indpole(j,i) = uinds(j,ii)
            inppole(j,i) = uinps(j,ii)
         end do
      end do
      call apbsinduce (indpole,pbeuind)
      call apbsnlinduce (inppole,pbeuinp)
c
c     perform deallocation of some local arrays
c
      deallocate (indpole)
      deallocate (inppole)
c
c     combine mutual induced dipole field and PB reaction field
c
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            fields(j,ii) = fields(j,ii) + pbeuind(j,i)
            fieldps(j,ii) = fieldps(j,ii) + pbeuinp(j,i)
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ulspred  --  induced dipole prediction coeffs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ulspred" uses standard extrapolation or a least squares fit
c     to set coefficients of an induced dipole predictor polynomial
c
c     literature references:
c
c     J. Kolafa, "Time-Reversible Always Stable Predictor-Corrector
c     Method for Molecular Dynamics of Polarizable Molecules", Journal
c     of Computational Chemistry, 25, 335-342 (2004)
c
c     W. Wang and R. D. Skeel, "Fast Evaluation of Polarizable Forces",
c     Journal of Chemical Physics, 123, 164107 (2005)
c
c
      subroutine ulspred
      use mpole
      use uprior
      implicit none
      integer i,j,k,m
      real*8 coeff,udk,upk
      real*8 amax,apmax
      real*8 b(maxualt)
      real*8 bp(maxualt)
      real*8 a(maxualt*(maxualt+1)/2)
      real*8 ap(maxualt*(maxualt+1)/2)
      real*8 c(maxualt,maxualt)
      real*8 cp(maxualt,maxualt)
c
c
c     set the Gear predictor binomial coefficients
c
      if (polpred .eq. 'GEAR') then
         do i = 1, nualt
            coeff = gear(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     set always stable predictor-corrector (ASPC) coefficients
c
      else if (polpred .eq. 'ASPC') then
         do i = 1, nualt
            coeff = aspc(i)
            bpred(i) = coeff
            bpredp(i) = coeff
            bpreds(i) = coeff
            bpredps(i) = coeff
         end do
c
c     derive normal equations corresponding to least squares fit
c
      else
         do k = 1, nualt
            b(k) = 0.0d0
            bp(k) = 0.0d0
            do m = k, nualt
               c(k,m) = 0.0d0
               cp(k,m) = 0.0d0
            end do
         end do
         do i = 1, npole
            do j = 1, 3
               do k = 1, nualt
                  udk = udalt(k,j,i)
                  upk = upalt(k,j,i)
                  do m = k, nualt
                     c(k,m) = c(k,m) + udk*udalt(m,j,i)
                     cp(k,m) = cp(k,m) + upk*upalt(m,j,i)
                  end do
               end do
            end do
         end do
         i = 0
         do k = 2, nualt
            b(k-1) = c(1,k)
            bp(k-1) = cp(1,k)
            do m = k, nualt
               i = i + 1
               a(i) = c(k,m)
               ap(i) = cp(k,m)
            end do
         end do
c
c     check for nonzero coefficients and solve normal equations
c
         k = nualt - 1
         amax = 0.0d0
         apmax = 0.0d0
         do i = 1, k*(k+1)/2
            amax = max(amax,a(i))
            apmax = max(apmax,ap(i))
         end do
         if (amax .ne. 0.0d0)  call cholesky (k,a,b)
         if (apmax .ne. 0.0d0)  call cholesky (k,ap,bp)
c
c     transfer the final solution to the coefficient vector
c
         do k = 1, nualt-1
            bpred(k) = b(k)
            bpredp(k) = bp(k)
            bpreds(k) = b(k)
            bpredps(k) = bp(k)
         end do
         bpred(nualt) = 0.0d0
         bpredp(nualt) = 0.0d0
         bpreds(nualt) = 0.0d0
         bpredps(nualt) = 0.0d0
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0a  --  dipole preconditioner via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0a" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a double loop
c
c
      subroutine uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
      use atoms
      use chgpen
      use couple
      use limits
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpcg
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 pgamma,off2
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do ii = 1, npole
            poli = udiag * max(polmin,polarity(ii))
            do j = 1, 3
               zrsd(j,ii) = poli * rsd(j,ii)
               zrsdp(j,ii) = poli * rsdp(j,ii)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
         off2 = usolvcut * usolvcut
         j = 0
         do ii = 1, npole-1
            i = ipole(ii)
            do kk = ii+1, npole
               k = ipole(kk)
               xr = x(k) - x(i)
               yr = y(k) - y(i)
               zr = z(k) - z(i)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  m1 = minv(j+1)
                  m2 = minv(j+2)
                  m3 = minv(j+3)
                  m4 = minv(j+4)
                  m5 = minv(j+5)
                  m6 = minv(j+6)
                  j = j + 6
                  zrsd(1,ii) = zrsd(1,ii) + m1*rsd(1,kk)
     &                            + m2*rsd(2,kk) + m3*rsd(3,kk)
                  zrsd(2,ii) = zrsd(2,ii) + m2*rsd(1,kk)
     &                            + m4*rsd(2,kk) + m5*rsd(3,kk)
                  zrsd(3,ii) = zrsd(3,ii) + m3*rsd(1,kk)
     &                            + m5*rsd(2,kk) + m6*rsd(3,kk)
                  zrsd(1,kk) = zrsd(1,kk) + m1*rsd(1,ii)
     &                            + m2*rsd(2,ii) + m3*rsd(3,ii)
                  zrsd(2,kk) = zrsd(2,kk) + m2*rsd(1,ii)
     &                            + m4*rsd(2,ii) + m5*rsd(3,ii)
                  zrsd(3,kk) = zrsd(3,kk) + m3*rsd(1,ii)
     &                            + m5*rsd(2,ii) + m6*rsd(3,ii)
                  zrsdp(1,ii) = zrsdp(1,ii) + m1*rsdp(1,kk)
     &                             + m2*rsdp(2,kk) + m3*rsdp(3,kk)
                  zrsdp(2,ii) = zrsdp(2,ii) + m2*rsdp(1,kk)
     &                             + m4*rsdp(2,kk) + m5*rsdp(3,kk)
                  zrsdp(3,ii) = zrsdp(3,ii) + m3*rsdp(1,kk)
     &                             + m5*rsdp(2,kk) + m6*rsdp(3,kk)
                  zrsdp(1,kk) = zrsdp(1,kk) + m1*rsdp(1,ii)
     &                             + m2*rsdp(2,ii) + m3*rsdp(3,ii)
                  zrsdp(2,kk) = zrsdp(2,kk) + m2*rsdp(1,ii)
     &                             + m4*rsdp(2,ii) + m5*rsdp(3,ii)
                  zrsdp(3,kk) = zrsdp(3,kk) + m3*rsdp(1,ii)
     &                             + m5*rsdp(2,ii) + m6*rsdp(3,ii)
               end if
            end do
         end do
c
c     construct off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (uscale(n))
         allocate (wscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            uscale(i) = 1.0d0
            wscale(i) = 1.0d0
         end do
c
c     determine the off-diagonal elements of the preconditioner
c
         off2 = usolvcut * usolvcut
         m = 0
         do ii = 1, npole-1
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            poli = polarity(ii)
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
c
c     evaluate all sites within the cutoff distance
c
            do kk = ii+1, npole
               k = ipole(kk)
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
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
                           scale5 = scale5 * (1.0d0-expdamp
     &                                           *(1.0d0-damp))
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
                  polik = poli * polarity(kk)
                  rr3 = scale3 * polik / (r*r2)
                  rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
                  minv(m+1) = rr5*xr*xr - rr3
                  minv(m+2) = rr5*xr*yr
                  minv(m+3) = rr5*xr*zr
                  minv(m+4) = rr5*yr*yr - rr3
                  minv(m+5) = rr5*yr*zr
                  minv(m+6) = rr5*zr*zr - rr3
                  m = m + 6
               end if
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
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0b  --  dipole preconditioner via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0b" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a neighbor pair list
c
c
      subroutine uscale0b (mode,rsd,rsdp,zrsd,zrsdp)
      use atoms
      use chgpen
      use couple
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpcg
      use polpot
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 pgamma
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8 dmpik(5)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8 rsd(3,*)
      real*8 rsdp(3,*)
      real*8 zrsd(3,*)
      real*8 zrsdp(3,*)
      real*8, allocatable :: zrsdt(:,:)
      real*8, allocatable :: zrsdtp(:,:)
      character*6 mode
c
c
c     apply the preconditioning matrix to the current residual
c
      if (mode .eq. 'APPLY') then
c
c     perform dynamic allocation of some local arrays
c
         allocate (zrsdt(3,npole))
         allocate (zrsdtp(3,npole))
c
c     use diagonal preconditioner elements as first approximation
c
         polmin = 0.00000001d0
         do ii = 1, npole
            poli = udiag * max(polmin,polarity(ii))
            do j = 1, 3
               zrsd(j,ii) = poli * rsd(j,ii)
               zrsdp(j,ii) = poli * rsdp(j,ii)
               zrsdt(j,ii) = 0.0d0
               zrsdtp(j,ii) = 0.0d0
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst,ulst,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do ii = 1, npole
            m = mindex(ii)
            do kkk = 1, nulst(ii)
               kk = ulst(kkk,ii)
               m1 = minv(m+1)
               m2 = minv(m+2)
               m3 = minv(m+3)
               m4 = minv(m+4)
               m5 = minv(m+5)
               m6 = minv(m+6)
               m = m + 6
               zrsdt(1,ii) = zrsdt(1,ii) + m1*rsd(1,kk)
     &                          + m2*rsd(2,kk) + m3*rsd(3,kk)
               zrsdt(2,ii) = zrsdt(2,ii) + m2*rsd(1,kk)
     &                          + m4*rsd(2,kk) + m5*rsd(3,kk)
               zrsdt(3,ii) = zrsdt(3,ii) + m3*rsd(1,kk)
     &                          + m5*rsd(2,kk) + m6*rsd(3,kk)
               zrsdt(1,kk) = zrsdt(1,kk) + m1*rsd(1,ii)
     &                          + m2*rsd(2,ii) + m3*rsd(3,ii)
               zrsdt(2,kk) = zrsdt(2,kk) + m2*rsd(1,ii)
     &                          + m4*rsd(2,ii) + m5*rsd(3,ii)
               zrsdt(3,kk) = zrsdt(3,kk) + m3*rsd(1,ii)
     &                          + m5*rsd(2,ii) + m6*rsd(3,ii)
               zrsdtp(1,ii) = zrsdtp(1,ii) + m1*rsdp(1,kk)
     &                           + m2*rsdp(2,kk) + m3*rsdp(3,kk)
               zrsdtp(2,ii) = zrsdtp(2,ii) + m2*rsdp(1,kk)
     &                           + m4*rsdp(2,kk) + m5*rsdp(3,kk)
               zrsdtp(3,ii) = zrsdtp(3,ii) + m3*rsdp(1,kk)
     &                           + m5*rsdp(2,kk) + m6*rsdp(3,kk)
               zrsdtp(1,kk) = zrsdtp(1,kk) + m1*rsdp(1,ii)
     &                           + m2*rsdp(2,ii) + m3*rsdp(3,ii)
               zrsdtp(2,kk) = zrsdtp(2,kk) + m2*rsdp(1,ii)
     &                           + m4*rsdp(2,ii) + m5*rsdp(3,ii)
               zrsdtp(3,kk) = zrsdtp(3,kk) + m3*rsdp(1,ii)
     &                           + m5*rsdp(2,ii) + m6*rsdp(3,ii)
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do ii = 1, npole
            do j = 1, 3
               zrsd(j,ii) = zrsd(j,ii) + zrsdt(j,ii)
               zrsdp(j,ii) = zrsdp(j,ii) + zrsdtp(j,ii)
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (zrsdt)
         deallocate (zrsdtp)
c
c     build the off-diagonal elements of preconditioning matrix
c
      else if (mode .eq. 'BUILD') then
         m = 0
         do ii = 1, npole
            mindex(ii) = m
            m = m + 6*nulst(ii)
         end do
c
c     perform dynamic allocation of some local arrays
c
         allocate (uscale(n))
         allocate (wscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            uscale(i) = 1.0d0
            wscale(i) = 1.0d0
         end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,polarity,
!$OMP& pdamp,thole,pcore,pval,palpha,u1scale,u2scale,u3scale,u4scale,
!$OMP& w2scale,w3scale,w4scale,w5scale,n12,i12,n13,i13,n14,i14,n15,i15,
!$OMP& np11,ip11,np12,ip12,np13,ip13,np14,ip14,use_thole,use_chgpen,
!$OMP& nulst,ulst,mindex,minv)
!$OMP& firstprivate (uscale,wscale)
c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)
         do ii = 1, npole
            i = ipole(ii)
            xi = x(i)
            yi = y(i)
            zi = z(i)
            poli = polarity(ii)
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
c
c     evaluate all sites within the cutoff distance
c
            m = mindex(ii)
            do kkk = 1, nulst(ii)
               kk = ulst(kkk,ii)
               k = ipole(kk)
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               r = sqrt(r2)
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
                        scale5 = scale5 * (1.0d0-expdamp
     &                                        *(1.0d0-damp))
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
               polik = poli * polarity(kk)
               rr3 = scale3 * polik / (r*r2)
               rr5 = 3.0d0 * scale5 * polik / (r*r2*r2)
               minv(m+1) = rr5*xr*xr - rr3
               minv(m+2) = rr5*xr*yr
               minv(m+3) = rr5*xr*zr
               minv(m+4) = rr5*yr*yr - rr3
               minv(m+5) = rr5*yr*zr
               minv(m+6) = rr5*zr*zr - rr3
               m = m + 6
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
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (uscale)
         deallocate (wscale)
      end if
      return
      end
