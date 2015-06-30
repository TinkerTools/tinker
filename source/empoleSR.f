c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine empoleSR  --  Ewald multipole derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "empoleSR" calculates the multipole and dipole polarization
c     energy and derivatives with respect to Cartesian coordinates
c     using particle mesh Ewald summation and a neighbor list
c
c
      subroutine empoleSR
      use sizes
      use atoms
      use boxes
      use chgpot
      use deriv
      use energi
      use ewald
      use inter
      use math
      use mpole
      use polar
      use polpot
      use virial
      implicit none
      integer i,j,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 xup,yup,zup
      real*8 xq,yq,zq
      real*8 xv,yv,zv,vterm
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 xdfield,xufield
      real*8 ydfield,yufield
      real*8 zdfield,zufield
      real*8 trq(3),trqi(3)
      real*8 frcx(3),frcy(3),frcz(3)
c
c
c     zero out multipole and polarization energy and derivatives
c
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         do j = 1, 3
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
         end do
      end do
      if (npole .eq. 0)  return
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c     should not need to change this
      call chkpole
c
c     rotate the multipole components into the global frame
c     should not need to change this
      call rotpole
c
c     compute the induced dipole moment at each atom
c     induce0a copied below -- subroutine for explicit solvent
      call induce0aSR
c
c     compute the real space part of the Ewald summation
c     ereal1d copied below -- uses lists
      call ereal1dSR (eintra)

c
c     compute the self-energy torque term due to induced dipole
c
      trq(1) = 0.0d0
      trq(2) = 0.0d0
      trq(3) = 0.0d0
      term = (4.0d0/3.0d0) * f * aewald**3 / sqrtpi
      do i = 1, npole
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         uix = 0.5d0 * (uind(1,i)+uinp(1,i))
         uiy = 0.5d0 * (uind(2,i)+uinp(2,i))
         uiz = 0.5d0 * (uind(3,i)+uinp(3,i))
         trqi(1) = term * (diy*uiz-diz*uiy)
         trqi(2) = term * (diz*uix-dix*uiz)
         trqi(3) = term * (dix*uiy-diy*uix)
         call torque (i,trq,trqi,frcx,frcy,frcz)
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
         xup = 0.0d0
         yup = 0.0d0
         zup = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i) + rpole(1,i)*x(ii)
            yd = yd + rpole(3,i) + rpole(1,i)*y(ii)
            zd = zd + rpole(4,i) + rpole(1,i)*z(ii)
            xu = xu + uind(1,i)
            yu = yu + uind(2,i)
            zu = zu + uind(3,i)
            xup = xup + uinp(1,i)
            yup = yup + uinp(2,i)
            zup = zup + uinp(3,i)
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
         do i = 1, npole
            ii = ipole(i)
            dem(1,ii) = dem(1,ii) + 2.0d0*term*rpole(1,i)*xd
            dem(2,ii) = dem(2,ii) + 2.0d0*term*rpole(1,i)*yd
            dem(3,ii) = dem(3,ii) + 2.0d0*term*rpole(1,i)*zd
            dep(1,ii) = dep(1,ii) + term*rpole(1,i)*(xu+xup)
            dep(2,ii) = dep(2,ii) + term*rpole(1,i)*(yu+yup)
            dep(3,ii) = dep(3,ii) + term*rpole(1,i)*(zu+zup)
         end do
         xdfield = -2.0d0 * term * xd
         ydfield = -2.0d0 * term * yd
         zdfield = -2.0d0 * term * zd
         xufield = -term * (xu+xup)
         yufield = -term * (yu+yup)
         zufield = -term * (zu+zup)
         do i = 1, npole
            trq(1) = rpole(3,i)*zdfield - rpole(4,i)*ydfield
            trq(2) = rpole(4,i)*xdfield - rpole(2,i)*zdfield
            trq(3) = rpole(2,i)*ydfield - rpole(3,i)*xdfield
            trqi(1) = rpole(3,i)*zufield - rpole(4,i)*yufield
            trqi(2) = rpole(4,i)*xufield - rpole(2,i)*zufield
            trqi(3) = rpole(2,i)*yufield - rpole(3,i)*xufield
            call torque (i,trq,trqi,frcx,frcy,frcz)
         end do
c
c     boundary correction to virial due to overall cell dipole
c
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xq = 0.0d0
         yq = 0.0d0
         zq = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            xd = xd + rpole(2,i)
            yd = yd + rpole(3,i)
            zd = zd + rpole(4,i)
            xq = xq + rpole(1,i)*x(ii)
            yq = yq + rpole(1,i)*y(ii)
            zq = zq + rpole(1,i)*z(ii)
         end do
         xv = xq * (xd+0.5d0*(xu+xup))
         yv = yq * (yd+0.5d0*(yu+yup))
         zv = zq * (zd+0.5d0*(zu+zup))
         vterm = term * (xq*xq + yq*yq + zq*zq + 2.0d0*(xv+yv+zv)
     &                      + xu*xup + yu*yup + zu*zup
     &                      + xd*(xd+xu+xup) + yd*(yd+yu+yup)
     &                      + zd*(zd+zu+zup))
         vir(1,1) = vir(1,1) + 2.0d0*term*(xq*xq+xv) + vterm
         vir(2,1) = vir(2,1) + 2.0d0*term*(xq*yq+xv)
         vir(3,1) = vir(3,1) + 2.0d0*term*(xq*zq+xv)
         vir(1,2) = vir(1,2) + 2.0d0*term*(yq*xq+yv)
         vir(2,2) = vir(2,2) + 2.0d0*term*(yq*yq+yv) + vterm
         vir(3,2) = vir(3,2) + 2.0d0*term*(yq*zq+yv)
         vir(1,3) = vir(1,3) + 2.0d0*term*(zq*xq+zv)
         vir(2,3) = vir(2,3) + 2.0d0*term*(zq*yq+zv)
         vir(3,3) = vir(3,3) + 2.0d0*term*(zq*zq+zv) + vterm
         if (poltyp .eq. 'DIRECT') then
            vterm = term * (xu*xup+yu*yup+zu*zup)
            vir(1,1) = vir(1,1) + vterm
            vir(2,2) = vir(2,2) + vterm
            vir(3,3) = vir(3,3) + vterm
         end if
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine induce0aSR  --  conjugate gradient dipole solver  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "induce0aSR" computes the induced dipole moments at polarizable
c     sites using a preconditioned conjugate gradient solver
c
c
      subroutine induce0aSR
      use sizes
      use atoms
      use inform
      use iounit
      use limits
      use mpole
      use polar
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
      real*8 sum,sump
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      real*8, allocatable :: udir(:,:)
      real*8, allocatable :: udirp(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: rsdp(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: zrsdp(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: conjp(:,:)
      real*8, allocatable :: vec(:,:)
      real*8, allocatable :: vecp(:,:)
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
      allocate (udir(3,npole))
      allocate (udirp(3,npole))
c
c     get the electrostatic field due to permanent multipoles
c
c      if (use_ewald) then
         call dfield0cSR (field,fieldp)
c      else if (use_mlist) then
c         call dfield0b (field,fieldp)
c      else
c         call dfield0a (field,fieldp)
c      end if
c
c     set induced dipoles to polarizability times direct field
c
      do i = 1, npole
         do j = 1, 3
            udir(j,i) = polarity(i) * field(j,i)
            udirp(j,i) = polarity(i) * fieldp(j,i)
            uind(j,i) = udir(j,i)
            uinp(j,i) = udirp(j,i)
         end do
      end do
c
c     set tolerances for computation of mutual induced dipoles
c
      if (poltyp .eq. 'MUTUAL') then
         done = .false.
         maxiter = 500
         iter = 0
         polmin = 0.00000001d0
         eps = 100.0d0
c
c     estimated induced dipoles from polynomial predictor
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
c         if (use_ewald) then
            call ufield0cSR (field,fieldp)
c         else if (use_mlist) then
c            call ufield0b (field,fieldp)
c         else
c            call ufield0a (field,fieldp)
c         end if
c
c     set initial conjugate gradient residual and conjugate vector
c
         do i = 1, npole
            poli(i) = max(polmin,polarity(i))
            do j = 1, 3
               rsd(j,i) = (udir(j,i)-uind(j,i))/poli(i)
     &                       + field(j,i)
               rsdp(j,i) = (udirp(j,i)-uinp(j,i))/poli(i)
     &                       + fieldp(j,i)
            end do
         end do
         mode = 'BUILD'
c         if (use_mlist) then
            call uscale0bSR (mode,rsd,rsdp,zrsd,zrsdp)
            mode = 'APPLY'
            call uscale0bSR (mode,rsd,rsdp,zrsd,zrsdp)
c         else
c            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c            mode = 'APPLY'
c            call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c         end if
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i)
               conjp(j,i) = zrsdp(j,i)
            end do
         end do
c
c     conjugate gradient iteration of the mutual induced dipoles
c
         do while (.not. done)
            iter = iter + 1
            do i = 1, npole
               do j = 1, 3
                  vec(j,i) = uind(j,i)
                  vecp(j,i) = uinp(j,i)
                  uind(j,i) = conj(j,i)
                  uinp(j,i) = conjp(j,i)
               end do
            end do
c            if (use_ewald) then
               call ufield0cSR (field,fieldp)
c            else if (use_mlist) then
c               call ufield0b (field,fieldp)
c            else
c               call ufield0a (field,fieldp)
c            end if
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = vec(j,i)
                  uinp(j,i) = vecp(j,i)
                  vec(j,i) = conj(j,i)/poli(i) - field(j,i)
                  vecp(j,i) = conjp(j,i)/poli(i) - fieldp(j,i)
               end do
            end do
            a = 0.0d0
            ap = 0.0d0
            sum = 0.0d0
            sump = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  a = a + conj(j,i)*vec(j,i)
                  ap = ap + conjp(j,i)*vecp(j,i)
                  sum = sum + rsd(j,i)*zrsd(j,i)
                  sump = sump + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (a .ne. 0.0d0)  a = sum / a
            if (ap .ne. 0.0d0)  ap = sump / ap
            do i = 1, npole
               do j = 1, 3
                  uind(j,i) = uind(j,i) + a*conj(j,i)
                  uinp(j,i) = uinp(j,i) + ap*conjp(j,i)
                  rsd(j,i) = rsd(j,i) - a*vec(j,i)
                  rsdp(j,i) = rsdp(j,i) - ap*vecp(j,i)
               end do
            end do
c            if (use_mlist) then
               call uscale0bSR (mode,rsd,rsdp,zrsd,zrsdp)
c            else
c               call uscale0a (mode,rsd,rsdp,zrsd,zrsdp)
c            end if
            b = 0.0d0
            bp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  b = b + rsd(j,i)*zrsd(j,i)
                  bp = bp + rsdp(j,i)*zrsdp(j,i)
               end do
            end do
            if (sum .ne. 0.0d0)  b = b / sum
            if (sump .ne. 0.0d0)  bp = bp / sump
            epsd = 0.0d0
            epsp = 0.0d0
            do i = 1, npole
               do j = 1, 3
                  conj(j,i) = zrsd(j,i) + b*conj(j,i)
                  conjp(j,i) = zrsdp(j,i) + bp*conjp(j,i)
                  epsd = epsd + rsd(j,i)*rsd(j,i)
                  epsp = epsp + rsdp(j,i)*rsdp(j,i)
               end do
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
   10             format (/,' Determination of Induced Dipole',
     &                       ' Moments :',
     &                    //,4x,'Iter',8x,'RMS Change (Debyes)',/)
               end if
               write (iout,20)  iter,eps
   20          format (i8,7x,f16.10)
            end if
            if (eps .lt. poleps)  done = .true.
            if (eps .gt. epsold)  done = .true.
            if (iter .ge. politer)  done = .true.
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
      deallocate (udir)
      deallocate (udirp)
      return
      end
c

c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ereal1dSR  --  ewald real space derivs via list  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ereal1dSR" evaluates the real space portion of the regular Ewald
c     summation energy and gradient due to atomic multipole interactions
c     and dipole polarizability
c
c
      subroutine ereal1dSR (eintra)
      use sizes
      use atoms
      use bound
      use boxes
      use cell
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use limits
      use math
      use mdstuf
      use molcul
      use mplpot
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use shunt
      use virial
      implicit none
      integer i,j,k
      integer ii,kk,kkk
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,bfac
      real*8 eintra,erfc
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale5
      real*8 scale7
      real*8 temp3,temp5,temp7
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 usc3,usc5
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 gfd,gfdr
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 erl,erli
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 emo,epo,eintrao
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ftm2r(3),ftm2ri(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 ttm2r(3),ttm3r(3)
      real*8 ttm2ri(3),ttm3ri(3)
      real*8 fdir(3),dixdk(3)
      real*8 dkxui(3),dixuk(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 bn(0:5)
      real*8 sc(10),gl(0:8)
      real*8 sci(8),scip(8)
      real*8 gli(7),glip(7)
      real*8 gf(7),gfi(6)
      real*8 gfr(7),gfri(6)
      real*8 gti(6),gtri(6)
      real*8 viro(3,3)
      real*8 REmpolecut2,REmpoletaper2
      real*8 res_lam,res_u,ru2,ru3,ru4,ru5,q_switch
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: demo1(:,:)
      real*8, allocatable :: demo2(:,:)
      real*8, allocatable :: depo1(:,:)
      real*8, allocatable :: depo2(:,:)
      logical dorl,dorli
      character*6 mode
      external erfc
      REmpolecut2 = REmpolecut**2
      REmpoletaper2 = REmpoletaper**2
      res_lam = REmpolecut - REmpoletaper
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (demo1(3,n))
      allocate (demo2(3,n))
      allocate (depo1(3,n))
      allocate (depo2(3,n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     initialize local variables for OpenMP calculation
c
      emo = 0.0d0
      epo = 0.0d0
      eintrao = eintra
      do i = 1, n
         do j = 1, 3
            demo1(j,i) = 0.0d0
            demo2(j,i) = 0.0d0
            depo1(j,i) = 0.0d0
            depo2(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            viro(j,i) = 0.0d0
         end do
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) firstprivate(f)
!$OMP& private(i,j,k,ii,kk,kkk,e,ei,bfac,damp,expdamp,
!$OMP& pdi,pti,pgamma,scale3,scale5,scale7,temp3,temp5,temp7,
!$OMP& dsc3,dsc5,dsc7,psc3,psc5,psc7,usc3,usc5,alsq2,alsq2n,
!$OMP& exp2a,ralpha,gfd,gfdr,xr,yr,zr,xix,yix,zix,
!$OMP& xiy,yiy,ziy,xiz,yiz,ziz,xkx,ykx,zkx,xky,yky,zky,
!$OMP& xkz,ykz,zkz,r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
!$OMP& erl,erli,iax,iay,iaz,kax,kay,kaz,vxx,vyy,vzz,vyx,vzx,vzy,
!$OMP& frcxi,frcyi,frczi,frcxk,frcyk,frczk,ci,di,qi,ck,dk,qk,
!$OMP& fridmp,findmp,ftm2,ftm2i,ftm2r,ftm2ri,ttm2,ttm3,
!$OMP& ttm2i,ttm3i,ttm2r,ttm3r,ttm2ri,ttm3ri,fdir,dixdk,
!$OMP& dkxui,dixuk,dixukp,dkxuip,uixqkr,ukxqir,uixqkrp,ukxqirp,
!$OMP& qiuk,qkui,qiukp,qkuip,rxqiuk,rxqkui,rxqiukp,rxqkuip,
!$OMP& qidk,qkdi,qir,qkr,qiqkr,qkqir,qixqk,rxqir,dixr,dkxr,
!$OMP& dixqkr,dkxqir,rxqkr,qkrxqir,rxqikr,rxqkir,rxqidk,rxqkdi,
!$OMP& ddsc3,ddsc5,ddsc7,bn,sc,gl,sci,scip,gli,glip,gf,gfi,
!$OMP& gfr,gfri,gti,gtri,dorl,dorli)
!$OMP& firstprivate(mscale,pscale,dscale,uscale)
!$OMP DO reduction(+:emo,epo,eintrao,demo1,demo2,depo1,depo2,viro)
!$OMP& schedule(guided)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, REnvlst(i)
            k = REvlst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. REmpolecut2) then
               r = sqrt(r2)
               if ((r2.gt.REmpoletaper2).and.(use_mpole_switch)) then
                  res_u = (r - REmpolecut + res_lam) / res_lam
                  ru2 = res_u * res_u
                  ru3 = ru2 * res_lam
                  ru4 = ru3 * res_lam
                  ru5 = ru4 * res_lam
                  q_switch = 1.0d0 + (15.0d0*ru4) - (6.0d0*ru5) 
     &                     - (10.0d0*ru3)
               else
                  q_switch = 1.0d0 
               end if
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     calculate the real space error function terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do j = 1, 5
                  bfac = dble(2*j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               dsc3 = 1.0d0 - scale3*dscale(kk)
               dsc5 = 1.0d0 - scale5*dscale(kk)
               dsc7 = 1.0d0 - scale7*dscale(kk)
               psc3 = 1.0d0 - scale3*pscale(kk)
               psc5 = 1.0d0 - scale5*pscale(kk)
               psc7 = 1.0d0 - scale7*pscale(kk)
               usc3 = 1.0d0 - scale3*uscale(kk)
               usc5 = 1.0d0 - scale5*uscale(kk)
c
c     construct necessary auxiliary vectors
c
               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,k) - di(3)*uind(2,k)
               dixuk(2) = di(3)*uind(1,k) - di(1)*uind(3,k)
               dixuk(3) = di(1)*uind(2,k) - di(2)*uind(1,k)
               dkxui(1) = dk(2)*uind(3,i) - dk(3)*uind(2,i)
               dkxui(2) = dk(3)*uind(1,i) - dk(1)*uind(3,i)
               dkxui(3) = dk(1)*uind(2,i) - dk(2)*uind(1,i)
               dixukp(1) = di(2)*uinp(3,k) - di(3)*uinp(2,k)
               dixukp(2) = di(3)*uinp(1,k) - di(1)*uinp(3,k)
               dixukp(3) = di(1)*uinp(2,k) - di(2)*uinp(1,k)
               dkxuip(1) = dk(2)*uinp(3,i) - dk(3)*uinp(2,i)
               dkxuip(2) = dk(3)*uinp(1,i) - dk(1)*uinp(3,i)
               dkxuip(3) = dk(1)*uinp(2,i) - dk(2)*uinp(1,i)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,k) + qi(4)*uind(2,k)
     &                      + qi(7)*uind(3,k)
               qiuk(2) = qi(2)*uind(1,k) + qi(5)*uind(2,k)
     &                      + qi(8)*uind(3,k)
               qiuk(3) = qi(3)*uind(1,k) + qi(6)*uind(2,k)
     &                      + qi(9)*uind(3,k)
               qkui(1) = qk(1)*uind(1,i) + qk(4)*uind(2,i)
     &                      + qk(7)*uind(3,i)
               qkui(2) = qk(2)*uind(1,i) + qk(5)*uind(2,i)
     &                      + qk(8)*uind(3,i)
               qkui(3) = qk(3)*uind(1,i) + qk(6)*uind(2,i)
     &                      + qk(9)*uind(3,i)
               qiukp(1) = qi(1)*uinp(1,k) + qi(4)*uinp(2,k)
     &                       + qi(7)*uinp(3,k)
               qiukp(2) = qi(2)*uinp(1,k) + qi(5)*uinp(2,k)
     &                       + qi(8)*uinp(3,k)
               qiukp(3) = qi(3)*uinp(1,k) + qi(6)*uinp(2,k)
     &                       + qi(9)*uinp(3,k)
               qkuip(1) = qk(1)*uinp(1,i) + qk(4)*uinp(2,i)
     &                       + qk(7)*uinp(3,i)
               qkuip(2) = qk(2)*uinp(1,i) + qk(5)*uinp(2,i)
     &                       + qk(8)*uinp(3,i)
               qkuip(3) = qk(3)*uinp(1,i) + qk(6)*uinp(2,i)
     &                       + qk(9)*uinp(3,i)
               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,i)*qkr(3) - uind(3,i)*qkr(2)
               uixqkr(2) = uind(3,i)*qkr(1) - uind(1,i)*qkr(3)
               uixqkr(3) = uind(1,i)*qkr(2) - uind(2,i)*qkr(1)
               ukxqir(1) = uind(2,k)*qir(3) - uind(3,k)*qir(2)
               ukxqir(2) = uind(3,k)*qir(1) - uind(1,k)*qir(3)
               ukxqir(3) = uind(1,k)*qir(2) - uind(2,k)*qir(1)
               uixqkrp(1) = uinp(2,i)*qkr(3) - uinp(3,i)*qkr(2)
               uixqkrp(2) = uinp(3,i)*qkr(1) - uinp(1,i)*qkr(3)
               uixqkrp(3) = uinp(1,i)*qkr(2) - uinp(2,i)*qkr(1)
               ukxqirp(1) = uinp(2,k)*qir(3) - uinp(3,k)*qir(2)
               ukxqirp(2) = uinp(3,k)*qir(1) - uinp(1,k)*qir(3)
               ukxqirp(3) = uinp(1,k)*qir(2) - uinp(2,k)*qir(1)
               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate the scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate the scalar products for induced components
c
               sci(1) = uind(1,i)*dk(1) + uind(2,i)*dk(2)
     &                     + uind(3,i)*dk(3) + di(1)*uind(1,k)
     &                     + di(2)*uind(2,k) + di(3)*uind(3,k)
               sci(2) = uind(1,i)*uind(1,k) + uind(2,i)*uind(2,k)
     &                     + uind(3,i)*uind(3,k)
               sci(3) = uind(1,i)*xr + uind(2,i)*yr + uind(3,i)*zr
               sci(4) = uind(1,k)*xr + uind(2,k)*yr + uind(3,k)*zr
               sci(7) = qir(1)*uind(1,k) + qir(2)*uind(2,k)
     &                     + qir(3)*uind(3,k)
               sci(8) = qkr(1)*uind(1,i) + qkr(2)*uind(2,i)
     &                     + qkr(3)*uind(3,i)
               scip(1) = uinp(1,i)*dk(1) + uinp(2,i)*dk(2)
     &                      + uinp(3,i)*dk(3) + di(1)*uinp(1,k)
     &                      + di(2)*uinp(2,k) + di(3)*uinp(3,k)
               scip(2) = uind(1,i)*uinp(1,k)+uind(2,i)*uinp(2,k)
     &                      + uind(3,i)*uinp(3,k)+uinp(1,i)*uind(1,k)
     &                      + uinp(2,i)*uind(2,k)+uinp(3,i)*uind(3,k)
               scip(3) = uinp(1,i)*xr + uinp(2,i)*yr + uinp(3,i)*zr
               scip(4) = uinp(1,k)*xr + uinp(2,k)*yr + uinp(3,k)*zr
               scip(7) = qir(1)*uinp(1,k) + qir(2)*uinp(2,k)
     &                      + qir(3)*uinp(3,k)
               scip(8) = qkr(1)*uinp(1,i) + qkr(2)*uinp(2,i)
     &                      + qkr(3)*uinp(3,i)
c
c     calculate the gl functions for permanent components
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
               gl(4) = sc(5)*sc(6)
               gl(5) = -4.0d0 * sc(9)
               gl(6) = sc(2)
               gl(7) = 2.0d0 * (sc(7)-sc(8))
               gl(8) = 2.0d0 * sc(10)
c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               e = bn(0)*gl(0) + bn(1)*(gl(1)+gl(6))
     &                + bn(2)*(gl(2)+gl(7)+gl(8))
     &                + bn(3)*(gl(3)+gl(5)) + bn(4)*gl(4)
               ei = 0.5d0 * (bn(1)*(gli(1)+gli(6))
     &                      + bn(2)*(gli(2)+gli(7)) + bn(3)*gli(3))
c
c     get the real energy without any screening function
c
               erl = rr1*gl(0) + rr3*(gl(1)+gl(6))
     &                  + rr5*(gl(2)+gl(7)+gl(8))
     &                  + rr7*(gl(3)+gl(5)) + rr9*gl(4)
               erli = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               e = e - erl * (1.0d0-mscale(kk))
               ei = ei - erli
               e = f * e
               ei = f * ei
               emo = emo + e
               epo = epo + ei
c
c     increment the total intramolecular energy; assumes
c     intramolecular distances are less than half of cell
c     length and less than the ewald cutoff
c
               if (molcule(ii) .eq. molcule(kk)) then
                  eintrao = eintrao + mscale(kk)*erl*f
                  eintrao = eintrao + 0.5d0*pscale(kk)
     &                         * (rr3*(gli(1)+gli(6))*scale3
     &                              + rr5*(gli(2)+gli(7))*scale5
     &                              + rr7*gli(3)*scale7)
               end if
c
c     set flags to compute components without screening
c
               dorl = .false.
               dorli = .false.
               if (mscale(kk) .ne. 1.0d0)  dorl = .true.
               if (psc3 .ne. 0.0d0)  dorli = .true.
               if (dsc3 .ne. 0.0d0)  dorli = .true.
               if (usc3 .ne. 0.0d0)  dorli = .true.
c
c     zero out force and torque components without screening
c
               do j = 1, 3
                  ftm2r(j) = 0.0d0
                  ftm2ri(j) = 0.0d0
                  ttm2r(j) = 0.0d0
                  ttm2ri(j) = 0.0d0
                  ttm3r(j) = 0.0d0
                  ttm3ri(j) = 0.0d0
               end do
c
c     get the permanent force with screening
c
               gf(1) = bn(1)*gl(0) + bn(2)*(gl(1)+gl(6))
     &                    + bn(3)*(gl(2)+gl(7)+gl(8))
     &                    + bn(4)*(gl(3)+gl(5)) + bn(5)*gl(4)
               gf(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gf(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gf(4) = 2.0d0 * bn(2)
               gf(5) = 2.0d0 * (-ck*bn(2)+sc(4)*bn(3)-sc(6)*bn(4))
               gf(6) = 2.0d0 * (-ci*bn(2)-sc(3)*bn(3)-sc(5)*bn(4))
               gf(7) = 4.0d0 * bn(3)
               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))
c
c     get the permanent force without screening
c
               if (dorl) then
                  gfr(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
     &                        + rr7*(gl(2)+gl(7)+gl(8))
     &                        + rr9*(gl(3)+gl(5)) + rr11*gl(4)
                  gfr(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
                  gfr(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
                  gfr(4) = 2.0d0 * rr5
                  gfr(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
                  gfr(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
                  gfr(7) = 4.0d0 * rr7
                  ftm2r(1) = gfr(1)*xr + gfr(2)*di(1) + gfr(3)*dk(1)
     &                          + gfr(4)*(qkdi(1)-qidk(1))
     &                          + gfr(5)*qir(1) + gfr(6)*qkr(1)
     &                          + gfr(7)*(qiqkr(1)+qkqir(1))
                  ftm2r(2) = gfr(1)*yr + gfr(2)*di(2) + gfr(3)*dk(2)
     &                          + gfr(4)*(qkdi(2)-qidk(2))
     &                          + gfr(5)*qir(2) + gfr(6)*qkr(2)
     &                          + gfr(7)*(qiqkr(2)+qkqir(2))
                  ftm2r(3) = gfr(1)*zr + gfr(2)*di(3) + gfr(3)*dk(3)
     &                          + gfr(4)*(qkdi(3)-qidk(3))
     &                          + gfr(5)*qir(3) + gfr(6)*qkr(3)
     &                          + gfr(7)*(qiqkr(3)+qkqir(3))
               end if
c
c     get the induced force with screening
c
               gfi(1) = 0.5d0*bn(2)*(gli(1)+glip(1)+gli(6)+glip(6))
     &                     + 0.5d0*bn(2)*scip(2)
     &                     + 0.5d0*bn(3)*(gli(2)+glip(2)+gli(7)+glip(7))
     &                     - 0.5d0*bn(3)*(sci(3)*scip(4)+scip(3)*sci(4))
     &                     + 0.5d0*bn(4)*(gli(3)+glip(3))
               gfi(2) = -ck*bn(1) + sc(4)*bn(2) - sc(6)*bn(3)
               gfi(3) = ci*bn(1) + sc(3)*bn(2) + sc(5)*bn(3)
               gfi(4) = 2.0d0 * bn(2)
               gfi(5) = bn(3) * (sci(4)+scip(4))
               gfi(6) = -bn(3) * (sci(3)+scip(3))
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &             (gfi(2)*(uind(1,i)+uinp(1,i))
     &            + bn(2)*(sci(4)*uinp(1,i)+scip(4)*uind(1,i))
     &            + gfi(3)*(uind(1,k)+uinp(1,k))
     &            + bn(2)*(sci(3)*uinp(1,k)+scip(3)*uind(1,k))
     &            + (sci(4)+scip(4))*bn(2)*di(1)
     &            + (sci(3)+scip(3))*bn(2)*dk(1)
     &            + gfi(4)*(qkui(1)+qkuip(1)-qiuk(1)-qiukp(1)))
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &             (gfi(2)*(uind(2,i)+uinp(2,i))
     &            + bn(2)*(sci(4)*uinp(2,i)+scip(4)*uind(2,i))
     &            + gfi(3)*(uind(2,k)+uinp(2,k))
     &            + bn(2)*(sci(3)*uinp(2,k)+scip(3)*uind(2,k))
     &            + (sci(4)+scip(4))*bn(2)*di(2)
     &            + (sci(3)+scip(3))*bn(2)*dk(2)
     &            + gfi(4)*(qkui(2)+qkuip(2)-qiuk(2)-qiukp(2)))
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr + 0.5d0*
     &             (gfi(2)*(uind(3,i)+uinp(3,i))
     &            + bn(2)*(sci(4)*uinp(3,i)+scip(4)*uind(3,i))
     &            + gfi(3)*(uind(3,k)+uinp(3,k))
     &            + bn(2)*(sci(3)*uinp(3,k)+scip(3)*uind(3,k))
     &            + (sci(4)+scip(4))*bn(2)*di(3)
     &            + (sci(3)+scip(3))*bn(2)*dk(3)
     &            + gfi(4)*(qkui(3)+qkuip(3)-qiuk(3)-qiukp(3)))
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     get the induced force without screening
c
               if (dorli) then
                  gfri(1) = 0.5d0*rr5*((gli(1)+gli(6))*psc3
     &                               + (glip(1)+glip(6))*dsc3
     &                               + scip(2)*usc3)
     &                    + 0.5d0*rr7*((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                        - (sci(3)*scip(4)+scip(3)*sci(4))*usc5)
     &                    + 0.5d0*rr9*(gli(3)*psc7+glip(3)*dsc7)
                  gfri(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
                  gfri(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
                  gfri(4) = 2.0d0 * rr5
                  gfri(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
                  gfri(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
                  ftm2ri(1) = gfri(1)*xr + 0.5d0*
     &              (- rr3*ck*(uind(1,i)*psc3+uinp(1,i)*dsc3)
     &               + rr5*sc(4)*(uind(1,i)*psc5+uinp(1,i)*dsc5)
     &               - rr7*sc(6)*(uind(1,i)*psc7+uinp(1,i)*dsc7))
     &               + (rr3*ci*(uind(1,k)*psc3+uinp(1,k)*dsc3)
     &               + rr5*sc(3)*(uind(1,k)*psc5+uinp(1,k)*dsc5)
     &               + rr7*sc(5)*(uind(1,k)*psc7+uinp(1,k)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &               + sci(3)*uinp(1,k)+scip(3)*uind(1,k))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &               + 0.5d0*gfri(4)*((qkui(1)-qiuk(1))*psc5
     &               + (qkuip(1)-qiukp(1))*dsc5)
     &               + gfri(5)*qir(1) + gfri(6)*qkr(1)
                  ftm2ri(2) = gfri(1)*yr + 0.5d0*
     &              (- rr3*ck*(uind(2,i)*psc3+uinp(2,i)*dsc3)
     &               + rr5*sc(4)*(uind(2,i)*psc5+uinp(2,i)*dsc5)
     &               - rr7*sc(6)*(uind(2,i)*psc7+uinp(2,i)*dsc7))
     &               + (rr3*ci*(uind(2,k)*psc3+uinp(2,k)*dsc3)
     &               + rr5*sc(3)*(uind(2,k)*psc5+uinp(2,k)*dsc5)
     &               + rr7*sc(5)*(uind(2,k)*psc7+uinp(2,k)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &               + sci(3)*uinp(2,k)+scip(3)*uind(2,k))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &               + 0.5d0*gfri(4)*((qkui(2)-qiuk(2))*psc5
     &               + (qkuip(2)-qiukp(2))*dsc5)
     &               + gfri(5)*qir(2) + gfri(6)*qkr(2)
                  ftm2ri(3) = gfri(1)*zr + 0.5d0*
     &              (- rr3*ck*(uind(3,i)*psc3+uinp(3,i)*dsc3)
     &               + rr5*sc(4)*(uind(3,i)*psc5+uinp(3,i)*dsc5)
     &               - rr7*sc(6)*(uind(3,i)*psc7+uinp(3,i)*dsc7))
     &               + (rr3*ci*(uind(3,k)*psc3+uinp(3,k)*dsc3)
     &               + rr5*sc(3)*(uind(3,k)*psc5+uinp(3,k)*dsc5)
     &               + rr7*sc(5)*(uind(3,k)*psc7+uinp(3,k)*dsc7))*0.5d0
     &               + rr5*usc5*(sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &               + sci(3)*uinp(3,k)+scip(3)*uind(3,k))*0.5d0
     &               + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &               + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &               + 0.5d0*gfri(4)*((qkui(3)-qiuk(3))*psc5
     &               + (qkuip(3)-qiukp(3))*dsc5)
     &               + gfri(5)*qir(3) + gfri(6)*qkr(3)
               end if
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify the forces for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c
               if (poltyp .eq. 'DIRECT') then
                  gfd = 0.5d0 * (bn(2)*scip(2)
     &                     - bn(3)*(scip(3)*sci(4)+sci(3)*scip(4)))
                  gfdr = 0.5d0 * (rr5*scip(2)*usc3
     &                     - rr7*(scip(3)*sci(4)
     &                           +sci(3)*scip(4))*usc5)
                  ftm2i(1) = ftm2i(1) - gfd*xr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                          +sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  ftm2i(2) = ftm2i(2) - gfd*yr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                          +sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  ftm2i(3) = ftm2i(3) - gfd*zr - 0.5d0*bn(2)*
     &                          (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                          +sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  fdir(1) = gfdr*xr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(1,i)+scip(4)*uind(1,i)
     &                        + sci(3)*uinp(1,k)+scip(3)*uind(1,k))
                  fdir(2) = gfdr*yr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(2,i)+scip(4)*uind(2,i)
     &                        + sci(3)*uinp(2,k)+scip(3)*uind(2,k))
                  fdir(3) = gfdr*zr + 0.5d0*usc5*rr5*
     &                         (sci(4)*uinp(3,i)+scip(4)*uind(3,i)
     &                        + sci(3)*uinp(3,k)+scip(3)*uind(3,k))
                  ftm2i(1) = ftm2i(1) + fdir(1) + findmp(1)
                  ftm2i(2) = ftm2i(2) + fdir(2) + findmp(2)
                  ftm2i(3) = ftm2i(3) + fdir(3) + findmp(3)
               end if
c
c     get the permanent torque with screening
c
               ttm2(1) = -bn(1)*dixdk(1) + gf(2)*dixr(1)
     &                      + gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqidk(1)-2.0d0*qixqk(1))
     &                      - gf(5)*rxqir(1)
     &                      - gf(7)*(rxqikr(1)+qkrxqir(1))
               ttm2(2) = -bn(1)*dixdk(2) + gf(2)*dixr(2)
     &                      + gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqidk(2)-2.0d0*qixqk(2))
     &                      - gf(5)*rxqir(2)
     &                      - gf(7)*(rxqikr(2)+qkrxqir(2))
               ttm2(3) = -bn(1)*dixdk(3) + gf(2)*dixr(3)
     &                      + gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqidk(3)-2.0d0*qixqk(3))
     &                      - gf(5)*rxqir(3)
     &                      - gf(7)*(rxqikr(3)+qkrxqir(3))
               ttm3(1) = bn(1)*dixdk(1) + gf(3)*dkxr(1)
     &                      - gf(4)*(dixqkr(1)+dkxqir(1)
     &                              +rxqkdi(1)-2.0d0*qixqk(1))
     &                      - gf(6)*rxqkr(1)
     &                      - gf(7)*(rxqkir(1)-qkrxqir(1))
               ttm3(2) = bn(1)*dixdk(2) + gf(3)*dkxr(2)
     &                      - gf(4)*(dixqkr(2)+dkxqir(2)
     &                              +rxqkdi(2)-2.0d0*qixqk(2))
     &                      - gf(6)*rxqkr(2)
     &                      - gf(7)*(rxqkir(2)-qkrxqir(2))
               ttm3(3) = bn(1)*dixdk(3) + gf(3)*dkxr(3)
     &                      - gf(4)*(dixqkr(3)+dkxqir(3)
     &                              +rxqkdi(3)-2.0d0*qixqk(3))
     &                      - gf(6)*rxqkr(3)
     &                      - gf(7)*(rxqkir(3)-qkrxqir(3))
c
c     get the permanent torque without screening
c
               if (dorl) then
                  ttm2r(1) = -rr3*dixdk(1) + gfr(2)*dixr(1)
     &                          + gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqidk(1)-2.0d0*qixqk(1))
     &                          - gfr(5)*rxqir(1)
     &                          - gfr(7)*(rxqikr(1)+qkrxqir(1))
                  ttm2r(2) = -rr3*dixdk(2) + gfr(2)*dixr(2)
     &                          + gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqidk(2)-2.0d0*qixqk(2))
     &                          - gfr(5)*rxqir(2)
     &                          - gfr(7)*(rxqikr(2)+qkrxqir(2))
                  ttm2r(3) = -rr3*dixdk(3) + gfr(2)*dixr(3)
     &                          + gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqidk(3)-2.0d0*qixqk(3))
     &                          - gfr(5)*rxqir(3)
     &                          - gfr(7)*(rxqikr(3)+qkrxqir(3))
                  ttm3r(1) = rr3*dixdk(1) + gfr(3)*dkxr(1)
     &                          - gfr(4)*(dixqkr(1)+dkxqir(1)
     &                                   +rxqkdi(1)-2.0d0*qixqk(1))
     &                          - gfr(6)*rxqkr(1)
     &                          - gfr(7)*(rxqkir(1)-qkrxqir(1))
                  ttm3r(2) = rr3*dixdk(2) + gfr(3)*dkxr(2)
     &                          - gfr(4)*(dixqkr(2)+dkxqir(2)
     &                                   +rxqkdi(2)-2.0d0*qixqk(2))
     &                          - gfr(6)*rxqkr(2)
     &                          - gfr(7)*(rxqkir(2)-qkrxqir(2))
                  ttm3r(3) = rr3*dixdk(3) + gfr(3)*dkxr(3)
     &                          - gfr(4)*(dixqkr(3)+dkxqir(3)
     &                                   +rxqkdi(3)-2.0d0*qixqk(3))
     &                          - gfr(6)*rxqkr(3)
     &                          - gfr(7)*(rxqkir(3)-qkrxqir(3))
               end if
c
c     get the induced torque with screening
c
               gti(2) = 0.5d0 * bn(2) * (sci(4)+scip(4))
               gti(3) = 0.5d0 * bn(2) * (sci(3)+scip(3))
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
               ttm2i(1) = -0.5d0*bn(1)*(dixuk(1)+dixukp(1))
     &                       + gti(2)*dixr(1) - gti(5)*rxqir(1)
     &                       + 0.5d0*gti(4)*(ukxqir(1)+rxqiuk(1)
     &                                      +ukxqirp(1)+rxqiukp(1))
               ttm2i(2) = -0.5d0*bn(1)*(dixuk(2)+dixukp(2))
     &                       + gti(2)*dixr(2) - gti(5)*rxqir(2)
     &                       + 0.5d0*gti(4)*(ukxqir(2)+rxqiuk(2)
     &                                      +ukxqirp(2)+rxqiukp(2))
               ttm2i(3) = -0.5d0*bn(1)*(dixuk(3)+dixukp(3))
     &                       + gti(2)*dixr(3) - gti(5)*rxqir(3)
     &                       + 0.5d0*gti(4)*(ukxqir(3)+rxqiuk(3)
     &                                      +ukxqirp(3)+rxqiukp(3))
               ttm3i(1) = -0.5d0*bn(1)*(dkxui(1)+dkxuip(1))
     &                       + gti(3)*dkxr(1) - gti(6)*rxqkr(1)
     &                       - 0.5d0*gti(4)*(uixqkr(1)+rxqkui(1)
     &                                      +uixqkrp(1)+rxqkuip(1))
               ttm3i(2) = -0.5d0*bn(1)*(dkxui(2)+dkxuip(2))
     &                       + gti(3)*dkxr(2) - gti(6)*rxqkr(2)
     &                       - 0.5d0*gti(4)*(uixqkr(2)+rxqkui(2)
     &                                       +uixqkrp(2)+rxqkuip(2))
               ttm3i(3) = -0.5d0*bn(1)*(dkxui(3)+dkxuip(3))
     &                       + gti(3)*dkxr(3) - gti(6)*rxqkr(3)
     &                       - 0.5d0*gti(4)*(uixqkr(3)+rxqkui(3)
     &                                      +uixqkrp(3)+rxqkuip(3))
c
c     get the induced torque without screening
c
               if (dorli) then
                  gtri(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
                  gtri(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
                  gtri(4) = gfri(4)
                  gtri(5) = gfri(5)
                  gtri(6) = gfri(6)
                  ttm2ri(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(1) - gtri(5)*rxqir(1)
     &                           + gtri(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &                             +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0
                  ttm2ri(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(2) - gtri(5)*rxqir(2)
     &                           + gtri(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &                             +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0
                  ttm2ri(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &                           + gtri(2)*dixr(3) - gtri(5)*rxqir(3)
     &                           + gtri(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &                             +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0
                  ttm3ri(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(1) - gtri(6)*rxqkr(1)
     &                           - gtri(4)*((uixqkr(1)+rxqkui(1))*psc5
     &                             +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0
                  ttm3ri(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(2) - gtri(6)*rxqkr(2)
     &                           - gtri(4)*((uixqkr(2)+rxqkui(2))*psc5
     &                             +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0
                  ttm3ri(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &                           + gtri(3)*dkxr(3) - gtri(6)*rxqkr(3)
     &                           - gtri(4)*((uixqkr(3)+rxqkui(3))*psc5
     &                             +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0
               end if
c
c     handle the case where scaling is used
c
               do j = 1, 3
                  ftm2(j) = f * (ftm2(j)-(1.0d0-mscale(kk))*ftm2r(j))
                  ftm2i(j) = f * (ftm2i(j)-ftm2ri(j))
                  ttm2(j) = f * (ttm2(j)-(1.0d0-mscale(kk))*ttm2r(j))
                  ttm2i(j) = f * (ttm2i(j)-ttm2ri(j))
                  ttm3(j) = f * (ttm3(j)-(1.0d0-mscale(kk))*ttm3r(j))
                  ttm3i(j) = f * (ttm3i(j)-ttm3ri(j))
               end do
c
c     increment gradient due to force and torque on first site
c
               demo1(1,ii) = demo1(1,ii) + ftm2(1)*q_switch
               demo1(2,ii) = demo1(2,ii) + ftm2(2)*q_switch
               demo1(3,ii) = demo1(3,ii) + ftm2(3)*q_switch
               depo1(1,ii) = depo1(1,ii) + ftm2i(1)*q_switch
               depo1(2,ii) = depo1(2,ii) + ftm2i(2)*q_switch
               depo1(3,ii) = depo1(3,ii) + ftm2i(3)*q_switch
               call torque3 (i,ttm2,ttm2i,frcxi,frcyi,frczi,demo1,depo1)
c
c     increment gradient due to force and torque on second site
c
               demo2(1,kk) = demo2(1,kk) - ftm2(1)*q_switch
               demo2(2,kk) = demo2(2,kk) - ftm2(2)*q_switch
               demo2(3,kk) = demo2(3,kk) - ftm2(3)*q_switch
               depo2(1,kk) = depo2(1,kk) - ftm2i(1)*q_switch
               depo2(2,kk) = depo2(2,kk) - ftm2i(2)*q_switch
               depo2(3,kk) = depo2(3,kk) - ftm2i(3)*q_switch
               call torque3 (k,ttm3,ttm3i,frcxk,frcyk,frczk,demo2,depo2)
c
c     increment the internal virial tensor components
c
               iaz = zaxis(i)
               iax = xaxis(i)
               iay = yaxis(i)
               kaz = zaxis(k)
               kax = xaxis(k)
               kay = yaxis(k)
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)
               viro(1,1) = viro(1,1) + vxx
               viro(2,1) = viro(2,1) + vyx
               viro(3,1) = viro(3,1) + vzx
               viro(1,2) = viro(1,2) + vyx
               viro(2,2) = viro(2,2) + vyy
               viro(3,2) = viro(3,2) + vzy
               viro(1,3) = viro(1,3) + vzx
               viro(2,3) = viro(2,3) + vzy
               viro(3,3) = viro(3,3) + vzz
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     end OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      em = em + emo
      ep = ep + epo
      eintra = eintrao
      do i = 1, n
         do j = 1, 3
            dem(j,i) = dem(j,i) + demo1(j,i) + demo2(j,i)
            dep(j,i) = dep(j,i) + depo1(j,i) + depo2(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + viro(j,i)
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (demo1)
      deallocate (demo2)
      deallocate (depo1)
      deallocate (depo2)
      return
      end

c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine dfield0cSR  --  direct induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "dfield0cSR" computes the mutual electrostatic field due to
c     permanent multipole moments via Ewald summation
c
c
      subroutine dfield0cSR (field,fieldp)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
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
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do

c
c     get the real space portion of the electrostatic field
c
c      if (use_mlist) then
         call udirect2bSR (field,fieldp)
c      else
c         call udirect2a (field,fieldp)
c      end if
c
c     get the self-energy portion of the electrostatic field
c
c      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c      do i = 1, npole
c         do j = 1, 3
c            field(j,i) = field(j,i) + term*rpole(j+1,i)
c            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
c         end do
c      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end     
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ufield0cSR  --  mutual induction via Ewald sum  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ufield0cSR" computes the mutual electrostatic field due to
c     induced dipole moments via Ewald summation
c
c
      subroutine ufield0cSR (field,fieldp)
      use sizes
      use atoms
      use boxes
      use ewald
      use limits
      use math
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 term
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c
c
c     zero out the electrostatic field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do

c
c     get the real space portion of the electrostatic field
c
c      if (use_mlist) then
         call umutual2bSR (field,fieldp)
c      else
c         call umutual2a (field,fieldp)
c      end if
c
c     get the self-energy portion of the electrostatic field
c
c      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c      do i = 1, npole
c         do j = 1, 3
c            field(j,i) = field(j,i) + term*uind(j,i)
c            fieldp(j,i) = fieldp(j,i) + term*uinp(j,i)
c         end do
c      end do
c
c     compute the cell dipole boundary correction to the field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
            ucellp(i) = 0.0d0
         end do
         do i = 1, npole
            do j = 1, 3
               ucell(j) = ucell(j) + uind(j,i)
               ucellp(j) = ucellp(j) + uinp(j,i)
            end do
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucellp(j)
            end do
         end do
      end if
      return
      end
c
c
c     
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine uscale0bSR  --  dipole preconditioner via list  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "uscale0bSR" builds and applies a preconditioner for the conjugate
c     gradient induced dipole solver using a neighbor pair list
c
c
      subroutine uscale0bSR (mode,rsd,rsdp,zrsd,zrsdp)
      use sizes
      use atoms
      use mpole
      use neigh
      use polar
      use polgrp
      use polpot
      use usolve
      implicit none
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 pdi,pti
      real*8 polmin
      real*8 poli,polik
      real*8 damp,expdamp
      real*8 pgamma
      real*8 scale3,scale5
      real*8 m1,m2,m3
      real*8 m4,m5,m6
      real*8, allocatable :: dscale(:)
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
         do i = 1, npole
            poli = udiag * max(polmin,polarity(i))
            do j = 1, 3
               zrsd(j,i) = 0.0d0
               zrsdp(j,i) = 0.0d0
               zrsdt(j,i) = poli * rsd(j,i)
               zrsdtp(j,i) = poli * rsdp(j,i)
            end do
         end do
c
c     use the off-diagonal preconditioner elements in second phase
c
!$OMP PARALLEL default(private) shared(npole,mindex,minv,nulst,ulst,
!$OMP& rsd,rsdp,zrsd,zrsdp,zrsdt,zrsdtp)
!$OMP DO reduction(+:zrsdt,zrsdtp) schedule(guided)
         do i = 1, npole
            m = mindex(i)
            do kk = 1, nulst(i)
               k = ulst(kk,i)
               m1 = minv(m+1)
               m2 = minv(m+2)
               m3 = minv(m+3)
               m4 = minv(m+4)
               m5 = minv(m+5)
               m6 = minv(m+6)
               m = m + 6
               zrsdt(1,i) = zrsdt(1,i) + m1*rsd(1,k) + m2*rsd(2,k)
     &                        + m3*rsd(3,k)
               zrsdt(2,i) = zrsdt(2,i) + m2*rsd(1,k) + m4*rsd(2,k)
     &                        + m5*rsd(3,k)
               zrsdt(3,i) = zrsdt(3,i) + m3*rsd(1,k) + m5*rsd(2,k)
     &                        + m6*rsd(3,k)
               zrsdt(1,k) = zrsdt(1,k) + m1*rsd(1,i) + m2*rsd(2,i)
     &                        + m3*rsd(3,i)
               zrsdt(2,k) = zrsdt(2,k) + m2*rsd(1,i) + m4*rsd(2,i)
     &                        + m5*rsd(3,i)
               zrsdt(3,k) = zrsdt(3,k) + m3*rsd(1,i) + m5*rsd(2,i)
     &                        + m6*rsd(3,i)
               zrsdtp(1,i) = zrsdtp(1,i) + m1*rsdp(1,k) + m2*rsdp(2,k)
     &                         + m3*rsdp(3,k)
               zrsdtp(2,i) = zrsdtp(2,i) + m2*rsdp(1,k) + m4*rsdp(2,k)
     &                         + m5*rsdp(3,k)
               zrsdtp(3,i) = zrsdtp(3,i) + m3*rsdp(1,k) + m5*rsdp(2,k)
     &                         + m6*rsdp(3,k)
               zrsdtp(1,k) = zrsdtp(1,k) + m1*rsdp(1,i) + m2*rsdp(2,i)
     &                         + m3*rsdp(3,i)
               zrsdtp(2,k) = zrsdtp(2,k) + m2*rsdp(1,i) + m4*rsdp(2,i)
     &                         + m5*rsdp(3,i)
               zrsdtp(3,k) = zrsdtp(3,k) + m3*rsdp(1,i) + m5*rsdp(2,i)
     &                         + m6*rsdp(3,i)
            end do
         end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
         do i = 1, npole
            do j = 1, 3
               zrsd(j,i) = zrsdt(j,i) + zrsd(j,i)
               zrsdp(j,i) = zrsdtp(j,i) + zrsdp(j,i)
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
         do i = 1, npole
            mindex(i) = m
            m = m + 6*nulst(i)
         end do
c
c     perform dynamic allocation of some local arrays
c
         allocate (dscale(n))
c
c     set array needed to scale connected atom interactions
c
         do i = 1, n
            dscale(i) = 1.0d0
         end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,
!$OMP& thole,polarity,u1scale,u2scale,u3scale,u4scale,np11,ip11,
!$OMP& np12,ip12,np13,ip13,np14,ip14,nulst,ulst,mindex,minv)
!$OMP& firstprivate (dscale)
c
c     determine the off-diagonal elements of the preconditioner
c
!$OMP DO schedule(guided)
         do i = 1, npole
            ii = ipole(i)
            xi = x(ii)
            yi = y(ii)
            zi = z(ii)
            pdi = pdamp(i)
            pti = thole(i)
            poli = polarity(i)
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = u4scale
            end do
            m = mindex(i)
            do kkk = 1, nulst(i)
               k = ulst(kkk,i)
               kk = ipole(k)
               xr = x(kk) - xi
               yr = y(kk) - yi
               zr = z(kk) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               r = sqrt(r2)
               scale3 = dscale(kk)
               scale5 = dscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-expdamp*(1.0d0-damp))
                  end if
               end if
               polik = poli * polarity(k)
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
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, np11(ii)
               dscale(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale(ip14(j,ii)) = 1.0d0
            end do
         end do
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
         deallocate (dscale)
      end if
      return
      end

c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine udirect2bSR  --  Ewald real direct field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "udirect2bSR" computes the real space contribution of the permanent
c     atomic multipole moments to the field via a neighbor list
c
c
      subroutine udirect2bSR (field,fieldp)
      use sizes
      use atoms
      use boxes
      use bound
      use couple
      use ewald
      use limits
      use math
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
      integer nlocal,maxlocal
      integer tid,toffset0
!$    integer omp_get_thread_num
      integer, allocatable :: toffset(:)
      integer, allocatable :: ilocal(:,:)
      real*8 xr,yr,zr,r,r2
      real*8 rr1,rr2,rr3
      real*8 rr5,rr7
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
      real*8 ralpha,aefac
      real*8 aesq2,aesq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 REmpoletaper2,REmpolecut2
      real*8 res_lam,res_u,ru2,ru3,ru4,ru5,q_switch
      real*8 bn(0:3),bcn(3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8, allocatable :: fieldt(:,:)
      real*8, allocatable :: fieldtp(:,:)
      real*8, allocatable :: dlocal(:,:)
      character*6 mode
      external erfc
      REmpoletaper2 = REmpoletaper**2
      REmpolecut2 = REmpolecut**2


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
      nlocal = 0
      toffset0 = 0
      maxlocal = int((npole*maxelst)/nthread)
c
c     perform dynamic allocation of some local arrays
c
      allocate (toffset(0:nthread-1))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (fieldt(3,npole))
      allocate (fieldtp(3,npole))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,npole,ipole,x,y,z,pdamp,thole,
!$OMP& rpole,p2scale,p3scale,p4scale,p41scale,p5scale,d1scale,d2scale,
!$OMP& d3scale,d4scale,u1scale,u2scale,u3scale,u4scale,n12,i12,n13,i13,
!$OMP& n14,i14,n15,i15,np11,ip11,np12,ip12,np13,ip13,np14,ip14,nelst,
!$OMP& elst,cut2,aewald,aesq2,aesq2n,poltyp,ntpair,tindex,tdipdip,
!$OMP& toffset,toffset0,field,fieldp,fieldt,fieldtp,maxlocal)
!$OMP& firstprivate(pscale,dscale,uscale,nlocal)
c
c     perform dynamic allocation of some local arrays
c
      if (poltyp .eq. 'MUTUAL') then
         allocate (ilocal(2,maxlocal))
         allocate (dlocal(6,maxlocal))
      end if
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     compute the real space portion of the Ewald summation
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do kkk = 1, REnvlst(i)
            k = REvlst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
            if (r2 .le. REmpolecut2) then
               r = sqrt(r2)
               rr1 = 1.0d0 / r
               rr2 = rr1 * rr1
               rr3 = rr2 * rr1
               rr5 = rr2 * rr3
               rr7 = rr2 * rr5
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
c
c     calculate the error function damping factors
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
c     compute the polarization damping scale factors
c
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
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
c
c     find the field terms for the current interaction
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
               bcn(1) = bn(1) - (1.0d0-scale3*dscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*dscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*dscale(kk))*rr7
               fimd(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimd(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimd(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmd(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmd(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmd(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
               bcn(1) = bn(1) - (1.0d0-scale3*pscale(kk))*rr3
               bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*pscale(kk))*rr5
               bcn(3) = bn(3) - 15.0d0*(1.0d0-scale7*pscale(kk))*rr7
               fimp(1) = -xr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkx + 2.0d0*bcn(2)*qkx
               fimp(2) = -yr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dky + 2.0d0*bcn(2)*qky
               fimp(3) = -zr*(bcn(1)*ck-bcn(2)*dkr+bcn(3)*qkr)
     &                     - bcn(1)*dkz + 2.0d0*bcn(2)*qkz
               fkmp(1) = xr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*dix - 2.0d0*bcn(2)*qix
               fkmp(2) = yr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diy - 2.0d0*bcn(2)*qiy
               fkmp(3) = zr*(bcn(1)*ci+bcn(2)*dir+bcn(3)*qir)
     &                     - bcn(1)*diz - 2.0d0*bcn(2)*qiz
c
c     find terms needed later to compute mutual polarization
c
               if (poltyp .eq. 'MUTUAL') then
                  bcn(1) = bn(1) - (1.0d0-scale3*uscale(kk))*rr3
                  bcn(2) = bn(2) - 3.0d0*(1.0d0-scale5*uscale(kk))*rr5
                  nlocal = nlocal + 1
                  ilocal(1,nlocal) = i
                  ilocal(2,nlocal) = k
                  dlocal(1,nlocal) = -bcn(1) + bcn(2)*xr*xr
                  dlocal(2,nlocal) = bcn(2)*xr*yr
                  dlocal(3,nlocal) = bcn(2)*xr*zr
                  dlocal(4,nlocal) = -bcn(1) + bcn(2)*yr*yr
                  dlocal(5,nlocal) = bcn(2)*yr*zr
                  dlocal(6,nlocal) = -bcn(1) + bcn(2)*zr*zr
               end if
c
c     increment the field at each site due to this interaction
c
c
c        SWITCH FUNCTION HERE I THINK
c
c
c
               do j = 1, 3
                  fieldt(j,i) = fieldt(j,i) + fimd(j)
                  fieldt(j,k) = fieldt(j,k) + fkmd(j)
                  fieldtp(j,i) = fieldtp(j,i) + fimp(j)
                  fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
               end do
            end if
         end do

c
c     reset interaction scaling coefficients for connected atoms
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
         do j = 1, np11(ii)
            uscale(ip11(j,ii)) = 1.0d0
            dscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            uscale(ip12(j,ii)) = 1.0d0
            dscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            uscale(ip13(j,ii)) = 1.0d0
            dscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            uscale(ip14(j,ii)) = 1.0d0
            dscale(ip14(j,ii)) = 1.0d0
         end do
      end do
!$OMP END DO
c
c     transfer the results from local to global arrays
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
         end do
      end do
!$OMP END DO
c
c     store terms needed later to compute mutual polarization
c
!$OMP CRITICAL
      tid = 0
!$    tid = omp_get_thread_num ()
      toffset(tid) = toffset0
      toffset0 = toffset0 + nlocal
      ntpair = toffset0
!$OMP END CRITICAL
      if (poltyp .eq. 'MUTUAL') then
         k = toffset(tid)
         do i = 1, nlocal
            m = k + i
            tindex(1,m) = ilocal(1,i)
            tindex(2,m) = ilocal(2,i)
            do j = 1, 6
               tdipdip(j,m) = dlocal(j,i)
            end do
         end do
         deallocate (ilocal)
         deallocate (dlocal)
      end if
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (toffset)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (fieldt)
      deallocate (fieldtp)
      return
      end
c

c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine umutual2bSR  --  Ewald real mutual field via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "umutual2bSR" computes the real space contribution of the induced
c     atomic dipole moments to the field via a neighbor list
c
c
      subroutine umutual2bSR (field,fieldp)
      use sizes
      use mpole
      use polar
      use tarray
      implicit none
      integer i,j,k,m
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
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
c     set OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(npole,uind,uinp,ntpair,tindex,
!$OMP& tdipdip,field,fieldp,fieldt,fieldtp)
c
c     initialize local variables for OpenMP calculation
c
!$OMP DO collapse(2)
      do i = 1, npole
         do j = 1, 3
            fieldt(j,i) = 0.0d0
            fieldtp(j,i) = 0.0d0
         end do
      end do
!$OMP END DO
c
c     find the field terms for each pairwise interaction
c
!$OMP DO reduction(+:fieldt,fieldtp) schedule(guided)
      do m = 1, ntpair
         i = tindex(1,m)
         k = tindex(2,m)
         fimd(1) = tdipdip(1,m)*uind(1,k) + tdipdip(2,m)*uind(2,k)
     &                + tdipdip(3,m)*uind(3,k)
         fimd(2) = tdipdip(2,m)*uind(1,k) + tdipdip(4,m)*uind(2,k)
     &                + tdipdip(5,m)*uind(3,k)
         fimd(3) = tdipdip(3,m)*uind(1,k) + tdipdip(5,m)*uind(2,k)
     &                + tdipdip(6,m)*uind(3,k)
         fkmd(1) = tdipdip(1,m)*uind(1,i) + tdipdip(2,m)*uind(2,i)
     &                + tdipdip(3,m)*uind(3,i)
         fkmd(2) = tdipdip(2,m)*uind(1,i) + tdipdip(4,m)*uind(2,i)
     &                + tdipdip(5,m)*uind(3,i)
         fkmd(3) = tdipdip(3,m)*uind(1,i) + tdipdip(5,m)*uind(2,i)
     &                + tdipdip(6,m)*uind(3,i)
         fimp(1) = tdipdip(1,m)*uinp(1,k) + tdipdip(2,m)*uinp(2,k)
     &                + tdipdip(3,m)*uinp(3,k)
         fimp(2) = tdipdip(2,m)*uinp(1,k) + tdipdip(4,m)*uinp(2,k)
     &                + tdipdip(5,m)*uinp(3,k)
         fimp(3) = tdipdip(3,m)*uinp(1,k) + tdipdip(5,m)*uinp(2,k)
     &                + tdipdip(6,m)*uinp(3,k)
         fkmp(1) = tdipdip(1,m)*uinp(1,i) + tdipdip(2,m)*uinp(2,i)
     &                + tdipdip(3,m)*uinp(3,i)
         fkmp(2) = tdipdip(2,m)*uinp(1,i) + tdipdip(4,m)*uinp(2,i)
     &                + tdipdip(5,m)*uinp(3,i)
         fkmp(3) = tdipdip(3,m)*uinp(1,i) + tdipdip(5,m)*uinp(2,i)
     &                + tdipdip(6,m)*uinp(3,i)
c
c     increment the field at each site due to this interaction
c
         do j = 1, 3
            fieldt(j,i) = fieldt(j,i) + fimd(j)
            fieldt(j,k) = fieldt(j,k) + fkmd(j)
            fieldtp(j,i) = fieldtp(j,i) + fimp(j)
            fieldtp(j,k) = fieldtp(j,k) + fkmp(j)
         end do
      end do
!$OMP END DO
c
c     end OpenMP directives for the major loop structure
c
!$OMP DO
      do i = 1, npole
         do j = 1, 3
            field(j,i) = fieldt(j,i) + field(j,i)
            fieldp(j,i) = fieldtp(j,i) + fieldp(j,i)
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
