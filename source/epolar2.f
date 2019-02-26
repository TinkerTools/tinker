c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine epolar2  --  induced dipole polarization Hessian  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "epolar2" calculates second derivatives of the dipole polarization
c     energy for a single atom at a time
c
c     it is incorrect to neglect interactions with atoms not directly
c     involved as the multipole site; to get better accuracy, "list"
c     should include all atoms by setting "biglist" to "true"
c
c     the "twosided" flag controls use of one-sided vs. two-sided
c     numerical derivatives; setting the flag to "true" gives a more
c     accurate Hessian at the expense of increased computation time
c
c     also, the "reinduce" flag controls whether the induced dipoles
c     are recomputed every time an atom is moved during computation
c     of the numerical Hessian; setting the flag to "true" produces a
c     much slower calculation, but aids convergence of minimizations,
c     accuracy of vibrational frequencies, etc.
c
c
      subroutine epolar2 (i)
      use atoms
      use deriv
      use hessn
      use limits
      use mpole
      use potent
      implicit none
      integer i,j,k
      integer nlist
      integer, allocatable :: list(:)
      real*8 eps,old
      real*8, allocatable :: d0(:,:)
      logical prior
      logical biglist
      logical twosided
      logical reinduce
      logical save_mpole
c
c
c     set the default stepsize and accuracy control flags
c
      eps = 1.0d-5
      biglist = .false.
      twosided = .false.
      reinduce = .false.
      if (n .le. 300) then
         biglist = .true.
         twosided = .true.
         reinduce = .true.
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(npole))
      allocate (d0(3,n))
c
c     perform dynamic allocation of some global arrays
c
      prior = .false.
      if (allocated(dep)) then
         prior = .true.
         if (size(dep) .lt. 3*n)  deallocate (dep)
      end if
      if (.not. allocated(dep))  allocate (dep(3,n))
c
c     find the multipole definition involving the current atom;
c     results in a faster but approximate Hessian calculation
c
      nlist = 0
      do k = 1, npole
         if (biglist .or. ipole(k).eq.i) then
            nlist = nlist + 1
            list(nlist) = k
         end if
      end do
c
c     turn off multipole term as needed for Ewald gradient calls
c
      save_mpole = use_mpole
      use_mpole = .false.
c
c     get multipole first derivatives for the base structure
c
      if (.not. twosided) then
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
c
c     find numerical x-components via perturbed structures
c
      old = x(i)
      if (twosided) then
         x(i) = x(i) - 0.5d0*eps
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      x(i) = x(i) + eps
      if (use_ewald) then
         if (use_mlist) then
            call epolar1d
         else
            call epolar1c
         end if
      else
         call epolar2a (nlist,list,reinduce)
      end if
      x(i) = old
      do k = 1, n
         do j = 1, 3
            hessx(j,k) = hessx(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical y-components via perturbed structures
c
      old = y(i)
      if (twosided) then
         y(i) = y(i) - 0.5d0*eps
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      y(i) = y(i) + eps
      if (use_ewald) then
         if (use_mlist) then
            call epolar1d
         else
            call epolar1c
         end if
      else
         call epolar2a (nlist,list,reinduce)
      end if
      y(i) = old
      do k = 1, n
         do j = 1, 3
            hessy(j,k) = hessy(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     find numerical z-components via perturbed structures
c
      old = z(i)
      if (twosided) then
         z(i) = z(i) - 0.5d0*eps
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            call epolar2a (nlist,list,reinduce)
         end if
         do k = 1, n
            do j = 1, 3
               d0(j,k) = dep(j,k)
            end do
         end do
      end if
      z(i) = z(i) + eps
      if (use_ewald) then
         if (use_mlist) then
            call epolar1d
         else
            call epolar1c
         end if
      else
         call epolar2a (nlist,list,reinduce)
      end if
      z(i) = old
      do k = 1, n
         do j = 1, 3
            hessz(j,k) = hessz(j,k) + (dep(j,k)-d0(j,k))/eps
         end do
      end do
c
c     restore permanent multipole term to its original status
c
      use_mpole = save_mpole
c
c     perform deallocation of some global arrays
c
      if (.not. prior) then
         deallocate (dep)
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (d0)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine epolar2a  --  numerical polarization Hessian  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "epolar2a" computes polarization first derivatives for a single
c     atom with respect to Cartesian coordinates; used to get finite
c     difference second derivatives
c
c
      subroutine epolar2a (nlist,list,reinduce)
      use atoms
      use bound
      use cell
      use chgpen
      use chgpot
      use couple
      use deriv
      use mplpot
      use mpole
      use polar
      use polgrp
      use polopt
      use polpot
      use poltcg
      use potent
      use shunt
      implicit none
      integer i,j,k,m
      integer ii,kk,iii
      integer nlist,jcell
      integer list(*)
      real*8 f,damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 temp3,temp5,temp7
      real*8 sc3,sc5,sc7
      real*8 sr3,sr5,sr7
      real*8 psr3,psr5,psr7
      real*8 dsr3,dsr5,dsr7
      real*8 psr3i,psr5i,psr7i
      real*8 dsr3i,dsr5i,dsr7i
      real*8 psr3k,psr5k,psr7k
      real*8 dsr3k,dsr5k,dsr7k
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 uix,uiy,uiz
      real*8 uixp,uiyp,uizp
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 ukx,uky,ukz
      real*8 ukxp,ukyp,ukzp
      real*8 dir,uir,uirp
      real*8 dkr,ukr,ukrp
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 uirm,ukrm
      real*8 uirt,ukrt
      real*8 tuir,tukr
      real*8 tixx,tiyy,tizz
      real*8 tixy,tixz,tiyz
      real*8 tkxx,tkyy,tkzz
      real*8 tkxy,tkxz,tkyz
      real*8 tix3,tiy3,tiz3
      real*8 tix5,tiy5,tiz5
      real*8 tkx3,tky3,tkz3
      real*8 tkx5,tky5,tkz5
      real*8 term1,term2,term3
      real*8 term4,term5,term6
      real*8 term7,term8
      real*8 term1core
      real*8 term1i,term2i,term3i
      real*8 term4i,term5i,term6i
      real*8 term7i,term8i
      real*8 term1k,term2k,term3k
      real*8 term4k,term5k,term6k
      real*8 term7k,term8k
      real*8 depx,depy,depz
      real*8 frcx,frcy,frcz
      real*8 rc3(3),rc5(3),rc7(3)
      real*8 tep(3),fix(3)
      real*8 fiy(3),fiz(3)
      real*8 uax(3),uay(3),uaz(3)
      real*8 ubx(3),uby(3),ubz(3)
      real*8 uaxp(3),uayp(3),uazp(3)
      real*8 ubxp(3),ubyp(3),ubzp(3)
      real*8 dmpi(9),dmpk(9)
      real*8 dmpik(9)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: ufld(:,:)
      real*8, allocatable :: dufld(:,:)
      logical reinduce
      character*6 mode
c
c
c     zero out the polarization derivative components
c
      do i = 1, n
         do j = 1, 3
            dep(j,i) = 0.0d0
         end do
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
c     compute the induced dipoles at each polarizable atom
c
      if (reinduce)  call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (ufld(3,n))
      allocate (dufld(6,n))
c
c     set exclusion coefficients and arrays to store fields
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         do j = 1, 3
            ufld(j,i) = 0.0d0
         end do
         do j = 1, 6
            dufld(j,i) = 0.0d0
         end do
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = 0.5d0 * electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     compute the dipole polarization gradient components
c
      do iii = 1, nlist
         ii = list(iii)
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            wscale(i12(j,i)) = w2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            wscale(i13(j,i)) = w3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            wscale(i14(j,i)) = w4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            if (ii .eq. kk)  goto 10
            k = ipole(kk)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            if (use_bounds)  call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     set initial values for tha damping scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                        sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                        temp3 = -damp * expdamp * rr5
                        temp5 = -3.0d0 * damp / r2
                        temp7 = -(1.0d0+3.0d0*damp) / r2
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                        rc7(1) = rc5(1) * temp7
                        rc7(2) = rc5(2) * temp7
                        rc7(3) = rc5(3) * temp7
                     end if
                  end if
                  sr3 = rr3 * sc3
                  sr5 = rr5 * sc5
                  sr7 = rr7 * sc7
                  dsr3 = sr3 * dscale(k)
                  dsr5 = sr5 * dscale(k)
                  dsr7 = sr7 * dscale(k)
                  psr3 = sr3 * pscale(k)
                  psr5 = sr5 * pscale(k)
                  psr7 = sr7 * pscale(k)
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  dsr3i = 2.0d0 * rr3 * dmpi(3) * dscale(k)
                  dsr5i = 2.0d0 * rr5 * dmpi(5) * dscale(k)
                  dsr7i = 2.0d0 * rr7 * dmpi(7) * dscale(k)
                  dsr3k = 2.0d0 * rr3 * dmpk(3) * dscale(k)
                  dsr5k = 2.0d0 * rr5 * dmpk(5) * dscale(k)
                  dsr7k = 2.0d0 * rr7 * dmpk(7) * dscale(k)
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = dsr3i*ukx
                  tiy3 = dsr3i*uky
                  tiz3 = dsr3i*ukz
                  tkx3 = dsr3k*uix
                  tky3 = dsr3k*uiy
                  tkz3 = dsr3k*uiz
                  tuir = -dsr5i*ukr
                  tukr = -dsr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 2.0d0 * (dsr5i*ukx)
                  tiy5 = 2.0d0 * (dsr5i*uky)
                  tiz5 = 2.0d0 * (dsr5i*ukz)
                  tkx5 = 2.0d0 * (dsr5k*uix)
                  tky5 = 2.0d0 * (dsr5k*uiy)
                  tkz5 = 2.0d0 * (dsr5k*uiz)
                  tuir = -dsr7i*ukr
                  tukr = -dsr7k*uir
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
               if (use_thole) then
                  term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
                  term2 = (sc3+sc5)*rr5*xr - rc3(1)
                  term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
                  term6 = xr * (sc7*rr9*xr-rc7(1))
                  tixx = ci*term1 + dix*term2 - dir*term3
     &                      - qixx*term4 + qix*term5 - qir*term6
     &                      + (qiy*yr+qiz*zr)*sc7*rr7
                  tkxx = ck*term1 - dkx*term2 + dkr*term3
     &                      - qkxx*term4 + qkx*term5 - qkr*term6
     &                      + (qky*yr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
                  term2 = (sc3+sc5)*rr5*yr - rc3(2)
                  term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
                  term6 = yr * (sc7*rr9*yr-rc7(2))
                  tiyy = ci*term1 + diy*term2 - dir*term3
     &                      - qiyy*term4 + qiy*term5 - qir*term6
     &                      + (qix*xr+qiz*zr)*sc7*rr7
                  tkyy = ck*term1 - dky*term2 + dkr*term3
     &                      - qkyy*term4 + qky*term5 - qkr*term6
     &                      + (qkx*xr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
                  term2 = (sc3+sc5)*rr5*zr - rc3(3)
                  term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
                  term6 = zr * (sc7*rr9*zr-rc7(3))
                  tizz = ci*term1 + diz*term2 - dir*term3
     &                      - qizz*term4 + qiz*term5 - qir*term6
     &                      + (qix*xr+qiy*yr)*sc7*rr7
                  tkzz = ck*term1 - dkz*term2 + dkr*term3
     &                      - qkzz*term4 + qkz*term5 - qkr*term6
     &                      + (qkx*xr+qky*yr)*sc7*rr7
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = yr * term2
                  term3 = sc5 * rr5 * yr
                  term4 = yr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * yr
                  term8 = yr * (sc7*rr9*xr-rc7(1))
                  tixy = -ci*term1 + diy*term2 + dix*term3
     &                      - dir*term4 - qixy*term5 + qiy*term6
     &                      + qix*term7 - qir*term8
                  tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                      + dkr*term4 - qkxy*term5 + qky*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*xr-rc7(1))
                  tixz = -ci*term1 + diz*term2 + dix*term3
     &                      - dir*term4 - qixz*term5 + qiz*term6
     &                      + qix*term7 - qir*term8
                  tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                      + dkr*term4 - qkxz*term5 + qkz*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*yr - rc3(2)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*yr-rc5(2))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*yr-rc7(2))
                  tiyz = -ci*term1 + diz*term2 + diy*term3
     &                      - dir*term4 - qiyz*term5 + qiz*term6
     &                      + qiy*term7 - qir*term8
                  tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                      + dkr*term4 - qkyz*term5 + qkz*term6
     &                      + qky*term7 - qkr*term8
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*xr*xr
                  term1core = rr3 - rr5*xr*xr
                  term2i = 2.0d0*rr5*dmpi(5)*xr 
                  term3i = rr7*dmpi(7)*xr*xr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*xr
                  term6i = rr9*dmpi(9)*xr*xr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*xr*xr
                  term2k = 2.0d0*rr5*dmpk(5)*xr
                  term3k = rr7*dmpk(7)*xr*xr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*xr
                  term6k = rr9*dmpk(9)*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7*dmpi(7)
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*yr*yr
                  term1core = rr3 - rr5*yr*yr
                  term2i = 2.0d0*rr5*dmpi(5)*yr
                  term3i = rr7*dmpi(7)*yr*yr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*yr
                  term6i = rr9*dmpi(9)*yr*yr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*yr*yr
                  term2k = 2.0d0*rr5*dmpk(5)*yr
                  term3k = rr7*dmpk(7)*yr*yr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*yr
                  term6k = rr9*dmpk(9)*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7*dmpi(7)
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*zr*zr
                  term1core = rr3 - rr5*zr*zr
                  term2i = 2.0d0*rr5*dmpi(5)*zr
                  term3i = rr7*dmpi(7)*zr*zr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*zr
                  term6i = rr9*dmpi(9)*zr*zr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*zr*zr
                  term2k = 2.0d0*rr5*dmpk(5)*zr
                  term3k = rr7*dmpk(7)*zr*zr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*zr
                  term6k = rr9*dmpk(9)*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7*dmpi(7)
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7*dmpk(7)
                  term2i = rr5*dmpi(5)*xr 
                  term1i = yr * term2i
                  term1core = rr5*xr*yr
                  term3i = rr5*dmpi(5)*yr
                  term4i = yr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*yr
                  term8i = yr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = yr * term2k
                  term3k = rr5*dmpk(5)*yr
                  term4k = yr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*yr
                  term8k = yr*rr9*dmpk(9)*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*xr
                  term1i = zr * term2i
                  term1core = rr5*xr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*yr
                  term1i = zr * term2i
                  term1core = rr5*yr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*yr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*yr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*yr
                  term2k = rr5*dmpk(5)*yr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*yr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*yr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
               end if
c
c     get the dEd/dR terms for Thole direct polarization force
c
               if (use_thole) then
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = dscale(k) * depx
                  frcy = dscale(k) * depy
                  frcz = dscale(k) * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + pscale(k)*depx
                  frcy = frcy + pscale(k)*depy
                  frcz = frcz + pscale(k)*depz
c
c     get the dEp/dR terms for chgpen direct polarization force
c
               else if (use_chgpen) then
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = 2.0d0*dscale(k)*depx
                  frcy = 2.0d0*dscale(k)*depy
                  frcz = 2.0d0*dscale(k)*depz
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + uscale(kk)*depx
                  frcy = frcy + uscale(kk)*depy
                  frcz = frcz + uscale(kk)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * dmpik(5) * rr5
                  term2 = term1*xr
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr 
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr 
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr 
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5*dmpik(5)*yr
                  term2 = rr5*dmpik(5)*xr 
                  term3 = yr * (rr7*dmpik(7)*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5 *dmpik(5) * zr
                  term3 = zr * (rr7*dmpik(7)*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5*dmpik(5)*yr 
                  term3 = zr * (rr7*dmpik(7)*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + wscale(kk)*depx
                  frcy = frcy + wscale(kk)*depy
                  frcz = frcz + wscale(kk)*depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp(1:3).eq.'OPT' .and. use_thole) then
                  do j = 0, coptmax-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, coptmax-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = (sc3+sc5) * rr5
                        term2 = term1*xr - rc3(1)
                        term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                        tixx = uopt(j,1,ii)*term2 + uirm*term3
                        tkxx = uopt(m,1,kk)*term2 + ukrm*term3
                        term2 = term1*yr - rc3(2)
                        term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                        tiyy = uopt(j,2,ii)*term2 + uirm*term3
                        tkyy = uopt(m,2,kk)*term2 + ukrm*term3
                        term2 = term1*zr - rc3(3)
                        term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                        tizz = uopt(j,3,ii)*term2 + uirm*term3
                        tkzz = uopt(m,3,kk)*term2 + ukrm*term3
                        term1 = sc5 * rr5 * yr
                        term2 = sc3*rr5*xr - rc3(1)
                        term3 = yr * (sc5*rr7*xr-rc5(1))
                        tixy = uopt(j,1,ii)*term1 + uopt(j,2,ii)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,kk)*term1 + uopt(m,2,kk)*term2
     &                            - ukrm*term3
                        term1 = sc5 * rr5 * zr
                        term3 = zr * (sc5*rr7*xr-rc5(1))
                        tixz = uopt(j,1,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        term2 = sc3*rr5*yr - rc3(2)
                        term3 = zr * (sc5*rr7*yr-rc5(2))
                        tiyz = uopt(j,2,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*uscale(k)*depx
                        frcy = frcy + copm(j+m+1)*uscale(k)*depy
                        frcz = frcz + copm(j+m+1)*uscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp(1:3).eq.'OPT' .and. use_chgpen) then
                  do j = 0, coptmax-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, coptmax-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * dmpik(5) * rr5
                        term2 = term1*xr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5*dmpik(5)*yr
                        term2 = rr5*dmpik(5)*xr
                        term3 = yr * (rr7*dmpik(7)*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2
     &                            - ukrm*term3
                        term1 = rr5 *dmpik(5) * zr
                        term3 = zr * (rr7*dmpik(7)*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5*dmpik(5)*yr
                        term3 = zr * (rr7*dmpik(7)*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx + copm(j+m+1)*wscale(k)*depx
                        frcy = frcy + copm(j+m+1)*wscale(k)*depy
                        frcz = frcz + copm(j+m+1)*wscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp(1:3) .eq. 'TCG') then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = uax(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uay(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = uaz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = uax(j)*term1 + uay(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = uax(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uay(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = ubx(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uby(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = ubz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = ubx(j)*term1 + uby(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = ubx(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uby(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                  end do
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) + frcx
               dep(2,i) = dep(2,i) + frcy
               dep(3,i) = dep(3,i) + frcz
               dep(1,k) = dep(1,k) - frcx
               dep(2,k) = dep(2,k) - frcy
               dep(3,k) = dep(3,k) - frcz
            end if
   10       continue
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c
      do iii = 1, nlist
         ii = list(iii)
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
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         uixp = uinp(1,ii)
         uiyp = uinp(2,ii)
         uizp = uinp(3,ii)
         do j = 1, tcgnab
            uax(j) = uad(1,ii,j)
            uay(j) = uad(2,ii,j)
            uaz(j) = uad(3,ii,j)
            uaxp(j) = uap(1,ii,j)
            uayp(j) = uap(2,ii,j)
            uazp(j) = uap(3,ii,j)
            ubx(j) = ubd(1,ii,j)
            uby(j) = ubd(2,ii,j)
            ubz(j) = ubd(3,ii,j)
            ubxp(j) = ubp(1,ii,j)
            ubyp(j) = ubp(2,ii,j)
            ubzp(j) = ubp(3,ii,j)
         end do
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         do j = 1, n12(i)
            pscale(i12(j,i)) = p2scale
            wscale(i12(j,i)) = w2scale
         end do
         do j = 1, n13(i)
            pscale(i13(j,i)) = p3scale
            wscale(i13(j,i)) = w3scale
         end do
         do j = 1, n14(i)
            pscale(i14(j,i)) = p4scale
            wscale(i14(j,i)) = w4scale
            do k = 1, np11(i)
                if (i14(j,i) .eq. ip11(k,i))
     &            pscale(i14(j,i)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(i)
            pscale(i15(j,i)) = p5scale
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
c
c     evaluate all sites within the cutoff distance
c
         do kk = 1, npole
            k = ipole(kk)
            do jcell = 2, ncell
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            if (use_bounds)  call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
            if (.not. (use_polymer .and. r2.le.polycut2)) then
               dscale(k) = 1.0d0
               pscale(k) = 1.0d0
               uscale(k) = 1.0d0
            end if
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
               ukx = uind(1,kk)
               uky = uind(2,kk)
               ukz = uind(3,kk)
               ukxp = uinp(1,kk)
               ukyp = uinp(2,kk)
               ukzp = uinp(3,kk)
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
               uir = uix*xr + uiy*yr + uiz*zr
               uirp = uixp*xr + uiyp*yr + uizp*zr
               ukr = ukx*xr + uky*yr + ukz*zr
               ukrp = ukxp*xr + ukyp*yr + ukzp*zr
c
c     get reciprocal distance terms for this interaction
c
               rr1 = f / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
c
c     set initial values for tha damping scale factors
c
               sc3 = 1.0d0
               sc5 = 1.0d0
               sc7 = 1.0d0
               do j = 1, 3
                  rc3(j) = 0.0d0
                  rc5(j) = 0.0d0
                  rc7(j) = 0.0d0
               end do
c
c     apply Thole polarization damping to scale factors
c
               if (use_thole) then
                  damp = pdi * pdamp(kk)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(kk))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        sc3 = 1.0d0 - expdamp
                        sc5 = 1.0d0 - (1.0d0-damp)*expdamp
                        sc7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                        temp3 = -damp * expdamp * rr5
                        temp5 = -3.0d0 * damp / r2
                        temp7 = -(1.0d0+3.0d0*damp) / r2
                        rc3(1) = xr * temp3
                        rc3(2) = yr * temp3
                        rc3(3) = zr * temp3
                        rc5(1) = rc3(1) * temp5
                        rc5(2) = rc3(2) * temp5
                        rc5(3) = rc3(3) * temp5
                        rc7(1) = rc5(1) * temp7
                        rc7(2) = rc5(2) * temp7
                        rc7(3) = rc5(3) * temp7
                     end if
                  end if
                  sr3 = rr3 * sc3
                  sr5 = rr5 * sc5
                  sr7 = rr7 * sc7
                  dsr3 = sr3 * dscale(k)
                  dsr5 = sr5 * dscale(k)
                  dsr7 = sr7 * dscale(k)
                  psr3 = sr3 * pscale(k)
                  psr5 = sr5 * pscale(k)
                  psr7 = sr7 * pscale(k)
c
c     apply charge penetration damping to scale factors
c
               else if (use_chgpen) then
                  corek = pcore(kk)
                  valk = pval(kk)
                  alphak = palpha(kk)
                  call damppole (r,9,alphai,alphak,dmpi,dmpk,dmpik)
                  dsr3i = 2.0d0 * rr3 * dmpi(3) * dscale(k)
                  dsr5i = 2.0d0 * rr5 * dmpi(5) * dscale(k)
                  dsr7i = 2.0d0 * rr7 * dmpi(7) * dscale(k)
                  dsr3k = 2.0d0 * rr3 * dmpk(3) * dscale(k)
                  dsr5k = 2.0d0 * rr5 * dmpk(5) * dscale(k)
                  dsr7k = 2.0d0 * rr7 * dmpk(7) * dscale(k)
               end if
c
c     get the induced dipole field used for dipole torques
c
               if (use_thole) then
                  tix3 = psr3*ukx + dsr3*ukxp
                  tiy3 = psr3*uky + dsr3*ukyp
                  tiz3 = psr3*ukz + dsr3*ukzp
                  tkx3 = psr3*uix + dsr3*uixp
                  tky3 = psr3*uiy + dsr3*uiyp
                  tkz3 = psr3*uiz + dsr3*uizp
                  tuir = -psr5*ukr - dsr5*ukrp
                  tukr = -psr5*uir - dsr5*uirp
               else if (use_chgpen) then
                  tix3 = dsr3i*ukx
                  tiy3 = dsr3i*uky
                  tiz3 = dsr3i*ukz
                  tkx3 = dsr3k*uix
                  tky3 = dsr3k*uiy
                  tkz3 = dsr3k*uiz
                  tuir = -dsr5i*ukr
                  tukr = -dsr5k*uir
               end if
               ufld(1,i) = ufld(1,i) + tix3 + xr*tuir
               ufld(2,i) = ufld(2,i) + tiy3 + yr*tuir
               ufld(3,i) = ufld(3,i) + tiz3 + zr*tuir
               ufld(1,k) = ufld(1,k) + tkx3 + xr*tukr
               ufld(2,k) = ufld(2,k) + tky3 + yr*tukr
               ufld(3,k) = ufld(3,k) + tkz3 + zr*tukr
c
c     get induced dipole field gradient used for quadrupole torques
c
               if (use_thole) then
                  tix5 = 2.0d0 * (psr5*ukx+dsr5*ukxp)
                  tiy5 = 2.0d0 * (psr5*uky+dsr5*ukyp)
                  tiz5 = 2.0d0 * (psr5*ukz+dsr5*ukzp)
                  tkx5 = 2.0d0 * (psr5*uix+dsr5*uixp)
                  tky5 = 2.0d0 * (psr5*uiy+dsr5*uiyp)
                  tkz5 = 2.0d0 * (psr5*uiz+dsr5*uizp)
                  tuir = -psr7*ukr - dsr7*ukrp
                  tukr = -psr7*uir - dsr7*uirp
               else if (use_chgpen) then
                  tix5 = 2.0d0 * (dsr5i*ukx)
                  tiy5 = 2.0d0 * (dsr5i*uky)
                  tiz5 = 2.0d0 * (dsr5i*ukz)
                  tkx5 = 2.0d0 * (dsr5k*uix)
                  tky5 = 2.0d0 * (dsr5k*uiy)
                  tkz5 = 2.0d0 * (dsr5k*uiz)
                  tuir = -dsr7i*ukr
                  tukr = -dsr7k*uir
               end if
               dufld(1,i) = dufld(1,i) + xr*tix5 + xr*xr*tuir
               dufld(2,i) = dufld(2,i) + xr*tiy5 + yr*tix5
     &                         + 2.0d0*xr*yr*tuir
               dufld(3,i) = dufld(3,i) + yr*tiy5 + yr*yr*tuir
               dufld(4,i) = dufld(4,i) + xr*tiz5 + zr*tix5
     &                         + 2.0d0*xr*zr*tuir
               dufld(5,i) = dufld(5,i) + yr*tiz5 + zr*tiy5
     &                         + 2.0d0*yr*zr*tuir
               dufld(6,i) = dufld(6,i) + zr*tiz5 + zr*zr*tuir
               dufld(1,k) = dufld(1,k) - xr*tkx5 - xr*xr*tukr
               dufld(2,k) = dufld(2,k) - xr*tky5 - yr*tkx5
     &                         - 2.0d0*xr*yr*tukr
               dufld(3,k) = dufld(3,k) - yr*tky5 - yr*yr*tukr
               dufld(4,k) = dufld(4,k) - xr*tkz5 - zr*tkx5
     &                         - 2.0d0*xr*zr*tukr
               dufld(5,k) = dufld(5,k) - yr*tkz5 - zr*tky5
     &                         - 2.0d0*yr*zr*tukr
               dufld(6,k) = dufld(6,k) - zr*tkz5 - zr*zr*tukr
c
c     get the field gradient for direct polarization force
c
               if (use_thole) then
                  term1 = sc3*(rr3-rr5*xr*xr) + rc3(1)*xr
                  term2 = (sc3+sc5)*rr5*xr - rc3(1)
                  term3 = sc5*(rr7*xr*xr-rr5) - rc5(1)*xr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*xr-rc5(1)+1.5d0*sc7*rr7*xr)
                  term6 = xr * (sc7*rr9*xr-rc7(1))
                  tixx = ci*term1 + dix*term2 - dir*term3
     &                      - qixx*term4 + qix*term5 - qir*term6
     &                      + (qiy*yr+qiz*zr)*sc7*rr7
                  tkxx = ck*term1 - dkx*term2 + dkr*term3
     &                      - qkxx*term4 + qkx*term5 - qkr*term6
     &                      + (qky*yr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*yr*yr) + rc3(2)*yr
                  term2 = (sc3+sc5)*rr5*yr - rc3(2)
                  term3 = sc5*(rr7*yr*yr-rr5) - rc5(2)*yr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*yr-rc5(2)+1.5d0*sc7*rr7*yr)
                  term6 = yr * (sc7*rr9*yr-rc7(2))
                  tiyy = ci*term1 + diy*term2 - dir*term3
     &                      - qiyy*term4 + qiy*term5 - qir*term6
     &                      + (qix*xr+qiz*zr)*sc7*rr7
                  tkyy = ck*term1 - dky*term2 + dkr*term3
     &                      - qkyy*term4 + qky*term5 - qkr*term6
     &                      + (qkx*xr+qkz*zr)*sc7*rr7
                  term1 = sc3*(rr3-rr5*zr*zr) + rc3(3)*zr
                  term2 = (sc3+sc5)*rr5*zr - rc3(3)
                  term3 = sc5*(rr7*zr*zr-rr5) - rc5(3)*zr
                  term4 = 2.0d0 * sc5 * rr5
                  term5 = 2.0d0 * (sc5*rr7*zr-rc5(3)+1.5d0*sc7*rr7*zr)
                  term6 = zr * (sc7*rr9*zr-rc7(3))
                  tizz = ci*term1 + diz*term2 - dir*term3
     &                      - qizz*term4 + qiz*term5 - qir*term6
     &                      + (qix*xr+qiy*yr)*sc7*rr7
                  tkzz = ck*term1 - dkz*term2 + dkr*term3
     &                      - qkzz*term4 + qkz*term5 - qkr*term6
     &                      + (qkx*xr+qky*yr)*sc7*rr7
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = yr * term2
                  term3 = sc5 * rr5 * yr
                  term4 = yr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * yr
                  term8 = yr * (sc7*rr9*xr-rc7(1))
                  tixy = -ci*term1 + diy*term2 + dix*term3
     &                      - dir*term4 - qixy*term5 + qiy*term6
     &                      + qix*term7 - qir*term8
                  tkxy = -ck*term1 - dky*term2 - dkx*term3
     &                      + dkr*term4 - qkxy*term5 + qky*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*xr - rc3(1)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*xr-rc5(1))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*xr-rc5(1))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*xr-rc7(1))
                  tixz = -ci*term1 + diz*term2 + dix*term3
     &                      - dir*term4 - qixz*term5 + qiz*term6
     &                      + qix*term7 - qir*term8
                  tkxz = -ck*term1 - dkz*term2 - dkx*term3
     &                      + dkr*term4 - qkxz*term5 + qkz*term6
     &                      + qkx*term7 - qkr*term8
                  term2 = sc3*rr5*yr - rc3(2)
                  term1 = zr * term2
                  term3 = sc5 * rr5 * zr
                  term4 = zr * (sc5*rr7*yr-rc5(2))
                  term5 = 2.0d0 * sc5 * rr5
                  term6 = 2.0d0 * (sc5*rr7*yr-rc5(2))
                  term7 = 2.0d0 * sc7 * rr7 * zr
                  term8 = zr * (sc7*rr9*yr-rc7(2))
                  tiyz = -ci*term1 + diz*term2 + diy*term3
     &                      - dir*term4 - qiyz*term5 + qiz*term6
     &                      + qiy*term7 - qir*term8
                  tkyz = -ck*term1 - dkz*term2 - dky*term3
     &                      + dkr*term4 - qkyz*term5 + qkz*term6
     &                      + qky*term7 - qkr*term8
c
c     get the field gradient for direct polarization force
c
               else if (use_chgpen) then
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*xr*xr
                  term1core = rr3 - rr5*xr*xr
                  term2i = 2.0d0*rr5*dmpi(5)*xr 
                  term3i = rr7*dmpi(7)*xr*xr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*xr
                  term6i = rr9*dmpi(9)*xr*xr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*xr*xr
                  term2k = 2.0d0*rr5*dmpk(5)*xr
                  term3k = rr7*dmpk(7)*xr*xr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*xr
                  term6k = rr9*dmpk(9)*xr*xr
                  tixx = vali*term1i + corei*term1core  
     &                      + dix*term2i - dir*term3i
     &                      - qixx*term4i + qix*term5i - qir*term6i
     &                      + (qiy*yr+qiz*zr)*rr7*dmpi(7)
                  tkxx = valk*term1k + corek*term1core
     &                      - dkx*term2k + dkr*term3k
     &                      - qkxx*term4k + qkx*term5k - qkr*term6k
     &                      + (qky*yr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*yr*yr
                  term1core = rr3 - rr5*yr*yr
                  term2i = 2.0d0*rr5*dmpi(5)*yr
                  term3i = rr7*dmpi(7)*yr*yr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*yr
                  term6i = rr9*dmpi(9)*yr*yr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*yr*yr
                  term2k = 2.0d0*rr5*dmpk(5)*yr
                  term3k = rr7*dmpk(7)*yr*yr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*yr
                  term6k = rr9*dmpk(9)*yr*yr
                  tiyy = vali*term1i + corei*term1core
     &                      + diy*term2i - dir*term3i
     &                      - qiyy*term4i + qiy*term5i - qir*term6i
     &                      + (qix*xr+qiz*zr)*rr7*dmpi(7)
                  tkyy = valk*term1k + corek*term1core
     &                      - dky*term2k + dkr*term3k
     &                      - qkyy*term4k + qky*term5k - qkr*term6k
     &                      + (qkx*xr+qkz*zr)*rr7*dmpk(7)
                  term1i = rr3*dmpi(3) - rr5*dmpi(5)*zr*zr
                  term1core = rr3 - rr5*zr*zr
                  term2i = 2.0d0*rr5*dmpi(5)*zr
                  term3i = rr7*dmpi(7)*zr*zr - rr5*dmpi(5)
                  term4i = 2.0d0*rr5*dmpi(5)
                  term5i = 5.0d0*rr7*dmpi(7)*zr
                  term6i = rr9*dmpi(9)*zr*zr
                  term1k = rr3*dmpk(3) - rr5*dmpk(5)*zr*zr
                  term2k = 2.0d0*rr5*dmpk(5)*zr
                  term3k = rr7*dmpk(7)*zr*zr - rr5*dmpk(5)
                  term4k = 2.0d0*rr5*dmpk(5)
                  term5k = 5.0d0*rr7*dmpk(7)*zr
                  term6k = rr9*dmpk(9)*zr*zr
                  tizz = vali*term1i + corei*term1core
     &                      + diz*term2i - dir*term3i
     &                      - qizz*term4i + qiz*term5i - qir*term6i
     &                      + (qix*xr+qiy*yr)*rr7*dmpi(7)
                  tkzz = valk*term1k + corek*term1core
     &                      - dkz*term2k + dkr*term3k
     &                      - qkzz*term4k + qkz*term5k - qkr*term6k
     &                      + (qkx*xr+qky*yr)*rr7*dmpk(7)
                  term2i = rr5*dmpi(5)*xr 
                  term1i = yr * term2i
                  term1core = rr5*xr*yr
                  term3i = rr5*dmpi(5)*yr
                  term4i = yr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*yr
                  term8i = yr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = yr * term2k
                  term3k = rr5*dmpk(5)*yr
                  term4k = yr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*yr
                  term8k = yr*rr9*dmpk(9)*xr
                  tixy = -vali*term1i - corei*term1core 
     &                      + diy*term2i + dix*term3i
     &                      - dir*term4i - qixy*term5i + qiy*term6i
     &                      + qix*term7i - qir*term8i
                  tkxy = -valk*term1k - corek*term1core 
     &                      - dky*term2k - dkx*term3k
     &                      + dkr*term4k - qkxy*term5k + qky*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*xr
                  term1i = zr * term2i
                  term1core = rr5*xr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*xr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*xr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*xr
                  term2k = rr5*dmpk(5)*xr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*xr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*xr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*xr
                  tixz = -vali*term1i - corei*term1core
     &                      + diz*term2i + dix*term3i
     &                      - dir*term4i - qixz*term5i + qiz*term6i
     &                      + qix*term7i - qir*term8i
                  tkxz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dkx*term3k
     &                      + dkr*term4k - qkxz*term5k + qkz*term6k
     &                      + qkx*term7k - qkr*term8k
                  term2i = rr5*dmpi(5)*yr
                  term1i = zr * term2i
                  term1core = rr5*yr*zr
                  term3i = rr5*dmpi(5)*zr
                  term4i = zr * (rr7*dmpi(7)*yr)
                  term5i = 2.0d0*rr5*dmpi(5)
                  term6i = 2.0d0*rr7*dmpi(7)*yr
                  term7i = 2.0d0*rr7*dmpi(7)*zr
                  term8i = zr*rr9*dmpi(9)*yr
                  term2k = rr5*dmpk(5)*yr
                  term1k = zr * term2k
                  term3k = rr5*dmpk(5)*zr
                  term4k = zr * (rr7*dmpk(7)*yr)
                  term5k = 2.0d0*rr5*dmpk(5)
                  term6k = 2.0d0*rr7*dmpk(7)*yr
                  term7k = 2.0d0*rr7*dmpk(7)*zr
                  term8k = zr*rr9*dmpk(9)*yr
                  tiyz = -vali*term1i - corei*term1core
     &                      + diz*term2i + diy*term3i
     &                      - dir*term4i - qiyz*term5i + qiz*term6i
     &                      + qiy*term7i - qir*term8i
                  tkyz = -valk*term1k - corek*term1core
     &                      - dkz*term2k - dky*term3k
     &                      + dkr*term4k - qkyz*term5k + qkz*term6k
     &                      + qky*term7k - qkr*term8k
               end if
c
c     get the dEd/dR terms for Thole direct polarization force
c
               if (use_thole) then
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      - tkxx*uixp - tkxy*uiyp - tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      - tkxy*uixp - tkyy*uiyp - tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      - tkxz*uixp - tkyz*uiyp - tkzz*uizp
                  frcx = dscale(k) * depx
                  frcy = dscale(k) * depy
                  frcz = dscale(k) * depz
c
c     get the dEp/dR terms for Thole direct polarization force
c
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = frcx + pscale(k)*depx
                  frcy = frcy + pscale(k)*depy
                  frcz = frcz + pscale(k)*depz
c
c     get the dEp/dR terms for chgpen direct polarization force
c
               else if (use_chgpen) then
                  depx = tixx*ukx + tixy*uky + tixz*ukz
     &                      - tkxx*uix - tkxy*uiy - tkxz*uiz
                  depy = tixy*ukx + tiyy*uky + tiyz*ukz
     &                      - tkxy*uix - tkyy*uiy - tkyz*uiz
                  depz = tixz*ukx + tiyz*uky + tizz*ukz
     &                      - tkxz*uix - tkyz*uiy - tkzz*uiz
                  frcx = 2.0d0*dscale(k)*depx
                  frcy = 2.0d0*dscale(k)*depy
                  frcz = 2.0d0*dscale(k)*depz
               end if
c
c     get the dtau/dr terms used for mutual polarization force
c
               if (poltyp.eq.'MUTUAL' .and. use_thole) then
                  term1 = (sc3+sc5) * rr5
                  term2 = term1*xr - rc3(1)
                  term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr - rc3(2)
                  term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr - rc3(3)
                  term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = sc5 * rr5 * yr
                  term2 = sc3*rr5*xr - rc3(1)
                  term3 = yr * (sc5*rr7*xr-rc5(1))
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = sc5 * rr5 * zr
                  term3 = zr * (sc5*rr7*xr-rc5(1))
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = sc3*rr5*yr - rc3(2)
                  term3 = zr * (sc5*rr7*yr-rc5(2))
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + uscale(kk)*depx
                  frcy = frcy + uscale(kk)*depy
                  frcz = frcz + uscale(kk)*depz
c
c     get the dtau/dr terms used for mutual polarization force
c
               else if (poltyp.eq.'MUTUAL' .and. use_chgpen) then
                  term1 = 2.0d0 * dmpik(5) * rr5
                  term2 = term1*xr
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr 
                  tixx = uix*term2 + uir*term3
                  tkxx = ukx*term2 + ukr*term3
                  term2 = term1*yr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr 
                  tiyy = uiy*term2 + uir*term3
                  tkyy = uky*term2 + ukr*term3
                  term2 = term1*zr 
                  term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr 
                  tizz = uiz*term2 + uir*term3
                  tkzz = ukz*term2 + ukr*term3
                  term1 = rr5*dmpik(5)*yr
                  term2 = rr5*dmpik(5)*xr 
                  term3 = yr * (rr7*dmpik(7)*xr)
                  tixy = uix*term1 + uiy*term2 - uir*term3
                  tkxy = ukx*term1 + uky*term2 - ukr*term3
                  term1 = rr5 *dmpik(5) * zr
                  term3 = zr * (rr7*dmpik(7)*xr)
                  tixz = uix*term1 + uiz*term2 - uir*term3
                  tkxz = ukx*term1 + ukz*term2 - ukr*term3
                  term2 = rr5*dmpik(5)*yr 
                  term3 = zr * (rr7*dmpik(7)*yr)
                  tiyz = uiy*term1 + uiz*term2 - uir*term3
                  tkyz = uky*term1 + ukz*term2 - ukr*term3
                  depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                      + tkxx*uixp + tkxy*uiyp + tkxz*uizp
                  depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                      + tkxy*uixp + tkyy*uiyp + tkyz*uizp
                  depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                      + tkxz*uixp + tkyz*uiyp + tkzz*uizp
                  frcx = frcx + wscale(kk)*depx
                  frcy = frcy + wscale(kk)*depy
                  frcz = frcz + wscale(kk)*depz
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp(1:3).eq.'OPT' .and. use_thole) then
                  do j = 0, coptmax-1
                     uirm = uopt(j,1,ii)*xr + uopt(j,2,ii)*yr
     &                          + uopt(j,3,ii)*zr
                     do m = 0, coptmax-j-1
                        ukrm = uopt(m,1,kk)*xr + uopt(m,2,kk)*yr
     &                             + uopt(m,3,kk)*zr
                        term1 = (sc3+sc5) * rr5
                        term2 = term1*xr - rc3(1)
                        term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                        tixx = uopt(j,1,ii)*term2 + uirm*term3
                        tkxx = uopt(m,1,kk)*term2 + ukrm*term3
                        term2 = term1*yr - rc3(2)
                        term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                        tiyy = uopt(j,2,ii)*term2 + uirm*term3
                        tkyy = uopt(m,2,kk)*term2 + ukrm*term3
                        term2 = term1*zr - rc3(3)
                        term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                        tizz = uopt(j,3,ii)*term2 + uirm*term3
                        tkzz = uopt(m,3,kk)*term2 + ukrm*term3
                        term1 = sc5 * rr5 * yr
                        term2 = sc3*rr5*xr - rc3(1)
                        term3 = yr * (sc5*rr7*xr-rc5(1))
                        tixy = uopt(j,1,ii)*term1 + uopt(j,2,ii)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,kk)*term1 + uopt(m,2,kk)*term2
     &                            - ukrm*term3
                        term1 = sc5 * rr5 * zr
                        term3 = zr * (sc5*rr7*xr-rc5(1))
                        tixz = uopt(j,1,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        term2 = sc3*rr5*yr - rc3(2)
                        term3 = zr * (sc5*rr7*yr-rc5(2))
                        tiyz = uopt(j,2,ii)*term1 + uopt(j,3,ii)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,kk)*term1 + uopt(m,3,kk)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,kk) + tkxx*uoptp(j,1,ii)
     &                       + tixy*uoptp(m,2,kk) + tkxy*uoptp(j,2,ii)
     &                       + tixz*uoptp(m,3,kk) + tkxz*uoptp(j,3,ii)
                        depy = tixy*uoptp(m,1,kk) + tkxy*uoptp(j,1,ii)
     &                       + tiyy*uoptp(m,2,kk) + tkyy*uoptp(j,2,ii)
     &                       + tiyz*uoptp(m,3,kk) + tkyz*uoptp(j,3,ii)
                        depz = tixz*uoptp(m,1,kk) + tkxz*uoptp(j,1,ii)
     &                       + tiyz*uoptp(m,2,kk) + tkyz*uoptp(j,2,ii)
     &                       + tizz*uoptp(m,3,kk) + tkzz*uoptp(j,3,ii)
                        frcx = frcx + copm(j+m+1)*uscale(k)*depx
                        frcy = frcy + copm(j+m+1)*uscale(k)*depy
                        frcz = frcz + copm(j+m+1)*uscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for OPT polarization force
c
               else if (poltyp(1:3).eq.'OPT' .and. use_chgpen) then
                  do j = 0, coptmax-1
                     uirm = uopt(j,1,i)*xr + uopt(j,2,i)*yr
     &                          + uopt(j,3,i)*zr
                     do m = 0, coptmax-j-1
                        ukrm = uopt(m,1,k)*xr + uopt(m,2,k)*yr
     &                             + uopt(m,3,k)*zr
                        term1 = 2.0d0 * dmpik(5) * rr5
                        term2 = term1*xr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*xr*xr
                        tixx = uopt(j,1,i)*term2 + uirm*term3
                        tkxx = uopt(m,1,k)*term2 + ukrm*term3
                        term2 = term1*yr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*yr*yr
                        tiyy = uopt(j,2,i)*term2 + uirm*term3
                        tkyy = uopt(m,2,k)*term2 + ukrm*term3
                        term2 = term1*zr
                        term3 = rr5*dmpik(5) - rr7*dmpik(7)*zr*zr
                        tizz = uopt(j,3,i)*term2 + uirm*term3
                        tkzz = uopt(m,3,k)*term2 + ukrm*term3
                        term1 = rr5*dmpik(5)*yr
                        term2 = rr5*dmpik(5)*xr
                        term3 = yr * (rr7*dmpik(7)*xr)
                        tixy = uopt(j,1,i)*term1 + uopt(j,2,i)*term2
     &                            - uirm*term3
                        tkxy = uopt(m,1,k)*term1 + uopt(m,2,k)*term2
     &                            - ukrm*term3
                        term1 = rr5 *dmpik(5) * zr
                        term3 = zr * (rr7*dmpik(7)*xr)
                        tixz = uopt(j,1,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkxz = uopt(m,1,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        term2 = rr5*dmpik(5)*yr
                        term3 = zr * (rr7*dmpik(7)*yr)
                        tiyz = uopt(j,2,i)*term1 + uopt(j,3,i)*term2
     &                            - uirm*term3
                        tkyz = uopt(m,2,k)*term1 + uopt(m,3,k)*term2
     &                            - ukrm*term3
                        depx = tixx*uoptp(m,1,k) + tkxx*uoptp(j,1,i)
     &                       + tixy*uoptp(m,2,k) + tkxy*uoptp(j,2,i)
     &                       + tixz*uoptp(m,3,k) + tkxz*uoptp(j,3,i)
                        depy = tixy*uoptp(m,1,k) + tkxy*uoptp(j,1,i)
     &                       + tiyy*uoptp(m,2,k) + tkyy*uoptp(j,2,i)
     &                       + tiyz*uoptp(m,3,k) + tkyz*uoptp(j,3,i)
                        depz = tixz*uoptp(m,1,k) + tkxz*uoptp(j,1,i)
     &                       + tiyz*uoptp(m,2,k) + tkyz*uoptp(j,2,i)
     &                       + tizz*uoptp(m,3,k) + tkzz*uoptp(j,3,i)
                        frcx = frcx + copm(j+m+1)*wscale(k)*depx
                        frcy = frcy + copm(j+m+1)*wscale(k)*depy
                        frcz = frcz + copm(j+m+1)*wscale(k)*depz
                     end do
                  end do
c
c     get the dtau/dr terms used for TCG polarization force
c
               else if (poltyp(1:3) .eq. 'TCG') then
                  do j = 1, tcgnab
                     ukx = ubd(1,kk,j)
                     uky = ubd(2,kk,j)
                     ukz = ubd(3,kk,j)
                     ukxp = ubp(1,kk,j)
                     ukyp = ubp(2,kk,j)
                     ukzp = ubp(3,kk,j)
                     uirt = uax(j)*xr + uay(j)*yr + uaz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = uax(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uay(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = uaz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = uax(j)*term1 + uay(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = uax(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uay(j)*term1 + uaz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*uaxp(j) + tkxy*uayp(j)
     &                         + tkxz*uazp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*uaxp(j) + tkyy*uayp(j)
     &                         + tkyz*uazp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*uaxp(j) + tkyz*uayp(j)
     &                         + tkzz*uazp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                     ukx = uad(1,kk,j)
                     uky = uad(2,kk,j)
                     ukz = uad(3,kk,j)
                     ukxp = uap(1,kk,j)
                     ukyp = uap(2,kk,j)
                     ukzp = uap(3,kk,j)
                     uirt = ubx(j)*xr + uby(j)*yr + ubz(j)*zr
                     ukrt = ukx*xr + uky*yr + ukz*zr
                     term1 = (sc3+sc5) * rr5
                     term2 = term1*xr - rc3(1)
                     term3 = sc5*(rr5-rr7*xr*xr) + rc5(1)*xr
                     tixx = ubx(j)*term2 + uirt*term3
                     tkxx = ukx*term2 + ukrt*term3
                     term2 = term1*yr - rc3(2)
                     term3 = sc5*(rr5-rr7*yr*yr) + rc5(2)*yr
                     tiyy = uby(j)*term2 + uirt*term3
                     tkyy = uky*term2 + ukrt*term3
                     term2 = term1*zr - rc3(3)
                     term3 = sc5*(rr5-rr7*zr*zr) + rc5(3)*zr
                     tizz = ubz(j)*term2 + uirt*term3
                     tkzz = ukz*term2 + ukrt*term3
                     term1 = sc5 * rr5 * yr
                     term2 = sc3*rr5*xr - rc3(1)
                     term3 = yr * (sc5*rr7*xr-rc5(1))
                     tixy = ubx(j)*term1 + uby(j)*term2 - uirt*term3
                     tkxy = ukx*term1 + uky*term2 - ukrt*term3
                     term1 = sc5 * rr5 * zr
                     term3 = zr * (sc5*rr7*xr-rc5(1))
                     tixz = ubx(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkxz = ukx*term1 + ukz*term2 - ukrt*term3
                     term2 = sc3*rr5*yr - rc3(2)
                     term3 = zr * (sc5*rr7*yr-rc5(2))
                     tiyz = uby(j)*term1 + ubz(j)*term2 - uirt*term3
                     tkyz = uky*term1 + ukz*term2 - ukrt*term3
                     depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
     &                         + tkxx*ubxp(j) + tkxy*ubyp(j)
     &                         + tkxz*ubzp(j)
                     depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
     &                         + tkxy*ubxp(j) + tkyy*ubyp(j)
     &                         + tkyz*ubzp(j)
                     depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
     &                         + tkxz*ubxp(j) + tkyz*ubyp(j)
     &                         + tkzz*ubzp(j)
                     frcx = frcx + uscale(k)*depx
                     frcy = frcy + uscale(k)*depy
                     frcz = frcz + uscale(k)*depz
                  end do
               end if
c
c     force and torque components scaled for self-interactions
c
               if (i .eq. k) then
                  frcx = 0.5d0 * frcx
                  frcy = 0.5d0 * frcy
                  frcz = 0.5d0 * frcz
                  psr3 = 0.5d0 * psr3
                  psr5 = 0.5d0 * psr5
                  psr7 = 0.5d0 * psr7
                  dsr3 = 0.5d0 * dsr3
                  dsr5 = 0.5d0 * dsr5
                  dsr7 = 0.5d0 * dsr7
               end if
c
c     increment force-based gradient on the interaction sites
c
               dep(1,i) = dep(1,i) + frcx
               dep(2,i) = dep(2,i) + frcy
               dep(3,i) = dep(3,i) + frcz
               dep(1,k) = dep(1,k) - frcx
               dep(2,k) = dep(2,k) - frcy
               dep(3,k) = dep(3,k) - frcz
            end if
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
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
      end do
      end if
c
c     torque is induced field and gradient cross permanent moments
c
      do ii = 1, npole
         i = ipole(ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         tep(1) = diz*ufld(2,i) - diy*ufld(3,i)
     &               + qixz*dufld(2,i) - qixy*dufld(4,i)
     &               + 2.0d0*qiyz*(dufld(3,i)-dufld(6,i))
     &               + (qizz-qiyy)*dufld(5,i)
         tep(2) = dix*ufld(3,i) - diz*ufld(1,i)
     &               - qiyz*dufld(2,i) + qixy*dufld(5,i)
     &               + 2.0d0*qixz*(dufld(6,i)-dufld(1,i))
     &               + (qixx-qizz)*dufld(4,i)
         tep(3) = diy*ufld(1,i) - dix*ufld(2,i)
     &               + qiyz*dufld(4,i) - qixz*dufld(5,i)
     &               + 2.0d0*qixy*(dufld(1,i)-dufld(3,i))
     &               + (qiyy-qixx)*dufld(2,i)
         call torque (ii,tep,fix,fiy,fiz,dep)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (ufld)
      deallocate (dufld)
      return
      end
