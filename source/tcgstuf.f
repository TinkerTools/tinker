c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2018 by Zhi Wang and Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine induce0b  --  truncated CG dipole solver  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "induce0b" computes and stores the induced dipoles via
c     the truncated conjugate gradient (TCG) method
c
c
      subroutine induce0b
      use poltcg
      implicit none
c
c
c     choose the options for computation of TCG induced dipoles
c
      if (tcgguess) then
         call indtcgb
      else
         call indtcga
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine indtcga  --  TCG zero guess and preconditioner  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "indtcga" computes the induced dipoles and intermediates used
c     in polarization force calculation for the TCG method with dp
c     cross terms = true, initial guess mu0 = 0 and using a diagonal
c     preconditioner
c
c
      subroutine indtcga
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      use potent
      implicit none
      integer i,j,order
      real*8 n0,np0,g0
      real*8 n1,np1,g1,beta1
      real*8 n2,np2,g2,beta2
      real*8 n3,beta3
      real*8 a100,a101,a102
      real*8 a103,b111
      real*8, allocatable :: rsd(:,:,:)
      real*8, allocatable :: r0(:,:,:)
      real*8, allocatable :: p1(:,:,:)
      real*8, allocatable :: p2(:,:,:)
      real*8, allocatable :: p3(:,:,:)
      real*8, allocatable :: tp(:,:,:)
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
c     set up nab based on tcgorder value
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some global arrays
c
      if (.not. allocated(uad))  allocate (uad(3,n,tcgnab))
      if (.not. allocated(uap))  allocate (uap(3,n,tcgnab))
      if (.not. allocated(ubd))  allocate (ubd(3,n,tcgnab))
      if (.not. allocated(ubp))  allocate (ubp(3,n,tcgnab))
      uad = 0.0d0
      uap = 0.0d0
      ubd = 0.0d0
      ubp = 0.0d0
c
c     perform dynamic allocation for some local arrays
c
      allocate (rsd(3,n,2))
      allocate (r0(3,n,2))
      allocate (p1(3,n,2))
      allocate (p2(3,n,2))
      allocate (p3(3,n,2))
      allocate (tp(3,n,2))
c
c     get the electrostaic field due to permanent multipoles
c     because mu0 = 0, r0 = field - T.mu0 = field
c
      if (use_ewald) then
         call dfield0c (r0(:,:,1),r0(:,:,2))
      else if (use_mlist) then
         call dfield0b (r0(:,:,1),r0(:,:,2))
      else
         call dfield0a (r0(:,:,1),r0(:,:,2))
      end if
c
c     udir = alpha.E = alpha.r0
c
      call tcg_alpha22 (r0(:,:,1),r0(:,:,2),udir,udirp)
c
c     compute the following tcg1 intermediates:
c     p0 = alpha*r0 = udir
c     n0 = r0*a*r0
c     np0 = p0*T*p0
c
      call tcg_alpha22 (r0(:,:,1),r0(:,:,2),udir,udirp)
      call tcg_alphaquad (n0,r0(:,:,1),r0(:,:,2))
      call tcg_t0 (udir,udirp,tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np0,3*npole,tp(:,:,1),udirp)
      g0 = n0 / np0
c
c     r1 = r0 - gamma0 * T*p0
c     n1 = r1*a*r1
c     p1 <- r1, p0
c
      rsd = r0 - g0 * tp
      call tcg_alphaquad (n1,rsd(:,:,1),rsd(:,:,2))
      beta1 = n1 / n0
      p1(:,:,1) = udir
      p1(:,:,2) = udirp
      call tcg_update (p1(:,:,1),rsd(:,:,1),beta1)
      call tcg_update (p1(:,:,2),rsd(:,:,2),beta1)
c
c     ua(1) = mu1 = g0 * p0
c     ub(1) <- p0
c     xde <- p0, p1
c
      uad(:,:,1) = g0*udir
      uap(:,:,1) = g0*udirp
      ubd(:,:,1) = ubd(:,:,1) + 0.5d0*g0*udir
      ubp(:,:,1) = ubp(:,:,1) + 0.5d0*g0*udirp
      uind = uind + g0*(1.0d0-beta1)*udir + g0*p1(:,:,1)
      uinp = uinp + g0*(1.0d0-beta1)*udirp + g0*p1(:,:,2)
c
c     the tcg1 energy and force are finished
c
      if (order .eq. 1)  goto 10
c
c     np1 = p1*T*p1
c     g1 = n1 / np1
c     r2 = r1 - g1 * T*p1
c     n2 = r2*a*r2
c     beta2 = n2 / n1
c
      call tcg_t0 (p1(:,:,1),p1(:,:,2),tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np1,3*npole,tp(:,:,1),p1(:,:,2))
      g1 = n1 / np1
      rsd = rsd - g1 * tp
      call tcg_alphaquad (n2,rsd(:,:,1),rsd(:,:,2))
      beta2 = n2 / n1
c
c     p2 <- r2, p1
c     np2 = p2*T*p2
c     g2 = n2 / np2
c
      p2 = p1
      call tcg_update (p2(:,:,1),rsd(:,:,1),beta2)
      call tcg_update (p2(:,:,2),rsd(:,:,2),beta2)
      call tcg_t0 (p2(:,:,1),p2(:,:,2),tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np2,3*npole,tp(:,:,1),p2(:,:,2))
      g2 = n2 / np2
c
c     r3 = r2 - g2 * T*p2
c     n3 = r3*a*r3
c     beta3 = n3 / n2
c
      rsd = rsd - g2 * tp
      call tcg_alphaquad (n3,rsd(:,:,1),rsd(:,:,2))
      beta3 = n3 / n2
c
c     p3 <- r3, p2
c
      p3 = p2
      call tcg_update (p3(:,:,1),rsd(:,:,1),beta3)
      call tcg_update (p3(:,:,2),rsd(:,:,2),beta3)
c
c     ua(2) = mu2 = g1 * p1
c     ub(1) <- p1, p2
c     ub(2) <- p1
c     xde <- p0, p1
c
      b111 = (1.0d0-beta2) * g1
      a103 = g0*g1/g2
      a102 = (1.0d0-beta2)*g0 + (1.0d0+beta1)*g1 - (1.0d0+beta3)*a103
      a101 = (beta2**2-1.0d0)*g0 + (1.0d0-beta2-beta1*beta2)*g1
     &          + beta2*a103
      a100 = (1.0d0-beta2)*g0*beta1
      uad(:,:,2) = g1*p1(:,:,1)
      uap(:,:,2) = g1*p1(:,:,2)
      ubd(:,:,1) = ubd(:,:,1) + b111*p1(:,:,1) + g1*p2(:,:,1)
      ubp(:,:,1) = ubp(:,:,1) + b111*p1(:,:,2) + g1*p2(:,:,2)
      ubd(:,:,2) = ubd(:,:,2) + 0.5d0*g1*p1(:,:,1)
      ubp(:,:,2) = ubp(:,:,2) + 0.5d0*g1*p1(:,:,2)
      uind = uind + a103*p3(:,:,1) + a102*p2(:,:,1) + a101*p1(:,:,1)
     &          + a100*udir
      uinp = uinp + a103*p3(:,:,2) + a102*p2(:,:,2) + a101*p1(:,:,2)
     &          + a100*udirp
c
c     perform deallocation for some local arrays
c
   10 continue
      deallocate (rsd)
      deallocate (r0)
      deallocate (p1)
      deallocate (p2)
      deallocate (p3)
      deallocate (tp)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine indtcgb  --  TCG direct guess and precondition  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "indtcgb" computes the induced dipoles and intermediates used
c     in polarization force calculation for the TCG method with dp
c     cross terms = true, initial guess mu0 = direct and using diagonal
c     preconditioner
c
c
      subroutine indtcgb
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      use potent
      implicit none
      integer i,j,order
      real*8 chi,xi0,xi1
      real*8 n0,np0,g0
      real*8 n1,np1,g1,beta1
      real*8 n2,np2,g2,beta2
      real*8 n3,beta3
      real*8 a100,a101,a102
      real*8 a103,b111
      real*8 c100,c101,c102
      real*8 c200,c201,c202
      real*8 c203,c204
      real*8 d111,d210,d211
      real*8 d212,d213,d222
      real*8, allocatable :: xdr0(:,:,:)
      real*8, allocatable :: rsd(:,:,:)
      real*8, allocatable :: p0(:,:,:)
      real*8, allocatable :: p1(:,:,:)
      real*8, allocatable :: p2(:,:,:)
      real*8, allocatable :: p3(:,:,:)
      real*8, allocatable :: tp(:,:,:)
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
c     set up nab based on tcgorder value
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some global arrays
c
      if (.not. allocated(uad))  allocate (uad(3,n,tcgnab))
      if (.not. allocated(uap))  allocate (uap(3,n,tcgnab))
      if (.not. allocated(ubd))  allocate (ubd(3,n,tcgnab))
      if (.not. allocated(ubp))  allocate (ubp(3,n,tcgnab))
      uad = 0.0d0
      uap = 0.0d0
      ubd = 0.0d0
      ubp = 0.0d0
c
c     perform dynamic allocation for some local arrays
c
      allocate (xdr0(3,n,2))
      allocate (rsd(3,n,2))
      allocate (p0(3,n,2))
      allocate (p1(3,n,2))
      allocate (p2(3,n,2))
      allocate (p3(3,n,2))
      allocate (tp(3,n,2))
      xdr0 = 0.0d0
c
c     chi = omega - 1
c
      chi = tcgpeek - 1.0d0
c
c     get the electrostatic field due to permanent multipoles
c     and mu0 = alpha*E; use tp to store the multipole field
c
      if (use_ewald) then
         call dfield0c (tp(:,:,1),tp(:,:,2))
      else if (use_mlist) then
         call dfield0b (tp(:,:,1),tp(:,:,2))
      else
         call dfield0a (tp(:,:,1),tp(:,:,2))
      end if
      call tcg_alpha22 (tp(:,:,1),tp(:,:,2),udir,udirp)
c
c     compute the following tcg1 intermediates:
c     r0 = -Tu*mu0
c     n0 = r0*a*r0
c     p0 = a*r0
c     xi0
c     np0 = p0*T*p0
c     g0
c
      call tcg_ufield (udir,udirp,rsd(:,:,1),rsd(:,:,2))
      call tcg_alphaquad (n0,rsd(:,:,1),rsd(:,:,2))
      call tcg_alpha22 (rsd(:,:,1),rsd(:,:,2),p0(:,:,1),p0(:,:,2))
      call tcg_dotprod (xi0,3*npole,rsd(:,:,1),udirp)
      xi0 = xi0 / n0
      call tcg_t0 (p0(:,:,1),p0(:,:,2),tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np0,3*npole,tp(:,:,1),p0(:,:,2))
      g0 = n0 / np0
c
c     r1 = r0 - g0 T*p0
c     n1, beta1
c
      rsd = rsd - g0 * tp
      call tcg_alphaquad (n1,rsd(:,:,1),rsd(:,:,2))
      beta1 = n1 / n0
c
c     p1 <- r1, p0
c
      p1 = p0
      call tcg_update (p1(:,:,1),rsd(:,:,1),beta1)
      call tcg_update (p1(:,:,2),rsd(:,:,2),beta1)
c
c     compute "Residual Mutual 1"
c     ua(1) <- mu0
c     ua(2) <- mu1 = g0 * p0
c     ub(2) <- p0
c     xdr0 <- p0, p1
c
      uad(:,:,1) = udir
      uap(:,:,1) = udirp
      ubd(:,:,1) = ubd(:,:,1) + 0.5d0*udir
      ubp(:,:,1) = ubp(:,:,1) + 0.5d0*udirp
      uad(:,:,2) = g0*p0(:,:,1)
      uap(:,:,2) = g0*p0(:,:,2)
      ubd(:,:,2) = ubd(:,:,2) + 0.5d0*g0*p0(:,:,1)
      ubp(:,:,2) = ubp(:,:,2) + 0.5d0*g0*p0(:,:,2)
      xdr0 = xdr0 + g0*(1.0d0-beta1)*p0 + g0*p1
c
c     get the tcg1 energy and force; tp array works as xde array
c
      if (order .eq. 1) then
         c100 = 0.5d0*(1.0d0-g0)
         c101 = (0.5d0 - beta1*(1.0d0-xi0))*g0
         c102 = g0*(1.0d0-xi0)
         d111 = 0.5d0*(1.0d0-xi0)*g0
         xdr0(:,:,1) = xdr0(:,:,1) + chi*(c100*udir
     &                  + c101*p0(:,:,1) + c102*p1(:,:,1))
         xdr0(:,:,2) = xdr0(:,:,2) + chi*(c100*udirp
     &                  + c101*p0(:,:,2) + c102*p1(:,:,2))
         ubd(:,:,1) = ubd(:,:,1) + xdr0(:,:,1)
         ubp(:,:,1) = ubp(:,:,1) + xdr0(:,:,2)
         ubd(:,:,2) = ubd(:,:,2) + chi*(d111*p0(:,:,1)+0.5d0*udir)
         ubp(:,:,2) = ubp(:,:,2) + chi*(d111*p0(:,:,2)+0.5d0*udirp)
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),tp(:,:,1),tp(:,:,2))
         call tcg_alpha12 (tp(:,:,1),tp(:,:,2))
         tp(:,:,1) = tp(:,:,1) + chi*0.5d0*p1(:,:,1)
     &             + (1.0d0-chi*beta1*0.5d0)*p0(:,:,1) + udir
         tp(:,:,2) = tp(:,:,2) + chi*0.5d0*p1(:,:,2)
     &             + (1.0d0-chi*beta1*0.5d0)*p0(:,:,2) + udirp
         goto 10
      end if
c
c     compute the tcg2 intermediates: xi1, np1 and g1
c
      call tcg_dotprod (xi1,3*npole,rsd(:,:,1),udirp)
      xi1 = xi1 / n1 + xi0
      call tcg_t0 (p1(:,:,1),p1(:,:,2),tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np1,3*npole,tp(:,:,1),p1(:,:,2))
      g1 = n1 / np1
c
c     r2 = r1 - g1 * T*p1
c     n2, beta2
c     p2 <- r2, p1
c     np2, g2
c
      rsd = rsd - g1 * tp
      call tcg_alphaquad (n2,rsd(:,:,1),rsd(:,:,2))
      beta2 = n2 / n1
      p2 = p1
      call tcg_update (p2(:,:,1),rsd(:,:,1),beta2)
      call tcg_update (p2(:,:,2),rsd(:,:,2),beta2)
      call tcg_t0 (p2(:,:,1),p2(:,:,2),tp(:,:,1),tp(:,:,2))
      call tcg_dotprod (np2,3*npole,tp(:,:,1),p2(:,:,2))
      g2 = n2/np2
c
c     r3 = r2 - g2 * T*p2
c     n3, beta3
c     p3 <- r3, p2
c
      rsd = rsd - g2 * tp
      call tcg_alphaquad (n3,rsd(:,:,1),rsd(:,:,2))
      beta3 = n3 / n2
      p3 = p2
      call tcg_update (p3(:,:,1),rsd(:,:,1),beta3)
      call tcg_update (p3(:,:,2),rsd(:,:,2),beta3)
c
c     compute "Residual Mutual 2"
c     ub(2) <- p1, p2
c     ua(3) <- mu2 = g1 * p1
c     ub(3) <- p1
c     xdr0 <- p0, p1, p2, p3
c
      b111 = (1.0d0-beta2) * g1
      a103 = g0*g1/g2
      a102 = (1.0d0-beta2)*g0 + (1.0d0+beta1)*g1 - (1.0d0+beta3)*a103
      a101 = (beta2**2-1.0d0)*g0 + (1.0d0-beta2-beta1*beta2)*g1
     &          + beta2*a103
      a100 = (1.0d0-beta2)*g0*beta1
      ubd(:,:,2) = ubd(:,:,2)+ b111*p1(:,:,1) + g1*p2(:,:,1)
      ubp(:,:,2) = ubp(:,:,2)+ b111*p1(:,:,2) + g1*p2(:,:,2)
      uad(:,:,3) = g1*p1(:,:,1)
      uap(:,:,3) = g1*p1(:,:,2)
      ubd(:,:,3) = ubd(:,:,3) + 0.5d0*g1*p1(:,:,1)
      ubp(:,:,3) = ubp(:,:,3) + 0.5d0*g1*p1(:,:,2)
      xdr0 = xdr0 + a100*p0 + a101*p1 + a102*p2 + a103*p3
c
c     get the tcg2 energy and force; tp array works as xde array
c
      if (order .eq. 2) then
         c200 = 0.5d0*((1.0d0-g0)*(1.0d0-g1)-beta1*g1)
         c201 = 0.5d0*(1.0d0-g1)*g0 + (xi0-xi1)*g1*beta1**2
     &             + (xi1-1.0d0)*beta1*beta2*g0
     &             + (xi0+g0-xi0*g0)*beta1*g1
         c202 = 0.5d0*g1 + (1.0d0-xi1)*beta2*g0*g1/g2
     &             + (1.0d0-xi1)*(beta2*g0-(1.0d0+beta1)*g1)*beta2
     &             + (xi0*g0-xi0-g0)*g1 + (beta2*g0-beta1*g1)*(xi0-xi1)
         c203 = (xi1-1.0d0)*(1.0d0+beta3)*g0*g1/g2
     &             + (g1+beta1*g1-beta2*g0)*(1.0d0-xi1) + (1.0d0-xi0)*g0
         c204 = (1.0d0-xi1)*g0*g1/g2
         d210 = 0.5d0*(1.0d0-g1)
         d211 = 0.5d0*(1.0d0-xi0)*(1.0d0-g1)*g0
         d212 = ((1.0d0-xi0)-(1.0d0-xi1)*beta2)*g1
         d213 = (1.0d0-xi1)*g1
         d222 = 0.5d0*d213
         xdr0(:,:,1) = xdr0(:,:,1) + chi*(c204*p3(:,:,1)
     &               + c203*p2(:,:,1) + c202*p1(:,:,1)
     &               + c201*p0(:,:,1) + c200*udir)
         xdr0(:,:,2) = xdr0(:,:,2) + chi*(c204*p3(:,:,2)
     &               + c203*p2(:,:,2) + c202*p1(:,:,2)
     &               + c201*p0(:,:,2) + c200*udirp)
         ubd(:,:,1) = ubd(:,:,1) + xdr0(:,:,1)
         ubp(:,:,1) = ubp(:,:,1) + xdr0(:,:,2)
         ubd(:,:,2) = ubd(:,:,2) + chi*(d213*p2(:,:,1) + d212*p1(:,:,1)
     &              + d211*p0(:,:,1) + d210*udir)
         ubp(:,:,2) = ubp(:,:,2) + chi*(d213*p2(:,:,2) + d212*p1(:,:,2)
     &              + d211*p0(:,:,2) + d210*udirp)
         ubd(:,:,3) = ubd(:,:,3) + chi*(d222*p1(:,:,1) + 0.5d0*udir)
         ubp(:,:,3) = ubp(:,:,3) + chi*(d222*p1(:,:,2) + 0.5d0*udirp)
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),tp(:,:,1),tp(:,:,2))
         call tcg_alpha12 (tp(:,:,1),tp(:,:,2))
         tp(:,:,1) = tp(:,:,1) + p0(:,:,1) + udir
     &             + chi*0.5d0*(p2(:,:,1)-beta2*p1(:,:,1))
         tp(:,:,2) = tp(:,:,2) + p0(:,:,2) + udirp
     &             + chi*0.5d0*(p2(:,:,2)-beta2*p1(:,:,2))
         goto 10
      end if
c
c     store induced dipoles from elements of the xde arrays
c
   10 continue
      uind = tp(:,:,1)
      uinp = tp(:,:,2)
c
c     perform deallocation for some local arrays
c
      deallocate (xdr0)
      deallocate (rsd)
      deallocate (p0)
      deallocate (p1)
      deallocate (p2)
      deallocate (p3)
      deallocate (tp)
      return
      end
c
c
c     ################################
c     ##                            ##
c     ##  subroutine tcg_alphaquad  ##
c     ##                            ##
c     ################################
c
c
c     "tcg_alphaquad" computes the quadratic form, <a*alpha*b>,
c     where alpha is the diagonal atomic polarizability matrix
c
c
      subroutine tcg_alphaquad (sum,a,b)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 sum
      real*8 a(3,*)
      real*8 b(3,*)
c
c
      sum = 0.0d0
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO reduction(+:sum) schedule(guided)
      do i = 1, npole
         do j = 1, 3
            sum = sum + a(j,i)*b(j,i)*polarity(i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_resource  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_resource" sets the number of mutual induced dipole
c     pairs based on the passed argument
c
c
      subroutine tcg_resource (order)
      use iounit
      use poltcg
      implicit none
      integer order
c
c
      if (order.lt.1 .or. order.gt.2) then
         write (iout,10)
   10    format (/,' TCG_RESOURCE -- Argument ORDER Is Out of Range')
         call fatal
      end if
      tcgnab = order
      if (tcgguess)  tcgnab = tcgnab + 1
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##  subroutine tcg_alpha12  ##
c     ##                          ##
c     ##############################
c
c
c     "tcg_alpha12" computes source1 = alpha*source1 and
c     source2 = alpha*source2
c
c
      subroutine tcg_alpha12 (source1,source2)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 source1(3,*)
      real*8 source2(3,*)
c
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, npole
         do j = 1, 3
            source1(j,i) = polarity(i) * source1(j,i)
            source2(j,i) = polarity(i) * source2(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##  subroutine tcg_alpha22  ##
c     ##                          ##
c     ##############################
c
c
c     "tcg_alpha22" computes result1 = alpha*source1 and
c     result2 = alpha*source2
c
c
      subroutine tcg_alpha22 (source1,source2,result1,result2)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 source1(3,*)
      real*8 source2(3,*)
      real*8 result1(3,*)
      real*8 result2(3,*)
c
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, npole
         do j = 1, 3
            result1(j,i) = polarity(i) * source1(j,i)
            result2(j,i) = polarity(i) * source2(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##############################
c     ##                          ##
c     ##  subroutine tcg_dotprod  ##
c     ##                          ##
c     ##############################
c
c
c     "tcg_dotprod" computes the dot product of two vectors
c     of length n elements
c
c
      subroutine tcg_dotprod (sum,n,a,b)
      implicit none
      integer i,n
      real*8 sum
      real*8 a(*)
      real*8 b(*)
c
c
c     find value of the scalar dot product
c
      sum = 0.0d0
!$OMP PARALLEL default(shared) private(i)
!$OMP DO reduction(+:sum) schedule(guided)
      do i = 1, n
         sum = sum + a(i)*b(i)
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #############################
c     ##                         ##
c     ##  subroutine tcg_ufield  ##
c     ##                         ##
c     #############################
c
c
c     "tcg_ufield" applies -Tu to ind/p and returns v3d/p
c
c
      subroutine tcg_ufield (ind,inp,v3d,v3p)
      use limits
      use mpole
      use polar
      implicit none
      real*8 ind(3,*)
      real*8 inp(3,*)
      real*8 v3d(3,*)
      real*8 v3p(3,*)
c
c
c     swap TCG components with induced dipoles
c
      call tcgswap (uind,uinp,ind,inp)
c
c     compute mutual field
c
      if (use_ewald) then
         call ufield0c (v3d,v3p)
      else if (use_mlist) then
         call ufield0b (v3d,v3p)
      else
         call ufield0a (v3d,v3p)
      end if
c
c     swap TCG components with induced dipoles
c
      call tcgswap (uind,uinp,ind,inp)
      return
      end
c
c
c     #########################
c     ##                     ##
c     ##  subroutine tcg_t0  ##
c     ##                     ##
c     #########################
c
c
c     "tcg_t0" applies T matrix to ind/p, and returns v3d/p
c     T = 1/alpha + Tu
c
c
      subroutine tcg_t0 (ind,inp,v3d,v3p)
      use limits
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 poli,polmin
      real*8 ind(3,*)
      real*8 inp(3,*)
      real*8 v3d(3,*)
      real*8 v3p(3,*)
c
c
c     apply -Tu to ind/p
c
      call tcg_ufield (ind,inp,v3d,v3p)
c
c     compute the 1/alpha contribution
c
      polmin = 0.00000001d0
!$OMP PARALLEL default(shared) private(i,j,poli)
!$OMP DO schedule(guided)
      do i = 1, npole
         if (douind(ipole(i))) then
            poli = max(polmin,polarity(i))
            do j = 1, 3
               v3d(j,i) = ind(j,i)/poli - v3d(j,i)
               v3p(j,i) = inp(j,i)/poli - v3p(j,i)
            end do
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine tcgswap  --  swap induced dipoles for TCG use  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "tcgswap" switches two sets of induced dipole quantities for
c     use with the TCG induced dipole solver
c
c
      subroutine tcgswap (uind1,uinp1,uind2,uinp2)
      use mpole
      implicit none
      integer i,j
      real*8 dterm,pterm
      real*8 uind1(3,*)
      real*8 uinp1(3,*)
      real*8 uind2(3,*)
      real*8 uinp2(3,*)
c
c
c     swap sets of induced dipoles for use with the TCG method
c
!$OMP PARALLEL default(shared) private(i,j,dterm,pterm)
!$OMP DO schedule(guided)
      do i = 1, npole
         do j = 1, 3
            dterm = uind1(j,i)
            pterm = uinp1(j,i)
            uind1(j,i) = uind2(j,i)
            uinp1(j,i) = uinp2(j,i)
            uind2(j,i) = dterm
            uinp2(j,i) = pterm
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine tcg_update  --  get an updated TCG p-vector  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "tcg_update" computes pvec = alpha*rvec + beta*pvec;
c     if the preconditioner is not used, then alpha = identity
c
c
      subroutine tcg_update (pvec,rvec,beta)
      use mpole
      use polar
      use poltcg
      implicit none
      integer i,j
      real*8 beta,alpha
      real*8 pvec(3,*)
      real*8 rvec(3,*)
c
c
c     computes an updated pvec from prior intermediates
c
!$OMP PARALLEL default(shared) private(i,j,alpha)
!$OMP DO schedule(guided)
      do i = 1, npole
         alpha = polarity(i)
         do j = 1, 3
            pvec(j,i) = alpha*rvec(j,i) + beta*pvec(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
