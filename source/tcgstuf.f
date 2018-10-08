c
c
c     ##############################
c     ##                          ##
c     ##  subroutine tcg_induce1  ##
c     ##                          ##
c     ##############################
c
c
c     "tcg_induce1" computes and stores the induced dipoles and
c     the derivatives of the induced dipoles used in polarization
c     force calculation for TCG method
c
c
      subroutine tcg_induce1
      use mpole
      use polar
      use poltcg
      use potent
      implicit none
      integer i,j
c
c
c     zero out the induced dipoles at each site
c
      if (.not. use_polar .or. npole .eq. 0)  return
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
c
c     select from alternative versions of the TCG algorithm
c
      if (tcgguess) then
         if (tcgprec) then
            call tcg_induce1a
         else
            call tcg_induce1b
         end if
      else
         if (tcgprec) then
            call tcg_induce1c
         else
            call tcg_induce1d
         end if
      end if
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_induce1a  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_induce1a" computes the induced dipoles and intermediates
c     used in polarization force calculation for the TCG method with
c     initial guess mu0 = direct and diagonal preconditioner = true
c
c
      subroutine tcg_induce1a
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      implicit none
      integer i,j
      integer order
      real*8 n0(2),t1(2)
      real*8 t4(2),nrsd1(2),n1(2)
      real*8 sp0(2),spp1(2)
      real*8 a110(2),a111(2),a112(2),a121(2)
      real*8 a1k10a(2),a1k11a(2),a1k11(2)
      real*8 a1k12(2),a1k20a(2),a1k21(2)
      real*8 t9(2),beta1(2),t2(2),np1(2)
      real*8 t8(2),t10(2),t3(2),gamma1(2)
      real*8 sp1(2),b1(2),b2(2),spp2(2)
      real*8 a210(2),a21n1(2),a211(2)
      real*8 a212(2),a213(2),a214(2)
      real*8 a220(2),a221(2),a222(2),a223(2)
      real*8 a231(2),a232(2),a241(2)
      real*8 a2k10a(2),a2k11a(2),a2k12a(2)
      real*8 a2k11(2),a2k12(2),a2k13(2),a2k14(2)
      real*8 a2k20a(2),a2k21a(2),a2k21(2)
      real*8 a2k22(2),a2k23(2)
      real*8 a2k30a(2),a2k31(2),a2k32(2),a2k41(2)
      real*8 a2kwt2(2),a2kwg1(2)
      real*8, allocatable :: field(:,:,:)
      real*8, allocatable :: xde(:,:,:)
      real*8, allocatable :: xdr0(:,:,:)
      real*8, allocatable :: r0(:,:,:)
      real*8, allocatable :: m0(:,:,:)
      real*8, allocatable :: P1(:,:,:)
      real*8, allocatable :: r1(:,:,:)
      real*8, allocatable :: t2m0(:,:,:)
      real*8, allocatable :: t3m0(:,:,:)
      logical converge
c
c
c     set up nab based on tcgorder
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some local arrays
c
      if (.not. allocated(uindt)) then
         allocate (uindt(3,n))
         allocate (uinpt(3,n))
         allocate (uad(3,n,tcgnab))
         allocate (uap(3,n,tcgnab))
         allocate (ubd(3,n,tcgnab))
         allocate (ubp(3,n,tcgnab))
      end if
      allocate (field(3,n,2))
      allocate (xde(3,n,2))
      allocate (xdr0(3,n,2))
      allocate (r0(3,n,2))
      allocate (m0(3,n,2))
      allocate (P1(3,n,2))
      allocate (r1(3,n,2))
      allocate (t2m0(3,n,2))
      allocate (t3m0(3,n,2))
c
c     get the electrostatic field due to permanent multipoles
c     and mu0 = alpha.E
c
      if (use_ewald) then
         call dfield0c (field(:,:,1),field(:,:,2))
      else if (use_mlist) then
         call dfield0b (field(:,:,1),field(:,:,2))
      else
         call dfield0a (field(:,:,1),field(:,:,2))
      end if
      call tcg_alpha22 (field(:,:,1),field(:,:,2),udir,udirp)
      uind = udir
      uinp = udirp
c
c     r0 = -Tu.mu0
c
      if (use_ewald) then
         call ufield0c (r0(:,:,1),r0(:,:,2))
      else if (use_mlist) then
         call ufield0b (r0(:,:,1),r0(:,:,2))
      else
         call ufield0a (r0(:,:,1),r0(:,:,2))
      end if
c
c     m0 = M.r0, M is the preconditioner matrix
c     n0 = r0.m0
c
      call tcg_alpha22 (r0(:,:,1),r0(:,:,2),m0(:,:,1),m0(:,:,2))
      call tcg_dotprod (n0(1),3*npole,r0(:,:,1),m0(:,:,1))
      call tcg_dotprod (n0(2),3*npole,r0(:,:,2),m0(:,:,2))
c
c     compute tcg1 intermediates
c     P1 = T.m0 = T.M.r0
c     t1 = m0.T.m0
c     (t4 or gamma0) = r0/t1
c
      call tcg_t0 (m0(:,:,1),m0(:,:,2),P1(:,:,1),P1(:,:,2))
      call tcg_dotprod (t1(1),3*npole,m0(:,:,1),P1(:,:,1))
      call tcg_dotprod (t1(2),3*npole,m0(:,:,2),P1(:,:,2))
      t4(1) = n0(1) / t1(1)
      t4(2) = n0(2) / t1(2)
c
c     mu1 = mu0 + gamma0 * p0 (or m0)
c
      uind = uind + t4(1) * m0(:,:,1)
      uinp = uinp + t4(2) * m0(:,:,2)
      if (order .gt. 1 .or. tcgpeek) then
c
c     r1 = r0 - gamma0 * T.p0 (or T.m0)
c
         r1(:,:,1) = r0(:,:,1) - t4(1)*P1(:,:,1)
         r1(:,:,2) = r0(:,:,2) - t4(2)*P1(:,:,2)
c
c     check convergence, stop at tcg1 level if n1 is small enough
c
         call tcg_dotprod (nrsd1(1),3*npole,r1(:,:,1),r1(:,:,1))
         call tcg_dotprod (nrsd1(2),3*npole,r1(:,:,2),r1(:,:,2))
         call tcg_converge (converge,nrsd1(1),nrsd1(2))
         if (converge) then
            order = 1
            call tcg_resource (order)
         end if
c
c     n1 = r1.M.r1
c
         call tcg_alphaquad (n1(1),r1(:,:,1),r1(:,:,1))
         call tcg_alphaquad (n1(2),r1(:,:,2),r1(:,:,2))
      end if
c
c     cross terms
c     sp0 = m0.E
c     spp1 = m0.T.alpha.E = m0.T.mu0
c
      call tcg_dotprod (sp0(1),3*npole,m0(:,:,1),field(:,:,2))
      call tcg_dotprod (sp0(2),3*npole,m0(:,:,2),field(:,:,1))
      if (tcgpeek) then
         call tcg_dotprod (spp1(1),3*npole,P1(:,:,1),udirp)
         call tcg_dotprod (spp1(2),3*npole,P1(:,:,2),udir)
      end if
c
c     tcg1 force and energy
c
      if (order .eq. 1) then
c
c     compute a(1) coefficients: a1...
c     and a(1k) coefficients: a1k...
c
         do i = 1, 2
            a110(i) = t4(i)
            a111(i) = 2.0d0 * sp0(i) / t1(i)
            a112(i) = -t4(i) * a111(i)
            a121(i) = 0.5d0 * a112(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a1k10a(i) = tcgomega
               a1k11a(i) = -tcgomega * t4(i)
               a1k11(i) = -2.0d0 * spp1(i) * tcgomega / t1(i)
               a1k12(i) = -t4(i) * a1k11(i)
               a1k20a(i) = a1k11a(i)
               a1k21(i) = 0.5d0 * a1k12(i)
            end do
c
c     mu1(peek) = mu1 + omega * alpha.r1
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = 0.5d0 * (-(a121(1)+a1k21(1))*m0(:,:,1)
     &                               - a1k20a(1)*udirp)
            ubd(:,:,1) = 0.5d0 * (-(a121(2)+a1k21(2))*m0(:,:,2)
     &                               - a1k20a(2)*udir)
            xdr0(:,:,1) = (a110(2)+a1k11a(2)+a1k10a(2))*field(:,:,1)
     &                  - a1k11a(2)*r0(:,:,1)
     &                  + (a111(2)+a1k11(2))*r0(:,:,2)
     &                  + (a112(2)+a1k12(2))*P1(:,:,2)
            xdr0(:,:,2) = (a110(1)+a1k11a(1)+a1k10a(1))*field(:,:,2)
     &                  - a1k11a(1)*r0(:,:,2)
     &                  + (a111(1)+a1k11(1))*r0(:,:,1)
     &                  + (a112(1)+a1k12(1))*P1(:,:,1)
         else
            ubp(:,:,1) = 0.5d0 * (-a121(1)*m0(:,:,1))
            ubd(:,:,1) = 0.5d0 * (-a121(2)*m0(:,:,2))
            xdr0(:,:,1) = a110(2)*field(:,:,1) + a111(2)*r0(:,:,2)
     &                  + a112(2)*P1(:,:,2)
            xdr0(:,:,2) = a110(1)*field(:,:,2) + a111(1)*r0(:,:,1)
     &                  + a112(1)*P1(:,:,1)
         end if
         call tcg_alpha12 (xdr0(:,:,1),xdr0(:,:,2))
c
c     xde: rhs array in <E' xde>
c     xde = mu + mu0 + alpha.(-Tu.x)
c
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),xde(:,:,1),xde(:,:,2))
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind + udir)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp + udirp)
         ubp(:,:,2) = 0.5d0 * xdr0(:,:,2)
         ubd(:,:,2) = 0.5d0 * xdr0(:,:,1)
         uad(:,:,1) = m0(:,:,1)
         uap(:,:,1) = m0(:,:,2)
         uad(:,:,2) = udir
         uap(:,:,2) = udirp
         goto 10
      end if
c
c     compute tcg2 intermediates, use "xde" as temporary storage
c     t2m0 = T.M.T.M.r0 = T.M.P1
c     t3m0 = T.(M.T)^2.M.r0 = T.M.t2m0
c     t9 = r0.M.T.M.T.M.T.M.r0 = t2m0.M.P1
c
      call tcg_alpha22 (P1(:,:,1),P1(:,:,2),xde(:,:,1),xde(:,:,2))
      call tcg_t0 (xde(:,:,1),xde(:,:,2),t2m0(:,:,1),t2m0(:,:,2))
      call tcg_alpha22 (t2m0(:,:,1),t2m0(:,:,2),xde(:,:,1),xde(:,:,2))
      call tcg_t0 (xde(:,:,1),xde(:,:,2),t3m0(:,:,1),t3m0(:,:,2))
      call tcg_alphaquad (t9(1),t2m0(:,:,1),P1(:,:,1))
      call tcg_alphaquad (t9(2),t2m0(:,:,2),P1(:,:,2))
c
c     beta1 = r1.r1/r0.r0 = n1/n0
c     t2 = 1 + beta1
c
      beta1(1) = n1(1) / n0(1)
      beta1(2) = n1(2) / n0(2)
      t2(1) = 1.0d0 + beta1(1)
      t2(2) = 1.0d0 + beta1(2)
c
c     np1 = P1.M.P1
c     t8 = t2*np1 - t4*t9
c     t10 = t1^2 - n0.|P1|^2
c     t3  = t1*t8
c     gamma1 = t10/t3
c
      call tcg_alphaquad (np1(1),P1(:,:,1),P1(:,:,1))
      call tcg_alphaquad (np1(2),P1(:,:,2),P1(:,:,2))
      t8(1) = t2(1)*np1(1) - t4(1)*t9(1)
      t8(2) = t2(2)*np1(2) - t4(2)*t9(2)
      t10(1) = t1(1)**2 - n0(1)*np1(1)
      t10(2) = t1(2)**2 - n0(2)*np1(2)
      t3(1) = t1(1)*t8(1)
      t3(2) = t1(2)*t8(2)
      gamma1(1) = t10(1) / t3(1)
      gamma1(2) = t10(2) / t3(2)
c
c     mu2 = mu1 + gamma1*p1 = mu1 + gamma1*(M.r1 + beta1*p0)
c         = mu1 + gamma1*(t2*p0 - t4*M.T.p0)
c         = mu1 + gamma1*(t2*M.r0 - t4*M.T.M.r0)
c         = mu1 + gamma1*M.(t2*r0 - t4*P1)
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = uind(j,i) + (t2(1)*r0(j,i,1) - t4(1)*P1(j,i,1))
     &                     * gamma1(1) * polarity(i)
            uinp(j,i) = uinp(j,i) + (t2(2)*r0(j,i,2) - t4(2)*P1(j,i,2))
     &                     * gamma1(2) * polarity(i)
         end do
      end do
      if (order .gt. 2 .or. tcgpeek) then
c
c     r2 = r1 - gamma1 * T.p1 = r1 - gamma1 * P2
c        = r1 - gamma1 * (t2*T.M.r0 - t4*T.M.T.M.r0)
c        = r1 - gamma1 * (t2*P1 - t4*t2m0)
c     reuse r1 as r2
c
         r1(:,:,1) = r1(:,:,1)
     &                  - gamma1(1)*(t2(1)*P1(:,:,1)-t4(1)*t2m0(:,:,1))
         r1(:,:,2) = r1(:,:,2)
     &                  - gamma1(2)*(t2(2)*P1(:,:,2)-t4(2)*t2m0(:,:,2))
      end if
c
c     cross terms
c     sp1 = m0.T.M.E = P1.M.E
c     b1 = sp0 - gamma1*sp1
c     b2 = sp0*t2 - t4*sp1
c     spp2 = r0.M.T.M.T.alpha.E = P1.M.(E-r0)
c          = P1.M.E - P1.M.r0 = sp1 - P1.M.r0
c
      call tcg_alphaquad (sp1(1),P1(:,:,1),field(:,:,2))
      call tcg_alphaquad (sp1(2),P1(:,:,2),field(:,:,1))
      b1(1) = sp0(1) - gamma1(1)*sp1(1)
      b1(2) = sp0(2) - gamma1(2)*sp1(2)
      b2(1) = sp0(1)*t2(1) - t4(1)*sp1(1)
      b2(2) = sp0(2)*t2(2) - t4(2)*sp1(2)
      if (tcgpeek) then
         call tcg_alphaquad (spp2(1),P1(:,:,1),r0(:,:,2))
         call tcg_alphaquad (spp2(2),P1(:,:,2),r0(:,:,1))
         spp2(1) = sp1(1) - spp2(1)
         spp2(2) = sp1(2) - spp2(2)
      end if
c
c     tcg2 force and energy
c
      if (order .eq. 2) then
c
c     compute a(2) coefficients: a2...
c     and a(2k) coefficients: a2k...
c
         do i = 1, 2
            a232(i) = t1(i)*t4(i)*gamma1(i)*b2(i)/t3(i)
            a241(i) = a232(i)
            a231(i) = -n0(i)*b2(i)/t3(i)
     &                -2.0d0*t1(i)*t2(i)*gamma1(i)*b2(i)/t3(i)
     &                +t4(i)*gamma1(i)*sp0(i)/t1(i)
            a223(i) = a232(i)
            a222(i) = a231(i)
            a221(i) = -t4(i)*b1(i)/t1(i) +2.0d0*t1(i)*b2(i)/t3(i)
     &                -t4(i)*t9(i)*gamma1(i)*b2(i)/t3(i)
     &                +2.0d0*t2(i)*np1(i)*gamma1(i)*b2(i)/t3(i)
     &                -t8(i)*gamma1(i)*b2(i)/t3(i)
     &                -2.0d0*t4(i)*np1(i)*sp0(i)*gamma1(i)/(t1(i)**2)
            a220(i) = -gamma1(i)*t4(i)
            a214(i) = 2.0d0*a232(i)
            a213(i) = 2.0d0*a231(i)
            a212(i) = 2.0d0*a221(i)
            a211(i) = 2.0d0*(b1(i)/t1(i) -np1(i)*b2(i)/t3(i)
     &                   -(np1(i)**2)*gamma1(i)*b2(i)/t3(i)/t1(i)
     &                   +t9(i)*gamma1(i)*b2(i)/t3(i)
     &                   +np1(i)*sp0(i)*gamma1(i)/(t1(i)**2))
            a21n1(i) = a220(i)
            a210(i) = t4(i) + gamma1(i)*t2(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a2kwt2(i) = tcgomega*(t2(i)*spp1(i)-t4(i)*spp2(i))
               a2kwg1(i) = tcgomega*(spp1(i)-gamma1(i)*spp2(i))
               a2k41(i) = -a2kwt2(i)*t1(i)*t4(i)*gamma1(i)/t3(i)
               a2k32(i) = a2k41(i)
               a2k31(i) = -tcgomega*t4(i)*gamma1(i)*spp1(i)/t1(i)
     &                    +a2kwt2(i)*(n0(i)/t3(i)
     &                       +2.0d0*t1(i)*t2(i)*gamma1(i)/t3(i))
               a2k30a(i) = tcgomega*gamma1(i)*t4(i)
               a2k23(i) = a2k41(i)
               a2k22(i) = a2k31(i)
               a2k21(i) = 2.0d0*t4(i)*np1(i)/(t1(i)**2)*
     &                       tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(-2.0d0*t1(i)+
     &                       (t4(i)*t9(i)-2.0d0*np1(i)*t2(i)+t8(i))*
     &                       gamma1(i))/t3(i)
     &                    +t4(i)*a2kwg1(i)/t1(i)
               a2k21a(i) = a2k30a(i)
               a2k20a(i) = -tcgomega*(gamma1(i)*t2(i)+t4(i))
               a2k14(i) = 2.0d0*a2k41(i)
               a2k13(i) = 2.0d0*a2k22(i)
               a2k12(i) = 2.0d0*a2k21(i)
               a2k11(i) = -np1(i)/(t1(i)**2)*tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(np1(i)
     &                       +(np1(i)**2)*gamma1(i)/t1(i)
     &                       -t9(i)*gamma1(i))/t3(i)
     &                    -a2kwg1(i)/t1(i)
               a2k11(i) = 2.0d0*a2k11(i)
               a2k12a(i) = a2k30a(i)
               a2k11a(i) = a2k20a(i)
               a2k10a(i) = tcgomega
            end do
c
c     mu2(peek) = mu2 + omega * alpha.r2
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = 0.5d0 * (-(a220(1)+a2k20a(1)
     &                               +a2k21a(1))*field(:,:,2)
     &                               +a2k21a(1)*r0(:,:,2)
     &                               -(a221(1)+a2k21(1))*r0(:,:,1)
     &                               -(a222(1)+a231(1)+a2k22(1)
     &                                    +a2k31(1))*P1(:,:,1)
     &                               -(a223(1)+a241(1)+a2k23(1)
     &                                    +a2k41(1))*t2m0(:,:,1))
            ubd(:,:,1) = 0.5d0 * (-(a220(2)+a2k20a(2)
     &                               +a2k21a(2))*field(:,:,1)
     &                               +a2k21a(2)*r0(:,:,1)
     &                               -(a221(2)+a2k21(2))*r0(:,:,2)
     &                               -(a222(2)+a231(2)+a2k22(2)
     &                                    +a2k31(2))*P1(:,:,2)
     &                               -(a223(2)+a241(2)+a2k23(2)
     &                                    +a2k41(2))*t2m0(:,:,2))
            ubp(:,:,2) = 0.5d0 * (-(a232(1)+a2k32(1))*P1(:,:,1)
     &                               -a2k30a(1)*field(:,:,2))
            ubd(:,:,2) = 0.5d0 * (-(a232(2)+a2k32(2))*P1(:,:,2)
     &                               -a2k30a(2)*field(:,:,1))
            xdr0(:,:,1) = (a210(2)+a21n1(2)+a2k10a(2)
     &                        +a2k11a(2)+a2k12a(2))*field(:,:,1)
     &                  - (a21n1(2)+a2k11a(2)+a2k12a(2))*r0(:,:,1)
     &                  + (a211(2)+a2k11(2))*r0(:,:,2)
     &                  - a2k12a(2)*P1(:,:,1)
     &                  + (a212(2)+a2k12(2))*P1(:,:,2)
     &                  + (a213(2)+a2k13(2))*t2m0(:,:,2)
     &                  + (a214(2)+a2k14(2))*t3m0(:,:,2)
            xdr0(:,:,2) = (a210(1)+a21n1(1)+a2k10a(1)
     &                        +a2k11a(1)+a2k12a(1))*field(:,:,2)
     &                  - (a21n1(1)+a2k11a(1)+a2k12a(1))*r0(:,:,2)
     &                  + (a211(1)+a2k11(1))*r0(:,:,1)
     &                  - a2k12a(1)*P1(:,:,2)
     &                  + a212(1)*P1(:,:,1)+a2k12(1)*P1(:,:,1)
     &                  + (a213(1)+a2k13(1))*t2m0(:,:,1)
     &                  + (a214(1)+a2k14(1))*t3m0(:,:,1)
         else
            ubp(:,:,1) = 0.5d0 * (-a220(1)*field(:,:,2)
     &                               -a221(1)*r0(:,:,1)
     &                               -(a222(1)+a231(1))*P1(:,:,1)
     &                               -(a223(1)+a241(1))*t2m0(:,:,1))
            ubd(:,:,1) = 0.5d0 * (-a220(2)*field(:,:,1)
     &                               -a221(2)*r0(:,:,2)
     &                               -(a222(2)+a231(2))*P1(:,:,2)
     &                               -(a223(2)+a241(2))*t2m0(:,:,2))
            ubp(:,:,2) = 0.5d0 * (-a232(1)*P1(:,:,1))
            ubd(:,:,2) = 0.5d0 * (-a232(2)*P1(:,:,2))
            xdr0(:,:,1) = (a210(2)+a21n1(2))*field(:,:,1)
     &                  - a21n1(2)*r0(:,:,1)
     &                  + a211(2)*r0(:,:,2) + a212(2)*P1(:,:,2)
     &                  + a213(2)*t2m0(:,:,2) + a214(2)*t3m0(:,:,2)
            xdr0(:,:,2) = (a210(1)+a21n1(1))*field(:,:,2)
     &                  - a21n1(1)*r0(:,:,2)
     &                  + a211(1)*r0(:,:,1) + a212(1)*P1(:,:,1)
     &                  + a213(1)*t2m0(:,:,1) + a214(1)*t3m0(:,:,1)
         end if
         call tcg_alpha12 (xdr0(:,:,1),xdr0(:,:,2))
c
c     xde = mu + mu0 + alpha.(-Tu.x)
c
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),xde(:,:,1),xde(:,:,2))
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind + udir)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp + udirp)
         call tcg_alpha12 (ubp(:,:,1),ubd(:,:,1))
         call tcg_alpha12 (ubp(:,:,2),ubd(:,:,2))
         ubp(:,:,3) = 0.5d0 * xdr0(:,:,2)
         ubd(:,:,3) = 0.5d0 * xdr0(:,:,1)
         uad(:,:,1) = m0(:,:,1)
         uap(:,:,1) = m0(:,:,2)
         call tcg_alpha22 (P1(:,:,1),P1(:,:,2),uad(:,:,2),uap(:,:,2))
         uad(:,:,3) = udir
         uap(:,:,3) = udirp
         goto 10
      end if
c
c     store induced dipoles and use uind/p to store xde arrays
c
   10 continue
      call tcg_store
      uind = xde(:,:,1)
      uinp = xde(:,:,2)
c
c     perform deallocation for some local arrays
c
      deallocate (field)
      deallocate (xde)
      deallocate (xdr0)
      deallocate (r0)
      deallocate (m0)
      deallocate (P1)
      deallocate (r1)
      deallocate (t2m0)
      deallocate (t3m0)
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_induce1b  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_induce1b" computes the induced dipoles and intermediates
c     used in polarization force calculation for the TCG method with
c     initial guess mu0 = direct and diagonal preconditioner = false
c
c
      subroutine tcg_induce1b
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      implicit none
      integer i,order
      real*8 n0(2),t1(2),t4(2),n1(2)
      real*8 sp0(2),spp1(2)
      real*8 a110(2),a111(2),a112(2),a121(2)
      real*8 a1k10a(2),a1k11a(2),a1k11(2)
      real*8 a1k12(2),a1k20a(2),a1k21(2)
      real*8 t9(2),beta1(2),t2(2),np1(2)
      real*8 t8(2),t10(2),t3(2),gamma1(2)
      real*8 sp1(2),b1(2),b2(2),spp2(2)
      real*8 a210(2),a21n1(2),a211(2)
      real*8 a212(2),a213(2),a214(2)
      real*8 a220(2),a221(2),a222(2),a223(2)
      real*8 a231(2),a232(2),a241(2)
      real*8 a2k10a(2),a2k11a(2),a2k12a(2)
      real*8 a2k11(2),a2k12(2),a2k13(2),a2k14(2)
      real*8 a2k20a(2),a2k21a(2),a2k21(2)
      real*8 a2k22(2),a2k23(2)
      real*8 a2k30a(2),a2k31(2),a2k32(2),a2k41(2)
      real*8 a2kwt2(2),a2kwg1(2)
      real*8, allocatable :: field(:,:,:)
      real*8, allocatable :: xde(:,:,:)
      real*8, allocatable :: xdr0(:,:,:)
      real*8, allocatable :: r0(:,:,:)
      real*8, allocatable :: P1(:,:,:)
      real*8, allocatable :: r1(:,:,:)
      real*8, allocatable :: t2r0(:,:,:)
      real*8, allocatable :: t3r0(:,:,:)
      real*8, allocatable :: te(:,:,:)
      logical converge
c
c
c     set up nab based on tcgorder
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some local arrays
c
      if (.not. allocated(uindt))  then
         allocate (uindt(3,n))
         allocate (uinpt(3,n))
         allocate (uad(3,n,tcgnab))
         allocate (uap(3,n,tcgnab))
         allocate (ubd(3,n,tcgnab))
         allocate (ubp(3,n,tcgnab))
      end if
      allocate (field(3,n,2))
      allocate (xde(3,n,2))
      allocate (xdr0(3,n,2))
      allocate (r0(3,n,2))
      allocate (P1(3,n,2))
      allocate (r1(3,n,2))
      allocate (t2r0(3,n,2))
      allocate (t3r0(3,n,2))
      allocate (te(3,n,2))
c
c     get the electrostatic field due to permanent multipoles
c
      if (use_ewald) then
         call dfield0c (field(:,:,1),field(:,:,2))
      else if (use_mlist) then
         call dfield0b (field(:,:,1),field(:,:,2))
      else
         call dfield0a (field(:,:,1),field(:,:,2))
      end if
c
c     mu0 = alpha.E
c
      call tcg_alpha22 (field(:,:,1),field(:,:,2),udir,udirp)
      uind = udir
      uinp = udirp
c
c     r0 = -Tu.mu0 = mutual field of mu0
c
      if (use_ewald) then
         call ufield0c (r0(:,:,1),r0(:,:,2))
      else if (use_mlist) then
         call ufield0b (r0(:,:,1),r0(:,:,2))
      else
         call ufield0a (r0(:,:,1),r0(:,:,2))
      end if
c
c     n0 = r0.r0
c
      call tcg_dotprod (n0(1),3*npole,r0(:,:,1),r0(:,:,1))
      call tcg_dotprod (n0(2),3*npole,r0(:,:,2),r0(:,:,2))
c
c     compute tcg1 intermediates
c     (P1 or t1r0) = T.r0
c     (t1 or rtr0) = r0.T.r0
c     (t4 or gamma0) = r0.r0/r0.T.r0
c
      call tcg_t0 (r0(:,:,1),r0(:,:,2),P1(:,:,1),P1(:,:,2))
      call tcg_dotprod (t1(1),3*npole,r0(:,:,1),P1(:,:,1))
      call tcg_dotprod (t1(2),3*npole,r0(:,:,2),P1(:,:,2))
      t4(1) = n0(1) / t1(1)
      t4(2) = n0(2) / t1(2)
c
c     mu1 = mu0 + gamma0 * p0 (or r0)
c
      uind = uind + t4(1) * r0(:,:,1)
      uinp = uinp + t4(2) * r0(:,:,2)
      if (order .gt. 1 .or. tcgpeek) then
c
c     r1 = r0 - gamma0 * T.p0 (or T.r0)
c
         r1(:,:,1) = r0(:,:,1) - t4(1)*P1(:,:,1)
         r1(:,:,2) = r0(:,:,2) - t4(2)*P1(:,:,2)
c
c     n1 = r1.r1
c     check convergence, stop at tcg1 level if n1 is small enough
c
         call tcg_dotprod (n1(1),3*npole,r1(:,:,1),r1(:,:,1))
         call tcg_dotprod (n1(2),3*npole,r1(:,:,2),r1(:,:,2))
         call tcg_converge (converge,n1(1),n1(2))
         if (converge) then
            order = 1
            call tcg_resource (order)
         end if
      end if
c
c     cross terms
c     sp0 = r0.E
c     spp1 = r0.T.alpha.E = r0.T.mu0
c
      call tcg_dotprod (sp0(1),3*npole,r0(:,:,1),field(:,:,2))
      call tcg_dotprod (sp0(2),3*npole,r0(:,:,2),field(:,:,1))
      if (tcgpeek) then
         call tcg_dotprod (spp1(1),3*npole,P1(:,:,1),udirp)
         call tcg_dotprod (spp1(2),3*npole,P1(:,:,2),udir)
      end if
c
c     tcg1 force and energy
c
      if (order .eq. 1) then
c
c     compute a(1) coefficients: a1...
c     and a(1k) coefficients: a1k...
c
         do i = 1, 2
            a110(i) = t4(i)
            a111(i) = 2.0d0 * sp0(i) / t1(i)
            a112(i) = -t4(i) * a111(i)
            a121(i) = 0.5d0 * a112(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a1k10a(i) = tcgomega
               a1k11a(i) = -tcgomega * t4(i)
               a1k11(i) = -2.0d0 * spp1(i) * tcgomega / t1(i)
               a1k12(i) = -t4(i) * a1k11(i)
               a1k20a(i) = a1k11a(i)
               a1k21(i) = 0.5d0 * a1k12(i)
            end do
c
c     mu1(peek) = mu1 + omega * alpha.r1
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = 0.5d0 * (-(a121(1)+a1k21(1))*r0(:,:,1)
     &                               - a1k20a(1)*udirp)
            ubd(:,:,1) = 0.5d0 * (-(a121(2)+a1k21(2))*r0(:,:,2)
     &                               - a1k20a(2)*udir)
            xdr0(:,:,1) = (a110(2)+a1k11a(2))*field(:,:,1)
     &                  + a1k10a(2)*udir - a1k11a(2)*r0(:,:,1)
     &                  + (a111(2)+a1k11(2))*r0(:,:,2)
     &                  + (a112(2)+a1k12(2))*P1(:,:,2)
            xdr0(:,:,2) = (a110(1)+a1k11a(1))*field(:,:,2)
     &                  + a1k10a(1)*udirp - a1k11a(1)*r0(:,:,2)
     &                  + (a111(1)+a1k11(1))*r0(:,:,1)
     &                  + (a112(1)+a1k12(1))*P1(:,:,1)
         else
            ubp(:,:,1) = 0.5d0 * (-a121(1)*r0(:,:,1))
            ubd(:,:,1) = 0.5d0 * (-a121(2)*r0(:,:,2))
            xdr0(:,:,1) = a110(2)*field(:,:,1) + a111(2)*r0(:,:,2)
     &                  + a112(2)*P1(:,:,2)
            xdr0(:,:,2) = a110(1)*field(:,:,2) + a111(1)*r0(:,:,1)
     &                  + a112(1)*P1(:,:,1)
         end if
c
c     xde: rhs array in <E' xde>
c     xde = mu + mu0 + alpha.(-Tu.x)
c
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),xde(:,:,1),xde(:,:,2))
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind + udir)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp + udirp)
         ubp(:,:,2) = 0.5d0 * xdr0(:,:,2)
         ubd(:,:,2) = 0.5d0 * xdr0(:,:,1)
         uad(:,:,1) = r0(:,:,1)
         uad(:,:,2) = udir
         uap(:,:,1) = r0(:,:,2)
         uap(:,:,2) = udirp
         goto 10
      end if
c
c     compute tcg2 intermediates
c     t2r0 = T^2.r0 = T.P1
c     t3r0 = T^3.r0
c     te = T.E
c     t9 = r0.T^3.r0 = r0.T^2.T.r0
c
      call tcg_t0 (P1(:,:,1),P1(:,:,2),t2r0(:,:,1),t2r0(:,:,2))
      call tcg_t0 (t2r0(:,:,1),t2r0(:,:,2),t3r0(:,:,1),t3r0(:,:,2))
      call tcg_t0 (field(:,:,1),field(:,:,2),te(:,:,1),te(:,:,2))
      call tcg_dotprod (t9(1),3*npole,t2r0(:,:,1),P1(:,:,1))
      call tcg_dotprod (t9(2),3*npole,t2r0(:,:,2),P1(:,:,2))
c
c     beta1 = r1.r1/r0.r0 = n1/n0
c     t2 = 1 + beta1
c
      beta1(1) = n1(1) / n0(1)
      beta1(2) = n1(2) / n0(2)
      t2(1) = 1.0d0 + beta1(1)
      t2(2) = 1.0d0 + beta1(2)
c
c     np1 = P1.P1
c     t8 = t5 = P1.T.p1 = P1.P2 = t2*np1 - t4*t9
c     t10 = t1^2 - n0.|P1|^2
c     t3 = t1*P1.P2 = t1*t8
c     gamma1 = t10/t3
c
      call tcg_dotprod (np1(1),3*npole,P1(:,:,1),P1(:,:,1))
      call tcg_dotprod (np1(2),3*npole,P1(:,:,2),P1(:,:,2))
      t8(1) = t2(1)*np1(1) - t4(1)*t9(1)
      t8(2) = t2(2)*np1(2) - t4(2)*t9(2)
      t10(1) = t1(1)**2 - n0(1)*np1(1)
      t10(2) = t1(2)**2 - n0(2)*np1(2)
      t3(1) = t1(1)*t8(1)
      t3(2) = t1(2)*t8(2)
      gamma1(1) = t10(1) / t3(1)
      gamma1(2) = t10(2) / t3(2)
c
c     mu2 = mu1 + gamma1*p1 = mu1 + gamma1*(r1 + beta1*p0)
c         = mu1 + gamma1*(r1 + beta1*r0)
c
      uind = uind + gamma1(1)*(r1(:,:,1) + beta1(1)*r0(:,:,1))
      uinp = uinp + gamma1(2)*(r1(:,:,2) + beta1(2)*r0(:,:,2))
      if (order .gt. 2 .or. tcgpeek) then
c
c     r2 = r1 - gamma1 * T.p1 = r1 - gamma1 * P2
c        = r1 - gamma1 * (t2*P1 - t4*t2r0)
c     reuse r1 as r2
c
         r1(:,:,1) = r1(:,:,1)
     &                  - gamma1(1)*(t2(1)*P1(:,:,1)-t4(1)*t2r0(:,:,1))
         r1(:,:,2) = r1(:,:,2)
     &                  - gamma1(2)*(t2(2)*P1(:,:,2)-t4(2)*t2r0(:,:,2))
      end if
c
c     cross terms
c     sp1 = r0.T.E = P1.E
c     b1 = sp0 - gamma1*sp1
c     b2 = sp0*t2 - t4*sp1
c     spp2 = r0.T.T.alpha.E
c
      call tcg_dotprod (sp1(1),3*npole,P1(:,:,1),field(:,:,2))
      call tcg_dotprod (sp1(2),3*npole,P1(:,:,2),field(:,:,1))
      b1(1) = sp0(1) - gamma1(1)*sp1(1)
      b1(2) = sp0(2) - gamma1(2)*sp1(2)
      b2(1) = sp0(1)*t2(1) - t4(1)*sp1(1)
      b2(2) = sp0(2)*t2(2) - t4(2)*sp1(2)
      if (tcgpeek) then
         call tcg_dotprod (spp2(1),3*npole,t2r0(:,:,1),udirp)
         call tcg_dotprod (spp2(2),3*npole,t2r0(:,:,2),udir)
      end if
c
c     tcg2 force and energy
c
      if (order .eq. 2) then
c
c     compute a(2) coefficients: a2...
c     and a(2k) coefficients: a2k...
c
         do i = 1, 2
            a232(i) = t1(i)*t4(i)*gamma1(i)*b2(i)/t3(i)
            a241(i) = a232(i)
            a231(i) = -n0(i)*b2(i)/t3(i)
     &                -2.0d0*t1(i)*t2(i)*gamma1(i)*b2(i)/t3(i)
     &                +t4(i)*gamma1(i)*sp0(i)/t1(i)
            a223(i) = a232(i)
            a222(i) = a231(i)
            a221(i) = -t4(i)*b1(i)/t1(i) +2.0d0*t1(i)*b2(i)/t3(i)
     &                -t4(i)*t9(i)*gamma1(i)*b2(i)/t3(i)
     &                +2.0d0*t2(i)*np1(i)*gamma1(i)*b2(i)/t3(i)
     &                -t8(i)*gamma1(i)*b2(i)/t3(i)
     &                -2.0d0*t4(i)*np1(i)*sp0(i)*gamma1(i)/(t1(i)**2)
            a220(i) = -gamma1(i)*t4(i)
            a214(i) = 2.0d0*a232(i)
            a213(i) = 2.0d0*a231(i)
            a212(i) = 2.0d0*a221(i)
            a211(i) = 2.0d0*(b1(i)/t1(i) -np1(i)*b2(i)/t3(i)
     &                   -(np1(i)**2)*gamma1(i)*b2(i)/t3(i)/t1(i)
     &                   +t9(i)*gamma1(i)*b2(i)/t3(i)
     &                   +np1(i)*sp0(i)*gamma1(i)/(t1(i)**2))
            a21n1(i) = a220(i)
            a210(i) = t4(i) + gamma1(i)*t2(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a2kwt2(i) = tcgomega*(t2(i)*spp1(i)-t4(i)*spp2(i))
               a2kwg1(i) = tcgomega*(spp1(i)-gamma1(i)*spp2(i))
               a2k41(i) = -a2kwt2(i)*t1(i)*t4(i)*gamma1(i)/t3(i)
               a2k32(i) = a2k41(i)
               a2k31(i) = -tcgomega*t4(i)*gamma1(i)*spp1(i)/t1(i)
     &                    +a2kwt2(i)*(n0(i)/t3(i)
     &                       +2.0d0*t1(i)*t2(i)*gamma1(i)/t3(i))
               a2k30a(i) = tcgomega*gamma1(i)*t4(i)
               a2k23(i) = a2k41(i)
               a2k22(i) = a2k31(i)
               a2k21(i) = 2.0d0*t4(i)*np1(i)/(t1(i)**2)*
     &                       tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(-2.0d0*t1(i)+
     &                       (t4(i)*t9(i)-2.0d0*np1(i)*t2(i)+t8(i))*
     &                       gamma1(i))/t3(i)
     &                    +t4(i)*a2kwg1(i)/t1(i)
               a2k21a(i) = a2k30a(i)
               a2k20a(i) = -tcgomega*(gamma1(i)*t2(i)+t4(i))
               a2k14(i) = 2.0d0*a2k41(i)
               a2k13(i) = 2.0d0*a2k22(i)
               a2k12(i) = 2.0d0*a2k21(i)
               a2k11(i) = -np1(i)/(t1(i)**2)*tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(np1(i)
     &                       +(np1(i)**2)*gamma1(i)/t1(i)
     &                       -t9(i)*gamma1(i))/t3(i)
     &                    -a2kwg1(i)/t1(i)
               a2k11(i) = 2.0d0*a2k11(i)
               a2k12a(i) = a2k30a(i)
               a2k11a(i) = a2k20a(i)
               a2k10a(i) = tcgomega
            end do
c
c     mu2(peek) = mu2 + omega * alpha.r2
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = 0.5d0 * (-(a220(1)+a2k21a(1))*field(:,:,2)
     &                   -(a221(1)+a2k21(1))*r0(:,:,1)
     &                   +a2k21a(1)*r0(:,:,2) -a2k20a(1)*udirp
     &                   -(a222(1)+a231(1)+a2k22(1)+a2k31(1))*P1(:,:,1)
     &                   -(a223(1)+a241(1)
     &                        +a2k23(1)+a2k41(1))*t2r0(:,:,1))
            ubp(:,:,2) = 0.5d0 * (-(a232(1)+a2k32(1))*P1(:,:,1)
     &                               -a2k30a(1)*udirp)
            ubd(:,:,1) = 0.5d0 * (-(a220(2)+a2k21a(2))*field(:,:,1)
     &                   -(a221(2)+a2k21(2))*r0(:,:,2)
     &                   +a2k21a(2)*r0(:,:,1) -a2k20a(2)*udir
     &                   -(a222(2)+a231(2)+a2k22(2)+a2k31(2))*P1(:,:,2)
     &                   -(a223(2)+a241(2)
     &                        +a2k23(2)+a2k41(2))*t2r0(:,:,2))
            ubd(:,:,2) = 0.5d0 * (-(a232(2)+a2k32(2))*P1(:,:,2)
     &                               -a2k30a(2)*udir)
            xdr0(:,:,1) = (a210(2)+a2k11a(2))*field(:,:,1)
     &                  + (a21n1(2)+a2k12a(2))*te(:,:,1)
     &                  + a2k10a(2)*udir - a2k11a(2)*r0(:,:,1)
     &                  - a2k12a(2)*P1(:,:,1)
     &                  + (a211(2)+a2k11(2))*r0(:,:,2)
     &                  + (a212(2)+a2k12(2))*P1(:,:,2)
     &                  + (a213(2)+a2k13(2))*t2r0(:,:,2)
     &                  + (a214(2)+a2k14(2))*t3r0(:,:,2)
            xdr0(:,:,2) = (a210(1)+a2k11a(1))*field(:,:,2)
     &                  + (a21n1(1)+a2k12a(1))*te(:,:,2)
     &                  + a2k10a(1)*udirp - a2k11a(1)*r0(:,:,2)
     &                  - a2k12a(1)*P1(:,:,2)
     &                  + (a211(1)+a2k11(1))*r0(:,:,1)
     &                  + (a212(1)+a2k12(1))*P1(:,:,1)
     &                  + (a213(1)+a2k13(1))*t2r0(:,:,1)
     &                  + (a214(1)+a2k14(1))*t3r0(:,:,1)
         else
            ubp(:,:,1) = 0.5d0 * (-a220(1)*field(:,:,2)
     &                               -a221(1)*r0(:,:,1)
     &                               -(a222(1)+a231(1))*P1(:,:,1)
     &                               -(a223(1)+a241(1))*t2r0(:,:,1))
            ubp(:,:,2) = 0.5d0 * (-a232(1)*P1(:,:,1))
            ubd(:,:,1) = 0.5d0 * (-a220(2)*field(:,:,1)
     &                               -a221(2)*r0(:,:,2)
     &                               -(a222(2)+a231(2))*P1(:,:,2)
     &                               -(a223(2)+a241(2))*t2r0(:,:,2))
            ubd(:,:,2) = 0.5d0 * (-a232(2)*P1(:,:,2))
            xdr0(:,:,1) = a210(2)*field(:,:,1) + a21n1(2)*te(:,:,1)
     &                  + a211(2)*r0(:,:,2) + a212(2)*P1(:,:,2)
     &                  + a213(2)*t2r0(:,:,2) + a214(2)*t3r0(:,:,2)
            xdr0(:,:,2) = a210(1)*field(:,:,2) + a21n1(1)*te(:,:,2)
     &                  + a211(1)*r0(:,:,1) + a212(1)*P1(:,:,1)
     &                  + a213(1)*t2r0(:,:,1) + a214(1)*t3r0(:,:,1)
         end if
c
c     xde = mu + mu0 + alpha.(-Tu.x)
c
         call tcg_ufield (xdr0(:,:,1),xdr0(:,:,2),xde(:,:,1),xde(:,:,2))
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind + udir)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp + udirp)
         ubp(:,:,3) = 0.5d0 * xdr0(:,:,2)
         ubd(:,:,3) = 0.5d0 * xdr0(:,:,1)
         uad(:,:,1) = r0(:,:,1)
         uad(:,:,2) = P1(:,:,1)
         uad(:,:,3) = udir 
         uap(:,:,1) = r0(:,:,2)
         uap(:,:,2) = P1(:,:,2)
         uap(:,:,3) = udirp
         goto 10
      end if
c
c     store induced dipoles and use uind/p to store xde arrays
c
   10 continue
      call tcg_store
      uind = xde(:,:,1)
      uinp = xde(:,:,2)
c
c     perform deallocation for some local arrays
c
      deallocate (field)
      deallocate (xde)
      deallocate (xdr0)
      deallocate (r0)
      deallocate (P1)
      deallocate (r1)
      deallocate (t2r0)
      deallocate (t3r0)
      deallocate (te)
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_induce1c  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_induce1c" computes the induced dipoles and intermediates
c     used in polarization force calculation for the TCG method with
c     initial guess mu0 = 0 and diagonal preconditioner = true
c
c
      subroutine tcg_induce1c
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      implicit none
      integer i,j
      integer order
      real*8 n0(2),t1(2),t4(2),n1(2),nrsd1(2)
      real*8 sp0,spp1
      real*8 a110(2),a111(2),a112(2),a121(2)
      real*8 a1k10a(2),a1k11a(2),a1k11(2)
      real*8 a1k12(2),a1k20a(2),a1k21(2)
      real*8 t9(2),beta1(2),t2(2),np1(2)
      real*8 t8(2),t10(2),t3(2),gamma1(2)
      real*8 sp1,b1(2),b2(2),spp2
      real*8 a210(2),a21n1(2),a211(2)
      real*8 a212(2),a213(2),a214(2)
      real*8 a220(2),a221(2),a222(2),a223(2)
      real*8 a231(2),a232(2),a241(2)
      real*8 a2k10a(2),a2k11a(2),a2k12a(2)
      real*8 a2k11(2),a2k12(2),a2k13(2),a2k14(2)
      real*8 a2k20a(2),a2k21a(2),a2k21(2)
      real*8 a2k22(2),a2k23(2)
      real*8 a2k30a(2),a2k31(2),a2k32(2),a2k41(2)
      real*8 a2kwt2(2),a2kwg1(2)
      real*8, allocatable :: r0(:,:,:)
      real*8, allocatable :: xde(:,:,:)
      real*8, allocatable :: P1(:,:,:)
      real*8, allocatable :: r1(:,:,:)
      real*8, allocatable :: t2m0(:,:,:)
      real*8, allocatable :: t3m0(:,:,:)
      logical converge
c
c
c     set up nab based on tcgorder
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some local arrays
c
      if (.not. allocated(uindt)) then
         allocate (uindt(3,n))
         allocate (uinpt(3,n))
         allocate (uad(3,n,tcgnab))
         allocate (uap(3,n,tcgnab))
         allocate (ubd(3,n,tcgnab))
         allocate (ubp(3,n,tcgnab))
      end if
      allocate (r0(3,n,2))
      allocate (xde(3,n,2))
      allocate (P1(3,n,2))
      allocate (r1(3,n,2))
      allocate (t2m0(3,n,2))
      allocate (t3m0(3,n,2))
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
c     m0 = M.r0 = alpha.E = udir
c
      call tcg_alpha22 (r0(:,:,1),r0(:,:,2),udir,udirp)
c
c     compute tcg1 intermediates
c     n0 = r0.M.r0 = E.alpha.E = r0.udir
c     P1 = T.m0 = tae
c     t1 = m0.T.m0 = P1.m0
c
      call tcg_dotprod (n0(1),3*npole,r0(:,:,1),udir)
      call tcg_dotprod (n0(2),3*npole,r0(:,:,2),udirp)
      call tcg_t0 (udir,udirp,P1(:,:,1),P1(:,:,2))
      call tcg_dotprod (t1(1),3*npole,P1(:,:,1),udir)
      call tcg_dotprod (t1(2),3*npole,P1(:,:,2),udirp)
      t4(1) = n0(1) / t1(1)
      t4(2) = n0(2) / t1(2)
c
c     mu1 = mu0 + gamma0 * p0 (or m0)
c
      uind = uind + t4(1) * udir
      uinp = uinp + t4(2) * udirp
      if (order .gt. 1 .or. tcgpeek) then
c
c     r1 = r0 - gamma0 * T.p0 (or T.m0)
c
         r1(:,:,1) = r0(:,:,1) - t4(1)*P1(:,:,1)
         r1(:,:,2) = r0(:,:,2) - t4(2)*P1(:,:,2)
c
c     check convergence, stop at tcg1 level if n1 is small enough
c
         call tcg_dotprod (nrsd1(1),3*npole,r1(:,:,1),r1(:,:,1))
         call tcg_dotprod (nrsd1(2),3*npole,r1(:,:,2),r1(:,:,2))
         call tcg_converge (converge,nrsd1(1),nrsd1(2))
         if (converge) then
            order = 1
            call tcg_resource (order)
         end if
c
c     n1 = r1.M.r1
c
         call tcg_alphaquad (n1(1),r1(:,:,1),r1(:,:,1))
         call tcg_alphaquad (n1(2),r1(:,:,2),r1(:,:,2))
      end if
c
c     cross terms
c     sp0 = r0.M.E = r0.udir
c     spp1 = m0.T.alpha.E = sp1, so spp1 is calculated
c     outside of if (tcgpeek) block
c
      call tcg_dotprod (sp0,3*npole,r0(:,:,1),udirp)
      call tcg_dotprod (spp1,3*npole,P1(:,:,1),udirp)
c
c     tcg1 force and energy
c
      if (order .eq. 1) then
c
c     compute a(1) coefficients: a1...
c     and a(1k) coefficients: a1k...
c
         do i = 1, 2
            a110(i) = t4(i)
            a111(i) = 2.0d0 * sp0 / t1(i)
            a112(i) = -t4(i) * a111(i)
            a121(i) = 0.5d0 * a112(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a1k10a(i) = tcgomega
               a1k11a(i) = -tcgomega * t4(i)
               a1k11(i) = -2.0d0 * spp1 * tcgomega / t1(i)
               a1k12(i) = -t4(i) * a1k11(i)
               a1k20a(i) = a1k11a(i)
               a1k21(i) = 0.5d0 * a1k12(i)
            end do
c
c     mu1(peek) = mu1 + omega * alpha.r1
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = -0.5d0*((a121(1)+a1k21(1))*udir
     &                              + a1k20a(1)*udirp)
            ubd(:,:,1) = -0.5d0*((a121(2)+a1k21(2))*udirp
     &                              + a1k20a(2)*udir)
            xde(:,:,1) = (a110(2)+a1k10a(2))*r0(:,:,1)
     &                 + (a111(2)+a1k11(2))*r0(:,:,2)
     &                 + (a112(2)+a1k12(2))*P1(:,:,2)
     &                 + a1k11a(2)*P1(:,:,1)
            xde(:,:,2) = (a110(1)+a1k10a(1))*r0(:,:,2)
     &                 + (a111(1)+a1k11(1))*r0(:,:,1)
     &                 + (a112(1)+a1k12(1))*P1(:,:,1)
     &                 + a1k11a(1)*P1(:,:,2)
         else
            ubp(:,:,1) = -0.5d0*a121(1)*udir
            ubd(:,:,1) = -0.5d0*a121(2)*udirp
            xde(:,:,1) = a110(2)*r0(:,:,1) + a111(2)*r0(:,:,2)
     &                      + a112(2)*P1(:,:,2)
            xde(:,:,2) = a110(1)*r0(:,:,2) + a111(1)*r0(:,:,1)
     &                      + a112(1)*P1(:,:,1)
         end if
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp)
         uad(:,:,1) = udir
         uap(:,:,1) = udirp
         goto 10
      end if
c
c     compute tcg2 intermediates, use "xde" as temporary storage
c     t2m0 = T.M.T.M.r0 = T.M.P1
c     t3m0 = T.(M.T)^2.M.r0 = T.M.t2m0
c     t9 = r0.M.T.M.T.M.T.M.r0 = t2m0.M.P1
c
      call tcg_alpha22 (P1(:,:,1),P1(:,:,2),xde(:,:,1),xde(:,:,2))
      call tcg_t0 (xde(:,:,1),xde(:,:,2),t2m0(:,:,1),t2m0(:,:,2))
      call tcg_alpha22 (t2m0(:,:,1),t2m0(:,:,2),xde(:,:,1),xde(:,:,2))
      call tcg_t0 (xde(:,:,1),xde(:,:,2),t3m0(:,:,1),t3m0(:,:,2))
      call tcg_alphaquad (t9(1),t2m0(:,:,1),P1(:,:,1))
      call tcg_alphaquad (t9(2),t2m0(:,:,2),P1(:,:,2))
c
c     beta1 = r1.r1/r0.r0 = n1/n0
c     t2 = 1 + beta1
c
      beta1(1) = n1(1) / n0(1)
      beta1(2) = n1(2) / n0(2)
      t2(1) = 1.0d0 + beta1(1)
      t2(2) = 1.0d0 + beta1(2)
c
c     np1 = P1.M.P1
c     t8 = t2*np1 - t4*t9
c     t10 = t1^2 - n0.|P1|^2
c     t3  = t1*t8
c     gamma1 = t10/t3
c
      call tcg_alphaquad (np1(1),P1(:,:,1),P1(:,:,1))
      call tcg_alphaquad (np1(2),P1(:,:,2),P1(:,:,2))
      t8(1) = t2(1)*np1(1) - t4(1)*t9(1)
      t8(2) = t2(2)*np1(2) - t4(2)*t9(2)
      t10(1) = t1(1)**2 - n0(1)*np1(1)
      t10(2) = t1(2)**2 - n0(2)*np1(2)
      t3(1) = t1(1)*t8(1)
      t3(2) = t1(2)*t8(2)
      gamma1(1) = t10(1) / t3(1)
      gamma1(2) = t10(2) / t3(2)
c
c     mu2 = mu1 + gamma1*p1 = mu1 + gamma1*(M.r1 + beta1*p0)
c         = mu1 + gamma1*(t2*p0 - t4*M.T.p0)
c         = mu1 + gamma1*(t2*M.r0 - t4*M.T.M.r0)
c         = mu1 + gamma1*M.(t2*r0 - t4*P1)
c
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = uind(j,i) + (t2(1)*r0(j,i,1) - t4(1)*P1(j,i,1))
     &                     * gamma1(1) * polarity(i)
            uinp(j,i) = uinp(j,i) + (t2(2)*r0(j,i,2) - t4(2)*P1(j,i,2))
     &                     * gamma1(2) * polarity(i)
         end do
      end do
      if (order .gt. 2 .or. tcgpeek) then
c
c     r2 = r1 - gamma1 * T.p1 = r1 - gamma1 * P2
c        = r1 - gamma1 * (t2*T.M.r0 - t4*T.M.T.M.r0)
c        = r1 - gamma1 * (t2*P1 - t4*t2m0)
c     reuse r1 as r2
c
         r1(:,:,1) = r1(:,:,1)
     &                  - gamma1(1)*(t2(1)*P1(:,:,1)-t4(1)*t2m0(:,:,1))
         r1(:,:,2) = r1(:,:,2)
     &                  - gamma1(2)*(t2(2)*P1(:,:,2)-t4(2)*t2m0(:,:,2))
      end if
c
c     cross terms
c     sp1 = P1.M.E = P1.udir = spp1
c     b1 = sp0 - gamma1*sp1
c     b2 = sp0*t2 - t4*sp1
c     spp2 = r0.M.T.M.T.alpha.E = P1.M.P1
c
      sp1 = spp1
      b1(1) = sp0 - gamma1(1)*sp1
      b1(2) = sp0 - gamma1(2)*sp1
      b2(1) = sp0*t2(1) - t4(1)*sp1
      b2(2) = sp0*t2(2) - t4(2)*sp1
      if (tcgpeek) then
         call tcg_alphaquad (spp2,P1(:,:,1),P1(:,:,2))
      end if
c
c     tcg2 force and energy
c
      if (order .eq. 2) then
c
c     compute a(2) coefficients: a2...
c     and a(2k) coefficients: a2k...
c
         do i = 1, 2
            a232(i) = t1(i)*t4(i)*gamma1(i)*b2(i)/t3(i)
            a241(i) = a232(i)
            a231(i) = -n0(i)*b2(i)/t3(i)
     &                -2.0d0*t1(i)*t2(i)*gamma1(i)*b2(i)/t3(i)
     &                +t4(i)*gamma1(i)*sp0/t1(i)
            a223(i) = a232(i)
            a222(i) = a231(i)
            a221(i) = -t4(i)*b1(i)/t1(i) +2.0d0*t1(i)*b2(i)/t3(i)
     &                -t4(i)*t9(i)*gamma1(i)*b2(i)/t3(i)
     &                +2.0d0*t2(i)*np1(i)*gamma1(i)*b2(i)/t3(i)
     &                -t8(i)*gamma1(i)*b2(i)/t3(i)
     &                -2.0d0*t4(i)*np1(i)*sp0*gamma1(i)/(t1(i)**2)
            a220(i) = -gamma1(i)*t4(i)
            a214(i) = 2.0d0*a232(i)
            a213(i) = 2.0d0*a231(i)
            a212(i) = 2.0d0*a221(i)
            a211(i) = 2.0d0*(b1(i)/t1(i) -np1(i)*b2(i)/t3(i)
     &                   -(np1(i)**2)*gamma1(i)*b2(i)/t3(i)/t1(i)
     &                   +t9(i)*gamma1(i)*b2(i)/t3(i)
     &                   +np1(i)*sp0*gamma1(i)/(t1(i)**2))
            a21n1(i) = a220(i)
            a210(i) = t4(i) + gamma1(i)*t2(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a2kwt2(i) = tcgomega*(t2(i)*spp1-t4(i)*spp2)
               a2kwg1(i) = tcgomega*(spp1-gamma1(i)*spp2)
               a2k41(i) = -a2kwt2(i)*t1(i)*t4(i)*gamma1(i)/t3(i)
               a2k32(i) = a2k41(i)
               a2k31(i) = -tcgomega*t4(i)*gamma1(i)*spp1/t1(i)
     &                    +a2kwt2(i)*(n0(i)/t3(i)
     &                       +2.0d0*t1(i)*t2(i)*gamma1(i)/t3(i))
               a2k30a(i) = tcgomega*gamma1(i)*t4(i)
               a2k23(i) = a2k41(i)
               a2k22(i) = a2k31(i)
               a2k21(i) = 2.0d0*t4(i)*np1(i)/(t1(i)**2)*
     &                       tcgomega*gamma1(i)*spp1
     &                    +a2kwt2(i)*(-2.0d0*t1(i)+
     &                       (t4(i)*t9(i)-2.0d0*np1(i)*t2(i)+t8(i))*
     &                       gamma1(i))/t3(i)
     &                    +t4(i)*a2kwg1(i)/t1(i)
               a2k21a(i) = a2k30a(i)
               a2k20a(i) = -tcgomega*(gamma1(i)*t2(i)+t4(i))
               a2k14(i) = 2.0d0*a2k41(i)
               a2k13(i) = 2.0d0*a2k22(i)
               a2k12(i) = 2.0d0*a2k21(i)
               a2k11(i) = -np1(i)/(t1(i)**2)*tcgomega*gamma1(i)*spp1
     &                    +a2kwt2(i)*(np1(i)
     &                       +(np1(i)**2)*gamma1(i)/t1(i)
     &                       -t9(i)*gamma1(i))/t3(i)
     &                    -a2kwg1(i)/t1(i)
               a2k11(i) = 2.0d0*a2k11(i)
               a2k12a(i) = a2k30a(i)
               a2k11a(i) = a2k20a(i)
               a2k10a(i) = tcgomega
            end do
c
c     mu2(peek) = mu2 + omega * alpha.r2
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = -0.5d0*((a220(1)+a2k20a(1))*r0(:,:,2)
     &                 + (a221(1)+a2k21(1))*r0(:,:,1)
     &                 + (a222(1)+a231(1)+a2k22(1)+a2k31(1))*P1(:,:,1)
     &                 + a2k21a(1)*P1(:,:,2)
     &                 + (a223(1)+a241(1)+a2k23(1)
     &                       +a2k41(1))*t2m0(:,:,1))
            ubd(:,:,1) = -0.5d0*((a220(2)+a2k20a(2))*r0(:,:,1)
     &                 + (a221(2)+a2k21(2))*r0(:,:,2)
     &                 + (a222(2)+a231(2)+a2k22(2)+a2k31(2))*P1(:,:,2)
     &                 + a2k21a(2)*P1(:,:,1)
     &                 + (a223(2)+a241(2)+a2k23(2)
     &                       +a2k41(2))*t2m0(:,:,2))
            ubp(:,:,2) = -0.5d0*((a232(1)+a2k32(1))*P1(:,:,1)
     &                              + a2k30a(1)*r0(:,:,2))
            ubd(:,:,2) = -0.5d0*((a232(2)+a2k32(2))*P1(:,:,2)
     &                              + a2k30a(2)*r0(:,:,1))
            xde(:,:,1) = (a210(2)+a2k10a(2))*r0(:,:,1)
     &                 + (a211(2)+a2k11(2))*r0(:,:,2)
     &                 + (a21n1(2)+a2k11a(2))*P1(:,:,1)
     &                 + (a212(2)+a2k12(2))*P1(:,:,2)
     &                 + a2k12a(2)*t2m0(:,:,1)
     &                 + (a213(2)+a2k13(2))*t2m0(:,:,2)
     &                 + (a214(2)+a2k14(2))*t3m0(:,:,2)
            xde(:,:,2) = (a210(1)+a2k10a(1))*r0(:,:,2)
     &                 + (a211(1)+a2k11(1))*r0(:,:,1)
     &                 + (a21n1(1)+a2k11a(1))*P1(:,:,2)
     &                 + (a212(1)+a2k12(1))*P1(:,:,1)
     &                 + a2k12a(1)*t2m0(:,:,2)
     &                 + (a213(1)+a2k13(1))*t2m0(:,:,1)
     &                 + (a214(1)+a2k14(1))*t3m0(:,:,1)
         else
            ubp(:,:,1) = -0.5d0*(a220(1)*r0(:,:,2) + a221(1)*r0(:,:,1)
     &                 + (a222(1)+a231(1))*P1(:,:,1)
     &                 + (a223(1)+a241(1))*t2m0(:,:,1))
            ubd(:,:,1) = -0.5d0*(a220(2)*r0(:,:,1) + a221(2)*r0(:,:,2)
     &                 + (a222(2)+a231(2))*P1(:,:,2)
     &                 + (a223(2)+a241(2))*t2m0(:,:,2))
            ubp(:,:,2) = -0.5d0*a232(1)*P1(:,:,1)
            ubd(:,:,2) = -0.5d0*a232(2)*P1(:,:,2)
            xde(:,:,1) = a210(2)*r0(:,:,1) + a211(2)*r0(:,:,2)
     &                 + a21n1(2)*P1(:,:,1) + a212(2)*P1(:,:,2)
     &                 + a213(2)*t2m0(:,:,2) + a214(2)*t3m0(:,:,2)
            xde(:,:,2) = a210(1)*r0(:,:,2) + a211(1)*r0(:,:,1)
     &                 + a21n1(1)*P1(:,:,2) + a212(1)*P1(:,:,1)
     &                 + a213(1)*t2m0(:,:,1) + a214(1)*t3m0(:,:,1)
         end if
         call tcg_alpha12 (xde(:,:,1),xde(:,:,2))
         xde(:,:,1) = 0.5d0 * (xde(:,:,1) + uind)
         xde(:,:,2) = 0.5d0 * (xde(:,:,2) + uinp)
         call tcg_alpha12 (ubp(:,:,1),ubd(:,:,1))
         call tcg_alpha12 (ubp(:,:,2),ubd(:,:,2))
         uad(:,:,1) = udir
         uap(:,:,1) = udirp
         call tcg_alpha22 (P1(:,:,1),P1(:,:,2),uad(:,:,2),uap(:,:,2))
         goto 10
      end if
c
c     store induced dipoles and use uind/p to store xde arrays
c
   10 continue
      call tcg_store
      uind = xde(:,:,1)
      uinp = xde(:,:,2)
c
c     perform deallocation for some local arrays
c
      deallocate (r0)
      deallocate (xde)
      deallocate (P1)
      deallocate (r1)
      deallocate (t2m0)
      deallocate (t3m0)
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_induce1d  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_induce1d" computes the induced dipoles and intermediates
c     used in polarization force calculation for the TCG method with
c     initial guess mu0 = 0 and diagonal preconditioner = false
c
c
      subroutine tcg_induce1d
      use atoms
      use limits
      use mpole
      use polar
      use poltcg
      implicit none
      integer i,order
      real*8 n0(2),t1(2),t4(2),n1(2)
      real*8 sp0,spp1(2)
      real*8 a110(2),a111(2),a112(2),a121(2)
      real*8 a1k10a(2),a1k11a(2),a1k11(2)
      real*8 a1k12(2),a1k20a(2),a1k21(2)
      real*8 t9(2),beta1(2),t2(2)
      real*8 np1(2),t8(2),t10(2),t3(2),gamma1(2)
      real*8 sp1,b1(2),b2(2),spp2(2)
      real*8 a210(2),a21n1(2),a211(2)
      real*8 a212(2),a213(2),a214(2)
      real*8 a220(2),a221(2),a222(2),a223(2)
      real*8 a231(2),a232(2),a241(2)
      real*8 a2k10a(2),a2k11a(2),a2k12a(2)
      real*8 a2k11(2),a2k12(2),a2k13(2),a2k14(2)
      real*8 a2k20a(2),a2k21a(2),a2k21(2)
      real*8 a2k22(2),a2k23(2)
      real*8 a2k30a(2),a2k31(2),a2k32(2),a2k41(2)
      real*8 a2kwt2(2),a2kwg1(2)
      real*8, allocatable :: r0(:,:,:)
      real*8, allocatable :: xde(:,:,:)
      real*8, allocatable :: P1(:,:,:)
      real*8, allocatable :: r1(:,:,:)
      real*8, allocatable :: tae(:,:,:)
      real*8, allocatable :: t2r0(:,:,:)
      real*8, allocatable :: t3r0(:,:,:)
      real*8, allocatable :: t2ae(:,:,:)
      logical converge
c
c
c     set up nab based on tcgorder
c
      order = tcgorder
      call tcg_resource (order)
c
c     perform dynamic allocation for some local arrays
c
      if (.not. allocated(uindt)) then
         allocate (uindt(3,n))
         allocate (uinpt(3,n))
         allocate (uad(3,n,tcgnab))
         allocate (uap(3,n,tcgnab))
         allocate (ubd(3,n,tcgnab))
         allocate (ubp(3,n,tcgnab))
      end if
      allocate (r0(3,n,2))
      allocate (xde(3,n,2))
      allocate (P1(3,n,2))
      allocate (r1(3,n,2))
      allocate (tae(3,n,2))
      allocate (t2r0(3,n,2))
      allocate (t3r0(3,n,2))
      allocate (t2ae(3,n,2))
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
c     compute tcg1 intermediates
c     n0 = r0.r0
c     P1 = T.r0 = T.E
c     t1 = r0.T.r0 = E.T.E
c     (t4 or gamma0) = r0.r0/r0.T.r0
c
      call tcg_dotprod (n0(1),3*npole,r0(:,:,1),r0(:,:,1))
      call tcg_dotprod (n0(2),3*npole,r0(:,:,2),r0(:,:,2))
      call tcg_t0 (r0(:,:,1),r0(:,:,2),P1(:,:,1),P1(:,:,2))
      call tcg_dotprod (t1(1),3*npole,r0(:,:,1),P1(:,:,1))
      call tcg_dotprod (t1(2),3*npole,r0(:,:,2),P1(:,:,2))
      t4(1) = n0(1) / t1(1)
      t4(2) = n0(2) / t1(2)
c
c     mu1 = mu0 + gamma0 * p0 (or r0)
c
      uind = uind + t4(1) * r0(:,:,1)
      uinp = uinp + t4(2) * r0(:,:,2)
      if (order .gt. 1 .or. tcgpeek) then
c
c     tae = T.alpha.E = T.udir
c
         call tcg_t0 (udir,udirp,tae(:,:,1),tae(:,:,2))
c
c     r1 = r0 - gamma0 * T.p0 (or T.r0, or P1)
c
         r1(:,:,1) = r0(:,:,1) - t4(1)*P1(:,:,1)
         r1(:,:,2) = r0(:,:,2) - t4(2)*P1(:,:,2)
c
c     n1 = r1.r1
c     check convergence, stop at tcg1 level if n1 is small enough
c
         call tcg_dotprod (n1(1),3*npole,r1(:,:,1),r1(:,:,1))
         call tcg_dotprod (n1(2),3*npole,r1(:,:,2),r1(:,:,2))
         call tcg_converge (converge,n1(1),n1(2))
         if (converge) then
            order = 1
            call tcg_resource (order)
         end if
      end if
c
c     cross terms
c     sp0 = r0.E
c     spp1 = r0.T.alpha.E = r0.T.mu0
c
      call tcg_dotprod (sp0,3*npole,r0(:,:,1),r0(:,:,2))
      if (tcgpeek) then
         call tcg_dotprod (spp1(1),3*npole,P1(:,:,1),udirp)
         call tcg_dotprod (spp1(2),3*npole,P1(:,:,2),udir)
      end if
c
c     tcg1 force and energy
c
      if (order .eq. 1) then
c
c     compute a(1) coefficients: a1...
c     and a(1k) coefficients: a1k...
c
         do i = 1, 2
            a110(i) = t4(i)
            a111(i) = 2.0d0 * sp0 / t1(i)
            a112(i) = -t4(i) * a111(i)
            a121(i) = 0.5d0 * a112(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a1k10a(i) = tcgomega
               a1k11a(i) = -tcgomega * t4(i)
               a1k11(i) = -2.0d0 * spp1(i) * tcgomega / t1(i)
               a1k12(i) = -t4(i) * a1k11(i)
               a1k20a(i) = a1k11a(i)
               a1k21(i) = 0.5d0 * a1k12(i)
            end do
c
c     mu1(peek) = mu1 + omega * alpha.r1
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = -0.5d0 * ((a121(1)+a1k21(1))*r0(:,:,1)
     &                                + a1k20a(1)*udirp)
            ubd(:,:,1) = -0.5d0 * ((a121(2)+a1k21(2))*r0(:,:,2)
     &                                + a1k20a(2)*udir)
            xde(:,:,1) = 0.5d0 * (uind + a110(2)*r0(:,:,1)
     &                 + (a111(2)+a1k11(2))*r0(:,:,2)
     &                 + (a112(2)+a1k12(2))*P1(:,:,2)
     &                 + a1k10a(2)*udir + a1k11a(2)*tae(:,:,1))
            xde(:,:,2) = 0.5d0 * (uinp + a110(1)*r0(:,:,2)
     &                 + (a111(1)+a1k11(1))*r0(:,:,1)
     &                 + (a112(1)+a1k12(1))*P1(:,:,1)
     &                 + a1k10a(1)*udirp + a1k11a(1)*tae(:,:,2))
         else
            ubp(:,:,1) = -0.5d0 * a121(1) * r0(:,:,1)
            ubd(:,:,1) = -0.5d0 * a121(2) * r0(:,:,2)
            xde(:,:,1) = ((t4(1)+a110(2))*r0(:,:,1) + a111(2)*r0(:,:,2)
     &                       + a112(2)*P1(:,:,2)) * 0.5d0
            xde(:,:,2) = ((t4(2)+a110(1))*r0(:,:,2) + a111(1)*r0(:,:,1)
     &                       + a112(1)*P1(:,:,1)) * 0.5d0
         end if
         uad(:,:,1) = r0(:,:,1)
         uap(:,:,1) = r0(:,:,2)
         goto 10
      end if
c
c     compute tcg2 intermediates
c     t2r0 = T^2.r0 = T.P1
c     t3r0 = T^3.r0
c     te = T.E = T.r0 = P1
c     t9 = r0.T^3.r0 = r0.T^2.T.r0
c
      call tcg_t0 (P1(:,:,1),P1(:,:,2),t2r0(:,:,1),t2r0(:,:,2))
      call tcg_t0 (t2r0(:,:,1),t2r0(:,:,2),t3r0(:,:,1),t3r0(:,:,2))
      call tcg_dotprod (t9(1),3*npole,t2r0(:,:,1),P1(:,:,1))
      call tcg_dotprod (t9(2),3*npole,t2r0(:,:,2),P1(:,:,2))
c
c     beta1 = r1.r1/r0.r0 = n1/n0
c     t2 = 1 + beta1
c
      beta1(1) = n1(1) / n0(1)
      beta1(2) = n1(2) / n0(2)
      t2(1) = 1.0d0 + beta1(1)
      t2(2) = 1.0d0 + beta1(2)
c
c     np1 = P1.P1
c     t8 = t5 = P1.T.p1 = P1.P2 = t2*np1 - t4*t9
c     t10 = t1^2 - n0.|P1|^2
c     t3 = t1*P1.P2 = t1*t8
c     gamma1 = t10/t3
c
      call tcg_dotprod (np1(1),3*npole,P1(:,:,1),P1(:,:,1))
      call tcg_dotprod (np1(2),3*npole,P1(:,:,2),P1(:,:,2))
      t8(1) = t2(1)*np1(1) - t4(1)*t9(1)
      t8(2) = t2(2)*np1(2) - t4(2)*t9(2)
      t10(1) = t1(1)**2 - n0(1)*np1(1)
      t10(2) = t1(2)**2 - n0(2)*np1(2)
      t3(1) = t1(1)*t8(1)
      t3(2) = t1(2)*t8(2)
      gamma1(1) = t10(1) / t3(1)
      gamma1(2) = t10(2) / t3(2)
c
c     mu2 = mu1 + gamma1*p1 = mu1 + gamma1*(r1 + beta1*p0)
c         = mu1 + gamma1*(r1 + beta1*r0)
c
      uind = uind + gamma1(1)*(r1(:,:,1) + beta1(1)*r0(:,:,1))
      uinp = uinp + gamma1(2)*(r1(:,:,2) + beta1(2)*r0(:,:,2))
      if (order .gt. 2 .or. tcgpeek) then
c
c     t2ae = T^2.alpha.E = T.tae
c
         call tcg_t0 (tae(:,:,1),tae(:,:,2),t2ae(:,:,1),t2ae(:,:,2))
c
c     r2 = r1 - gamma1 * T.p1 = r1 - gamma1 * P2
c        = r1 - gamma1 * (t2*P1 - t4*t2r0)
c     reuse r1 as r2
c
         r1(:,:,1) = r1(:,:,1)
     &                  - gamma1(1)*(t2(1)*P1(:,:,1)-t4(1)*t2r0(:,:,1))
         r1(:,:,2) = r1(:,:,2)
     &                  - gamma1(2)*(t2(2)*P1(:,:,2)-t4(2)*t2r0(:,:,2))
      end if
c
c     cross terms
c     sp1 = r0.T.E = P1.E
c     b1 = sp0 - gamma1*sp1
c     b2 = sp0*t2 - t4*sp1
c     spp2 = r0.T.T.alpha.E
c
      call tcg_dotprod (sp1,3*npole,P1(:,:,1),r0(:,:,2))
      b1(1) = sp0 - gamma1(1)*sp1
      b1(2) = sp0 - gamma1(2)*sp1
      b2(1) = sp0*t2(1) - t4(1)*sp1
      b2(2) = sp0*t2(2) - t4(2)*sp1
      if (tcgpeek) then
         call tcg_dotprod (spp2(1),3*npole,t2r0(:,:,1),udirp)
         call tcg_dotprod (spp2(2),3*npole,t2r0(:,:,2),udir)
      end if
c
c     tcg2 force and energy
c
      if (order .eq. 2) then
         if (tcgpeek) then
         end if
c
c     compute a(2) coefficients: a2...
c     and a(2k) coefficients: a2k...
c
         do i = 1, 2
            a232(i) = t1(i)*t4(i)*gamma1(i)*b2(i)/t3(i)
            a241(i) = a232(i)
            a231(i) = -n0(i)*b2(i)/t3(i)
     &                -2.0d0*t1(i)*t2(i)*gamma1(i)*b2(i)/t3(i)
     &                +t4(i)*gamma1(i)*sp0/t1(i)
            a223(i) = a232(i)
            a222(i) = a231(i)
            a221(i) = -t4(i)*b1(i)/t1(i) +2.0d0*t1(i)*b2(i)/t3(i)
     &                -t4(i)*t9(i)*gamma1(i)*b2(i)/t3(i)
     &                +2.0d0*t2(i)*np1(i)*gamma1(i)*b2(i)/t3(i)
     &                -t8(i)*gamma1(i)*b2(i)/t3(i)
     &                -2.0d0*t4(i)*np1(i)*sp0*gamma1(i)/(t1(i)**2)
            a220(i) = -gamma1(i)*t4(i)
            a214(i) = 2.0d0*a232(i)
            a213(i) = 2.0d0*a231(i)
            a212(i) = 2.0d0*a221(i)
            a211(i) = 2.0d0*(b1(i)/t1(i) -np1(i)*b2(i)/t3(i)
     &                   -(np1(i)**2)*gamma1(i)*b2(i)/t3(i)/t1(i)
     &                   +t9(i)*gamma1(i)*b2(i)/t3(i)
     &                   +np1(i)*sp0*gamma1(i)/(t1(i)**2))
            a21n1(i) = a220(i)
            a210(i) = t4(i) + gamma1(i)*t2(i)
         end do
         if (tcgpeek) then
            do i = 1, 2
               a2kwt2(i) = tcgomega*(t2(i)*spp1(i)-t4(i)*spp2(i))
               a2kwg1(i) = tcgomega*(spp1(i)-gamma1(i)*spp2(i))
               a2k41(i) = -a2kwt2(i)*t1(i)*t4(i)*gamma1(i)/t3(i)
               a2k32(i) = a2k41(i)
               a2k31(i) = -tcgomega*t4(i)*gamma1(i)*spp1(i)/t1(i)
     &                    +a2kwt2(i)*(n0(i)/t3(i)
     &                       +2.0d0*t1(i)*t2(i)*gamma1(i)/t3(i))
               a2k30a(i) = tcgomega*gamma1(i)*t4(i)
               a2k23(i) = a2k41(i)
               a2k22(i) = a2k31(i)
               a2k21(i) = 2.0d0*t4(i)*np1(i)/(t1(i)**2)*
     &                       tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(-2.0d0*t1(i)+
     &                       (t4(i)*t9(i)-2.0d0*np1(i)*t2(i)+t8(i))*
     &                       gamma1(i))/t3(i)
     &                    +t4(i)*a2kwg1(i)/t1(i)
               a2k21a(i) = a2k30a(i)
               a2k20a(i) = -tcgomega*(gamma1(i)*t2(i)+t4(i))
               a2k14(i) = 2.0d0*a2k41(i)
               a2k13(i) = 2.0d0*a2k22(i)
               a2k12(i) = 2.0d0*a2k21(i)
               a2k11(i) = -np1(i)/(t1(i)**2)*tcgomega*gamma1(i)*spp1(i)
     &                    +a2kwt2(i)*(np1(i)
     &                       +(np1(i)**2)*gamma1(i)/t1(i)
     &                       -t9(i)*gamma1(i))/t3(i)
     &                    -a2kwg1(i)/t1(i)
               a2k11(i) = 2.0d0*a2k11(i)
               a2k12a(i) = a2k30a(i)
               a2k11a(i) = a2k20a(i)
               a2k10a(i) = tcgomega
            end do
c
c     mu2(peek) = mu2 + omega * alpha.r2
c
            call tcg_add_omega_alpha2 (tcgomega,uind,uinp,
     &                                    r1(:,:,1),r1(:,:,2))
c
c     mutual and direct induced dipole components
c
            ubp(:,:,1) = -0.5d0*(a220(1)*r0(:,:,2)
     &                 + (a221(1)+a2k21(1))*r0(:,:,1)
     &                 + (a222(1)+a231(1)+a2k22(1)+a2k31(1))*P1(:,:,1)
     &                 + (a223(1)+a241(1)+a2k23(1)+a2k41(1))*t2r0(:,:,1)
     &                 + a2k21a(1)*tae(:,:,2) + a2k20a(1)*udirp)
            ubp(:,:,2) = -0.5d0*((a232(1)+a2k32(1))*P1(:,:,1)
     &                               + a2k30a(1)*udirp)
            ubd(:,:,1) = -0.5d0*(a220(2)*r0(:,:,1)
     &                 + (a221(2)+a2k21(2))*r0(:,:,2)
     &                 + (a222(2)+a231(2)+a2k22(2)+a2k31(2))*P1(:,:,2)
     &                 + (a223(2)+a241(2)+a2k23(2)+a2k41(2))*t2r0(:,:,2)
     &                 + a2k21a(2)*tae(:,:,1) + a2k20a(2)*udir)
            ubd(:,:,2) = -0.5d0*((a232(2)+a2k32(2))*P1(:,:,2)
     &                              + a2k30a(2)*udir)
            xde(:,:,1) = (a210(2)*r0(:,:,1) + a21n1(2)*P1(:,:,1)
     &                 + (a211(2)+a2k11(2))*r0(:,:,2)
     &                 + (a212(2)+a2k12(2))*P1(:,:,2)
     &                 + (a213(2)+a2k13(2))*t2r0(:,:,2)
     &                 + (a214(2)+a2k14(2))*t3r0(:,:,2)
     &                 + a2k11a(2)*tae(:,:,1) + a2k12a(2)*t2ae(:,:,1)
     &                 + a2k10a(2)*udir + uind) * 0.5d0
            xde(:,:,2) = (a210(1)*r0(:,:,2) + a21n1(1)*P1(:,:,2)
     &                 + (a211(1)+a2k11(1))*r0(:,:,1)
     &                 + (a212(1)+a2k12(1))*P1(:,:,1)
     &                 + (a213(1)+a2k13(1))*t2r0(:,:,1)
     &                 + (a214(1)+a2k14(1))*t3r0(:,:,1)
     &                 + a2k11a(1)*tae(:,:,2) + a2k12a(1)*t2ae(:,:,2)
     &                 + a2k10a(1)*udirp + uinp) * 0.5d0
         else
            ubp(:,:,1) = -0.5d0*(a220(1)*r0(:,:,2)
     &                              + a221(1)*r0(:,:,1)
     &                              + (a222(1)+a231(1))*P1(:,:,1)
     &                              + (a223(1)+a241(1))*t2r0(:,:,1))
            ubp(:,:,2) = -0.5d0*a232(1)*P1(:,:,1)
            ubd(:,:,1) = -0.5d0*(a220(2)*r0(:,:,1)
     &                              + a221(2)*r0(:,:,2)
     &                              + (a222(2)+a231(2))*P1(:,:,2)
     &                              + (a223(2)+a241(2))*t2r0(:,:,2))
            ubd(:,:,2) = -0.5d0*a232(2)*P1(:,:,2)
            xde(:,:,1) = (a210(2)*r0(:,:,1) + a21n1(2)*P1(:,:,1)
     &                 + a211(2)*r0(:,:,2) + a212(2)*P1(:,:,2)
     &                 + a213(2)*t2r0(:,:,2) + a214(2)*t3r0(:,:,2)
     &                 + uind) * 0.5d0
            xde(:,:,2) = (a210(1)*r0(:,:,2) + a21n1(1)*P1(:,:,2)
     &                 + a211(1)*r0(:,:,1) + a212(1)*P1(:,:,1)
     &                 + a213(1)*t2r0(:,:,1) + a214(1)*t3r0(:,:,1)
     &                 + uinp) * 0.5d0
         end if
         uad(:,:,1) = r0(:,:,1)
         uad(:,:,2) = P1(:,:,1)
         uap(:,:,1) = r0(:,:,2)
         uap(:,:,2) = P1(:,:,2)
         goto 10
      end if
c
c     store induced dipoles and use uind/p to store xde arrays
c
   10 continue
      call tcg_store
      uind = xde(:,:,1)
      uinp = xde(:,:,2)
c
c     perform deallocation for some local arrays
c
      deallocate (r0)
      deallocate (xde)
      deallocate (P1)
      deallocate (r1)
      deallocate (tae)
      deallocate (t2r0)
      deallocate (t3r0)
      deallocate (t2ae)
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
c     "tcg_alphaquad" computes the quadratic form, <a ALPHA b>,
c     where ALPHA is the diagonal polarization matrix
c
c
      subroutine tcg_alphaquad (scalar,a,b)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 scalar,a(3,*),b(3,*)
c
c
      scalar = 0.0d0
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO reduction(+:scalar) schedule(guided)
      do i = 1, npole
         do j = 1, 3
         scalar = scalar + a(j,i)*b(j,i)*polarity(i)
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
c     "tcg_resource" sets the number of mutual induced dipole pairs
c     based on the argument passed in
c
c
      subroutine tcg_resource (order)
      use iounit
      use poltcg
      implicit none
      integer order
c
c
      if (order .lt. 1 .or. order .gt. 2) then
         write (iout,10)
   10    format (/,' TCG_RESOURCE -- Argument ORDER Is Out of Range')
         call fatal
      end if
      tcgnab = order
      if (tcgguess) then
         tcgnab = tcgnab + 1
      end if
      return
      end
c
c
c     ###############################
c     ##                           ##
c     ##  subroutine tcg_converge  ##
c     ##                           ##
c     ###############################
c
c
c     "tcg_converge" checks the convergence of the residuals and writes
c     the result to the logical variable passed in
c
c
      subroutine tcg_converge (ifconv,rsdsq1,rsdsq2)
      use polar
      use polpot
      use units
      implicit none
      logical ifconv
      real*8 rsdsq1,rsdsq2
      real*8 eps
c
c
      ifconv = .false.
      eps = max(rsdsq1,rsdsq2)
      eps = debye * sqrt(eps/dble(npolar))
      if (eps .lt. poleps) then
         ifconv = .true.
      end if
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
c     "tcg_alpha12" computes source1 = alpha.source1 and
c     source2 = alpha.source2
c
c
      subroutine tcg_alpha12 (source1,source2)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 source1(3,*),source2(3,*)
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
c     "tcg_alpha22" computes result1 = alpha.source1 and
c     result2 = alpha.source2
c
c
      subroutine tcg_alpha22 (source1,source2,result1,result2)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 source1(3,*),source2(3,*)
      real*8 result1(3,*),result2(3,*)
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
c     #######################################
c     ##                                   ##
c     ##  subroutine tcg_add_omega_alpha2  ##
c     ##                                   ##
c     #######################################
c
c
c     "tcg_add_omega_alpha2" computes
c     u1 = u1 + omega * alpha.r1 and
c     u2 = u2 + omega * alpha.r2
c
c
      subroutine tcg_add_omega_alpha2 (omega,u1,u2,r1,r2)
      use mpole
      use polar
      implicit none
      integer i,j
      real*8 omega
      real*8 u1(3,*),u2(3,*)
      real*8 r1(3,*),r2(3,*)
c
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, npole
         do j = 1, 3
            u1(j,i) = u1(j,i) + omega*polarity(i)*r1(j,i)
            u2(j,i) = u2(j,i) + omega*polarity(i)*r2(j,i)
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
c     "tcg_dotprod" computes the dot product of two 1D arrays,
c     both of which are of n elements
c
c
      subroutine tcg_dotprod (scalar,n,a,b)
      implicit none
      integer i,n
      real*8 scalar,a(*),b(*)
c
c
      scalar = 0.0d0
!$OMP PARALLEL default(shared) private(i)
!$OMP DO reduction(+:scalar) schedule(guided)
      do i = 1, n
         scalar = scalar + a(i)*b(i)
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
      integer i,j
      real*8 ind(3,*),inp(3,*)
      real*8 v3d(3,*),v3p(3,*)
c
c
c     copy induced dipoles to storage arrays
c
      call tcg_store
c
c     move TCG components to induced dipoles
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = ind(j,i)
            uinp(j,i) = inp(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
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
c     restore the global induced dipoles
c
      call tcg_restore
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
      real*8 ind(3,*),inp(3,*)
      real*8 v3d(3,*),v3p(3,*)
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
c     ############################
c     ##                        ##
c     ##  subroutine tcg_store  ##
c     ##                        ##
c     ############################
c
c
c     "tcg_store" moves the value in global variables uind/uinp to
c     uindt/uinpt
c
c
      subroutine tcg_store
      use atoms
      use polar
      use poltcg
      implicit none
      integer i,j
c
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, n
         do j = 1, 3
            uindt(j,i) = uind(j,i)
            uinpt(j,i) = uinp(j,i)
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
c     ##  subroutine tcg_restore  ##
c     ##                          ##
c     ##############################
c
c
c     "tcg_restore" moves the value in global variables uindt/uindp
c     to uind/uinp
c
c
      subroutine tcg_restore
      use atoms
      use polar
      use poltcg
      implicit none
      integer i,j
c
c
!$OMP PARALLEL default(shared) private(i,j)
!$OMP DO schedule(guided)
      do i = 1, n
         do j = 1, 3
            uind(j,i) = uindt(j,i)
            uinp(j,i) = uinpt(j,i)
         end do
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar1d0  --  single-loop polarization energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar1d0" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field, assuming
c     the induced dipoles have already been calculated
c
c
      subroutine epolar1d0
      use atoms
      use boxes
      use chgpot
      use energi
      use ewald
      use limits
      use math
      use mpole
      use polar
      use polpot
      use poltcg
      implicit none
      integer i,j,ii
      real*8 e,f,fi,term
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
c
c
c     zero out the total polarization energy
c
      ep = 0.0d0
      if (npole .eq. 0)  return
c
c     restore the induced dipoles if using the TCG method
c
      if (poltyp(1:3) .eq. 'TCG') then
!$OMP    PARALLEL default(shared) private(i,j)
!$OMP    DO schedule(guided)
         do i = 1, npole
            do j = 1, 3
               uind(j,i) = uindt(j,i)
               uinp(j,i) = uinpt(j,i)
            end do
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     set the energy conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(shared) private(i,j,fi,e)
!$OMP DO reduction(+:ep) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (douind(i)) then
            fi = 0.5d0 * f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*(uind(j,i)*udirp(j,i)+uinp(j,i)*udir(j,i))
            end do
            ep = ep + e
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
            term = (2.0d0/3.0d0) * f * (pi/volbox)
            e = term * (xd*xu+yd*yu+zd*zu)
            ep = ep + e
         end if
      end if
      return
      end
