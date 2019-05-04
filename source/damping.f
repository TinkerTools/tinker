c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ## 
c     ############################################################ 
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampewald  --  find Ewald damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampewald" generates coefficients for Ewald error function
c     damping for powers of the interatomic distance
c
c
      subroutine dampewald (rorder,r,r2,scale)
      use ewald
      use math
      implicit none
      integer i,maxi
      integer rorder
      real*8 r,r2,bfac
      real*8 alsq2,alsq2n
      real*8 exp2a,ralpha
      real*8 scale(*)
      real*8, allocatable :: bn(:)
c
c
c     initialize Ewald damping factors and set storage size
c
      do i = 1, rorder
         scale(i) = 1.0d0
      end do
      maxi = (rorder-1) / 2
c
c     perform dynamic allocation of some local arrays
c
      allocate (bn(0:maxi))
c     
c     compute the successive Ewald damping factors
c
      ralpha = aewald * r
      bn(0) = erfc(ralpha) / r
      scale(1) = bn(0)
      alsq2 = 2.0d0 * aewald**2
      alsq2n = 0.0d0
      if (aewald .gt. 0.0d0)  alsq2n = 1.0d0 / (sqrtpi*aewald)
      exp2a = exp(-ralpha**2)
      do i = 1, maxi
         bfac = dble(2*i-1)
         alsq2n = alsq2 * alsq2n
         bn(i) = (bfac*bn(i-1)+alsq2n*exp2a) / r2
         scale(2*i+1) = bn(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (bn)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dampthole  --  find Thole damping coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dampthole" generates coefficients for the Thole damping
c     function for powers of the interatomic distance
c
c     literature reference:
c
c     B. T. Thole, "Molecular Polarizabilities Calculated with a
c     Modified Dipole Interaction", Chemical Physics, 59, 341-350 (1981)
c
c
      subroutine dampthole (i,k,rorder,r,scale)
      use polar
      implicit none
      integer i,j,k
      integer rorder
      real*8 r,damp
      real*8 damp2
      real*8 damp3
      real*8 expdamp
      real*8 pdi,pti
      real*8 pdk,ptk
      real*8 pgamma
      real*8 scale(*)
c
c
c     initialize the Thole damping factors to a value of one
c
      do j = 1, rorder
         scale(j) = 1.0d0
      end do
c
c     get the Thole damping parameters for the two sites
c
      pdi = pdamp(i)
      pti = thole(i)
      pdk = pdamp(k)
      ptk = thole(k)
c
c     assign the Thole polarization model damping factors
c
      damp = pdi * pdk
      if (damp .ne. 0.0d0) then
         pgamma = min(pti,ptk)
         damp = pgamma * (r/damp)**3
         if (damp .lt. 50.0d0) then
            expdamp = exp(-damp)
            scale(3) = 1.0d0 - expdamp
            scale(5) = 1.0d0 - (1.0d0+damp)*expdamp
            if (rorder .ge. 7) then
               damp2 = damp * damp
               scale(7) = 1.0d0 - (1.0d0+damp+0.6d0*damp2)*expdamp
               if (rorder .ge. 9) then
                  damp3 = damp * damp2
                  scale(9) = 1.0d0 - (1.0d0 + damp
     &                          + (18.0d0/35.0d0)*damp2
     &                          + (9.0d0/35.0d0)*damp3) * expdamp
               end if
            end if
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine damppole  --  penetration damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "damppole" generates coefficients for the charge penetration
c     damping function for powers of the interatomic distance
c
c     literature references:
c
c     L. V. Slipchenko and M. S. Gordon, "Electrostatic Energy in the
c     Effective Fragment Potential Method: Theory and Application to
c     the Benzene Dimer", Journal of Computational Chemistry, 28,
c     276-291 (2007)  [Gordon f1 and f2 models]
c
c     J. A. Rackers, Q. Wang, C. Liu, J.-P. Piquemal, P. Ren and
c     J. W. Ponder, "An Optimized Charge Penetration Model for Use with
c     the AMOEBA Force Field", Physical Chemistry Chemical Physics, 19,
c     276-291 (2017)
c
c
      subroutine damppole (r,rorder,alphai,alphak,dmpi,dmpk,dmpik)
      implicit none
      integer rorder
      real*8 r,termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 alphai2,alphak2
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampi3
      real*8 dampi4,dampi5
      real*8 dampi6,dampi7
      real*8 dampi8
      real*8 dampk2,dampk3
      real*8 dampk4,dampk5
      real*8 dampk6
      real*8 dmpi(*)
      real*8 dmpk(*)
      real*8 dmpik(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      dampi2 = dampi * dampi
      dampi3 = dampi * dampi2
      dampi4 = dampi2 * dampi2
      dampi5 = dampi2 * dampi3
      dmpi(1) = 1.0d0 - (1.0d0 + 0.5d0*dampi)*expi
      dmpi(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
      dmpi(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &             + dampi3/6.0d0)*expi
      dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &             + dampi3/6.0d0 + dampi4/30.0d0)*expi
      dmpi(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &             + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
     &             + dampi5/210.0d0)*expi
      if (diff .lt. eps) then
         dmpk(1) = dmpi(1)
         dmpk(3) = dmpi(3)
         dmpk(5) = dmpi(5)
         dmpk(7) = dmpi(7)
         dmpk(9) = dmpi(9)
      else
         dampk2 = dampk * dampk
         dampk3 = dampk * dampk2
         dampk4 = dampk2 * dampk2
         dampk5 = dampk2 * dampk3
         dmpk(1) = 1.0d0 - (1.0d0 + 0.5d0*dampk)*expk
         dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
         dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0)*expk
         dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0 + dampk4/30.0d0)*expk
         dmpk(9) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0 + 4.0d0*dampk4/105.0d0
     &                + dampk5/210.0d0)*expk
      end if
c
c     valence-valence charge penetration damping for Gordon f1
c
      if (diff .lt. eps) then
         dampi6 = dampi3 * dampi3
         dampi7 = dampi3 * dampi4
         dmpik(1) = 1.0d0 - (1.0d0 + 11.0d0*dampi/16.0d0
     &                 + 3.0d0*dampi2/16.0d0
     &                 + dampi3/48.0d0)*expi
         dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + 7.0d0*dampi3/48.0d0
     &                 + dampi4/48.0d0)*expi
         dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/24.0d0
     &                 + dampi5/144.0d0)*expi
         dmpik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/24.0d0
     &                 + dampi5/120.0d0 + dampi6/720.0d0)*expi
         dmpik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/24.0d0
     &                 + dampi5/120.0d0 + dampi6/720.0d0
     &                 + dampi7/5040.0d0)*expi
         if (rorder .ge. 11) then
            dampi8 = dampi4 * dampi4
            dmpik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                     + dampi3/6.0d0 + dampi4/24.0d0
     &                     + dampi5/120.0d0 + dampi6/720.0d0
     &                     + dampi7/5040.0d0 + dampi8/45360.0d0)*expi
         end if
      else
         alphai2 = alphai * alphai
         alphak2 = alphak * alphak
         termi = alphak2 / (alphak2-alphai2)
         termk = alphai2 / (alphai2-alphak2)
         termi2 = termi * termi
         termk2 = termk * termk
         dmpik(1) = 1.0d0 - termi2*(1.0d0 + 2.0d0*termk
     &                 + 0.5d0*dampi)*expi
     &              - termk2*(1.0d0 + 2.0d0*termi
     &                   + 0.5d0*dampk)*expk
         dmpik(3) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2)*expi
     &                 - termk2*(1.0d0+dampk+0.5d0*dampk2)*expk
     &                 - 2.0d0*termi2*termk*(1.0d0+dampi)*expi
     &                 - 2.0d0*termk2*termi*(1.0d0+dampk)*expk
         dmpik(5) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0)*expi
     &              - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0)*expk
     &              - 2.0d0*termi2*termk
     &                   *(1.0d0 + dampi + dampi2/3.0d0)*expi
     &              - 2.0d0*termk2*termi
     &                   *(1.0d0 + dampk + dampk2/3.0d0)*expk
         dmpik(7) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/30.0d0)*expi
     &              - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0 + dampk4/30.0d0)*expk
     &              - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                   + 2.0d0*dampi2/5.0d0 + dampi3/15.0d0)*expi
     &              - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                   + 2.0d0*dampk2/5.0d0 + dampk3/15.0d0)*expk
         dmpik(9) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
     &                 + dampi5/210.0d0)*expi
     &              - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                   + dampk3/6.0d0 + 4.0d0*dampk4/105.0d0
     &                   + dampk5/210.0d0)*expk
     &              - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                   + 3.0d0*dampi2/7.0d0
     &                   + 2.0d0*dampi3/21.0d0
     &                   + dampi4/105.0d0)*expi 
     &              - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                   + 3.0d0*dampk2/7.0d0
     &                   + 2.0d0*dampk3/21.0d0
     &                   + dampk4/105.0d0)*expk
         if (rorder .ge. 11) then
            dampi6 = dampi3 * dampi3
            dampk6 = dampk3 * dampk3
            dmpik(11) = 1.0d0 - termi2*(1.0d0 + dampi
     &                     + 0.5d0*dampi2 + dampi3/6.0d0
     &                     + 5.0d0*dampi4/126.0d0
     &                     + 2.0d0*dampi5/315.0d0
     &                     + dampi6/1890.0d0)*expi
     &                  - termk2*(1.0d0 + dampk
     &                       + 0.5d0*dampk2 + dampk3/6.0d0
     &                       + 5.0d0*dampk4/126.0d0
     &                       + 2.0d0*dampk5/315.0d0
     &                       + dampk6/1890.0d0)*expk
     &                  - 2.0d0*termi2*termk*(1.0d0 + dampi
     &                       + 4.0d0*dampi2/9.0d0 + dampi3/9.0d0
     &                       + dampi4/63.0d0 + dampi5/945.0d0)*expi
     &                  - 2.0d0*termk2*termi*(1.0d0 + dampk
     &                       + 4.0d0*dampk2/9.0d0 + dampk3/9.0d0
     &                       + dampk4/63.0d0 + dampk5/945.0d0)*expk 
         end if
      end if
c
c     core-valence charge penetration damping for Gordon f2
c
c     dampi2 = dampi * dampi
c     dampi3 = dampi * dampi2
c     dmpi(1) = 1.0d0 - dampi
c     dmpi(3) = 1.0d0 - (1.0d0 + dampi)*expi
c     dmpi(5) = 1.0d0 - (1.0d0 + dampi + dampi2/3.0d0)*expi
c     dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.4d0*dampi2
c    &             + dampi3/15.0d0)*expi
c     if (diff .lt. eps) then
c        dmpk(1) = dmpi(1)
c        dmpk(3) = dmpi(3)
c        dmpk(5) = dmpi(5)
c        dmpk(7) = dmpi(7)
c     else
c        dampk2 = dampk * dampk
c        dampk3 = dampk * dampk2
c        dmpk(1) = 1.0d0 - dampk
c        dmpk(3) = 1.0d0 - (1.0d0 + dampk)*expk
c        dmpk(5) = 1.0d0 - (1.0d0 + dampk + dampk2/3.0d0)*expk
c        dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.4d0*dampk2
c    &                + dampk3/15.0d0)*expk
c     end if
c
c     valence-valence charge penetration damping for Gordon f2
c
c     dampi4 = dampi2 * dampi2
c     dampi5 = dampi2 * dampi3
c     if (diff .lt. eps) then
c        dampk6 = dampk3 * dampk3
c        dmpik(1) = 1.0d0 - (1.0d0 + 0.5d0*dampi)*expi
c        dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
c        dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
c    &                 + dampi3/6.0d0)*expi
c        dmpik(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
c    &                 + dampi3/6.0d0 + dampi4/30.0d0)*expi
c        dmpik(9) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
c    &                 + dampi3/6.0d0 + 4.0d0*dampi4/105.0d0
c    &                 + dampi5/210.0d0)*expi
c        if (rorder .ge. 11) then
c           dmpik(11) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
c    &                     + dampi3/6.0d0 + 5.0d0*dampi4/126.0d0
c    &                     + 2.0d0*dampi5/315.0d0
c    &                     + dampi6/1890.0d0)*expi
c        end if
c     else
c        dampk4 = dampk2 * dampk2
c        dampk5 = dampk2 * dampk3
c        alphai2 = alphai * alphai
c        alphak2 = alphak * alphak
c        termi = alphak2 / (alphak2-alphai2)
c        termk = alphai2 / (alphai2-alphak2)
c        dmpik(1) = 1.0d0 - termi*expi - termk*expk
c        dmpik(3) = 1.0d0 - termi*(1.0d0 + dampi)*expi
c    &                 - termk*(1.0d0 + dampk)*expk
c        dmpik(5) = 1.0d0 - termi*(1.0d0 + dampi + dampi2/3.0d0)*expi
c    &                 - termk*(1.0d0 + dampk + dampk2/3.0d0)*expk
c        dmpik(7) = 1.0d0 - termi*(1.0d0 + dampi + 0.4d0*dampi2
c    &                 + dampi3/15.0d0)*expi
c    &              - termk*(1.0d0 + dampk + 0.4d0*dampk2
c    &                   + dampk3/15.0d0)*expk
c        dmpik(9) = 1.0d0 - termi*(1.0d0 + dampi + 3.0d0*dampi2/7.0d0
c    &                 + 2.0d0*dampi3/21.0d0 + dampi4/105.0d0)*expi
c    &              - termk*(1.0d0 + dampk + 3.0d0*dampk2/7.0d0
c    &                   + 2.0d0*dampk3/21.0d0 + dampk4/105.0d0)*expk
c        if (rorder .ge. 11) then
c           dmpik(11) = 1.0d0 - termi*(1.0d0 + dampi
c    &                     + 4.0d0*dampi2/9.0d0 + dampi3/9.0d0
c    &                     + dampi4/63.0d0 + dampi5/945.0d0)*expi
c    &                  - termk*(1.0d0 + dampk
c    &                       + 4.0d0*dampk2/9.0d0 + dampk3/9.0d0
c    &                       + dampk4/63.0d0 + dampk5/945.0d0)*expk
c        end if
c     end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine damppolar  --  polarization damping coefficents  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "damppolar" generates coefficients for the charge penetration
c     damping function used for polarization interactions
c
c
      subroutine damppolar (r,alphai,alphak,dmpi,dmpk,dmpik)
      implicit none
      real*8 r,termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 alphai2,alphak2
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampi3
      real*8 dampi4,dampi5
      real*8 dampk2,dampk3
      real*8 dmpi(*)
      real*8 dmpk(*)
      real*8 dmpik(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      dampi2 = dampi * dampi
      dampi3 = dampi * dampi2
      dmpi(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
      dmpi(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2 
     &             + dampi3/6.0d0)*expi
      if (diff .lt. eps) then
         dmpk(3) = dmpi(3)
         dmpk(5) = dmpi(5)
      else
         dampk2 = dampk * dampk
         dampk3 = dampk * dampk2
         dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
         dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0)*expk
      end if
c
c     valence-valence charge penetration damping for Gordon f1
c
      if (diff .lt. eps) then
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + 7.0d0*dampi3/48.0d0
     &                 + dampi4/48.0d0)*expi
         dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/24.0d0
     &                 + dampi5/144.0d0)*expi
      else
         alphai2 = alphai * alphai
         alphak2 = alphak * alphak
         termi = alphak2 / (alphak2-alphai2)
         termk = alphai2 / (alphai2-alphak2)
         termi2 = termi * termi
         termk2 = termk * termk
         dmpik(3) = 1.0d0 - termi2*(1.0d0 + dampi
     &                         + 0.5d0*dampi2)*expi
     &                 - termk2*(1.0d0 + dampk
     &                      + 0.5d0*dampk2)*expk
     &                 - 2.0d0*termi2*termk*(1.0d0+dampi)*expi
     &                 - 2.0d0*termk2*termi*(1.0d0+dampk)*expk
         dmpik(5) = 1.0d0 - termi2*(1.0d0 + dampi + 0.5d0*dampi2
     &                                 + dampi3/6.0d0)*expi
     &                 - termk2*(1.0d0 + dampk + 0.5d0*dampk2
     &                              + dampk3/6.00)*expk
     &                 - 2.0d0*termi2*termk
     &                      *(1.0d0 + dampi + dampi2/3.0d0)*expi
     &                 - 2.0d0*termk2*termi
     &                      *(1.0d0 + dampk + dampk2/3.0d0)*expk
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampdir  --  direct field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampdir" generates coefficients for the direct field damping
c     function for powers of the interatomic distance
c
c
      subroutine dampdir (r,alphai,alphak,dmpi,dmpk)
      implicit none
      real*8 r,termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampk2
      real*8 dampi3,dampk3
      real*8 dampi4,dampk4
      real*8 dmpi(*)
      real*8 dmpk(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      dampi2 = dampi * dampi
      dampi3 = dampi * dampi2
      dampi4 = dampi2 * dampi2
      dmpi(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2)*expi
      dmpi(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2 
     &             + dampi3/6.0d0)*expi
      dmpi(7) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &             + dampi3/6.0d0 + dampi4/30.0d0)*expi
      if (diff .lt. eps) then
         dmpk(3) = dmpi(3)
         dmpk(5) = dmpi(5)
         dmpk(7) = dmpi(7)
      else
         dampk2 = dampk * dampk
         dampk3 = dampk * dampk2
         dampk4 = dampk2 * dampk2
         dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
         dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0)*expk
         dmpk(7) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &                + dampk3/6.0d0 + dampk4/30.0d0)*expk
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dampmut  --  mutual field damping coefficents  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dampmut" generates coefficients for the mutual field damping
c     function for powers of the interatomic distance
c
c
      subroutine dampmut (r,alphai,alphak,dmpik)
      implicit none
      real*8 r,termi,termk
      real*8 termi2,termk2
      real*8 alphai,alphak
      real*8 alphai2,alphak2
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 dampi2,dampi3
      real*8 dampi4,dampi5
      real*8 dampk2,dampk3
      real*8 dmpik(*)
c
c
c     compute tolerance and exponential damping factors
c
      eps = 0.001d0
      diff = abs(alphai-alphak)
      dampi = alphai * r
      dampk = alphak * r
      expi = exp(-dampi)
      expk = exp(-dampk)
c
c     valence-valence charge penetration damping for Gordon f1
c
      dampi2 = dampi * dampi
      dampi3 = dampi * dampi2
      if (diff .lt. eps) then
         dampi4 = dampi2 * dampi2
         dampi5 = dampi2 * dampi3
         dmpik(3) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + 7.0d0*dampi3/48.0d0
     &                 + dampi4/48.0d0)*expi
         dmpik(5) = 1.0d0 - (1.0d0 + dampi + 0.5d0*dampi2
     &                 + dampi3/6.0d0 + dampi4/24.0d0
     &                 + dampi5/144.0d0)*expi
      else
         dampk2 = dampk * dampk
         dampk3 = dampk * dampk2
         alphai2 = alphai * alphai
         alphak2 = alphak * alphak
         termi = alphak2 / (alphak2-alphai2)
         termk = alphai2 / (alphai2-alphak2)
         termi2 = termi * termi
         termk2 = termk * termk
         dmpik(3) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2)*expi
     &                 - termk2*(1.0d0+dampk+0.5d0*dampk2)*expk
     &                 - 2.0d0*termi2*termk*(1.0d0+dampi)*expi
     &                 - 2.0d0*termk2*termi*(1.0d0+dampk)*expk
         dmpik(5) = 1.0d0 - termi2*(1.0d0+dampi+0.5d0*dampi2
     &                         +dampi3/6.0d0)*expi
     &                 - termk2*(1.0d0+dampk+0.5d0*dampk2
     &                      +dampk3/6.00)*expk
     &                 - 2.0d0*termi2*termk
     &                      *(1.0+dampi+dampi2/3.0d0)*expi
     &                 - 2.0d0*termk2*termi
     &                      *(1.0+dampk+dampk2/3.0d0)*expk
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine damppot  --  electrostatic potential damping  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "damppot" generates coefficients for the charge penetration
c     damping function used for the electrostatic potential
c
c
      subroutine damppot (r,alphak,dmpk)
      implicit none
      real*8 r,alphak
      real*8 expk,dampk
      real*8 dampk2,dampk3
      real*8 dmpk(*)
c
c
c     compute common exponential factors for damping
c
      dampk = alphak * r
      expk = exp(-dampk)
c
c     core-valence charge penetration damping for Gordon f1
c
      dampk2 = dampk * dampk
      dampk3 = dampk * dampk2
      dmpk(1) = 1.0d0 - (1.0d0 + 0.5d0*dampk)*expk
      dmpk(3) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2)*expk
      dmpk(5) = 1.0d0 - (1.0d0 + dampk + 0.5d0*dampk2
     &             + dampk3/6.0d0)*expk
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine damprep  --  find repulsion damping coefficents  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "damprep" generates coefficients for the Pauli repulsion
c     damping function for powers of the interatomic distance
c
c     literature reference:
c
c     J. A. Rackers and J. W. Ponder, "Classical Pauli Repulsion: An
c     Anisotropic, Atomic Multipole Model", Journal of Chemical Physics,
c     150, 084104 (2019)
c
c
      subroutine damprep (r,r2,rr1,rr3,rr5,rr7,rr9,rr11,
     &                       rorder,dmpi,dmpk,dmpik)
      implicit none
      integer rorder
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7,r8
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9,rr11
      real*8 s,ds,d2s
      real*8 d3s,d4s,d5s
      real*8 dmpi,dmpk
      real*8 dmpi2,dmpk2
      real*8 dmpi22,dmpi23
      real*8 dmpi24,dmpi25
      real*8 dmpi26,dmpi27
      real*8 dmpk22,dmpk23
      real*8 dmpk24,dmpk25
      real*8 dmpk26
      real*8 eps,diff
      real*8 expi,expk
      real*8 dampi,dampk
      real*8 pre,term,tmp
      real*8 dmpik(*)
c
c
c     compute tolerance value for damping exponents
c
      eps = 0.001d0
      diff = abs(dmpi-dmpk)
c
c     treat the case where alpha damping exponents are equal
c
      if (diff .lt. eps) then
         r3 = r2 * r
         r4 = r3 * r
         r5 = r4 * r
         r6 = r5 * r
         r7 = r6 * r
         dmpi2 = 0.5d0 * dmpi
         dampi = dmpi2 * r
         expi = exp(-dampi)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpi25 = dmpi24 * dmpi2
         dmpi26 = dmpi25 * dmpi2
         pre = 128.0d0        
         s = (r + dmpi2*r2 + dmpi22*r3/3.0d0) * expi
         ds = (dmpi22*r3 + dmpi23*r4) * expi / 3.0d0
         d2s = dmpi24 * expi * r5 / 9.0d0
         d3s = dmpi25 * expi * r6 / 45.0d0
         d4s = (dmpi25*r6 + dmpi26*r7) * expi / 315.0d0
         if (rorder .ge. 11) then
            r8 = r7 * r
            dmpi27 = dmpi2 * dmpi26
            d5s = (dmpi25*r6 + dmpi26*r7 + dmpi27*r8/3.0d0)
     &                * expi / 945.0d0
         end if
c
c     treat the case where alpha damping exponents are unequal
c
      else
         r3 = r2 * r
         r4 = r3 * r
         r5 = r4 * r
         dmpi2 = 0.5d0 * dmpi
         dmpk2 = 0.5d0 * dmpk
         dampi = dmpi2 * r
         dampk = dmpk2 * r
         expi = exp(-dampi)
         expk = exp(-dampk)
         dmpi22 = dmpi2 * dmpi2
         dmpi23 = dmpi22 * dmpi2
         dmpi24 = dmpi23 * dmpi2
         dmpi25 = dmpi24 * dmpi2
         dmpk22 = dmpk2 * dmpk2
         dmpk23 = dmpk22 * dmpk2
         dmpk24 = dmpk23 * dmpk2
         dmpk25 = dmpk24 * dmpk2
         term = dmpi22 - dmpk22
         pre = 8192.0d0 * dmpi23 * dmpk23 / term**4
         tmp = 4.0d0 * dmpi2 * dmpk2 / term
         s = (dampi-tmp)*expk + (dampk+tmp)*expi
         ds = (dmpi2*dmpk2*r2 - 4.0d0*dmpi2*dmpk22*r/term
     &            - 4.0d0*dmpi2*dmpk2/term) * expk
     &      + (dmpi2*dmpk2*r2 + 4.0d0*dmpi22*dmpk2*r/term
     &            + 4.0d0*dmpi2*dmpk2/term) * expi
         d2s = (dmpi2*dmpk2*r2/3.0d0
     &             + dmpi2*dmpk22*r3/3.0d0
     &             - (4.0d0/3.0d0)*dmpi2*dmpk23*r2/term
     &             - 4.0d0*dmpi2*dmpk22*r/term
     &             - 4.0d0*dmpi2*dmpk2/term) * expk
     &       + (dmpi2*dmpk2*r2/3.0d0
     &             + dmpi22*dmpk2*r3/3.0d0
     &             + (4.0d0/3.0d0)*dmpi23*dmpk2*r2/term
     &             + 4.0d0*dmpi22*dmpk2*r/term
     &             + 4.0d0*dmpi2*dmpk2/term) * expi
         d3s = (dmpi2*dmpk23*r4/15.0d0
     &             + dmpi2*dmpk22*r3/5.0d0
     &             + dmpi2*dmpk2*r2/5.0d0
     &             - (4.0d0/15.0d0)*dmpi2*dmpk24*r3/term
     &             - (8.0d0/5.0d0)*dmpi2*dmpk23*r2/term
     &             - 4.0d0*dmpi2*dmpk22*r/term
     &             - 4.0d0/term*dmpi2*dmpk2) * expk
     &       + (dmpi23*dmpk2*r4/15.0d0 
     &             + dmpi22*dmpk2*r3/5.0d0
     &             + dmpi2*dmpk2*r2/5.0d0 
     &             + (4.0d0/15.0d0)*dmpi24*dmpk2*r3/term
     &             + (8.0d0/5.0d0)*dmpi23*dmpk2*r2/term
     &             + 4.0d0*dmpi22*dmpk2*r/term
     &             + 4.0d0/term*dmpi2*dmpk2) * expi
         d4s = (dmpi2*dmpk24*r5/105.0d0
     &             + (2.0d0/35.0d0)*dmpi2*dmpk23*r4
     &             + dmpi2*dmpk22*r3/7.0d0
     &             + dmpi2*dmpk2*r2/7.0d0
     &             - (4.0d0/105.0d0)*dmpi2*dmpk25*r4/term
     &             - (8.0d0/21.0d0)*dmpi2*dmpk24*r3/term
     &             - (12.0d0/7.0d0)*dmpi2*dmpk23*r2/term
     &             - 4.0d0*dmpi2*dmpk22*r/term
     &             - 4.0d0*dmpi2*dmpk2/term) * expk
     &       + (dmpi24*dmpk2*r5/105.0d0
     &             + (2.0d0/35.0d0)*dmpi23*dmpk2*r4
     &             + dmpi22*dmpk2*r3/7.0d0
     &             + dmpi2*dmpk2*r2/7.0d0
     &             + (4.0d0/105.0d0)*dmpi25*dmpk2*r4/term
     &             + (8.0d0/21.0d0)*dmpi24*dmpk2*r3/term
     &             + (12.0d0/7.0d0)*dmpi23*dmpk2*r2/term
     &             + 4.0d0*dmpi22*dmpk2*r/term
     &             + 4.0d0*dmpi2*dmpk2/term) * expi
         if (rorder .ge. 11) then
            r6 = r5 * r
            dmpi26 = dmpi25 * dmpi2
            dmpk26 = dmpk25 * dmpk2
            d5s = (dmpi2*dmpk25*r6/945.0d0
     &                + (2.0d0/189.0d0)*dmpi2*dmpk24*r5
     &                + dmpi2*dmpk23*r4/21.0d0
     &                + dmpi2*dmpk22*r3/9.0d0
     &                + dmpi2*dmpk2*r2/9.0d0
     &                - (4.0d0/945.0d0)*dmpi2*dmpk26*r5/term
     &                - (4.0d0/63.0d0)*dmpi2*dmpk25*r4/term
     &                - (4.0d0/9.0d0)*dmpi2*dmpk24*r3/term
     &                - (16.0d0/9.0d0)*dmpi2*dmpk23*r2/term
     &                - 4.0d0*dmpi2*dmpk22*r/term
     &                - 4.0d0*dmpi2*dmpk2/term) * expk
     &          + (dmpi25*dmpk2*r6/945.0d0
     &                + (2.0d0/189.0d0)*dmpi24*dmpk2*r5
     &                + dmpi23*dmpk2*r4/21.0d0
     &                + dmpi22*dmpk2*r3/9.0d0
     &                + dmpi2*dmpk2*r2/9.0d0
     &                + (4.0d0/945.0d0)*dmpi26*dmpk2*r5/term
     &                + (4.0d0/63.0d0)*dmpi25*dmpk2*r4/term
     &                + (4.0d0/9.0d0)*dmpi24*dmpk2*r3/term
     &                + (16.0d0/9.0d0)*dmpi23*dmpk2*r2/term
     &                + 4.0d0*dmpi22*dmpk2*r/term
     &                + 4.0d0*dmpi2*dmpk2/term) * expi
         end if
      end if
c
c     convert partial derivatives into full derivatives
c
      s = s * rr1
      ds = ds * rr3
      d2s = d2s * rr5
      d3s = d3s * rr7
      d4s = d4s * rr9
      d5s = d5s * rr11
      dmpik(1) = 0.5d0 * pre * s * s
      dmpik(3) = pre * s * ds
      dmpik(5) = pre * (s*d2s + ds*ds)
      dmpik(7) = pre * (s*d3s + 3.0d0*ds*d2s)
      dmpik(9) = pre * (s*d4s + 4.0d0*ds*d3s + 3.0d0*d2s*d2s)
      if (rorder .ge. 11) then
         dmpik(11) = pre * (s*d5s + 5.0d0*ds*d4s + 10.0d0*d2s*d3s)
      end if
      return
      end
