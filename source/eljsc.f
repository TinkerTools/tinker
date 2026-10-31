c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine eljsc  --  soft core Lennard-Jones pair  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "eljsc" finds the soft core Lennard-Jones energy "e" and, if
c     "dograd" is set, its derivative "de" with respect to distance
c     for a pair at distance "rik", whose square is "rik2", with
c     minimum energy radius "rv" and well depth "eps", where the pair
c     takes the soft core scale "lsc" and shift "shf" of ligand group
c     "ig"; energy only callers may leave "rik" unset
c
c     the pair energy is eps * lambda^n * T(s) with T = 4*(1-s)/s^2
c     and s = 2*(rik/rv)^6 plus the shift, which is the form of
c     Beutler et al. for the default linear scale and shift
c
c     literature reference:
c
c     T. C. Beutler, A. E. Mark, R. C. van Schaik, P. R. Gerber and
c     W. F. van Gunsteren, "Avoiding Singularities and Numerical
c     Instabilities in Free Energy Calculations Based on Molecular
c     Simulations", Chemical Physics Letters, 222, 529-539 (1994)
c
c
      subroutine eljsc (rik,rik2,rv,eps,ig,lsc,shf,dograd,e,de)
      implicit none
      integer ig
      real*8 rik,rik2,rv,eps
      real*8 e,de
      real*8 p6,sc
      real*8 term,dterm
      real*8 lsc(2)
      real*8 shf(2)
      logical dograd
c
c
c     get the soft core energy and its distance derivative; atoms
c     at the same point have no direction and take zero force
c
      p6 = 2.0d0 * rik2**3 / rv**6
      sc = p6 + shf(ig)
      term = 4.0d0 * lsc(ig) * eps / (sc*sc)
      e = term * (1.0d0-sc)
      de = 0.0d0
      if (dograd .and. rik.ne.0.0d0) then
         dterm = -6.0d0 * p6 * term / rik
         de = dterm * (1.0d0+2.0d0*(1.0d0-sc)/sc)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine eljsc4  --  soft core LJ pair lambda derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "eljsc4" finds the soft core Lennard-Jones energy "e" and its
c     distance derivative "de" for a pair at distance "rik", whose
c     square is "rik2", with radius "rv" and well depth "eps" in ligand
c     group "ig", along with the first lambda derivative "dlambda" and,
c     if second derivatives are requested, the second lambda derivative
c     "dlambda2" and the lambda derivative "dlde" of the distance
c     derivative
c
c     with the pair energy eps * L * T(s), the scale L and the shift
c     h in s change with lambda as given by "grpsoftd", so that
c
c        dE/dl   = eps * (L'T + L T'h')
c        d2E/dl2 = eps * (L''T + 2L'T'h' + L(T''h'^2 + T'h''))
c
c     and the distance derivative eps * L * T' * ds/dr changes with
c     lambda as eps * (L'T' + L T''h') * ds/dr
c
c
      subroutine eljsc4 (rik,rik2,rv,eps,ig,lsc,dlsc,d2lsc,shf,dshf,
     &                   d2shf,e,de,dlambda,dlambda2,dlde)
      use dlmda
      implicit none
      integer ig
      real*8 rik,rik2,rv,eps
      real*8 e,de,dlambda
      real*8 dlambda2,dlde
      real*8 p6,sc,sc2
      real*8 dsdr
      real*8 l,dl,d2l
      real*8 h1,h2
      real*8 tt,ts,tss
      real*8 lsc(2),dlsc(2)
      real*8 d2lsc(2),shf(2)
      real*8 dshf(2),d2shf(2)
c
c
c     get the soft core energy and its distance derivative, then
c     the soft core variable and its distance derivative
c
      call eljsc (rik,rik2,rv,eps,ig,lsc,shf,.true.,e,de)
      p6 = 2.0d0 * rik2**3 / rv**6
      sc = p6 + shf(ig)
      dsdr = 0.0d0
      if (rik .ne. 0.0d0)  dsdr = 6.0d0 * p6 / rik
c
c     get T and its first two derivatives with respect to s
c
      sc2 = sc * sc
      tt = 4.0d0 * (1.0d0-sc) / sc2
      ts = 4.0d0 * (sc-2.0d0) / (sc2*sc)
      tss = 8.0d0 * (3.0d0-sc) / (sc2*sc2)
c
c     get the first lambda derivative of the energy
c
      l = lsc(ig)
      dl = dlsc(ig)
      h1 = dshf(ig)
      dlambda = eps * (dl*tt+l*ts*h1)
c
c     get the second lambda derivative and the lambda derivative of
c     the distance derivative
c
      dlambda2 = 0.0d0
      dlde = 0.0d0
      if (use_d2lmda) then
         d2l = d2lsc(ig)
         h2 = d2shf(ig)
         dlambda2 = eps * (d2l*tt + 2.0d0*dl*ts*h1
     &                        + l*(tss*h1*h1+ts*h2))
         dlde = eps * (dl*ts+l*tss*h1) * dsdr
      end if
      return
      end
