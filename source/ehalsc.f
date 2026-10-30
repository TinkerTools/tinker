c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine ehalsc  --  soft core buffered 14-7 pair  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "ehalsc" finds the soft core buffered 14-7 energy "e" and its
c     derivative "de" with respect to distance for a pair at distance
c     "rik" with minimum energy radius "rv" and well depth "eps", where
c     the pair takes the soft core terms of ligand group "ig"
c
c
      subroutine ehalsc (rik,rv,eps,ig,vlsc,vscal,e,de)
      use vdwpot
      implicit none
      integer ig
      real*8 rik,rv,eps
      real*8 e,de
      real*8 epsl,scal
      real*8 rho,rho6,rho7
      real*8 s1,s2,t1,t2
      real*8 dt1drho,dt2drho
      real*8 vlsc(2)
      real*8 vscal(2)
c
c
c     get the soft core energy and its distance derivative
c
      rho = rik / rv
      rho6 = rho**6
      rho7 = rho6 * rho
      epsl = eps * vlsc(ig)
      scal = vscal(ig)
      s1 = 1.0d0 / (scal+(rho+dhal)**7)
      s2 = 1.0d0 / (scal+rho7+ghal)
      t1 = (1.0d0+dhal)**7 * s1
      t2 = (1.0d0+ghal) * s2
      dt1drho = -7.0d0*(rho+dhal)**6 * t1 * s1
      dt2drho = -7.0d0*rho6 * t2 * s2
      e = epsl * t1 * (t2-2.0d0)
      de = epsl * (dt1drho*(t2-2.0d0)+t1*dt2drho) / rv
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ehalsc4  --  soft core 14-7 pair lambda derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ehalsc4" finds the soft core buffered 14-7 energy "e" and its
c     distance derivative "de" for a pair at distance "rik" with radius
c     "rv" and well depth "eps" in ligand group "ig", along with the
c     first lambda derivative "dlambda" and, if second derivatives are
c     requested, the second lambda derivative "dlambda2" and the lambda
c     derivative "dlde" of the distance derivative
c
c     the pair energy is eps * lambda^n * T(scal) with T = t1*(t2-2),
c     so each lambda derivative is a power of lambda times derivatives
c     of T with respect to the softening "scal", weighted by the scaled
c     softening derivatives from "grpsoft4", and no power of lambda goes
c     negative
c
c
      subroutine ehalsc4 (rik,rv,eps,ig,vlsc,vlsc1,vlsc2,vscal,vscal1,
     &                    vscal2,e,de,dlambda,dlambda2,dlde)
      use dlmda
      use mutant
      use vdwpot
      implicit none
      integer ig
      real*8 rik,rv,eps
      real*8 e,de,dlambda
      real*8 dlambda2,dlde
      real*8 epsl,scal
      real*8 sc1,sc2
      real*8 dhal17,ghal1
      real*8 rho,rho6,rho7
      real*8 rhopdhal,rhopdhal6,rhopdhal7
      real*8 s1,s2,t1,t2,t2m2
      real*8 dt1drho,dt2drho
      real*8 tt,ts,tss
      real*8 tr,trs
      real*8 vlsc(2),vlsc1(2)
      real*8 vlsc2(2),vscal(2)
      real*8 vscal1(2),vscal2(2)
c
c
c     get the soft core energy and its distance derivative, where
c     the scaled lambda derivative from "grpsoft4" already carries
c     the sign of the ligand group
c
      dhal17 = (1.0d0+dhal)**7
      ghal1 = 1.0d0 + ghal
      rho = rik / rv
      rho6 = rho**6
      rho7 = rho6 * rho
      rhopdhal = rho + dhal
      rhopdhal6 = rhopdhal**6
      rhopdhal7 = rhopdhal6 * rhopdhal
      epsl = eps * vlsc(ig)
      scal = vscal(ig)
      s1 = 1.0d0 / (scal+rhopdhal7)
      s2 = 1.0d0 / (scal+rho7+ghal)
      t1 = dhal17 * s1
      t2 = ghal1 * s2
      t2m2 = t2 - 2.0d0
      dt1drho = -7.0d0*rhopdhal6 * t1 * s1
      dt2drho = -7.0d0*rho6 * t2 * s2
      e = epsl * t1 * t2m2
      de = epsl * (dt1drho*t2m2+t1*dt2drho) / rv
c
c     get T and its derivative with respect to the softening, where
c     t1 and t2 change as minus t1*s1 and minus t2*s2
c
      sc1 = vscal1(ig)
      tt = t1 * t2m2
      ts = -t1*s1*t2m2 - t1*t2*s2
c
c     get the first lambda derivative of the energy
c
      dlambda = eps * vlsc1(ig) * (scexp*tt+sc1*ts)
c
c     get the second lambda derivative, then the lambda derivative of
c     the distance derivative from the distance derivative of T and its
c     own softening derivative
c
      dlambda2 = 0.0d0
      dlde = 0.0d0
      if (use_d2lmda) then
         sc2 = vscal2(ig)
         tss = 2.0d0*t1*s1*s1*t2m2 + 2.0d0*t1*t2*s1*s2
     &            + 2.0d0*t1*t2*s2*s2
         dlambda2 = eps * vlsc2(ig) * (scexp*(scexp-1.0d0)*tt
     &                 + 2.0d0*scexp*sc1*ts + sc1*sc1*tss + sc2*ts)
         tr = dt1drho*t2m2 + t1*dt2drho
         trs = -dt1drho*(2.0d0*s1*t2m2+t2*s2)
     &            - t1*dt2drho*(s1+2.0d0*s2)
         dlde = eps * vlsc1(ig) * (scexp*tr+sc1*trs) / rv
      end if
      return
      end
