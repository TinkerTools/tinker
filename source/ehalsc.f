c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine halsc  --  soft core terms per ligand group  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "halsc" sets the soft core scale "vlsc" and softening "vscal"
c     of the buffered 14-7 potential for each ligand group, where a
c     second ligand group couples as the complement of vlambda
c
c
      subroutine halsc (vlsc,vscal)
      use mutant
      implicit none
      integer ig
      real*8 vlmd
      real*8 vlsc(2)
      real*8 vscal(2)
c
c
c     set the soft core lambda terms for each ligand group
c
      do ig = 1, 2
         vlmd = vlambda
         if (ig .eq. 2)  vlmd = 1.0d0 - vlambda
         vlsc(ig) = vlmd**scexp
         vscal(ig) = scalphav * (1.0d0-vlmd)**2
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine halsc4  --  soft core terms and lambda derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "halsc4" sets the soft core terms of the buffered 14-7 potential
c     for each ligand group along with the factors needed for their
c     first and second derivatives with respect to the group lambda
c
c
      subroutine halsc4 (vlsc,vscal,vlsc1,vlsc2,dvscal)
      use mutant
      implicit none
      integer ig
      real*8 vlmd
      real*8 vlsc(2)
      real*8 vscal(2)
      real*8 vlsc1(2)
      real*8 vlsc2(2)
      real*8 dvscal(2)
c
c
c     set the soft core terms and their lambda derivative factors
c
      call halsc (vlsc,vscal)
      do ig = 1, 2
         vlmd = vlambda
         if (ig .eq. 2)  vlmd = 1.0d0 - vlambda
         vlsc1(ig) = vlmd**(scexp-1)
         vlsc2(ig) = vlmd**(scexp-2)
         dvscal(ig) = 2.0d0 * scalphav * (1.0d0-vlmd)
      end do
      return
      end
c
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
c
      subroutine ehalsc4 (rik,rv,eps,ig,vlsc,vlsc1,vlsc2,vscal,dvscal,
     &                    e,de,dlambda,dlambda2,dlde)
      use dlmda
      use mutant
      use vdwpot
      implicit none
      integer ig
      real*8 rik,rv,eps
      real*8 e,de,dlambda
      real*8 dlambda2,dlde
      real*8 vsgn,epsl
      real*8 scal,dscal
      real*8 dhal17,ghal1
      real*8 rho,rho6,rho7
      real*8 rhopdhal,rhopdhal6,rhopdhal7
      real*8 s1,s2,t1,t2,t2m2
      real*8 dt1drho,dt2drho
      real*8 dt0dl,dt1dl,dt2dl
      real*8 ds1dl,ds2dl
      real*8 d2t0dl2,d2t1dl2,d2t2dl2
      real*8 d2t1dldrho,d2t2dldrho
      real*8 vlsc(2),vlsc1(2)
      real*8 vlsc2(2),vscal(2)
      real*8 dvscal(2)
c
c
c     the second ligand group couples as the complement of vlambda,
c     which flips the sign of its odd lambda derivatives
c
      vsgn = 1.0d0
      if (ig .eq. 2)  vsgn = -1.0d0
      dhal17 = (1.0d0+dhal)**7
      ghal1 = 1.0d0 + ghal
c
c     get the soft core energy and its distance derivative
c
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
c     get the first lambda derivative of the energy
c
      dt0dl = eps * scexp * vlsc1(ig)
      dscal = dvscal(ig)
      ds1dl = dscal * s1 * s1
      ds2dl = dscal * s2 * s2
      dt1dl = dhal17 * ds1dl
      dt2dl = ghal1 * ds2dl
      dlambda = dt0dl * t1 * t2m2
     &             + epsl * dt1dl * t2m2
     &             + epsl * t1 * dt2dl
      dlambda = vsgn * dlambda
c
c     get the second lambda derivative of the energy and the lambda
c     derivative of its distance derivative
c
      dlambda2 = 0.0d0
      dlde = 0.0d0
      if (use_d2lmda) then
         d2t0dl2 = eps*scexp*(scexp-1) * vlsc2(ig)
         d2t1dl2 = dhal17 * (-2.0d0*scalphav*s1*s1
     &                + 2.0d0*dscal*s1*ds1dl)
         d2t2dl2 = ghal1 * (-2.0d0*scalphav*s2*s2
     &                + 2.0d0*dscal*s2*ds2dl)
         dlambda2 = d2t0dl2*t1*t2m2
     &                 + epsl*d2t1dl2*t2m2
     &                 + epsl*t1*d2t2dl2
     &                 + 2.0d0*dt0dl*dt1dl*t2m2
     &                 + 2.0d0*dt0dl*t1*dt2dl
     &                 + 2.0d0*epsl*dt1dl*dt2dl
         d2t1dldrho = -14.0d0*dhal17*s1*ds1dl*rhopdhal6
         d2t2dldrho = -14.0d0*ghal1*s2*ds2dl*rho6
         dlde = scexp*vlsc1(ig)
     &             * (dt1drho*t2m2 + t1*dt2drho)
     &             + vlsc(ig)
     &             * (d2t1dldrho*t2m2 + t1*d2t2dldrho
     &             + dt1dl*dt2drho + dt1drho*dt2dl)
         dlde = eps / rv * dlde
         dlde = vsgn * dlde
      end if
      return
      end
