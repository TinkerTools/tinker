c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine erepsc  --  soft core Pauli repulsion scaling  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erepsc" finds the soft core scaling "soft" of a mutated Pauli
c     repulsion pair at distance "r" in ligand group "ig" and its
c     derivative "dsoft" with respect to distance
c
c
      subroutine erepsc (r,ig,rlsc,rshf,soft,dsoft)
      implicit none
      integer ig
      real*8 r,soft,dsoft
      real*8 q,q3
      real*8 rlsc(2)
      real*8 rshf(2)
c
c
c     get the soft core scaling and its distance derivative
c
      q = 1.0d0 / sqrt(rshf(ig)+r*r)
      q3 = q * q * q
      soft = rlsc(ig) * r * q
      dsoft = rlsc(ig) * rshf(ig) * q3
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine erepsc4  --  soft core repulsion lambda derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "erepsc4" finds the soft core scaling "soft" of a mutated Pauli
c     repulsion pair at distance "r" in ligand group "ig" and its
c     distance derivative "dsoft", along with the first lambda
c     derivative "dsdl" and, if second derivatives are requested, the
c     second lambda derivative "d2sdl2" and the lambda derivative
c     "d2sdldr" of "dsoft"
c
c
      subroutine erepsc4 (r,ig,rlsc,rshf,rlsc1,rlsc2,rshf1,rshf2,
     &                    soft,dsoft,dsdl,d2sdl2,d2sdldr)
      use dlmda
      use mutant
      implicit none
      integer ig
      real*8 r
      real*8 lsc,shf
      real*8 lsc1,lsc2
      real*8 shf1,shf2
      real*8 soft,dsoft,dsdl
      real*8 d2sdl2,d2sdldr
      real*8 q,q2,q3
      real*8 rlsc(2),rshf(2)
      real*8 rlsc1(2),rlsc2(2)
      real*8 rshf1(2),rshf2(2)
c
c
c     get the soft core terms of the ligand group, whose scaled
c     lambda derivative already carries the sign of the group
c
      lsc = rlsc(ig)
      shf = rshf(ig)
      lsc1 = rlsc1(ig)
      shf1 = rshf1(ig)
c
c     get the soft core scaling and its distance derivative with the
c     same expressions as "erepsc", so both give equal values
c
      q = 1.0d0 / sqrt(shf+r*r)
      q2 = q * q
      q3 = q2 * q
      soft = lsc * r * q
      dsoft = lsc * shf * q3
c
c     get the first lambda derivative of the scaling
c
      dsdl = lsc1 * r * q * (scexpr-0.5d0*shf1*q2)
c
c     get the second lambda derivative of the scaling and the lambda
c     derivative of its distance derivative
c
      d2sdl2 = 0.0d0
      d2sdldr = 0.0d0
      if (use_d2lmda) then
         lsc2 = rlsc2(ig)
         shf2 = rshf2(ig)
         d2sdl2 = lsc2 * r * q * (scexpr*(scexpr-1.0d0)
     &               - scexpr*shf1*q2 - 0.5d0*shf2*q2
     &               + 0.75d0*shf1*shf1*q2*q2)
         d2sdldr = lsc1 * q3 * (scexpr*shf + shf1
     &                - 1.5d0*shf*shf1*q2)
      end if
      return
      end
