c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine repsc  --  soft core Pauli repulsion lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "repsc" sets the soft core scale "rlsc" and shift "rshf" of the
c     Pauli repulsion, which multiplies each mutated pair energy by
c     s = lambda^n * r / sqrt(alpha*lambda^m*(1-lambda)^p + r^2)
c
c
      subroutine repsc (rlsc,rshf)
      use mutant
      implicit none
      real*8 rlsc,rshf
c
c
c     set the soft core lambda terms of the Pauli repulsion
c
      rlsc = vlambda**scexpr
      rshf = scalphar * vlambda**scexprm * (1.0d0-vlambda)**scexprp
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine repsc4  --  soft core repulsion lambda derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "repsc4" sets the soft core terms of the Pauli repulsion along
c     with the factors needed for their first and second derivatives
c     with respect to lambda, where "rshf1" and "rshf2" hold lambda
c     and lambda squared times the shift derivatives so that no power
c     of lambda goes negative, and the second derivative factors are
c     only set when requested
c
c
      subroutine repsc4 (rlsc,rshf,rlsc1,rlsc2,rshf1,rshf2)
      use dlmda
      use mutant
      implicit none
      real*8 rlsc,rshf
      real*8 rlsc1,rlsc2
      real*8 rshf1,rshf2
      real*8 vlm,omvl
      real*8 alm,alm1,alm2
c
c
c     set the soft core terms and their lambda derivative factors,
c     skipping any shift term whose coefficient vanishes, so that no
c     power of one minus lambda goes negative
c
      call repsc (rlsc,rshf)
      vlm = vlambda
      omvl = 1.0d0 - vlambda
      alm = scalphar * vlm**scexprm
      alm1 = alm * vlm
      alm2 = alm1 * vlm
      rlsc1 = vlm**(scexpr-1.0d0)
      rshf1 = 0.0d0
      if (scexprm .ne. 0.0d0) then
         rshf1 = scexprm * alm * omvl**scexprp
      end if
      rshf1 = rshf1 - scexprp*alm1*omvl**(scexprp-1.0d0)
      rlsc2 = 0.0d0
      rshf2 = 0.0d0
      if (use_d2lmda) then
         rlsc2 = vlm**(scexpr-2.0d0)
         if (scexprm .ne. 0.0d0) then
            rshf2 = scexprm*(scexprm-1.0d0) * alm * omvl**scexprp
     &                 - 2.0d0*scexprm*scexprp * alm1
     &                      * omvl**(scexprp-1.0d0)
         end if
         if (scexprp .ne. 1.0d0) then
            rshf2 = rshf2 + scexprp*(scexprp-1.0d0) * alm2
     &                         * omvl**(scexprp-2.0d0)
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine erepsc  --  soft core Pauli repulsion scaling  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "erepsc" finds the soft core scaling "soft" of a mutated Pauli
c     repulsion pair at distance "r" and its derivative "dsoft" with
c     respect to distance
c
c
      subroutine erepsc (r,rlsc,rshf,soft,dsoft)
      implicit none
      real*8 r,rlsc,rshf
      real*8 soft,dsoft
      real*8 q,q3
c
c
c     get the soft core scaling and its distance derivative
c
      q = 1.0d0 / sqrt(rshf+r*r)
      q3 = q * q * q
      soft = rlsc * r * q
      dsoft = rlsc * rshf * q3
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
c     repulsion pair at distance "r" and its distance derivative
c     "dsoft", along with the first lambda derivative "dsdl" and, if
c     second derivatives are requested, the second lambda derivative
c     "d2sdl2" and the lambda derivative "d2sdldr" of "dsoft"
c
c
      subroutine erepsc4 (r,rlsc,rshf,rlsc1,rlsc2,rshf1,rshf2,
     &                    soft,dsoft,dsdl,d2sdl2,d2sdldr)
      use dlmda
      use mutant
      implicit none
      real*8 r,rlsc,rshf
      real*8 rlsc1,rlsc2
      real*8 rshf1,rshf2
      real*8 soft,dsoft,dsdl
      real*8 d2sdl2,d2sdldr
      real*8 q,q2,q3
c
c
c     get the soft core scaling and its distance derivative with the
c     same expressions as "erepsc", so both give equal values
c
      q = 1.0d0 / sqrt(rshf+r*r)
      q2 = q * q
      q3 = q2 * q
      soft = rlsc * r * q
      dsoft = rlsc * rshf * q3
c
c     get the first lambda derivative of the scaling
c
      dsdl = rlsc1 * r * q * (scexpr-0.5d0*rshf1*q2)
c
c     get the second lambda derivative of the scaling and the lambda
c     derivative of its distance derivative
c
      d2sdl2 = 0.0d0
      d2sdldr = 0.0d0
      if (use_d2lmda) then
         d2sdl2 = rlsc2 * r * q * (scexpr*(scexpr-1.0d0)
     &               - scexpr*rshf1*q2 - 0.5d0*rshf2*q2
     &               + 0.75d0*rshf1*rshf1*q2*q2)
         d2sdldr = rlsc1 * q3 * (scexpr*rshf + rshf1
     &                - 1.5d0*rshf*rshf1*q2)
      end if
      return
      end
