c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dspsc  --  soft core damped dispersion lambda  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dspsc" sets the lambda power "dlsc" and the lambda power at
c     contact "dlsq" of the damped dispersion, which multiplies each
c     mutated pair energy by s = lambda^n * (1-S) + lambda^(n+q) * S
c     with switch S = (1-r^2/rs^2)^3 inside the radius rs, else zero,
c     so that beyond rs every mutated pair carries the lambda power
c     that also weights the reciprocal space sum
c
c
      subroutine dspsc (dlsc,dlsq)
      use mutant
      implicit none
      real*8 dlsc,dlsq
c
c
c     set the soft core lambda terms of the damped dispersion
c
      dlsc = vlambda**scexpd
      dlsq = vlambda**(scexpd+scexpdq)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine dspsc4  --  soft core dispersion lambda derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dspsc4" sets the soft core terms of the damped dispersion along
c     with their first and second derivatives with respect to lambda,
c     where the second derivatives are only set when requested
c
c
      subroutine dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
      use dlmda
      use mutant
      implicit none
      real*8 dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
c
c
c     set the soft core terms and their lambda derivatives
c
      call lmdapwr (vlambda,scexpd,use_d2lmda,dlsc,dlsc1,dlsc2)
      call lmdapwr (vlambda,scexpd+scexpdq,use_d2lmda,dlsq,dlsq1,dlsq2)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edspsc  --  soft core damped dispersion scaling  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edspsc" finds the soft core scaling "soft" of a mutated damped
c     dispersion pair at distance "r" and its derivative "dsoft" with
c     respect to distance
c
c
      subroutine edspsc (r,dlsc,dlsq,soft,dsoft)
      use mutant
      implicit none
      real*8 r,dlsc,dlsq
      real*8 soft,dsoft
      real*8 rs2,t,sw
c
c
c     get the soft core scaling and its distance derivative, which
c     reduce to the plain lambda power outside the switch radius
c
      soft = dlsc
      dsoft = 0.0d0
      rs2 = scrdsp * scrdsp
      if (r*r .lt. rs2) then
         t = 1.0d0 - r*r/rs2
         sw = t * t * t
         soft = dlsc - (dlsc-dlsq)*sw
         dsoft = 6.0d0 * (dlsc-dlsq) * r * t * t / rs2
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine edspsc4  --  soft core dispersion lambda derivs  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "edspsc4" finds the soft core scaling "soft" of a mutated damped
c     dispersion pair at distance "r" and its distance derivative
c     "dsoft", along with the first lambda derivative "dsdl" and, if
c     second derivatives are requested, the second lambda derivative
c     "d2sdl2" and the lambda derivative "d2sdldr" of "dsoft"
c
c
      subroutine edspsc4 (r,dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2,
     &                    soft,dsoft,dsdl,d2sdl2,d2sdldr)
      use dlmda
      use mutant
      implicit none
      real*8 r,dlsc,dlsq
      real*8 dlsc1,dlsq1
      real*8 dlsc2,dlsq2
      real*8 soft,dsoft,dsdl
      real*8 d2sdl2,d2sdldr
      real*8 rs2,t,sw
c
c
c     get the soft core scaling, its distance derivative, the lambda
c     derivatives of the scaling and the lambda derivative of its
c     distance derivative, with the same expressions as "edspsc" so
c     both routines give equal values
c
      soft = dlsc
      dsoft = 0.0d0
      dsdl = dlsc1
      d2sdl2 = 0.0d0
      d2sdldr = 0.0d0
      if (use_d2lmda)  d2sdl2 = dlsc2
      rs2 = scrdsp * scrdsp
      if (r*r .lt. rs2) then
         t = 1.0d0 - r*r/rs2
         sw = t * t * t
         soft = dlsc - (dlsc-dlsq)*sw
         dsoft = 6.0d0 * (dlsc-dlsq) * r * t * t / rs2
         dsdl = dlsc1 - (dlsc1-dlsq1)*sw
         if (use_d2lmda) then
            d2sdl2 = dlsc2 - (dlsc2-dlsq2)*sw
            d2sdldr = 6.0d0 * (dlsc1-dlsq1) * r * t * t / rs2
         end if
      end if
      return
      end
