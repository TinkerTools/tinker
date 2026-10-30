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
c     contact "dlsq" of the damped dispersion for each ligand group,
c     which multiplies each mutated pair energy by
c     s = lambda^n * (1-S) + lambda^(n+q) * S
c     with switch S = (1-r^2/rs^2)^3 inside the radius rs, else zero,
c     so that beyond rs every mutated pair carries the lambda power
c     that also weights the reciprocal space sum, and where a second
c     ligand group couples as the complement of vlambda
c
c
      subroutine dspsc (dlsc,dlsq)
      use mutant
      implicit none
      integer ig
      real*8 vlmd,vsgn
      real*8 dlsc(2)
      real*8 dlsq(2)
c
c
c     set the soft core lambda terms for each ligand group
c
      do ig = 1, 2
         call grpvlmda (ig,vlmd,vsgn)
         dlsc(ig) = vlmd**scexpd
         dlsq(ig) = vlmd**(scexpd+scexpdq)
      end do
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
c     "dspsc4" sets the soft core terms of the damped dispersion for
c     each ligand group along with their first and second derivatives
c     with respect to vlambda, where the first derivatives of a second
c     ligand group change sign as it couples as one minus vlambda, and
c     the second derivatives are only set when requested
c
c
      subroutine dspsc4 (dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2)
      use dlmda
      use mutant
      implicit none
      integer ig
      real*8 vlmd,vsgn
      real*8 dlsc(2),dlsq(2)
      real*8 dlsc1(2),dlsq1(2)
      real*8 dlsc2(2),dlsq2(2)
c
c
c     set the soft core terms and their lambda derivatives
c
      do ig = 1, 2
         call grpvlmda (ig,vlmd,vsgn)
         call lmdapwr (vlmd,scexpd,use_d2lmda,dlsc(ig),dlsc1(ig),
     &                 dlsc2(ig))
         call lmdapwr (vlmd,scexpd+scexpdq,use_d2lmda,dlsq(ig),
     &                 dlsq1(ig),dlsq2(ig))
         dlsc1(ig) = vsgn * dlsc1(ig)
         dlsq1(ig) = vsgn * dlsq1(ig)
      end do
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
c     dispersion pair at distance "r" in ligand group "ig" and its
c     derivative "dsoft" with respect to distance
c
c
      subroutine edspsc (r,ig,dlsc,dlsq,soft,dsoft)
      use mutant
      implicit none
      integer ig
      real*8 r,soft,dsoft
      real*8 rs2,t,sw
      real*8 dlsc(2)
      real*8 dlsq(2)
c
c
c     get the soft core scaling and its distance derivative, which
c     reduce to the plain lambda power outside the switch radius
c
      soft = dlsc(ig)
      dsoft = 0.0d0
      rs2 = scrdsp * scrdsp
      if (r*r .lt. rs2) then
         t = 1.0d0 - r*r/rs2
         sw = t * t * t
         soft = dlsc(ig) - (dlsc(ig)-dlsq(ig))*sw
         dsoft = 6.0d0 * (dlsc(ig)-dlsq(ig)) * r * t * t / rs2
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
c     dispersion pair at distance "r" in ligand group "ig" and its
c     distance derivative "dsoft", along with the first lambda
c     derivative "dsdl" and, if second derivatives are requested, the
c     second lambda derivative "d2sdl2" and the lambda derivative
c     "d2sdldr" of "dsoft"
c
c
      subroutine edspsc4 (r,ig,dlsc,dlsq,dlsc1,dlsq1,dlsc2,dlsq2,
     &                    soft,dsoft,dsdl,d2sdl2,d2sdldr)
      use dlmda
      use mutant
      implicit none
      integer ig
      real*8 r
      real*8 lsc,lsq
      real*8 lsc1,lsq1
      real*8 lsc2,lsq2
      real*8 soft,dsoft,dsdl
      real*8 d2sdl2,d2sdldr
      real*8 rs2,t,sw
      real*8 dlsc(2),dlsq(2)
      real*8 dlsc1(2),dlsq1(2)
      real*8 dlsc2(2),dlsq2(2)
c
c
c     get the soft core terms of the ligand group, whose lambda
c     derivatives already carry the sign of the group
c
      lsc = dlsc(ig)
      lsq = dlsq(ig)
      lsc1 = dlsc1(ig)
      lsq1 = dlsq1(ig)
      lsc2 = dlsc2(ig)
      lsq2 = dlsq2(ig)
c
c     get the soft core scaling, its distance derivative, the lambda
c     derivatives of the scaling and the lambda derivative of its
c     distance derivative, with the same expressions as "edspsc" so
c     both routines give equal values
c
      soft = lsc
      dsoft = 0.0d0
      dsdl = lsc1
      d2sdl2 = 0.0d0
      d2sdldr = 0.0d0
      if (use_d2lmda)  d2sdl2 = lsc2
      rs2 = scrdsp * scrdsp
      if (r*r .lt. rs2) then
         t = 1.0d0 - r*r/rs2
         sw = t * t * t
         soft = lsc - (lsc-lsq)*sw
         dsoft = 6.0d0 * (lsc-lsq) * r * t * t / rs2
         dsdl = lsc1 - (lsc1-lsq1)*sw
         if (use_d2lmda) then
            d2sdl2 = lsc2 - (lsc2-lsq2)*sw
            d2sdldr = 6.0d0 * (lsc1-lsq1) * r * t * t / rs2
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine dhalwt  --  Ewald dispersion weight under hal  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "dhalwt" finds the weight "w" of the real space Ewald dispersion
c     pair term when a dispersion tail takes over from a
c     buffered 14-7 vdw across its switching window, as one minus the
c     vdw switch, along with its derivative with respect to distance
c
c
      subroutine dhalwt (r,w,dw)
      use dsppot
      implicit none
      real*8 r,w,dw
      real*8 taper,dtaper
      real*8 d2taper
c
c
c     the tail is absent inside the switch and whole beyond it, and
c     the switch is evaluated in reduced form to keep its precision
c     over a narrow window
c
      call quintaper (r,dhcut,dhoff,taper,dtaper,d2taper)
      w = 1.0d0 - taper
      dw = -dtaper
      return
      end
