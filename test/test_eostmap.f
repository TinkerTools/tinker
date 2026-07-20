c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine test_eostmap  --  EOST map tests  ##
c     ##                                               ##
c     ###################################################
c
c
c     "test_eostmap" checks lambda map and taper helper
c     routines used by orthogonal space tempering
c
c
      subroutine test_eostmap
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_eostmap')
c
c
      if (skiptest(tname,'ost'))  return
      call initial
      call test_eostmap_mapsub
      call test_eostmap_sublmda
      call test_eostmap_taper
      call test_eostmap_lmdachain
      call final
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_eostmap_mapsub  --  map tests  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "test_eostmap_mapsub" checks mapsublmda for exponential,
c     inverse-power and quintic lambda maps
c
c
      subroutine test_eostmap_mapsub
      use bound
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 tref,dtref,d2tref
c
c
c     test exponential sublambda maps and chain rule derivatives
c
      call resetost (5,5,1)
      ostlambda = 0.25d0
      ostpmap = 'EXP'
      ostemap = 'EXP'
      ostvmap = 'EXP'
      ostepexp = 2
      ostemexp = 3
      ostevexp = 4
      call mapsublmda
      call assert_real (plambda,0.0625d0,1.0d-12,
     &                  'mapsublmda exponential plambda')
      call assert_real (dpldlmda,0.5d0,1.0d-12,
     &                  'mapsublmda exponential dpldlmda')
      call assert_real (d2pldlmda2,2.0d0,1.0d-12,
     &                  'mapsublmda exponential d2pldlmda2')
      call assert_real (elambda,0.015625d0,1.0d-12,
     &                  'mapsublmda exponential elambda')
      call assert_real (deldlmda,0.1875d0,1.0d-12,
     &                  'mapsublmda exponential deldlmda')
      call assert_real (d2eldlmda2,1.5d0,1.0d-12,
     &                  'mapsublmda exponential d2eldlmda2')
      call assert_real (vlambda,0.00390625d0,1.0d-12,
     &                  'mapsublmda exponential vlambda')
      call assert_real (dvldlmda,0.0625d0,1.0d-12,
     &                  'mapsublmda exponential dvldlmda')
      call assert_real (d2vldlmda2,0.75d0,1.0d-12,
     &                  'mapsublmda exponential d2vldlmda2')
c
c     test shifted inverse-power sublambda maps and derivatives
c
      call resetost (5,5,1)
      ostlambda = 0.25d0
      ostpmap = 'INV'
      ostemap = 'INV'
      ostvmap = 'INV'
      ostinvepn = 2
      ostinvemn = 3
      ostinvevn = 4
      ostinvepeps = 0.01d0
      ostinvemeps = 0.02d0
      ostinveveps = 0.03d0
      call mapsublmda
      call assert_real (plambda,0.452936557937477d0,1.0d-12,
     &                  'mapsublmda invpower plambda')
      call assert_real (dpldlmda,1.08352945028593d0,1.0d-12,
     &                  'mapsublmda invpower dpldlmda')
      call assert_real (d2pldlmda2,-2.08371048131911d0,1.0d-12,
     &                  'mapsublmda invpower d2pldlmda2')
      call assert_real (elambda,0.509927040983045d0,1.0d-12,
     &                  'mapsublmda invpower elambda')
      call assert_real (deldlmda,1.08536378202464d0,1.0d-12,
     &                  'mapsublmda invpower deldlmda')
      call assert_real (d2eldlmda2,-2.67991057290036d0,1.0d-12,
     &                  'mapsublmda invpower d2eldlmda2')
      call assert_real (vlambda,0.526434441030379d0,1.0d-12,
     &                  'mapsublmda invpower vlambda')
      call assert_real (dvldlmda,1.09852311505228d0,1.0d-12,
     &                  'mapsublmda invpower dvldlmda')
      call assert_real (d2vldlmda2,-2.94247262960432d0,1.0d-12,
     &                  'mapsublmda invpower d2vldlmda2')
c
c     any map other than EXP or INV falls back to the quintic taper,
c     where the sublambda is the complement of the taper
c
      use_bounds = .false.
      call resetost (5,5,1)
      ostplmda0 = 0.2d0
      ostplmda1 = 0.8d0
      ostelmda0 = 0.3d0
      ostelmda1 = 0.7d0
      ostvlmda0 = 0.1d0
      ostvlmda1 = 0.9d0
      ostpmap = 'QNT'
      ostemap = 'QNT'
      ostvmap = 'QNT'
      ostlambda = 0.5d0
c
c     each sublambda uses its own window, so the midpoint value is
c     one half for all three but the slopes differ
c
      call mapsublmda
      call assert_real (plambda,0.5d0,1.0d-12,
     &                  'mapsublmda taper plambda')
      call assert_real (elambda,0.5d0,1.0d-12,
     &                  'mapsublmda taper elambda')
      call assert_real (vlambda,0.5d0,1.0d-12,
     &                  'mapsublmda taper vlambda')
      call assert_real (dpldlmda,3.125d0,1.0d-12,
     &                  'mapsublmda taper dpldlmda')
      call assert_real (deldlmda,4.6875d0,1.0d-12,
     &                  'mapsublmda taper deldlmda')
      call assert_real (dvldlmda,2.34375d0,1.0d-12,
     &                  'mapsublmda taper dvldlmda')
      call assert_real (d2pldlmda2,0.0d0,1.0d-12,
     &                  'mapsublmda taper d2pldlmda2')
      call assert_real (d2eldlmda2,0.0d0,1.0d-12,
     &                  'mapsublmda taper d2eldlmda2')
      call assert_real (d2vldlmda2,0.0d0,1.0d-12,
     &                  'mapsublmda taper d2vldlmda2')
c
c     the taper branch must negate the taper derivatives, checked
c     off center where the second derivative is nonzero
c
      ostlambda = 0.35d0
      call refquintic (0.35d0,0.2d0,0.8d0,tref,dtref,d2tref)
      call mapsublmda
      call assert_real (plambda,1.0d0-tref,1.0d-12,
     &                  'mapsublmda taper offcenter plambda')
      call assert_real (dpldlmda,-dtref,1.0d-12,
     &                  'mapsublmda taper offcenter dpldlmda')
      call assert_real (d2pldlmda2,-d2tref,1.0d-12,
     &                  'mapsublmda taper offcenter d2pldlmda2')
c
c     below the polarization window the sublambda is fully off and
c     above it the sublambda is fully on
c
      ostlambda = 0.1d0
      call mapsublmda
      call assert_real (plambda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window dpldlmda')
      ostlambda = 0.9d0
      call mapsublmda
      call assert_real (plambda,1.0d0,1.0d-12,
     &                  'mapsublmda taper above window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper above window dpldlmda')
c
c     a QNT polarization map sets the initial and final polarization
c     flags by comparing lambda against the polarization window
c
      ostpmap = 'QNT'
      ostlambda = 0.1d0
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i below window')
      call assert_logical (use_pol4f,.false.,
     &                     'mapsublmda qnt pol4f below window')
      ostlambda = 0.5d0
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i mid window')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f mid window')
      ostlambda = 0.9d0
      call mapsublmda
      call assert_logical (use_pol4i,.false.,
     &                     'mapsublmda qnt pol4i above window')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f above window')
c
c     anywhere inside the window both endpoint states still
c     contribute, since plambda has not yet pinned to zero or one
c
      ostlambda = 0.65d0
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i upper ramp')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f upper ramp')
      ostlambda = 0.30d0
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i lower ramp')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f lower ramp')
c
c     the window bounds are inclusive, so a lambda sitting exactly
c     on either edge must keep both endpoint states
c
      ostlambda = ostplmda1
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i at upper bound')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f at upper bound')
      ostlambda = ostplmda0
      call mapsublmda
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i at lower bound')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f at lower bound')
c
c     a non QNT polarization map must leave the flags untouched
c
      ostpmap = 'EXP'
      ostepexp = 2
      use_pol4i = .false.
      use_pol4f = .false.
      ostlambda = 0.5d0
      call mapsublmda
      call assert_logical (use_pol4i,.false.,
     &                     'mapsublmda exp leaves pol4i')
      call assert_logical (use_pol4f,.false.,
     &                     'mapsublmda exp leaves pol4f')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eostmap_sublmda  --  sublambda tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eostmap_sublmda" checks endpoint, clamp and
c     shift behavior for elementary sublambda maps
c
c
      subroutine test_eostmap_sublmda
      implicit none
      real*8 lmda,dlmda,d2lmda
      real*8 lmda0,dlmda0,d2lmda0
c
c
c     at and below x=0 the exponential map is zero, but the first
c     and second derivatives survive for exponents of one and two
c
      call sublmdaexp (0.0d0,1,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=1 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=1 d2lmda')
      call sublmdaexp (0.0d0,2,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=2 lmda')
      call assert_real (dlmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=2 dlmda')
      call assert_real (d2lmda,2.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=2 d2lmda')
      call sublmdaexp (0.0d0,3,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=3 lmda')
      call assert_real (dlmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=3 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x=0 n=3 d2lmda')
      call sublmdaexp (-0.5d0,2,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x below 0 lmda')
      call assert_real (d2lmda,2.0d0,1.0d-12,
     &                  'sublmdaexp x below 0 d2lmda')
c
c     at and above x=1 the exponential map saturates with the
c     analytic power law derivatives
c
      call sublmdaexp (1.0d0,3,lmda,dlmda,d2lmda)
      call assert_real (lmda,1.0d0,1.0d-12,
     &                  'sublmdaexp x=1 n=3 lmda')
      call assert_real (dlmda,3.0d0,1.0d-12,
     &                  'sublmdaexp x=1 n=3 dlmda')
      call assert_real (d2lmda,6.0d0,1.0d-12,
     &                  'sublmdaexp x=1 n=3 d2lmda')
      call sublmdaexp (1.5d0,1,lmda,dlmda,d2lmda)
      call assert_real (lmda,1.0d0,1.0d-12,
     &                  'sublmdaexp x above 1 n=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaexp x above 1 n=1 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp x above 1 n=1 d2lmda')
c
c     the interior exponential map is continuous with the x=1
c     endpoint for an exponent of one
c
      call sublmdaexp (0.5d0,1,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.5d0,1.0d-12,
     &                  'sublmdaexp interior n=1 lmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaexp interior n=1 d2lmda')
c
c     an inverse power of one or less is the identity map
c
      call sublmdainvpower (0.4d0,1,0.1d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.4d0,1.0d-12,
     &                  'sublmdainvpower n=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdainvpower n=1 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdainvpower n=1 d2lmda')
      call sublmdainvpower (0.4d0,0,0.1d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.4d0,1.0d-12,
     &                  'sublmdainvpower n=0 lmda')
c
c     the identity map still clamps out of range coordinates
c
      call sublmdainvpower (1.5d0,1,0.1d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,1.0d0,1.0d-12,
     &                  'sublmdainvpower n=1 clamp high')
      call sublmdainvpower (-0.5d0,1,0.1d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdainvpower n=1 clamp low')
c
c     out of range coordinates clamp to the endpoint values of the
c     shifted inverse power map, which is normalized to zero and one
c
      call sublmdainvpower (1.5d0,3,0.02d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,1.0d0,1.0d-12,
     &                  'sublmdainvpower clamp high')
      call sublmdainvpower (-0.5d0,3,0.02d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdainvpower clamp low')
c
c     a nonpositive shift falls back to a default shift of 0.1
c
      call sublmdainvpower (0.25d0,2,0.1d0,lmda0,dlmda0,d2lmda0)
      call sublmdainvpower (0.25d0,2,0.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,lmda0,1.0d-12,
     &                  'sublmdainvpower zero eps lmda')
      call assert_real (dlmda,dlmda0,1.0d-12,
     &                  'sublmdainvpower zero eps dlmda')
      call assert_real (d2lmda,d2lmda0,1.0d-12,
     &                  'sublmdainvpower zero eps d2lmda')
      call sublmdainvpower (0.25d0,2,-1.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,lmda0,1.0d-12,
     &                  'sublmdainvpower negative eps lmda')
      call assert_real (dlmda,dlmda0,1.0d-12,
     &                  'sublmdainvpower negative eps dlmda')
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_eostmap_taper  --  taper tests  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "test_eostmap_taper" checks quintic taper values,
c     derivatives and per-mode lambda windows
c
c
      subroutine test_eostmap_taper
      use bound
      use ost
      implicit none
      real*8 taper,dtaper,d2taper
      real*8 tref,dtref,d2tref
      character*6 mode
c
c
c     avoid the periodic image search inside the switch routine
c
      use_bounds = .false.
c
c     use a distinct taper window for each sublambda so a mode
c     that reads the wrong window cannot pass
c
      ostplmda0 = 0.2d0
      ostplmda1 = 0.8d0
      ostelmda0 = 0.3d0
      ostelmda1 = 0.7d0
      ostvlmda0 = 0.1d0
      ostvlmda1 = 0.9d0
      mode = 'OSTPOL'
c
c     below and at the lower bound the taper is fully on and flat
c
      call sublmdataper (mode,0.05d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'sublmdataper below cut taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'sublmdataper below cut dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'sublmdataper below cut d2taper')
      call sublmdataper (mode,0.2d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'sublmdataper at cut taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'sublmdataper at cut dtaper')
c
c     above and at the upper bound the taper is fully off and flat
c
      call sublmdataper (mode,0.95d0,taper,dtaper,d2taper)
      call assert_real (taper,0.0d0,1.0d-12,
     &                  'sublmdataper above off taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'sublmdataper above off dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'sublmdataper above off d2taper')
      call sublmdataper (mode,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.0d0,1.0d-12,
     &                  'sublmdataper at off taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'sublmdataper at off dtaper')
c
c     the quintic is the unique polynomial that is one at cut and
c     zero at off with vanishing first and second derivatives at
c     both ends, so it must equal the analytic smoothstep form
c
      call sublmdataper (mode,0.5d0,taper,dtaper,d2taper)
      call refquintic (0.5d0,0.2d0,0.8d0,tref,dtref,d2tref)
      call assert_real (taper,0.5d0,1.0d-12,
     &                  'sublmdataper midpoint taper')
      call assert_real (taper,tref,1.0d-12,
     &                  'sublmdataper midpoint reference')
      call assert_real (dtaper,dtref,1.0d-12,
     &                  'sublmdataper midpoint dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'sublmdataper midpoint d2taper')
      call sublmdataper (mode,0.35d0,taper,dtaper,d2taper)
      call refquintic (0.35d0,0.2d0,0.8d0,tref,dtref,d2tref)
      call assert_real (taper,tref,1.0d-12,
     &                  'sublmdataper offcenter taper')
      call assert_real (dtaper,dtref,1.0d-12,
     &                  'sublmdataper offcenter dtaper')
      call assert_real (d2taper,d2tref,1.0d-12,
     &                  'sublmdataper offcenter d2taper')
      call sublmdataper (mode,0.70d0,taper,dtaper,d2taper)
      call refquintic (0.70d0,0.2d0,0.8d0,tref,dtref,d2tref)
      call assert_real (taper,tref,1.0d-12,
     &                  'sublmdataper upper half taper')
      call assert_real (dtaper,dtref,1.0d-12,
     &                  'sublmdataper upper half dtaper')
      call assert_real (d2taper,d2tref,1.0d-12,
     &                  'sublmdataper upper half d2taper')
c
c     the same lambda gives different results per mode because each
c     mode selects its own taper window
c
      mode = 'OSTELE'
      call sublmdataper (mode,0.25d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'sublmdataper ostele below own cut')
      mode = 'OSTVDW'
      call sublmdataper (mode,0.25d0,taper,dtaper,d2taper)
      call refquintic (0.25d0,0.1d0,0.9d0,tref,dtref,d2tref)
      call assert_real (taper,tref,1.0d-12,
     &                  'sublmdataper ostvdw own window taper')
      call assert_real (dtaper,dtref,1.0d-12,
     &                  'sublmdataper ostvdw own window dtaper')
      mode = 'OSTPOL'
      call sublmdataper (mode,0.25d0,taper,dtaper,d2taper)
      call refquintic (0.25d0,0.2d0,0.8d0,tref,dtref,d2tref)
      call assert_real (taper,tref,1.0d-12,
     &                  'sublmdataper ostpol own window taper')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine refquintic  --  analytic taper reference  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "refquintic" evaluates the quintic taper and its
c     first two derivatives from the smoothstep formula
c
c
      subroutine refquintic (x,cutval,offval,taper,dtaper,d2taper)
      implicit none
      real*8 x,cutval,offval
      real*8 taper,dtaper,d2taper
      real*8 t,t2,t3,t4,t5
      real*8 w
c
c
      w = offval - cutval
      t = (x-cutval) / w
      t2 = t * t
      t3 = t2 * t
      t4 = t2 * t2
      t5 = t2 * t3
      taper = 1.0d0 - (6.0d0*t5-15.0d0*t4+10.0d0*t3)
      dtaper = -(30.0d0*t4-60.0d0*t3+30.0d0*t2) / w
      d2taper = -(120.0d0*t3-180.0d0*t2+60.0d0*t) / (w*w)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eostmap_lmdachain  --  chain tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eostmap_lmdachain" checks lambda chain-rule
c     scaling for energies, forces and virials
c
c
      subroutine test_eostmap_lmdachain
      use atoms
      use dlmda
      implicit none
      integer i,j
      real*8 fpref(3,2)
      real*8 fmref(3,2)
      real*8 fvref(3,2)
      real*8 vpref(3,3)
      real*8 vmref(3,3)
      real*8 vvref(3,3)
c
c
c     use distinct sublambda derivatives so a swapped polarization,
c     multipole or van der Waals factor cannot pass
c
      n = 2
      if (allocated(dfpdl))  deallocate (dfpdl)
      if (allocated(dfmdl))  deallocate (dfmdl)
      if (allocated(dfvdl))  deallocate (dfvdl)
      allocate (dfpdl(3,n))
      allocate (dfmdl(3,n))
      allocate (dfvdl(3,n))
      dpldlmda = 2.0d0
      d2pldlmda2 = 3.0d0
      deldlmda = 11.0d0
      d2eldlmda2 = 13.0d0
      dvldlmda = 5.0d0
      d2vldlmda2 = 7.0d0
c
c     set the sublambda energy derivatives and expected results
c
      depdl = 1.5d0
      d2epdl2 = 0.5d0
      demdl = 0.25d0
      d2emdl2 = 4.0d0
      devdl = 2.5d0
      d2evdl2 = 1.5d0
c
c     set force and virial derivatives that vary by component
c
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = dble(j+3*(i-1))
            dfmdl(j,i) = dble(j+3*(i-1)) + 0.5d0
            dfvdl(j,i) = dble(j+3*(i-1)) - 0.5d0
            fpref(j,i) = dfpdl(j,i) * dpldlmda
            fmref(j,i) = dfmdl(j,i) * deldlmda
            fvref(j,i) = dfvdl(j,i) * dvldlmda
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = dble(j+3*(i-1))
            demvirdl(j,i) = dble(j+3*(i-1)) + 0.5d0
            devvirdl(j,i) = dble(j+3*(i-1)) - 0.5d0
            vpref(j,i) = depvirdl(j,i) * dpldlmda
            vmref(j,i) = demvirdl(j,i) * deldlmda
            vvref(j,i) = devvirdl(j,i) * dvldlmda
         end do
      end do
      call lmdachain
c
c     the second derivatives must use the old first derivatives,
c     so d2epdl2 = 0.5*2^2 + 1.5*3 and depdl = 1.5*2
c
      call assert_real (depdl,3.0d0,1.0d-12,
     &                  'lmdachain depdl')
      call assert_real (d2epdl2,6.5d0,1.0d-12,
     &                  'lmdachain d2epdl2')
      call assert_real (demdl,2.75d0,1.0d-12,
     &                  'lmdachain demdl')
      call assert_real (d2emdl2,487.25d0,1.0d-12,
     &                  'lmdachain d2emdl2')
      call assert_real (devdl,12.5d0,1.0d-12,
     &                  'lmdachain devdl')
      call assert_real (d2evdl2,55.0d0,1.0d-12,
     &                  'lmdachain d2evdl2')
c
c     forces and virials scale by the first derivative only
c
      call assert_array2 (dfpdl,fpref,3,n,1.0d-12,
     &                    'lmdachain dfpdl')
      call assert_array2 (dfmdl,fmref,3,n,1.0d-12,
     &                    'lmdachain dfmdl')
      call assert_array2 (dfvdl,fvref,3,n,1.0d-12,
     &                    'lmdachain dfvdl')
      call assert_array2 (depvirdl,vpref,3,3,1.0d-12,
     &                    'lmdachain depvirdl')
      call assert_array2 (demvirdl,vmref,3,3,1.0d-12,
     &                    'lmdachain demvirdl')
      call assert_array2 (devvirdl,vvref,3,3,1.0d-12,
     &                    'lmdachain devvirdl')
c
c     an identity sublambda map must leave everything unchanged
c
      dpldlmda = 1.0d0
      d2pldlmda2 = 0.0d0
      deldlmda = 1.0d0
      d2eldlmda2 = 0.0d0
      dvldlmda = 1.0d0
      d2vldlmda2 = 0.0d0
      depdl = 1.5d0
      d2epdl2 = 0.5d0
      call lmdachain
      call assert_real (depdl,1.5d0,1.0d-12,
     &                  'lmdachain identity depdl')
      call assert_real (d2epdl2,0.5d0,1.0d-12,
     &                  'lmdachain identity d2epdl2')
      deallocate (dfpdl)
      deallocate (dfmdl)
      deallocate (dfvdl)
      return
      end
