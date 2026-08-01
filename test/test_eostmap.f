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
      call test_eostmap_relstage
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
      plmdamap = 'EXP'
      elmdamap = 'EXP'
      vlmdamap = 'EXP'
      plmdaexp = 2
      elmdaexp = 3
      vlmdaexp = 4
      call mapsublmda (ostlambda)
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
      plmdamap = 'INV'
      elmdamap = 'INV'
      vlmdamap = 'INV'
      plmdainvn = 2
      elmdainvn = 3
      vlmdainvn = 4
      plmdainveps = 0.01d0
      elmdainveps = 0.02d0
      vlmdainveps = 0.03d0
      call mapsublmda (ostlambda)
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
      qntplmda0 = 0.2d0
      qntplmda1 = 0.8d0
      qntelmda0 = 0.3d0
      qntelmda1 = 0.7d0
      qntvlmda0 = 0.1d0
      qntvlmda1 = 0.9d0
      plmdamap = 'QNT'
      elmdamap = 'QNT'
      vlmdamap = 'QNT'
      ostlambda = 0.5d0
c
c     each sublambda uses its own window, so the midpoint value is
c     one half for all three but the slopes differ
c
      call mapsublmda (ostlambda)
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
      tref = 0.896484375d0
      dtref = -1.7578125d0
      d2tref = -15.625d0
      call mapsublmda (ostlambda)
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
      call mapsublmda (ostlambda)
      call assert_real (plambda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window dpldlmda')
      ostlambda = 0.9d0
      call mapsublmda (ostlambda)
      call assert_real (plambda,1.0d0,1.0d-12,
     &                  'mapsublmda taper above window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper above window dpldlmda')
c
c     a QNT polarization map sets the initial and final polarization
c     flags by comparing lambda against the polarization window
c
      plmdamap = 'QNT'
      ostlambda = 0.1d0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i below window')
      call assert_logical (use_pol4f,.false.,
     &                     'mapsublmda qnt pol4f below window')
      ostlambda = 0.5d0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i mid window')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f mid window')
      ostlambda = 0.9d0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.false.,
     &                     'mapsublmda qnt pol4i above window')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f above window')
c
c     anywhere inside the window both endpoint states still
c     contribute, since plambda has not yet pinned to zero or one
c
      ostlambda = 0.65d0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i upper ramp')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f upper ramp')
      ostlambda = 0.30d0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i lower ramp')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f lower ramp')
c
c     the window bounds are inclusive, so a lambda sitting exactly
c     on either edge must keep both endpoint states
c
      ostlambda = qntplmda1
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i at upper bound')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f at upper bound')
      ostlambda = qntplmda0
      call mapsublmda (ostlambda)
      call assert_logical (use_pol4i,.true.,
     &                     'mapsublmda qnt pol4i at lower bound')
      call assert_logical (use_pol4f,.true.,
     &                     'mapsublmda qnt pol4f at lower bound')
c
c     QNT electrostatic and van der Waals maps use their own windows
c     to select the required endpoint states
c
      ostlambda = 0.1d0
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.true.,
     &                     'mapsublmda qnt ele4i below window')
      call assert_logical (use_ele4f,.false.,
     &                     'mapsublmda qnt ele4f below window')
      ostlambda = 0.5d0
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.true.,
     &                     'mapsublmda qnt ele4i mid window')
      call assert_logical (use_ele4f,.true.,
     &                     'mapsublmda qnt ele4f mid window')
      ostlambda = 0.9d0
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.false.,
     &                     'mapsublmda qnt ele4i above window')
      call assert_logical (use_ele4f,.true.,
     &                     'mapsublmda qnt ele4f above window')
      ostlambda = qntelmda0
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.true.,
     &                     'mapsublmda qnt ele4i at lower bound')
      call assert_logical (use_ele4f,.true.,
     &                     'mapsublmda qnt ele4f at lower bound')
      ostlambda = qntelmda1
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.true.,
     &                     'mapsublmda qnt ele4i at upper bound')
      call assert_logical (use_ele4f,.true.,
     &                     'mapsublmda qnt ele4f at upper bound')
      ostlambda = 0.0d0
      call mapsublmda (ostlambda)
      call assert_logical (use_vdw4i,.true.,
     &                     'mapsublmda qnt vdw4i below window')
      call assert_logical (use_vdw4f,.false.,
     &                     'mapsublmda qnt vdw4f below window')
      ostlambda = 0.5d0
      call mapsublmda (ostlambda)
      call assert_logical (use_vdw4i,.true.,
     &                     'mapsublmda qnt vdw4i mid window')
      call assert_logical (use_vdw4f,.true.,
     &                     'mapsublmda qnt vdw4f mid window')
      ostlambda = 1.0d0
      call mapsublmda (ostlambda)
      call assert_logical (use_vdw4i,.false.,
     &                     'mapsublmda qnt vdw4i above window')
      call assert_logical (use_vdw4f,.true.,
     &                     'mapsublmda qnt vdw4f above window')
      ostlambda = qntvlmda0
      call mapsublmda (ostlambda)
      call assert_logical (use_vdw4i,.true.,
     &                     'mapsublmda qnt vdw4i at lower bound')
      call assert_logical (use_vdw4f,.true.,
     &                     'mapsublmda qnt vdw4f at lower bound')
      ostlambda = qntvlmda1
      call mapsublmda (ostlambda)
      call assert_logical (use_vdw4i,.true.,
     &                     'mapsublmda qnt vdw4i at upper bound')
      call assert_logical (use_vdw4f,.true.,
     &                     'mapsublmda qnt vdw4f at upper bound')
c
c     non-QNT maps must leave all endpoint flags untouched
c
      plmdamap = 'EXP'
      elmdamap = 'EXP'
      vlmdamap = 'EXP'
      plmdaexp = 2
      elmdaexp = 2
      vlmdaexp = 2
      use_ele4i = .false.
      use_ele4f = .false.
      use_pol4i = .false.
      use_pol4f = .false.
      use_vdw4i = .false.
      use_vdw4f = .false.
      ostlambda = 0.5d0
      call mapsublmda (ostlambda)
      call assert_logical (use_ele4i,.false.,
     &                     'mapsublmda exp leaves ele4i')
      call assert_logical (use_ele4f,.false.,
     &                     'mapsublmda exp leaves ele4f')
      call assert_logical (use_pol4i,.false.,
     &                     'mapsublmda exp leaves pol4i')
      call assert_logical (use_pol4f,.false.,
     &                     'mapsublmda exp leaves pol4f')
      call assert_logical (use_vdw4i,.false.,
     &                     'mapsublmda exp leaves vdw4i')
      call assert_logical (use_vdw4f,.false.,
     &                     'mapsublmda exp leaves vdw4f')
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
      implicit none
      real*8 taper,dtaper,d2taper
c
c
c     below and at the lower bound the taper is fully on and flat
c
      call quintaper (0.05d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'quintaper below cut taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'quintaper below cut dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'quintaper below cut d2taper')
      call quintaper (0.2d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'quintaper at cut taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'quintaper at cut dtaper')
c
c     above and at the upper bound the taper is fully off and flat
c
      call quintaper (0.95d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.0d0,1.0d-12,
     &                  'quintaper above off taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'quintaper above off dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'quintaper above off d2taper')
      call quintaper (0.8d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.0d0,1.0d-12,
     &                  'quintaper at off taper')
      call assert_real (dtaper,0.0d0,1.0d-12,
     &                  'quintaper at off dtaper')
c
c     the quintic is the unique polynomial that is one at cut and
c     zero at off with vanishing first and second derivatives at
c     both ends, so it must equal the analytic smoothstep form; the
c     reference values below are that form evaluated exactly
c
      call quintaper (0.5d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.5d0,1.0d-12,
     &                  'quintaper midpoint taper')
      call assert_real (dtaper,-3.125d0,1.0d-12,
     &                  'quintaper midpoint dtaper')
      call assert_real (d2taper,0.0d0,1.0d-12,
     &                  'quintaper midpoint d2taper')
      call quintaper (0.35d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.896484375d0,1.0d-12,
     &                  'quintaper offcenter taper')
      call assert_real (dtaper,-1.7578125d0,1.0d-12,
     &                  'quintaper offcenter dtaper')
      call assert_real (d2taper,-15.625d0,1.0d-12,
     &                  'quintaper offcenter d2taper')
      call quintaper (0.70d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.035493827160493825d0,1.0d-12,
     &                  'quintaper upper half taper')
      call assert_real (dtaper,-0.96450617283950613d0,1.0d-12,
     &                  'quintaper upper half dtaper')
      call assert_real (d2taper,15.432098765432098d0,1.0d-12,
     &                  'quintaper upper half d2taper')
c
c     the same lambda gives different results for different windows,
c     since the bounds are passed in rather than read from a mode
c
      call quintaper (0.25d0,0.3d0,0.7d0,taper,dtaper,d2taper)
      call assert_real (taper,1.0d0,1.0d-12,
     &                  'quintaper narrow window below cut')
      call quintaper (0.25d0,0.1d0,0.9d0,taper,dtaper,d2taper)
      call assert_real (taper,0.95123100280761719d0,1.0d-12,
     &                  'quintaper wide window taper')
      call assert_real (dtaper,-0.87032318115234375d0,1.0d-12,
     &                  'quintaper wide window dtaper')
      call quintaper (0.25d0,0.2d0,0.8d0,taper,dtaper,d2taper)
      call assert_real (taper,0.99491222993827155d0,1.0d-12,
     &                  'quintaper middle window taper')
c
c     the reduced-coordinate form stays accurate on a window narrow
c     enough that the "switch" coefficients lose most of their digits;
c     the bounds are exact binary fractions half a step either side of
c     one half, so the reference values below are exact as well
c
      call quintaper (0.5d0,0.4990234375d0,0.5009765625d0,
     &                taper,dtaper,d2taper)
      call assert_real (taper,0.5d0,1.0d-12,
     &                  'quintaper narrow midpoint taper')
      call assert_real (dtaper,-960.0d0,1.0d-9,
     &                  'quintaper narrow midpoint dtaper')
      call assert_real (d2taper,0.0d0,1.0d-9,
     &                  'quintaper narrow midpoint d2taper')
      call quintaper (0.50048828125d0,0.4990234375d0,0.5009765625d0,
     &                taper,dtaper,d2taper)
      call assert_real (taper,0.103515625d0,1.0d-12,
     &                  'quintaper narrow offcenter taper')
      call assert_real (dtaper,-540.0d0,1.0d-9,
     &                  'quintaper narrow offcenter dtaper')
      call assert_real (d2taper,1474560.0d0,1.0d-6,
     &                  'quintaper narrow offcenter d2taper')
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
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eostmap_relstage  --  staged lambda  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eostmap_relstage" checks the staged relative free energy
c     schedule, in which the electrostatic sublambda discharges one
c     ligand over a low lambda window, stays zero while van der Waals
c     morphs between the two ligands, then recharges the other ligand
c     over a high lambda window
c
c
      subroutine test_eostmap_relstage
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 edge,eps
      real*8 dbound,window
      real*8 whi,wlo,dhi,dlo
      real*8 tref(3),dtref(3),d2tref(3)
      real*8 probe(3)
c
c
c     drive the schedule directly; only the scalar mapping is
c     exercised here, the endpoint mixing needs a real system
c
      use_relstage = .true.
      relstg1lmda0 = 0.0d0
      relstg1lmda1 = 0.3d0
      relstg2lmda0 = 0.7d0
      relstg2lmda1 = 1.0d0
      vlmdamap = 'QNT'
      qntvlmda0 = 0.3d0
      qntvlmda1 = 0.7d0
c
c     lambda of one, ligand 1 fully coupled and van der Waals with it
c
      call mapsublmda (1.0d0)
      call assert_logical (relstage.eq.'LIG1',.true.,
     &                     'maprelstage relstage at lambda one')
      call assert_real (elambda,1.0d0,1.0d-14,
     &                  'maprelstage elambda at lambda one')
      call assert_logical (relstagemix,.false.,
     &                     'maprelstage relstagemix at lambda one')
      call assert_real (vlambda,1.0d0,1.0d-14,
     &                  'maprelstage vlambda at lambda one')
      call assert_real (deldlmda,0.0d0,1.0d-14,
     &                  'maprelstage deldlmda at lambda one')
c
c     interior of the ligand 1 discharge leg, weight grows with lambda
c
      call mapsublmda (0.85d0)
      call assert_logical (relstage.eq.'LIG1',.true.,
     &                     'maprelstage relstage on ligand 1 leg')
      call assert_real (elambda,0.5d0,1.0d-12,
     &                  'maprelstage elambda on ligand 1 leg')
      call assert_logical (relstagemix,.true.,
     &                     'maprelstage relstagemix on ligand 1 leg')
      call assert_logical (deldlmda.gt.0.0d0,.true.,
     &                     'maprelstage deldlmda sign on ligand 1 leg')
      call assert_real (vlambda,1.0d0,1.0d-14,
     &                  'maprelstage vlambda on ligand 1 leg')
c
c     the van der Waals morph window, both ligands decoupled
c
      probe(1) = 0.7d0
      probe(2) = 0.5d0
      probe(3) = 0.3d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_logical (relstage.eq.'VDWM',.true.,
     &                        'maprelstage relstage in morph window')
         call assert_real (elambda,0.0d0,1.0d-14,
     &                     'maprelstage elambda in morph window')
         call assert_logical (relstagemix,.false.,
     &                        'maprelstage relstagemix in morph window')
         call assert_real (deldlmda,0.0d0,1.0d-14,
     &                     'maprelstage deldlmda in morph window')
         call assert_real (d2eldlmda2,0.0d0,1.0d-14,
     &                     'maprelstage d2eldlmda2 in morph window')
      end do
      call mapsublmda (0.5d0)
      call assert_real (vlambda,0.5d0,1.0d-12,
     &                  'maprelstage vlambda in morph window')
c
c     interior of the ligand 0 recharge leg, weight grows as lambda falls
c
      call mapsublmda (0.15d0)
      call assert_logical (relstage.eq.'LIG0',.true.,
     &                     'maprelstage relstage on ligand 0 leg')
      call assert_real (elambda,0.5d0,1.0d-12,
     &                  'maprelstage elambda on ligand 0 leg')
      call assert_logical (relstagemix,.true.,
     &                     'maprelstage relstagemix on ligand 0 leg')
      call assert_logical (deldlmda.lt.0.0d0,.true.,
     &                     'maprelstage deldlmda sign on ligand 0 leg')
      call assert_real (vlambda,0.0d0,1.0d-14,
     &                  'maprelstage vlambda on ligand 0 leg')
c
c     lambda of zero, ligand 0 fully coupled
c
      call mapsublmda (0.0d0)
      call assert_logical (relstage.eq.'LIG0',.true.,
     &                     'maprelstage relstage at lambda zero')
      call assert_real (elambda,1.0d0,1.0d-14,
     &                  'maprelstage elambda at lambda zero')
      call assert_logical (relstagemix,.false.,
     &                     'maprelstage relstagemix at lambda zero')
      call assert_real (vlambda,0.0d0,1.0d-14,
     &                  'maprelstage vlambda at lambda zero')
      call assert_real (deldlmda,0.0d0,1.0d-14,
     &                  'maprelstage deldlmda at lambda zero')
c
c     just inside a leg the weight is built by cancellation and
c     collapses onto zero, or a little past it, for about 7e-7 of
c     main lambda past the decoupled edge; that has to clamp to
c     zero and report the morph window, since a leg with no mix
c     means a weight of one to the energy routines and would switch
c     a whole ligand on inside a window lambda dynamics can visit
c
      probe(1) = 1.0d-7
      probe(2) = 5.0d-7
      do i = 1, 2
         call mapsublmda (0.7d0+probe(i))
         call assert_real (elambda,0.0d0,0.0d0,
     &                     'maprelstage collapsed weight above morph')
         call assert_logical (relstage.eq.'VDWM',.true.,
     &                        'maprelstage collapsed leg above morph')
         call assert_logical (relstagemix,.false.,
     &                        'maprelstage collapsed mix above morph')
         call mapsublmda (0.3d0-probe(i))
         call assert_real (elambda,0.0d0,0.0d0,
     &                     'maprelstage collapsed weight below morph')
         call assert_logical (relstage.eq.'VDWM',.true.,
     &                        'maprelstage collapsed leg below morph')
         call assert_logical (relstagemix,.false.,
     &                        'maprelstage collapsed mix below morph')
      end do
c
c     past the collapse the weight survives and the leg resumes
c     with a mix rather than a bare endpoint
c
      call mapsublmda (0.7d0+1.0d-5)
      call assert_logical (elambda.gt.0.0d0,.true.,
     &                     'maprelstage revived weight above morph')
      call assert_logical (relstage.eq.'VDWM',.false.,
     &                     'maprelstage revived leg above morph')
      call assert_logical (relstagemix,.true.,
     &                     'maprelstage revived mix above morph')
      call mapsublmda (0.3d0-1.0d-5)
      call assert_logical (elambda.gt.0.0d0,.true.,
     &                     'maprelstage revived weight below morph')
      call assert_logical (relstage.eq.'VDWM',.false.,
     &                     'maprelstage revived leg below morph')
      call assert_logical (relstagemix,.true.,
     &                     'maprelstage revived mix below morph')
c
c     polarization tracks the multipoles exactly across the whole
c     schedule, so the fused multipole plus polarization path stays
c     usable and the two terms never see different states
c
      probe(1) = 0.95d0
      probe(2) = 0.72d0
      probe(3) = 0.05d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_real (plambda,elambda,1.0d-15,
     &                     'maprelstage plambda tracks elambda')
         call assert_real (dpldlmda,deldlmda,1.0d-15,
     &                     'maprelstage dpldlmda tracks deldlmda')
         call assert_real (d2pldlmda2,d2eldlmda2,1.0d-15,
     &                     'maprelstage d2pldlmda2 tracks d2eldlmda2')
      end do
c
c     the weight and its derivative approach each leg boundary
c     continuously from the inside, so dU/dlambda has no step where
c     one leg hands over to the next; the weight is flat to machine
c     precision and the derivative vanishes quadratically, so it is
c     checked against the analytic bound 30*(eps/w)**2/w for a window
c     of width w, doubled to leave room for round-off
c
      eps = 1.0d-7
      window = 0.3d0
      dbound = 60.0d0 * (eps/window) * (eps/window) / window
      do i = 1, 2
         if (i .eq. 1) then
            edge = 0.7d0
         else
            edge = 0.3d0
         end if
         call mapsublmda (edge+eps)
         whi = elambda
         dhi = deldlmda
         call mapsublmda (edge-eps)
         wlo = elambda
         dlo = deldlmda
         call assert_real (whi,wlo,1.0d-14,
     &                     'maprelstage weight across leg boundary')
         call assert_logical (abs(dhi-dlo).le.dbound,.true.,
     &                        'maprelstage deriv step at leg boundary')
         call assert_logical (abs(dhi).le.dbound,.true.,
     &                        'maprelstage deriv above leg boundary')
         call assert_logical (abs(dlo).le.dbound,.true.,
     &                        'maprelstage deriv below leg boundary')
      end do
c
c     each leg is the quintic taper of its own window, so the staged
c     weight is smooth in the main lambda all the way across
c
      probe(1) = 0.75d0
      probe(2) = 0.85d0
      probe(3) = 0.95d0
      tref(1) = 0.96450617283950613d0
      tref(2) = 0.5d0
      tref(3) = 0.035493827160493825d0
      dtref(1) = -1.9290123456790123d0
      dtref(2) = -6.25d0
      dtref(3) = -1.9290123456790123d0
      d2tref(1) = -61.728395061728392d0
      d2tref(2) = 0.0d0
      d2tref(3) = 61.728395061728392d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_real (elambda,1.0d0-tref(i),1.0d-12,
     &                     'maprelstage ligand 1 leg weight')
         call assert_real (deldlmda,-dtref(i),1.0d-12,
     &                     'maprelstage ligand 1 leg deldlmda')
         call assert_real (d2eldlmda2,-d2tref(i),1.0d-12,
     &                     'maprelstage ligand 1 leg d2eldlmda2')
      end do
      probe(1) = 0.05d0
      probe(2) = 0.15d0
      probe(3) = 0.25d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_real (elambda,tref(i),1.0d-12,
     &                     'maprelstage ligand 0 leg weight')
         call assert_real (deldlmda,dtref(i),1.0d-12,
     &                     'maprelstage ligand 0 leg deldlmda')
         call assert_real (d2eldlmda2,d2tref(i),1.0d-12,
     &                     'maprelstage ligand 0 leg d2eldlmda2')
      end do
c
c     the ordinary maps must be untouched when staging is off
c
      use_relstage = .false.
      return
      end
