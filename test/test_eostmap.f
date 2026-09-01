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
      call test_eostmap_refresh
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
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eostmap_refresh  --  map owner routing  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eostmap_refresh" checks that refreshsublmda maps the one
c     main lambda whichever sampling method is active, and leaves
c     explicit sublambdas alone when no main lambda is present
c
c
      subroutine test_eostmap_refresh
      use dlmda
      use mutant
      implicit none
c
c
c     use identity maps so a sublambda tracks the main lambda exactly
c
      use_dlmda = .true.
      use_elmdamap = .true.
      use_plmdamap = .true.
      use_vlmdamap = .true.
      use_relstage = .false.
      elmdamap = 'EXP'
      plmdamap = 'EXP'
      vlmdamap = 'EXP'
      elmdaexp = 1
      plmdaexp = 1
      vlmdaexp = 1
c
c     a static main lambda drives the map with no sampling method
c
      use_ost = .false.
      use_meta = .false.
      use_ti = .false.
      use_mainlmda = .true.
      lambda = 0.25d0
      call refreshsublmda
      call assert_real (elambda,0.25d0,0.0d0,
     &                  'refreshsublmda static lambda')
c
c     thermodynamic integration reads the same main lambda
c
      use_ti = .true.
      lambda = 0.50d0
      call refreshsublmda
      call assert_real (elambda,0.50d0,0.0d0,
     &                  'refreshsublmda TI lambda')
c
c     so do OST and metadynamics
c
      use_ti = .false.
      use_ost = .true.
      lambda = 0.75d0
      call refreshsublmda
      call assert_real (elambda,0.75d0,0.0d0,
     &                  'refreshsublmda OST lambda')
      use_ost = .false.
      use_meta = .true.
      lambda = 0.40d0
      call refreshsublmda
      call assert_real (elambda,0.40d0,0.0d0,
     &                  'refreshsublmda metadynamics lambda')
c
c     explicit sublambdas remain unchanged without a main-lambda owner
c
      use_meta = .false.
      use_mainlmda = .false.
      elambda = 0.125d0
      call refreshsublmda
      call assert_real (elambda,0.125d0,0.0d0,
     &                  'refreshsublmda no owner')
      use_dlmda = .false.
      use_elmdamap = .false.
      use_plmdamap = .false.
      use_vlmdamap = .false.
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
c     inverse-power, asymmetric power and quintic lambda maps
c
c
      subroutine test_eostmap_mapsub
      use bound
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 tref,dtref,d2tref
      real*8 weight1,dweight1,d2weight1
      logical need0,need1
c
c
c     the dual topology exponents weight the endpoint states, and a real
c     calculation takes them from "mutate_dlmda" while parsing keywords;
c     these cases drive the lambda maps on their own, so they carry the
c     exponents themselves as they do every other map parameter
c
      emdtexp = 1
      epdtexp = 1
      evdtexp = 1
c
c     test exponential sublambda maps and chain rule derivatives
c
      call resetost (5,5,1)
      lambda = 0.25d0
      use_plmdamap = .true.
      use_elmdamap = .true.
      use_vlmdamap = .true.
      plmdamap = 'EXP'
      elmdamap = 'EXP'
      vlmdamap = 'EXP'
      plmdaexp = 2
      elmdaexp = 3
      vlmdaexp = 4
      call mapsublmda (lambda)
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
      lambda = 0.25d0
      use_plmdamap = .true.
      use_elmdamap = .true.
      use_vlmdamap = .true.
      plmdamap = 'INV'
      elmdamap = 'INV'
      vlmdamap = 'INV'
      plmdainvn = 2
      elmdainvn = 3
      vlmdainvn = 4
      plmdainveps = 0.01d0
      elmdainveps = 0.02d0
      vlmdainveps = 0.03d0
      call mapsublmda (lambda)
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
c     test asymmetric power sublambda maps and derivatives, the three
c     terms taking different powers and slope ratios so that a crossed
c     argument would land the wrong shape on the wrong sublambda
c
      call resetost (5,5,1)
      lambda = 0.25d0
      use_plmdamap = .true.
      use_elmdamap = .true.
      use_vlmdamap = .true.
      plmdamap = 'APM'
      elmdamap = 'APM'
      vlmdamap = 'APM'
      plmdaapmn = 2
      elmdaapmn = 3
      vlmdaapmn = 4
      plmdaapmrho = 1.5d0
      elmdaapmrho = 2.0d0
      vlmdaapmrho = 2.5d0
      call mapsublmda (lambda)
      call assert_real (plambda,0.3203125d0,1.0d-12,
     &                  'mapsublmda asympower plambda')
      call assert_real (dpldlmda,1.09375d0,1.0d-12,
     &                  'mapsublmda asympower dpldlmda')
      call assert_real (d2pldlmda2,-1.25d0,1.0d-12,
     &                  'mapsublmda asympower d2pldlmda2')
      call assert_real (elambda,0.3818359375d0,1.0d-12,
     &                  'mapsublmda asympower elambda')
      call assert_real (deldlmda,1.140625d0,1.0d-12,
     &                  'mapsublmda asympower deldlmda')
      call assert_real (d2eldlmda2,-2.4375d0,1.0d-12,
     &                  'mapsublmda asympower d2eldlmda2')
      call assert_real (vlambda,0.43017578125d0,1.0d-12,
     &                  'mapsublmda asympower vlambda')
      call assert_real (dvldlmda,1.134765625d0,1.0d-12,
     &                  'mapsublmda asympower dvldlmda')
      call assert_real (d2vldlmda2,-3.34375d0,1.0d-12,
     &                  'mapsublmda asympower d2vldlmda2')
c
c     any map other than EXP, INV or APM falls back to the quintic
c     taper, where the sublambda is the complement of the taper
c
      use_bounds = .false.
      call resetost (5,5,1)
      qntplmda0 = 0.2d0
      qntplmda1 = 0.8d0
      qntelmda0 = 0.3d0
      qntelmda1 = 0.7d0
      qntvlmda0 = 0.1d0
      qntvlmda1 = 0.9d0
      use_plmdamap = .true.
      use_elmdamap = .true.
      use_vlmdamap = .true.
      plmdamap = 'QNT'
      elmdamap = 'QNT'
      vlmdamap = 'QNT'
      lambda = 0.5d0
c
c     each sublambda uses its own window, so the midpoint value is
c     one half for all three but the slopes differ
c
      call mapsublmda (lambda)
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
      lambda = 0.35d0
      tref = 0.896484375d0
      dtref = -1.7578125d0
      d2tref = -15.625d0
      call mapsublmda (lambda)
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
      lambda = 0.1d0
      call mapsublmda (lambda)
      call assert_real (plambda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper below window dpldlmda')
      lambda = 0.9d0
      call mapsublmda (lambda)
      call assert_real (plambda,1.0d0,1.0d-12,
     &                  'mapsublmda taper above window plambda')
      call assert_real (dpldlmda,0.0d0,1.0d-12,
     &                  'mapsublmda taper above window dpldlmda')
c
c     a QNT polarization map sets the initial and final polarization
c     flags by comparing lambda against the polarization window
c
      plmdamap = 'QNT'
      lambda = 0.1d0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt pol need0 below window')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt pol need1 below window')
      lambda = 0.5d0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt pol need0 mid window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt pol need1 mid window')
      lambda = 0.9d0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt pol need0 above window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt pol need1 above window')
c
c     anywhere inside the window both endpoint states still
c     contribute, since plambda has not yet pinned to zero or one
c
      lambda = 0.65d0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt pol need0 upper ramp')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt pol need1 upper ramp')
      lambda = 0.30d0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt pol need0 lower ramp')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt pol need1 lower ramp')
c
c     the window bounds are inclusive, so a lambda sitting exactly
c     on either edge must keep both endpoint states
c
      lambda = qntplmda1
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt pol need0 at upper bound')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt pol need1 at upper bound')
      lambda = qntplmda0
      call mapsublmda (lambda)
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt pol need0 at lower bound')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt pol need1 at lower bound')
c
c     QNT electrostatic and van der Waals maps use their own windows
c     to select the required endpoint states
c
      lambda = 0.1d0
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt ele need0 below window')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt ele need1 below window')
      lambda = 0.5d0
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt ele need0 mid window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt ele need1 mid window')
      lambda = 0.9d0
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt ele need0 above window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt ele need1 above window')
      lambda = qntelmda0
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt ele need0 at lower bound')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt ele need1 at lower bound')
      lambda = qntelmda1
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt ele need0 at upper bound')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt ele need1 at upper bound')
      lambda = 0.0d0
      call mapsublmda (lambda)
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt vdw need0 below window')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt vdw need1 below window')
      lambda = 0.5d0
      call mapsublmda (lambda)
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt vdw need0 mid window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt vdw need1 mid window')
      lambda = 1.0d0
      call mapsublmda (lambda)
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt vdw need0 above window')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt vdw need1 above window')
      lambda = qntvlmda0
      call mapsublmda (lambda)
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda qnt vdw need0 at lower bound')
      call assert_logical (need1,.false.,
     &                     'mapsublmda qnt vdw need1 at lower bound')
      lambda = qntvlmda1
      call mapsublmda (lambda)
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'mapsublmda qnt vdw need0 at upper bound')
      call assert_logical (need1,.true.,
     &                     'mapsublmda qnt vdw need1 at upper bound')
c
c     a non-QNT map has no window to fall flat outside of, so away from
c     the ends of the main lambda both endpoint states stay live
c
      plmdamap = 'EXP'
      elmdamap = 'EXP'
      vlmdamap = 'EXP'
      plmdaexp = 2
      elmdaexp = 2
      vlmdaexp = 2
      lambda = 0.5d0
      call mapsublmda (lambda)
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda exp ele need0')
      call assert_logical (need1,.true.,
     &                     'mapsublmda exp ele need1')
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dpldlmda,d2pldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda exp pol need0')
      call assert_logical (need1,.true.,
     &                     'mapsublmda exp pol need1')
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'mapsublmda exp vdw need0')
      call assert_logical (need1,.true.,
     &                     'mapsublmda exp vdw need1')
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
c
c     the asymmetric power map pins both endpoint values, and carries
c     the slope ratio at the decoupled end against a unit slope at the
c     coupled end, which is what the normalization exists to enforce
c
      call sublmdaapm (0.0d0,12,4.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.0d0,1.0d-12,
     &                  'sublmdaapm x=0 lmda')
      call assert_real (dlmda,4.0d0,1.0d-12,
     &                  'sublmdaapm x=0 dlmda')
      call assert_real (d2lmda,-39.272727272727273d0,1.0d-12,
     &                  'sublmdaapm x=0 d2lmda')
      call sublmdaapm (1.0d0,12,4.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,1.0d0,1.0d-12,
     &                  'sublmdaapm x=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaapm x=1 dlmda')
      call assert_real (d2lmda,3.272727272727273d0,1.0d-12,
     &                  'sublmdaapm x=1 d2lmda')
c
c     the interior map runs ahead of the linear schedule, its slope
c     having already fallen below one by the midpoint
c
      call sublmdaapm (0.5d0,12,4.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.6153564453125d0,1.0d-12,
     &                  'sublmdaapm interior lmda')
      call assert_real (dlmda,0.728138316761364d0,1.0d-12,
     &                  'sublmdaapm interior dlmda')
      call assert_real (d2lmda,-0.017578125d0,1.0d-12,
     &                  'sublmdaapm interior d2lmda')
c
c     a power below two or a slope ratio of one or less leaves no
c     asymmetry to impose, so the map degenerates to the identity
c
      call sublmdaapm (0.4d0,1,4.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.4d0,1.0d-12,
     &                  'sublmdaapm n=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaapm n=1 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaapm n=1 d2lmda')
      call sublmdaapm (0.4d0,12,1.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.4d0,1.0d-12,
     &                  'sublmdaapm rho=1 lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaapm rho=1 dlmda')
      call assert_real (d2lmda,0.0d0,1.0d-12,
     &                  'sublmdaapm rho=1 d2lmda')
c
c     a slope ratio at or above the power sits on the pole of the
c     normalization, and falls back to the identity rather than it
c
      call sublmdaapm (0.4d0,4,4.0d0,lmda,dlmda,d2lmda)
      call assert_real (lmda,0.4d0,1.0d-12,
     &                  'sublmdaapm rho=n lmda')
      call assert_real (dlmda,1.0d0,1.0d-12,
     &                  'sublmdaapm rho=n dlmda')
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
      if (allocated(dfsumdl))  deallocate (dfsumdl)
      allocate (dfpdl(3,n))
      allocate (dfmdl(3,n))
      allocate (dfvdl(3,n))
      allocate (dfsumdl(3,n))
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
      deallocate (dfsumdl)
      return
      end
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eostmap_relstage  --  staged lambda  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eostmap_relstage" checks the staged relative free energy
c     schedule, in which each run walks one declared leg, the LIG2 leg
c     discharging ligand 2 over its electrostatic window, the VDWM leg
c     morphing van der Waals between the two ligands while both stay
c     electrostatically decoupled, and the LIG1 leg charging ligand 1
c
c
      subroutine test_eostmap_relstage
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 tref(3),dtref(3),d2tref(3)
      real*8 probe(3)
      real*8 wlo,dlo,vlo,dvlo
      real*8 weight1,dweight1,d2weight1
      logical need0,need1
c
c
c     drive the schedule directly; only the scalar mapping is
c     exercised here, the endpoint mixing needs a real system
c
      use_relstage = .true.
      elmdamap = 'QNT'
      vlmdamap = 'QNT'
      qntvlmda0 = 0.3d0
      qntvlmda1 = 0.7d0
c
c     a staged leg spends its exponents on the taper, so the endpoint
c     weight is linear in each sublambda
c
      emdtexp = 1
      epdtexp = 1
      evdtexp = 1
c
c     the ligand 1 leg charges ligand 1 over the upper window against
c     the decoupled reference, with van der Waals already morphed on
c
      relstage = 'LIG1'
      qntelmda0 = 0.7d0
      qntelmda1 = 1.0d0
c
c     lambda of one, ligand 1 fully coupled and van der Waals with it,
c     so only the coupled endpoint has to be built
c
      call mapsublmda (1.0d0)
      call assert_logical (erelst0.eq.relnone,.true.,
     &                     'maprelstage lig1 lower state at lambda one')
      call assert_logical (erelst1.eq.rellig1,.true.,
     &                     'maprelstage lig1 upper state at lambda one')
      call assert_real (elambda,1.0d0,1.0d-14,
     &                  'maprelstage lig1 elambda at lambda one')
      call assert_real (vlambda,1.0d0,1.0d-14,
     &                  'maprelstage lig1 vlambda at lambda one')
      call assert_real (deldlmda,0.0d0,1.0d-14,
     &                  'maprelstage lig1 deldlmda at lambda one')
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'maprelstage lig1 need0 at lambda one')
      call assert_logical (need1,.true.,
     &                     'maprelstage lig1 need1 at lambda one')
c
c     interior of the ligand 1 leg, the weight growing with lambda and
c     both endpoint states live
c
      call mapsublmda (0.85d0)
      call assert_real (elambda,0.5d0,1.0d-12,
     &                  'maprelstage lig1 elambda on leg')
      call assert_logical (deldlmda.gt.0.0d0,.true.,
     &                     'maprelstage lig1 deldlmda sign on leg')
      call assert_real (vlambda,1.0d0,1.0d-14,
     &                  'maprelstage lig1 vlambda on leg')
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need0,.true.,
     &                     'maprelstage lig1 need0 on leg')
      call assert_logical (need1,.true.,
     &                     'maprelstage lig1 need1 on leg')
c
c     the van der Waals endpoints are the two coupled states on every
c     leg, pinned onto one of them by a flat weight here
c
      call assert_logical (vrelst0.eq.rellig2,.true.,
     &                     'maprelstage lig1 vdw lower state')
      call assert_logical (vrelst1.eq.rellig1,.true.,
     &                     'maprelstage lig1 vdw upper state')
      call relpowerwt (vlambda,evdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              dvldlmda,d2vldlmda2,need0,need1)
      call assert_logical (need0,.false.,
     &                     'maprelstage lig1 vdw need0 on leg')
c
c     the flat end of the ligand 1 leg, where the weight has fallen to
c     zero and only the decoupled reference has to be built
c
      call mapsublmda (0.7d0)
      call assert_real (elambda,0.0d0,1.0d-14,
     &                  'maprelstage lig1 elambda at leg end')
      call assert_real (deldlmda,0.0d0,1.0d-14,
     &                  'maprelstage lig1 deldlmda at leg end')
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need1,.false.,
     &                     'maprelstage lig1 need1 at leg end')
c
c     the ligand 2 leg discharges ligand 2 over the lower window, so
c     its weight is the taper itself and falls as lambda rises
c
      relstage = 'LIG2'
      qntelmda0 = 0.0d0
      qntelmda1 = 0.3d0
      call mapsublmda (0.0d0)
      call assert_logical (erelst1.eq.rellig2,.true.,
     &                     'maprelstage lig2 upper state at zero')
      call assert_real (elambda,1.0d0,1.0d-14,
     &                  'maprelstage lig2 elambda at lambda zero')
      call assert_real (vlambda,0.0d0,1.0d-14,
     &                  'maprelstage lig2 vlambda at lambda zero')
      call assert_real (deldlmda,0.0d0,1.0d-14,
     &                  'maprelstage lig2 deldlmda at lambda zero')
      call mapsublmda (0.15d0)
      call assert_real (elambda,0.5d0,1.0d-12,
     &                  'maprelstage lig2 elambda on leg')
      call assert_logical (deldlmda.lt.0.0d0,.true.,
     &                     'maprelstage lig2 deldlmda sign on leg')
      call assert_real (vlambda,0.0d0,1.0d-14,
     &                  'maprelstage lig2 vlambda on leg')
c
c     the middle leg holds both ligands decoupled, so electrostatics
c     sit at the reference state and van der Waals morphs across
c
      relstage = 'VDWM'
      probe(1) = 0.3d0
      probe(2) = 0.5d0
      probe(3) = 0.7d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_logical (erelst0.eq.relnone,.true.,
     &                        'maprelstage vdwm lower state')
         call assert_logical (erelst1.eq.relnone,.true.,
     &                        'maprelstage vdwm upper state')
         call assert_real (elambda,0.0d0,1.0d-14,
     &                     'maprelstage vdwm elambda')
         call assert_real (deldlmda,0.0d0,1.0d-14,
     &                     'maprelstage vdwm deldlmda')
         call assert_real (d2eldlmda2,0.0d0,1.0d-14,
     &                     'maprelstage vdwm d2eldlmda2')
         call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
         call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
         call assert_logical (need1,.false.,
     &                       'maprelstage vdwm need1')
      end do
      call mapsublmda (0.5d0)
      call assert_real (vlambda,0.5d0,1.0d-12,
     &                  'maprelstage vdwm vlambda in morph window')
c
c     each leg hands its end state to the next, so the three legs of a
c     transformation compose; the shared boundary must map to the same
c     sublambdas and the same flat derivatives from either side, which
c     is what lets the legs be run as separate calculations
c
      relstage = 'LIG2'
      qntelmda0 = 0.0d0
      qntelmda1 = 0.3d0
      call mapsublmda (0.3d0)
      wlo = elambda
      dlo = deldlmda
      vlo = vlambda
      dvlo = dvldlmda
      relstage = 'VDWM'
      call mapsublmda (0.3d0)
      call assert_real (elambda,wlo,0.0d0,
     &                  'maprelstage lig2 to vdwm elambda')
      call assert_real (deldlmda,dlo,0.0d0,
     &                  'maprelstage lig2 to vdwm deldlmda')
      call assert_real (vlambda,vlo,0.0d0,
     &                  'maprelstage lig2 to vdwm vlambda')
      call assert_real (dvldlmda,dvlo,0.0d0,
     &                  'maprelstage lig2 to vdwm dvldlmda')
      call mapsublmda (0.7d0)
      wlo = elambda
      dlo = deldlmda
      vlo = vlambda
      dvlo = dvldlmda
      relstage = 'LIG1'
      qntelmda0 = 0.7d0
      qntelmda1 = 1.0d0
      call mapsublmda (0.7d0)
      call assert_real (elambda,wlo,0.0d0,
     &                  'maprelstage vdwm to lig1 elambda')
      call assert_real (deldlmda,dlo,0.0d0,
     &                  'maprelstage vdwm to lig1 deldlmda')
      call assert_real (vlambda,vlo,0.0d0,
     &                  'maprelstage vdwm to lig1 vlambda')
      call assert_real (dvldlmda,dvlo,0.0d0,
     &                  'maprelstage vdwm to lig1 dvldlmda')
c
c     just inside a leg the weight is built by cancellation and
c     collapses onto zero for about 7e-7 of main lambda past the
c     decoupled edge; that has to clamp to zero so the energy is the
c     bare reference there, since a weight of one would switch a whole
c     ligand on inside a window lambda dynamics can visit
c
c     the taper slope survives the collapse, though, so the coupled
c     endpoint stays live and dU/dlambda keeps sampling it; the weight
c     alone going flat is not enough to drop an endpoint
c
      probe(1) = 1.0d-7
      probe(2) = 5.0d-7
      do i = 1, 2
         call mapsublmda (0.7d0+probe(i))
         call assert_real (elambda,0.0d0,0.0d0,
     &                     'maprelstage lig1 collapsed weight')
         call assert_logical (deldlmda.gt.0.0d0,.true.,
     &                        'maprelstage lig1 collapsed slope')
         call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
         call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
         call assert_logical (need1,.true.,
     &                        'maprelstage lig1 collapsed need1')
      end do
c
c     past the collapse the weight survives and the leg resumes with a
c     mix rather than a bare endpoint
c
      call mapsublmda (0.7d0+1.0d-5)
      call assert_logical (elambda.gt.0.0d0,.true.,
     &                     'maprelstage lig1 revived weight')
      call relpowerwt (elambda,emdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &              deldlmda,d2eldlmda2,need0,need1)
      call assert_logical (need1,.true.,
     &                     'maprelstage lig1 revived need1')
c
c     polarization tracks the multipoles exactly across the whole
c     schedule, so the fused multipole plus polarization path stays
c     usable and the two terms never see different states
c
      probe(1) = 0.95d0
      probe(2) = 0.72d0
      probe(3) = 0.85d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_real (plambda,elambda,1.0d-15,
     &                     'maprelstage plambda tracks elambda')
         call assert_real (dpldlmda,deldlmda,1.0d-15,
     &                     'maprelstage dpldlmda tracks deldlmda')
         call assert_real (d2pldlmda2,d2eldlmda2,1.0d-15,
     &                     'maprelstage d2pldlmda2 tracks d2eldlmda2')
         call assert_logical (prelst0.eq.erelst0,.true.,
     &                        'maprelstage pol lower state tracks ele')
         call assert_logical (prelst1.eq.erelst1,.true.,
     &                        'maprelstage pol upper state tracks ele')
      end do
c
c     the staged weight of each leg is the quintic taper of its own
c     window, rising with lambda on the ligand 1 leg
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
     &                     'maprelstage lig1 leg weight')
         call assert_real (deldlmda,-dtref(i),1.0d-12,
     &                     'maprelstage lig1 leg deldlmda')
         call assert_real (d2eldlmda2,-d2tref(i),1.0d-12,
     &                     'maprelstage lig1 leg d2eldlmda2')
      end do
c
c     and falling with lambda on the ligand 2 leg, the same taper of a
c     window of the same width read the other way
c
      relstage = 'LIG2'
      qntelmda0 = 0.0d0
      qntelmda1 = 0.3d0
      probe(1) = 0.05d0
      probe(2) = 0.15d0
      probe(3) = 0.25d0
      do i = 1, 3
         call mapsublmda (probe(i))
         call assert_real (elambda,tref(i),1.0d-12,
     &                     'maprelstage lig2 leg weight')
         call assert_real (deldlmda,dtref(i),1.0d-12,
     &                     'maprelstage lig2 leg deldlmda')
         call assert_real (d2eldlmda2,d2tref(i),1.0d-12,
     &                     'maprelstage lig2 leg d2eldlmda2')
      end do
c
c     the ordinary maps must be untouched when staging is off
c
      use_relstage = .false.
      relstage = 'VDWM'
      return
      end
