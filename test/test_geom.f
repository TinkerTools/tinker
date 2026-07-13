c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_geom  --  geometric restraint tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_geom" checks the geometric restraint cases exercised by
c     the tinker-gpu geom.cpp tests
c
c
      subroutine test_geom
      implicit none
c
c
      call test_geom_group
      call test_geom_distance
      call test_geom_angle
      call test_geom_torsion
      call test_geom_position
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_geom_group  --  group restraint test  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_geom_group
      implicit none
c
c
      call test_geom_case ('group','geom.1.txt','group',
     &                     'geom_group',1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_geom_distance  --  distance restraint test  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine test_geom_distance
      implicit none
c
c
      call test_geom_case ('distance','geom.2.txt','distance',
     &                     'geom_distance',1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_geom_angle  --  angle restraint test  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_geom_angle
      implicit none
c
c
      call test_geom_case ('angle','geom.3.txt','angle',
     &                     'geom_angle',1.0d-4,1.0d-4,1.0d-3)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_geom_torsion  --  torsion restraint test  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_geom_torsion
      implicit none
c
c
      call test_geom_case ('torsion','geom.4.txt','torsion',
     &                     'geom_torsion',5.0d-4,1.0d-1,1.5d-3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_geom_position  --  position restraint test  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine test_geom_position
      implicit none
c
c
      call test_geom_case ('position','geom.5.txt','position',
     &                     'geom_position',5.0d-4,1.0d-1,1.5d-3)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_geom_case  --  run one restraint case  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_geom_case" loads the shared local frame fixture with one
c     restraint key, evaluates the
c     geometric restraint term through energy, gradient and analysis
c     paths, and compares energy, count, virial and gradient to the
c     stored reference
c
c
      subroutine test_geom_case (key,reffile,kind,tname,
     &                           eps_e,eps_g,eps_v)
      use action
      use atoms
      use energi
      use restrn
      use virial
      implicit none
      integer nat,refcnt,cnt
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest
      character*(*) key,reffile,kind,tname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
c
c     load the restraint-only fixture and reference values
c
      call pushdir ('file/geom')
      call loadfix ('test_local_frame2',key//'.key')
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('geom',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      call load_engcnt (rpath,'Geometric Restraints',refeng,refcnt)
c
c     setup and level 0  --  total geometric restraint energy
c
      cnt = 0
      if (kind .eq. 'group')  cnt = ngfix
      if (kind .eq. 'distance')  cnt = ndfix
      if (kind .eq. 'angle')  cnt = nafix
      if (kind .eq. 'torsion')  cnt = ntfix
      if (kind .eq. 'position')  cnt = npfix
      call assert_int (cnt,refcnt,tpre//tname//' count')
      e = energy ()
      call assert_real (esum,ref_e,eps_e,tpre//tname//' energy (v0)')
      call assert_real (eg,refeng,eps_e,tpre//tname//' geom (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,tpre//tname//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,tpre//tname//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,tpre//tname//' virial (v1)')
c
c     level 3  --  energy and per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,
     &                  tpre//tname//' analysis (v3)')
      call check_engcnt (rpath,'Geometric Restraints',eg,neg,eps_e,
     &                   tpre//tname//' geom (v3)')
c
c     clean up
c
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
