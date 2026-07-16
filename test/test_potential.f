c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_potential  --  potential program tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     this file collects Tinker system tests that run the
c     "test_potential" program as a subprocess and compare its output
c     against a stored reference with a floating-point tolerance
c
c
      subroutine test_potential
      implicit none
c
c
      call test_potential_tip3p
      call test_potential_amoeba
      call test_potential_hippo
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_potential_tip3p  --  TIP3P potential fit  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_potential_tip3p" reproduces the TIP3P run of the manual
c     water potential-fitting example and checks the full program output
c     against the reference, allowing floats to differ within tolerance
c
c
      subroutine test_potential_tip3p
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='potential_tip3p')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'tip3p'))  return
      call pushdir ('file/potential')
      call run_prog ('potential',
     &               '-k water.tip3p 6 water water n 0.001',
     &               'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('potential',tname//'.txt',rpath)
      call assert_files (rpath,'out.txt',1.0d-3,tpre//tname)
      call execute_command_line ('rm -f water.key* out.txt')
      call popdir
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_potential_amoeba  --  AMOEBA potential fit  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_potential_amoeba" reproduces the AMOEBA run of the manual
c     water potential-fitting example and checks the full program output
c     against the reference, allowing floats to differ within tolerance
c
c
      subroutine test_potential_amoeba
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='potential_amoeba')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'amoeba'))  return
      call pushdir ('file/potential')
      call run_prog ('potential',
     &               '-k water.amoeba 6 water water n 0.001',
     &               'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('potential',tname//'.txt',rpath)
      call assert_files (rpath,'out.txt',1.0d-3,tpre//tname)
      call execute_command_line ('rm -f water.key* out.txt')
      call popdir
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_potential_hippo  --  HIPPO potential fit  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_potential_hippo" reproduces the HIPPO run of the manual
c     water potential-fitting example and checks the full program output
c     against the reference, allowing floats to differ within tolerance
c
c
      subroutine test_potential_hippo
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      character*(*) tpre
      parameter (tname='potential_hippo')
      parameter (tpre='test_')
c
c
      if (skiptest(tpre//tname,'hippo'))  return
      call pushdir ('file/potential')
      call run_prog ('potential',
     &               '-k water.hippo 6 water water n 0.001',
     &               'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('potential',tname//'.txt',rpath)
      call assert_files (rpath,'out.txt',1.0d-3,tpre//tname)
      call execute_command_line ('rm -f water.key* out.txt')
      call popdir
      return
      end
