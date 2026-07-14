c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_localframe2  --  local frame vdW  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_localframe2" checks the local-frame periodic vdW cases
c     exercised by the tinker-gpu localframe2.cpp tests
c
c
      subroutine test_localframe2
      implicit none
c
c
      call test_localframe2_triclinic
      call test_localframe2_monoclinic
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_localframe2_triclinic  --  triclinic  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine test_localframe2_triclinic
      implicit none
c
c
      call test_localframe2_case ('localframe2_triclinic.key',
     &   'localframe2_triclinic','localframe2.1.txt',.true.)
      call test_localframe2_case ('localframe2_triclinic.key',
     &   'localframe2_triclinic','localframe2.1.txt',.false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_localframe2_monoclinic  --  monoclinic  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_localframe2_monoclinic
      implicit none
c
c
      call test_localframe2_case ('localframe2_monoclinic.key',
     &   'localframe2_monoclinic','localframe2.2.txt',.true.)
      call test_localframe2_case ('localframe2_monoclinic.key',
     &   'localframe2_monoclinic','localframe2.2.txt',.false.)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_localframe2_case  --  run vdW  ##
c     ##                                                 ##
c     #####################################################
c
c
      subroutine test_localframe2_case (key,tname,rfile,uselist)
      use action
      use atoms
      use energi
      implicit none
      integer refcnt
      real*8 e,refeng
      logical skiptest,uselist
      character*(*) key,tname,rfile
      character*72 cname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (uselist) then
         cname = tpre//trim(tname)//'_list'
      else
         cname = tpre//trim(tname)//'_nolist'
      end if
      if (skiptest(cname,'amoeba'))  return
      call pushdir ('file/localframe')
      if (uselist) then
         call loadfix_keyadd ('localframe2',key,'neighbor-list')
      else
         call loadfix ('localframe2',key)
      end if
      call refpath ('localframe',rfile,rpath)
      call load_engcnt (rpath,'Van der Waals',refeng,refcnt)
      call analysis (e)
      call assert_real (ev,refeng,1.0d-4,trim(cname)//
     &                  ' vdw (v3)')
      call assert_int (nev,refcnt,trim(cname)//' count (v3)')
      call popdir
      call final
      return
      end
