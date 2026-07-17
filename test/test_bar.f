c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_bar  --  BAR perturbed energy file tests  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_bar" runs "bar" in mode 1 on two pairs of water trajectory
c     archives, creating a Tinker BAR file of perturbed potential
c     energies for each, and checks each generated file against its
c     stored reference with the tolerant "assert_files" comparison
c
c
      subroutine test_bar
      implicit none
      logical skiptest
c
c
      if (skiptest('test_bar','amoeba'))  return
      call test_bar_case ('water-00-0000-0950','water-00-0000-0900')
      call test_bar_case ('water-00-0950-1000','water-00-0900-1000')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_bar_case  --  one BAR file regression  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_bar_case" runs "bar 1 <A>.arc 298 <B>.arc 298 N" from the
c     fixture directory, which writes "<A>.bar", and compares it to the
c     reference of the same name
c
c
      subroutine test_bar_case (basea,baseb)
      implicit none
      real*8 eps
      integer ist
      character*(*) basea,baseb
      character*240 rpath
      character*512 args
c
c
      eps = 1.0d-3
      call pushdir ('file/bar')
      call execute_command_line
     &   ('rm -f '//trim(basea)//'.bar out.txt')
      args = '1 '//trim(basea)//'.arc 298 '//trim(baseb)//'.arc 298 N'
      call run_prog ('bar',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('bar',trim(basea)//'.bar',rpath)
      call assert_files (rpath,trim(basea)//'.bar',eps,
     &                   'test_bar '//trim(basea))
      call execute_command_line
     &   ('rm -f '//trim(basea)//'.bar out.txt')
      call popdir
      return
      end
