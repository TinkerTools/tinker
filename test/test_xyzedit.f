c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_xyzedit  --  every XYZEDIT menu option  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_xyzedit" runs each of the XYZEDIT menu options once on a
c     small AMBER99SB fixture, with answers given as arguments or on
c     standard input, and compares the edited coordinates against a
c     reference; options taking answers both ways are run both ways
c     against one reference; the template and solvent file prompts
c     must stop after five missing files, ion placement must stop
c     when too few solvent molecules can be replaced, and a simple
c     XYZ input missing its last atom line must stop with an error
c
c     the references are XYZEDIT output, each checked independently
c     when made; option 18 depends on the signs chosen by the inertia
c     tensor diagonalization, which are stable for the asymmetric
c     alanine dipeptide
c
c
      subroutine test_xyzedit
      implicit none
c
c
      call test_xyzedit_run ('test_xyzedit_m01','simple.xyz 1',
     &                       '/dev/null','simple','m01.xyz')
      call test_xyzedit_run ('test_xyzedit_m02','ala.xyz 2',
     &                       '/dev/null','ala','m02.xyz')
      call test_xyzedit_run ('test_xyzedit_m03','ala.xyz 3 100',
     &                       '/dev/null','ala','m03.xyz')
      call test_xyzedit_run ('test_xyzedit_m04','ala.xyz 4',
     &                       'm04.in','ala','m04.xyz')
      call test_xyzedit_run ('test_xyzedit_m05','ala.xyz 5',
     &                       'm05.in','ala','m05.xyz')
      call test_xyzedit_run ('test_xyzedit_m06',
     &                       'ala.xyz -k m06.key 6',
     &                       '/dev/null','ala','m06.xyz')
      call test_xyzedit_run ('test_xyzedit_m07','ala.xyz 7',
     &                       'm07.in','ala','m07.xyz')
      call test_xyzedit_run ('test_xyzedit_m07_range','ala.xyz 7',
     &                       'm07range.in','ala','m07.xyz')
      call test_xyzedit_run ('test_xyzedit_m08','ala.xyz 8',
     &                       'm08.in','ala','m08.xyz')
      call test_xyzedit_run ('test_xyzedit_m08_args',
     &                       'ala.xyz 8 341 14',
     &                       '/dev/null','ala','m08.xyz')
      call test_xyzedit_run ('test_xyzedit_m09','simple.xyz 9',
     &                       '/dev/null','simple','m09.xyz')
      call test_xyzedit_run ('test_xyzedit_m10','simple.xyz 10',
     &                       '/dev/null','simple','m10.xyz')
      call test_xyzedit_run ('test_xyzedit_m11','ala.xyz 11',
     &                       '/dev/null','ala','m11.xyz')
      call test_xyzedit_run ('test_xyzedit_m12','ala.xyz 12 tmpl.xyz',
     &                       '/dev/null','ala','m12.xyz')
      call test_xyzedit_run ('test_xyzedit_m13','ala.xyz 13',
     &                       '/dev/null','ala','m13.xyz')
      call test_xyzedit_run ('test_xyzedit_m14','ala.xyz 14',
     &                       '/dev/null','ala','m14.xyz')
      call test_xyzedit_run ('test_xyzedit_m15','ala.xyz 15 1 2 3',
     &                       '/dev/null','ala','m15.xyz')
      call test_xyzedit_run ('test_xyzedit_m16','ala.xyz 16',
     &                       '/dev/null','ala','m16.xyz')
      call test_xyzedit_run ('test_xyzedit_m17','ala.xyz 17 5',
     &                       '/dev/null','ala','m17.xyz')
      call test_xyzedit_run ('test_xyzedit_m17_range','ala.xyz 17',
     &                       'm17.in','ala','m17.xyz')
      call test_xyzedit_run ('test_xyzedit_m18','ala.xyz 18',
     &                       '/dev/null','ala','m18.xyz')
      call test_xyzedit_run ('test_xyzedit_m19','ala.xyz 19',
     &                       'm19.in','ala','m19.xyz')
      call test_xyzedit_run ('test_xyzedit_m19_args',
     &                       'ala.xyz 19 z 90',
     &                       '/dev/null','ala','m19.xyz')
      call test_xyzedit_run ('test_xyzedit_m19_x','ala.xyz 19',
     &                       'm19x.in','ala','m19x.xyz')
      call test_xyzedit_run ('test_xyzedit_m19_y','ala.xyz 19',
     &                       'm19y.in','ala','m19y.xyz')
      call test_xyzedit_run ('test_xyzedit_m20',
     &                       'ala.xyz 20 1 2 3 10 20 30',
     &                       '/dev/null','ala','m20.xyz')
      call test_xyzedit_run ('test_xyzedit_m21','stray.xyz 21',
     &                       '/dev/null','stray','m21.xyz')
      call test_xyzedit_run ('test_xyzedit_m22','wbox.xyz 22 8 8 8',
     &                       '/dev/null','wbox','m22.xyz')
      call test_xyzedit_run ('test_xyzedit_m23','wbox.xyz 23',
     &                       '/dev/null','wbox','m23.xyz')
      call test_xyzedit_run ('test_xyzedit_m24','wbox.xyz 24',
     &                       '/dev/null','wbox','m24.xyz')
      call test_xyzedit_run ('test_xyzedit_m25','ala.xyz 25 water.xyz',
     &                       '/dev/null','ala','m25.xyz')
      call test_xyzedit_run ('test_xyzedit_m26',
     &                       'water.xyz 26 20 10 10 10 N',
     &                       '/dev/null','water','m26.xyz')
      call test_xyzedit_run ('test_xyzedit_m27','ala.xyz 27 wbox.xyz',
     &                       '/dev/null','ala','m27.xyz')
      call test_xyzedit_run ('test_xyzedit_m28','wbox.xyz 28',
     &                       'm28.in','wbox','m28.xyz')
      call test_xyzedit_run ('test_xyzedit_m28_args',
     &                       'wbox.xyz 28 1 2 3',
     &                       'm28ion.in','wbox','m28.xyz')
      call test_xyzedit_fatal ('test_xyzedit_m12_missing','ala.xyz 12',
     &                         'missing.in','ala')
      call test_xyzedit_fatal ('test_xyzedit_m27_missing','ala.xyz 27',
     &                         'missing.in','ala')
      call test_xyzedit_fatal ('test_xyzedit_m28_toomany','wbox.xyz 28',
     &                         'm28many.in','wbox')
      call test_xyzedit_fatal ('test_xyzedit_truncated','trunc.xyz 1',
     &                         '/dev/null','trunc')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_xyzedit_run  --  run and compare output  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_xyzedit_run" runs XYZEDIT with the given arguments and
c     standard input in the fixture directory, then compares the new
c     coordinate file against the reference
c
c
      subroutine test_xyzedit_run (tname,args,infile,base,ref)
      implicit none
      integer ist
      logical skiptest
      character*(*) tname,args,infile,base,ref
      character*240 rpath
c
c
      if (skiptest(tname,'xyzedit'))  return
      call pushdir ('file/xyzedit')
      call execute_command_line ('rm -f '//base//'.xyz_* out.txt')
      call run_prog_input ('xyzedit',args,infile,'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call assert_int (ist,0,tname//' exit status')
      call refpath ('xyzedit',ref,rpath)
      call assert_files (rpath,base//'.xyz_2',1.0d-5,tname)
      call execute_command_line ('rm -f '//base//'.xyz_* out.txt')
      call popdir
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_xyzedit_fatal  --  stop on missing files  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_xyzedit_fatal" runs XYZEDIT with names of missing files on
c     standard input and checks that it stops through "fatal" rather
c     than prompting until the input runs out
c
c
      subroutine test_xyzedit_fatal (tname,args,infile,base)
      implicit none
      integer ist
      logical skiptest
      character*(*) tname,args,infile,base
c
c
      if (skiptest(tname,'xyzedit'))  return
      call pushdir ('file/xyzedit')
      call execute_command_line ('rm -f '//base//'.xyz_* out.txt')
      call run_prog_input ('xyzedit',args,infile,'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call execute_command_line
     &   ('grep -q "Unable to Continue" out.txt',exitstat=ist)
      call assert_int (ist,0,tname//' stops in fatal')
      call execute_command_line ('rm -f '//base//'.xyz_* out.txt')
      call popdir
      return
      end
