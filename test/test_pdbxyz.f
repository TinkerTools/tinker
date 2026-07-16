c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_pdbxyz  --  PDB/CIF to Tinker XYZ           ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_pdbxyz" converts scorpion.cif to a Tinker XYZ file and a
c     sequence file with pdbxyz and checks both against the reference;
c     because the fixture already holds scorpion.xyz and scorpion.seq
c     (inputs for the xyzpdb test), pdbxyz writes versioned "_2" copies
c
c
      subroutine test_pdbxyz
      implicit none
      real*8 eps
      integer ist
      character*240 rpath
      logical skiptest
c
c
      eps = 1.0d-4
      if (skiptest('test_pdbxyz','amber'))  return
      call pushdir ('file/scorpion')
      call execute_command_line
     &   ('rm -f scorpion.xyz_* scorpion.seq_* out.txt')
      call run_prog ('pdbxyz','scorpion.cif 0','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('scorpion','scorpion.seq',rpath)
      call assert_files (rpath,'scorpion.seq_2',eps,'test_pdbxyz seq')
      call refpath ('scorpion','scorpion.xyz',rpath)
      call assert_files (rpath,'scorpion.xyz_2',eps,'test_pdbxyz xyz')
      call execute_command_line
     &   ('rm -f scorpion.xyz_* scorpion.seq_* out.txt')
      call popdir
      return
      end
