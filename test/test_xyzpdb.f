c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_xyzpdb  --  Tinker XYZ to PDB               ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_xyzpdb" converts scorpion.xyz to PDB format with xyzpdb and
c     checks the resulting scorpion.pdb against the reference
c
c
      subroutine test_xyzpdb
      implicit none
      integer ist
      character*240 rpath
      logical skiptest
c
c
      if (skiptest('test_xyzpdb','amber'))  return
      call pushdir ('file/scorpion')
      call execute_command_line ('rm -f scorpion.pdb* out.txt')
      call run_prog ('xyzpdb','scorpion PDB','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('scorpion','scorpion.pdb',rpath)
      call assert_files (rpath,'scorpion.pdb',1.0d-4,'test_xyzpdb')
      call execute_command_line ('rm -f scorpion.pdb* out.txt')
      call popdir
      return
      end
