c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_kinetic  --  kinetic energy test  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_kinetic" checks the arbox kinetic energy case exercised
c     by the tinker-gpu kinetic.cpp test
c
c
      subroutine test_kinetic
      use atoms
      use mdstuf
      use moldyn
      implicit none
      real*8 eksum,temp
      real*8 ekin(3,3)
      logical skiptest
      character*(*) tname
      parameter (tname='test_kinetic_arbox')
c
c
      if (skiptest(tname,'amoeba'))  return
      call pushdir ('file/kinetic')
      call loadfix ('arbox','arbox.key')
      if (allocated(v))  deallocate (v)
      allocate (v(3,n))
      call test_kinetic_readvel ('arbox.dyn_2')
      nfree = 3 * n
      integrate = 'VERLET'
      call kinetic (eksum,ekin,temp)
      call assert_real (eksum,100446.40376d0,1.0d-4,
     &                  tname//' eksum')
      call assert_real (temp,156008.001336d0,1.0d-4,
     &                  tname//' temp')
      call popdir
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_kinetic_readvel  --  read dyn file  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_kinetic_readvel (fname)
      use atoms
      use moldyn
      implicit none
      integer i,idyn,freeunit
      character*(*) fname
      character*240 record
      logical found
c
c
      found = .false.
      idyn = freeunit ()
      open (unit=idyn,file=fname,status='old')
      do while (.not. found)
         read (idyn,10,end=20)  record
   10    format (a240)
         if (index(record,'Current Atomic Velocities') .gt. 0)
     &      found = .true.
      end do
      do i = 1, n
         read (idyn,*)  v(1,i),v(2,i),v(3,i)
      end do
      close (unit=idyn)
      return
   20 continue
      close (unit=idyn)
      call assert_int (0,1,'test_kinetic_arbox velocities found')
      return
      end
