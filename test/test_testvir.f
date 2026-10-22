c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_testvir  --  testvir program tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     this file collects Tinker system tests that run the "testvir"
c     program as a subprocess and check that its analytical virial
c     agrees with the numerical virial from the lattice derivatives
c
c     the systems carry a net charge, so the Ewald uniform background
c     correction contributes to each diagonal element of the virial;
c     the lambda cases decouple the chloride of a +2 environment on
c     a quintic map, so the background term goes through "empole4"
c     with a lambda scaled net charge
c
c     cases 08 and 09 scale the charge flux of an AMOEBA+ water
c     hydrogen whose flux bond reaches into the environment, so the
c     net charge changes with geometry and the background term adds
c     a charge flux virial, through "empole4" with lambda derivatives
c     and through "empole1" without them
c
c     case 10 annihilates a neutral ligand with Ewald dispersion,
c     whose reciprocal space virial combines two rescaled sums
c
c     case 11 scales HIPPO water charge transfer by a squared dual
c     topology weight through "echgtrn1", in a box small enough for
c     the cutoff to reach replicate cells
c
c
      subroutine test_testvir
      implicit none
c
c
      call test_testvir_case ('01_charge_ewald','ions')
      call test_testvir_case ('02_charge_ewald_lights','ions')
      call test_testvir_case ('03_charge_ewald_nlist','ions')
      call test_testvir_case ('04_mpole_ewald','ionwat')
      call test_testvir_case ('05_mpole_ewald_nlist','ionwat')
      call test_testvir_case ('06_mpole_lambda_l05','ionwat')
      call test_testvir_case ('07_mpole_lambda_nlist_l05','ionwat')
      call test_testvir_case ('08_mpole_cflux_lambda_l04','apgas')
      call test_testvir_case ('09_mpole_cflux_l04','apgas')
      call test_testvir_case ('10_disp_dewald_ann_l05',
     &                        '../ephippo/dmso')
      call test_testvir_case ('11_ct_lambda_l04','../testlmda/hipwat')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_testvir_case  --  one testvir case  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_testvir_case" runs "testvir" on the given coordinates with
c     the key file that shares the given base name, and checks that the
c     analytical and numerical virial tensors in its output agree; the
c     program prints three decimals, and the tolerance also absorbs the
c     finite difference noise in the numerical virial
c
c
      subroutine test_testvir_case (base,xyz)
      use assert
      implicit none
      integer ist
      real*8 vira(3,3),virn(3,3)
      logical skiptest,founda,foundn
      character*(*) base,xyz
      character*512 args
c
c
      if (skiptest('test_testvir_'//base,'testvir'))  return
      call pushdir ('file/testvir')
      call execute_command_line ('rm -f out.txt')
      args = '-k '//base//'.key '//xyz
      call run_prog ('testvir',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call read_tensor ('out.txt','Analytical Virial Tensor',
     &                     vira,founda)
      call read_tensor ('out.txt','Numerical Virial Tensor',
     &                     virn,foundn)
      if (.not.(founda .and. foundn)) then
         nfail = nfail + 1
         write (*,10)  base
   10    format (1x,'FAIL  test_testvir ',a,' (virial not in output)')
         call popdir
         return
      end if
      call assert_grad (vira,virn,3,1.0d-2,'test_testvir '//base)
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine read_tensor  --  read labeled 3x3 tensor  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "read_tensor" reads a 3x3 tensor written as one labeled row and
c     two continuation rows of three components each, as printed by the
c     "testvir" program, and reports whether the label was found
c
c
      subroutine read_tensor (fname,label,t,found)
      implicit none
      integer i,ios,next
      real*8 t(3,3)
      logical found
      character*(*) fname,label
      character*256 line
c
c
      found = .false.
      open (unit=15,file=fname,status='old',iostat=ios)
      if (ios .ne. 0)  return
   20 continue
      read (15,'(a)',iostat=ios) line
      if (ios .ne. 0)  go to 30
      if (index(line,label) .eq. 0)  go to 20
      next = index(line,':') + 1
      do i = 1, 3
         if (i .gt. 1) then
            read (15,'(a)',iostat=ios) line
            if (ios .ne. 0)  go to 30
            next = 1
         end if
         call getfloat (line,t(1,i),next)
         call getfloat (line,t(2,i),next)
         call getfloat (line,t(3,i),next)
      end do
      found = .true.
   30 continue
      close (unit=15)
      return
      end
