c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar4  --  polarization energy & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar4" calculates the induced dipole polarization energy
c     and first derivatives with respect to Cartesian coordinates
c
c
      subroutine epolar4
      use dlmda
      use iounit
      use limits
      use mplpot
      use mutant
      use polpot
      use virial
      implicit none
      integer i,j
c
c
c     check for use of TCG polarization with charge penetration
c
      if (poltyp.eq.'TCG' .and. use_chgpen) then
         write (iout,10)
   10    format (/,' EPOLAR4  --  TCG Polarization not Available',
     &              ' with Charge Penetration')
         call fatal
      end if
c
c     compute polarization interactions
c
      call epolar4f
c
c     modify the gradient and virial for exchange polarization
c
      if (use_expol) then
         call dexpol
      end if
c
c     add the polarization virial to main virial
c
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = vir(j,i) + epvir(j,i)
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar4f  --  dual topology lambda derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar4f" calculates the polarization energy, derivatives
c     with respect to Cartesian coordinates, and lambda derivatives
c     with dual topology method
c
c
      subroutine epolar4f
      use atoms
      use energi
      use deriv
      use dlmda
      use limits
      use mutant
      use ost
      use polar
      use polpot
      use virial
      implicit none
      integer i,j
      real*8 ep1,ep0
      real*8 plambdaorig
      real*8 plambdaexp
      real*8 dplambdaexp
      real*8 d2plambdaexp
      real*8 epvir1(3,3)
      real*8 epvir0(3,3)
      real*8, allocatable :: dep1(:,:)
      real*8, allocatable :: dep0(:,:)
      character*6 mode
c
c
c     copy original plambda
c
      plambdaorig = plambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (dep0(3,n))
      allocate (dep1(3,n))
c
c     compute energy, force, and virial of the lambda = 0 state
c
      if (use_pol4i) then
         plambda = 0.0d0
         call altpolr
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            if (use_mlist) then
               call epolar1b
            else
               call epolar1a
            end if
         end if
         if (use_expol) then
            call dexpol
         end if
c
c     copy energy, force, and virial of the lambda = 0 state
c
         ep0 = ep
         do i = 1, n
            do j = 1, 3
               dep0(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir0(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     compute energy of the lambda = 1 state
c
      if (use_pol4f) then
         plambda = 1.0d0
         call altpolr
         if (use_ewald) then
            if (use_mlist) then
               call epolar1d
            else
               call epolar1c
            end if
         else
            if (use_mlist) then
               call epolar1b
            else
               call epolar1a
            end if
         end if
         if (use_expol) then
            call dexpol
         end if
c
c     copy energy, force, and virial of the lambda = 1 state
c
         ep1 = ep
         do i = 1, n
            do j = 1, 3
               dep1(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir1(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     copy energy, force, and virial if only one state is computed
c
      if (use_pol4i .and. .not.use_pol4f) then
         ep1 = ep0
         do i = 1, n
            do j = 1, 3
               dep1(j,i) = dep0(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir1(j,i) = epvir0(j,i)
            end do
         end do
      else if (.not.use_pol4i .and. use_pol4f) then
         ep0 = ep1
         do i = 1, n
            do j = 1, 3
               dep0(j,i) = dep1(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir0(j,i) = epvir1(j,i)
            end do
         end do
      end if
c
c     set original plambda
c
      plambda = plambdaorig
      call altelec
c
c     interpolate energy, force, and virial
c
      plambdaexp = plambda**epdtexp
      ep = plambdaexp * ep1 + (1.0d0 - plambdaexp) * ep0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = plambdaexp * dep1(j,i)
     &                 + (1.0d0 - plambdaexp) * dep0(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            epvir(j,i) = plambdaexp * epvir1(j,i)
     &                 + (1.0d0 - plambdaexp) * epvir0(j,i)
         end do
      end do
c
c     compute lambda derivative
c
      dplambdaexp = epdtexp * plambda**(epdtexp-1)
      if (epdtexp .gt. 1) then
         d2plambdaexp = dble(epdtexp*(epdtexp-1))*plambda**(epdtexp-2)
      else
         d2plambdaexp = 0.0d0
      end if
      depdl = dplambdaexp * (ep1 - ep0)
      depdl2 = d2plambdaexp * (ep1 - ep0)
      do i = 1, n
         do j = 1, 3
            dep(j,i) = plambdaexp * dep1(j,i)
     &                 + (1.0d0 - plambdaexp) * dep0(j,i)
            dfpdl(j,i) = dplambdaexp * (dep1(j,i) - dep0(j,i))
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = dplambdaexp * (epvir1(j,i) - epvir0(j,i))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dep0)
      deallocate (dep1)
      return
      end
