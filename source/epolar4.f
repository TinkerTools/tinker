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
c     "epolar4" calculates the induced dipole polarization energy,
c     first derivatives with respect to Cartesian coordinates, and
c     derivatives with respect to the polarization lambda
c
c
      subroutine epolar4
      use dlmda
      use iounit
      use mplpot
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
c     compute polarization interactions, where exchange polarization
c     is already applied to each state by "epolar1calc"
c
      if (use_prst) then
         call epolar4s
      else
         call epolar4f
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
c     ##  subroutine epolar4s  --  single topology polar dU/dl  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar4s" calculates the single topology polarization energy
c     and Cartesian derivatives together with the first energy
c     derivative with respect to polarization lambda; second lambda,
c     force lambda and virial lambda derivatives are not computed
c
c     the lambda derivative is the explicit plambda derivative of the
c     stationary mutual polarization functional; the permanent field
c     derivative is exact because "altpolr" scales the permanent
c     electrostatic parameters of the charging ligand linearly, so it is
c     just the field of those multipoles held at their unscaled values,
c     and the other ligand of a staged leg is annihilated and adds
c     nothing; the inverse-polarizability derivative of the charging
c     ligand is evaluated from the local total fields, which remains
c     finite at plambda equal to zero
c
c
      subroutine epolar4s
      use atoms
      use chgpot
      use dlmda
      use extfld
      use limits
      use mpole
      use mutant
      use polar
      implicit none
      integer i,j,ii
      real*8 plmda
      real*8 fldd,fldp
      real*8 term
      real*8, allocatable :: dfldd(:,:)
      real*8, allocatable :: dfldp(:,:)
      real*8, allocatable :: field0(:,:)
      real*8, allocatable :: fieldp0(:,:)
      real*8, allocatable :: ufield(:,:)
      real*8, allocatable :: ufieldp(:,:)
      real*8 sc(0:2)
      real*8 dsc(0:2)
      logical atzero,same
c
c
c     compute the polarization energy and gradient at plambda; the
c     "gradient" routine only calls "epolar4" when polarization
c     follows the main lambda
c
      plmda = plambda
      call altepset (same)
      call epolar1calc
c
c     only the charging ligand carries a polarizability derivative
c
      call grpscale (plmda,sc,dsc)
c
c     derivative is the field of the mutated multipoles at the endpoint
c
      atzero = (plmda .eq. 0.0d0)
      allocate (dfldd(3,n))
      allocate (dfldp(3,n))
      allocate (ufield(3,n))
      allocate (ufieldp(3,n))
      call altepdt (1.0d0)
      mutfield = .true.
      if (use_ewald) then
         call dfield0c (dfldd,dfldp)
      else if (use_mlist) then
         call dfield0b (dfldd,dfldp)
      else
         call dfield0a (dfldd,dfldp)
      end if
      mutfield = .false.
      call altepdt (plmda)
c
c     get the mutual fields from the converged dipoles
c
      if (use_ewald) then
         call ufield0c (ufield,ufieldp)
      else if (use_mlist) then
         call ufield0b (ufield,ufieldp)
      else
         call ufield0a (ufield,ufieldp)
      end if
c
c     pol vanishes at plmda=0, so compute the permanent field
c
      if (atzero) then
         allocate (field0(3,n))
         allocate (fieldp0(3,n))
         if (use_ewald) then
            call dfield0c (field0,fieldp0)
         else if (use_mlist) then
            call dfield0b (field0,fieldp0)
         else
            call dfield0a (field0,fieldp0)
         end if
         if (use_exfld) then
            do ii = 1, npole
               i = ipole(ii)
               do j = 1, 3
                  field0(j,i) = field0(j,i) + texfld(j)
                  fieldp0(j,i) = fieldp0(j,i) + texfld(j)
               end do
            end do
         end if
      end if
c
c     differentiate the polarization function
c
      term = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            term = term + uinp(j,i)*dfldd(j,i) + uind(j,i)*dfldp(j,i)
            if (dsc(mutg(i)).ne.0.0d0 .and. douindorig(i)) then
               if (atzero) then
                  fldd = field0(j,i) + ufield(j,i)
                  fldp = fieldp0(j,i) + ufieldp(j,i)
               else
                  fldd = udir(j,i)/polarity(i) + ufield(j,i)
                  fldp = udirp(j,i)/polarity(i) + ufieldp(j,i)
               end if
               term = term + polarityorig(i)*fldd*fldp
            end if
         end do
      end do
      depdl = -0.5d0 * (electric/dielec) * term
c
c     perform deallocation of local arrays
c
      deallocate (dfldd)
      deallocate (dfldp)
      deallocate (ufield)
      deallocate (ufieldp)
      if (atzero) then
         deallocate (field0)
         deallocate (fieldp0)
      end if
c
c     restore the electrostatic parameter state for subsequent terms,
c     which is left at plambda above
c
      if (.not. same)  call alteprst
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
      use deriv
      use dlmda
      use energi
      use mutant
      use virial
      implicit none
      integer i,j
      real*8 weight1,dweight1,d2weight1
      real*8 ep0,plambdaorig
      real*8 epvir0(3,3)
      real*8, allocatable :: dep0(:,:)
      logical need0,need1
c
c
c     an endpoint is live when it carries weight or a lambda derivative
c
      plambdaorig = plambda
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &                 dpldlmda,d2pldlmda2,need0,need1)
c
c     build the needed endpoint states, then restore plambda
c
      allocate (dep0(3,n))
      call epolar1dt (need0,need1,ep0,dep0,epvir0)
      plambda = plambdaorig
      call alteprst
c
c     a single state carries all of the weight and no lambda
c     derivative, and is already in place
c
      depdl = 0.0d0
      if (need0 .and. need1) then
c
c     compute lambda derivative, along with the second, force and
c     virial lambda derivatives only when they are requested
c
         depdl = dweight1 * (ep - ep0)
         if (use_d2lmda) then
            d2epdl2 = d2weight1 * (ep - ep0)
            do i = 1, n
               do j = 1, 3
                  dfpdl(j,i) = dweight1 * (dep(j,i) - dep0(j,i))
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  depvirdl(j,i) = dweight1 * (epvir(j,i)-epvir0(j,i))
               end do
            end do
         end if
c
c     interpolate energy, force, and virial between the two states
c
         ep = weight1 * ep + (1.0d0 - weight1) * ep0
         do i = 1, n
            do j = 1, 3
               dep(j,i) = weight1 * dep(j,i)
     &                       + (1.0d0 - weight1) * dep0(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir(j,i) = weight1 * epvir(j,i)
     &                         + (1.0d0 - weight1) * epvir0(j,i)
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dep0)
      return
      end
