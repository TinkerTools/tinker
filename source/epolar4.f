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
      if (use_rel) then
         call epolar4fr
      else if (use_past) then
         call epolar4s
      else
         call epolar4f
      end if
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
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar4s  --  single topology polar dU/dl  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar4s" calculates the absolute single topology polarization
c     energy and Cartesian derivatives together with the first energy
c     derivative with respect to polarization lambda; second lambda,
c     force lambda and virial lambda derivatives are not computed
c
c
      subroutine epolar4s
      use bound
      use iounit
      use mutant
      implicit none
      real*8 plambdaorig
c
c
c     add the scalar derivative with respect to plambda
c
      plambdaorig = plambda
      call altepdt (plambdaorig)
      call epolar1calc
      call epolar4d (plambdaorig)
c
c     restore the electrostatic parameter state for subsequent terms
c
      call alteprst
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine epolar4d  --  single topology polar dU/dl  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "epolar4d" evaluates the explicit plambda derivative of the
c     stationary mutual polarization functional.  The permanent field
c     derivative is exact because "altpolr" scales every mutated
c     permanent electrostatic parameter linearly, so it is just the
c     field of the mutated multipoles held at their unscaled values.
c     The inverse-polarizability derivative is evaluated from the local
c     total fields, which remains finite at plambda equal to zero.
c
c
      subroutine epolar4d (plmda)
      use atoms
      use chgpot
      use dlmda
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
      logical atzero
c
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
      call dfield0a (dfldd,dfldp)
      mutfield = .false.
      call altepdt (plmda)
c
c     get the mutual fields from the converged dipoles
c
      call ufield0a (ufield,ufieldp)
c
c     pol vanishes at plmda=0, so compute the permanent field
c
      if (atzero) then
         allocate (field0(3,n))
         allocate (fieldp0(3,n))
         call dfield0a (field0,fieldp0)
      end if
c
c     differentiate the polarization function
c
      term = 0.0d0
      do ii = 1, npole
         i = ipole(ii)
         do j = 1, 3
            term = term + uinp(j,i)*dfldd(j,i) + uind(j,i)*dfldp(j,i)
            if (mut(i) .and. douindorig(i)) then
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
      use polar
      use polpot
      use potent
      use virial
      implicit none
      real*8 weight1,dweight1,d2weight1
      logical need0,need1
      integer i,j
      real*8 ep1,ep0
      real*8 plambdaorig
      real*8 elambdaorig
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
      elambdaorig = elambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (dep0(3,n))
      allocate (dep1(3,n))
c
c     compute energy, force, and virial of the lambda = 0 state
c
c
c     an endpoint is live when it carries weight or a lambda derivative
c
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &                 dpldlmda,d2pldlmda2,need0,need1)
      if (need0) then
         call altepdt (0.0d0)
         call epolar1calc
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
      if (need1) then
         call altepdt (1.0d0)
         call epolar1calc
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
      if (need0 .and. .not.need1) then
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
      else if (.not.need0 .and. need1) then
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
      if (use_mpole) then
         call altemdt (elambdaorig)
      else
         call altepdt (plambdaorig)
      end if
c
c     interpolate energy, force, and virial
c
      ep = weight1 * ep1 + (1.0d0 - weight1) * ep0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = weight1 * dep1(j,i)
     &                 + (1.0d0 - weight1) * dep0(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            epvir(j,i) = weight1 * epvir1(j,i)
     &                 + (1.0d0 - weight1) * epvir0(j,i)
         end do
      end do
c
c     compute lambda derivative
c
      depdl = dweight1 * (ep1 - ep0)
      d2epdl2 = d2weight1 * (ep1 - ep0)
      do i = 1, n
         do j = 1, 3
            dep(j,i) = weight1 * dep1(j,i)
     &                 + (1.0d0 - weight1) * dep0(j,i)
            dfpdl(j,i) = dweight1 * (dep1(j,i) - dep0(j,i))
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = dweight1 * (epvir1(j,i) - epvir0(j,i))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dep0)
      deallocate (dep1)
      return
      end
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar4fr  --  relative dual topo polar dlmda  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar4fr" interpolates between the two coupling states of a
c     two-ligand relative dual topology calculation, each state a sum
c     of parameter-zeroed subsystem energies,
c
c        E = weight1*E(prelst1) + (1-weight1)*E(prelst0)
c
c
      subroutine epolar4fr
      use atoms
      use deriv
      use dlmda
      use energi
      use mutant
      use virial
      implicit none
      real*8 weight1,dweight1,d2weight1
      integer i,j,k
      real*8 ep0,ep1
      real*8 epvir0(3,3),epvir1(3,3)
      logical la,lb,le
      logical in0,in1
      logical need0,need1
      real*8, allocatable :: dep0(:,:)
      real*8, allocatable :: dep1(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dep0(3,n))
      allocate (dep1(3,n))
c
c     an endpoint is live when it carries weight or a lambda derivative
c
      call relpowerwt (plambda,epdtexp,weight1,dweight1,d2weight1)
      call relneed (weight1,dweight1,d2weight1,
     &                 dpldlmda,d2pldlmda2,need0,need1)
c
c     zero out the two endpoint accumulators
c
      ep0 = 0.0d0
      ep1 = 0.0d0
      do i = 1, n
         do j = 1, 3
            dep0(j,i) = 0.0d0
            dep1(j,i) = 0.0d0
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            epvir0(j,i) = 0.0d0
            epvir1(j,i) = 0.0d0
         end do
      end do
c
c     build each subsystem once, add to the endpoints
c
      do k = 1, nrelsub
         call relslot (k,prelst0,prelst1,la,lb,le,in0,in1)
         in0 = in0 .and. need0
         in1 = in1 .and. need1
         if (.not. (in0 .or. in1))  cycle
         call altpolrsub (la,lb,le)
         call epolar1calc
         if (in0) then
            ep0 = ep0 + ep
            do i = 1, n
               do j = 1, 3
                  dep0(j,i) = dep0(j,i) + dep(j,i)
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  epvir0(j,i) = epvir0(j,i) + epvir(j,i)
               end do
            end do
         end if
         if (in1) then
            ep1 = ep1 + ep
            do i = 1, n
               do j = 1, 3
                  dep1(j,i) = dep1(j,i) + dep(j,i)
               end do
            end do
            do i = 1, 3
               do j = 1, 3
                  epvir1(j,i) = epvir1(j,i) + epvir(j,i)
               end do
            end do
         end if
      end do
c
c     restore the original full system parameters
c
      call altpolrsub (.true.,.true.,.true.)
c
c     copy energy if only one endpoint state is computed
c
      if (.not. need0) then
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
      else if (.not. need1) then
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
      end if
c
c     interpolate between the two endpoint states
c
      ep = weight1*ep1 + (1.0d0-weight1)*ep0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = weight1*dep1(j,i) + (1.0d0-weight1)*dep0(j,i)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            epvir(j,i) = weight1*epvir1(j,i)
     &                + (1.0d0-weight1)*epvir0(j,i)
         end do
      end do
c
c     interpolate the lambda derivative
c
      depdl = dweight1 * (ep1-ep0)
      d2epdl2 = d2weight1 * (ep1-ep0)
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = dweight1 * (epvir1(j,i)-epvir0(j,i))
         end do
      end do
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = dweight1 * (dep1(j,i)-dep0(j,i))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dep0)
      deallocate (dep1)
      return
      end
