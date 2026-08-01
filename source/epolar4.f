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
         if (use_relstage) then
            call epolar4frs
         else
            call epolar4fr
         end if
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
      integer i,j
      real*8 ep1,ep0
      real*8 plambdaorig
      real*8 elambdaorig
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
      elambdaorig = elambda
c
c     perform dynamic allocation of some local arrays
c
      allocate (dep0(3,n))
      allocate (dep1(3,n))
c
c     compute energy, force, and virial of the lambda = 0 state
c
      if (use_pol4i) then
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
      if (use_pol4f) then
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
      if (use_mpole) then
         call altemdt (elambdaorig)
      else
         call altepdt (plambdaorig)
      end if
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
      d2epdl2 = d2plambdaexp * (ep1 - ep0)
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
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine epolar4fr  --  relative dual topology pol deriv  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "epolar4fr" calculates the polarization energy, Cartesian gradient
c     and lambda derivatives for a two-ligand relative dual topology
c     calculation by combining four self-consistent subsystem states,
c     E1 = E(A+env) + E(B) and E0 = E(B+env) + E(A), which are
c     independent of plambda so the lambda derivatives follow directly
c     from the interpolation weight as in "epolar4f"
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
      integer i,j
      real*8 epae,epbe
      real*8 epa,epb
      real*8 ep1,ep0
      real*8 plambdaexp
      real*8 dplambdaexp,d2plambdaexp
      real*8 epvirae(3,3),epvirbe(3,3)
      real*8 epvira(3,3),epvirb(3,3)
      real*8, allocatable :: depae(:,:)
      real*8, allocatable :: depbe(:,:)
      real*8, allocatable :: depa(:,:)
      real*8, allocatable :: depb(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (depae(3,n))
      allocate (depbe(3,n))
      allocate (depa(3,n))
      allocate (depb(3,n))
c
c     ligand A coupled to environment, group B fully decoupled
c
      if (use_pol4f) then
         call altpolrsub (.true.,.false.,.true.)
         call epolar1calc
         epae = ep
         do i = 1, n
            do j = 1, 3
               depae(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvirae(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     ligand B coupled to environment, group A fully decoupled
c
      if (use_pol4i) then
         call altpolrsub (.false.,.true.,.true.)
         call epolar1calc
         epbe = ep
         do i = 1, n
            do j = 1, 3
               depbe(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvirbe(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     ligand A alone, giving its intramolecular polarization energy
c
      if (use_pol4i) then
         call altpolrsub (.true.,.false.,.false.)
         call epolar1calc
         epa = ep
         do i = 1, n
            do j = 1, 3
               depa(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvira(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     ligand B alone, giving its intramolecular polarization energy
c
      if (use_pol4f) then
         call altpolrsub (.false.,.true.,.false.)
         call epolar1calc
         epb = ep
         do i = 1, n
            do j = 1, 3
               depb(j,i) = dep(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvirb(j,i) = epvir(j,i)
            end do
         end do
      end if
c
c     alias the omitted composite endpoint to the computed endpoint
c
      if (use_pol4i .and. .not.use_pol4f) then
         epae = epbe
         epb = epa
         do i = 1, n
            do j = 1, 3
               depae(j,i) = depbe(j,i)
               depb(j,i) = depa(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvirae(j,i) = epvirbe(j,i)
               epvirb(j,i) = epvira(j,i)
            end do
         end do
      else if (.not.use_pol4i .and. use_pol4f) then
         epbe = epae
         epa = epb
         do i = 1, n
            do j = 1, 3
               depbe(j,i) = depae(j,i)
               depa(j,i) = depb(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvirbe(j,i) = epvirae(j,i)
               epvira(j,i) = epvirb(j,i)
            end do
         end do
      end if
c
c     restore full system and interpolate the dual topology result
c
      call altpolrsub (.true.,.true.,.true.)
      plambdaexp = plambda**epdtexp
      ep1 = epae + epb
      ep0 = epbe + epa
      ep = plambdaexp*ep1 + (1.0d0-plambdaexp)*ep0
      do i = 1, n
         do j = 1, 3
            dep(j,i) = plambdaexp*(depae(j,i)+depb(j,i))
     &               + (1.0d0-plambdaexp)*(depbe(j,i)+depa(j,i))
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            epvir(j,i) = plambdaexp*(epvirae(j,i)+epvirb(j,i))
     &                 + (1.0d0-plambdaexp)*(epvirbe(j,i)+epvira(j,i))
         end do
      end do
c
c     lambda derivatives of energy, gradient and virial
c
      dplambdaexp = epdtexp * plambda**(epdtexp-1)
      if (epdtexp .gt. 1) then
         d2plambdaexp = dble(epdtexp*(epdtexp-1))*plambda**(epdtexp-2)
      else
         d2plambdaexp = 0.0d0
      end if
      depdl = dplambdaexp * (ep1 - ep0)
      d2epdl2 = d2plambdaexp * (ep1 - ep0)
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = dplambdaexp
     &         * ((depae(j,i)+depb(j,i))-(depbe(j,i)+depa(j,i)))
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = dplambdaexp
     &         * ((epvirae(j,i)+epvirb(j,i))-(epvirbe(j,i)+epvira(j,i)))
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (depae)
      deallocate (depbe)
      deallocate (depa)
      deallocate (depb)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine epolar4frs  --  staged rel dual topo pol deriv  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "epolar4frs" calculates the polarization energy, Cartesian
c     gradient and lambda derivatives for a two-ligand relative dual
c     topology calculation run on the staged schedule, combining the
c     same subsystem states as "epolar1frs", which are themselves
c     independent of plambda, so the sublambda derivatives follow
c     directly from the interpolation weight
c
c
      subroutine epolar4frs
      use atoms
      use deriv
      use dlmda
      use energi
      use mutant
      use virial
      implicit none
      integer i,j
      real*8 ep0,ep1
      real*8 weight1,weight0
      real*8 epvir0(3,3),epvir1(3,3)
      logical lig1,domix,dovdwm,needref
      real*8, allocatable :: dep0(:,:)
      real*8, allocatable :: dep1(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (dep0(3,n))
      allocate (dep1(3,n))
c
c     decide which leg of the staged schedule is active
c
      lig1 = (relstage .eq. 'LIG1')
      dovdwm = (relstage .eq. 'VDWM')
      domix = relstagemix
      needref = (dovdwm .or. domix)
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
c     environment alone, part of the decoupled reference
c
      if (needref) then
         call altpolrsub (.false.,.false.,.true.)
         call epolar1calc
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
c
c     ligand A alone, in the reference and in the ligand 0 endpoint
c
      if (needref .or. (.not.dovdwm .and. .not.lig1)) then
         call altpolrsub (.true.,.false.,.false.)
         call epolar1calc
         if (needref) then
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
         if (.not.dovdwm .and. .not.lig1) then
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
      end if
c
c     ligand B alone, in the reference and in the ligand 1 endpoint
c
      if (needref .or. (.not.dovdwm .and. lig1)) then
         call altpolrsub (.false.,.true.,.false.)
         call epolar1calc
         if (needref) then
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
         if (.not.dovdwm .and. lig1) then
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
      end if
c
c     the active ligand coupled to the environment
c
      if (.not. dovdwm) then
         if (lig1) then
            call altpolrsub (.true.,.false.,.true.)
         else
            call altpolrsub (.false.,.true.,.true.)
         end if
         call epolar1calc
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
c
c     restore the original full system parameters
c
      call altpolrsub (.true.,.true.,.true.)
c
c     interpolate the active leg, or take the reference in the middle
c
      if (dovdwm) then
         ep = ep0
         do i = 1, n
            do j = 1, 3
               dep(j,i) = dep0(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir(j,i) = epvir0(j,i)
            end do
         end do
      else if (domix) then
         weight1 = plambda
         weight0 = 1.0d0 - weight1
         ep = weight1*ep1 + weight0*ep0
         do i = 1, n
            do j = 1, 3
               dep(j,i) = weight1*dep1(j,i) + weight0*dep0(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir(j,i) = weight1*epvir1(j,i) + weight0*epvir0(j,i)
            end do
         end do
      else
         ep = ep1
         do i = 1, n
            do j = 1, 3
               dep(j,i) = dep1(j,i)
            end do
         end do
         do i = 1, 3
            do j = 1, 3
               epvir(j,i) = epvir1(j,i)
            end do
         end do
      end if
c
c     the interpolation weight is linear in plambda, so its first
c     derivative is one and its second derivative vanishes; on the
c     middle window and at a flat leg end there is nothing to mix and
c     every lambda derivative is zero, as "dpldlmda" is zero there too
c
      depdl = 0.0d0
      d2epdl2 = 0.0d0
      do i = 1, 3
         do j = 1, 3
            depvirdl(j,i) = 0.0d0
         end do
      end do
      do i = 1, n
         do j = 1, 3
            dfpdl(j,i) = 0.0d0
         end do
      end do
      if (domix) then
         depdl = ep1 - ep0
         do i = 1, 3
            do j = 1, 3
               depvirdl(j,i) = epvir1(j,i) - epvir0(j,i)
            end do
         end do
         do i = 1, n
            do j = 1, 3
               dfpdl(j,i) = dep1(j,i) - dep0(j,i)
            end do
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dep0)
      deallocate (dep1)
      return
      end
