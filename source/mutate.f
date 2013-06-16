c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c
      subroutine mutate
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      include 'mutant.i'
      integer i,j,k
      integer ihyb,ia
      integer it0,it1,next
      integer list(20)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     zero number of hybrid atoms and hybrid atom list
c
      nmut = 0
      do i = 1, n
         mut(i) = .false.
      end do
      lambda = 1.0d0
      vlambda = 1.0d0
      elambda = 1.0d0
      scexp = 5.0d0
      scalpha = 0.7d0
      do i = 1, 20
         list(i) = 0
      end do
c
c     search keywords for free energy perturbation options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  lambda
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  vlambda
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:120) 
            read (string,*,err=20)  elambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=20)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:120)
            read (string,*,err=10,end=10)  (list(k),k=1,20)
   10       continue
            k = 1
            do while (list(k) .ne. 0) 
               if (list(k) .gt. 0) then
                  ia = list(k)
                  nmut = nmut + 1
                  imut(nmut) = ia
                  mut(ia) = .true.
                  type0(nmut) = 0
                  type1(nmut) = type(ia)
                  class0(nmut) = 0
                  class1(nmut) = class(ia)
                  k = k + 1
               else
                  do j = abs(list(k)), abs(list(k+1))
                     nmut = nmut + 1
                     imut(nmut) = i
                     mut(j) = .true.
                     type0(nmut) = 0
                     type1(nmut) = type(i)
                     class0(nmut) = 0
                     class1(nmut) = class(i)
                  end do
                  k = k + 2
               end if
            end do
         end if
   20    continue
      end do
c
c     scale electrostatic parameter values based on lambda
c
      if (elambda.ge.0.0d0 .and. elambda.lt.1.0d0)  call altelec
c
c     write the status of the current free energy perturbation step
c
      if (nmut .ne. 0) then
         write (iout,30)  vlambda
   30    format (/,' Free Energy Perturbation :',f15.3,
     &              ' Lambda for van der Waals')
         write (iout,40)  elambda
   40    format (' Free Energy Perturbation :',f15.3,
     &              ' Lambda for Electrostatics')
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altelec  --  mutated electrostatic parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altelec" constructs the mutated electrostatic parameters
c     based on the lambda mutation parameter "elmd"
c
c
      subroutine altelec
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'mpole.i'
      include 'mutant.i'
      include 'polar.i'
      include 'potent.i'
      integer i,j,k
c
c
c     set electrostatic parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            if (mut(i)) then
               pchg(i) = pchg(i) * elambda
            end if
         end do
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole (j,i) = pole(j,i) * elambda
               end do
            end if
         end do
         do i = 1, npolar
            if (mut(i)) then
               polarity(i) = polarity(i) * elambda
            end if
         end do
      end if
      return
      end
