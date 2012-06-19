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
      include 'atoms.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'keys.i'
      include 'mutant.i'
      integer i,j,k,ihyb
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
         alter(i) = .false.
      end do
      elmd = 1.0d0
      vlmd = 1.0d0
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
            nmut = 0
            lambda = 0.0d0
            string = record(next:120)
            read (string,*,err=20)  lambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=20)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
            alter(ihyb) = .true.
         else if (keyword(1:7) .eq. 'ELAMDA ') then
            string = record(next:120) 
            read (string,*,err=20)  elmd
         else if (keyword(1:7) .eq. 'VLAMDA ') then
            string = record(next:120)
            read (string,*,err=20)  vlmd
         else if (keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:120)
            read (string,*,err=10,end=10)  (list(k),k=1,20)
   10       continue
            k = 1
            do while (list(k) .ne. 0) 
               if (list(k) .gt. 0) then
                  alter(list(k)) = .true. 
                  k = k + 1
               else
                  do j = abs(list(k)), abs(list(k+1))
                     alter(j) = .true.
                  end do
                  k = k + 2
               end if
            end do
         end if
   20    continue
      end do
      if (elmd.ge.0.0d0 .and. elmd.lt.1.0d0)  call altelec
      if (elmd.ne.1.0d0 .and. vlmd.ne.1.0d0) then
         write (*,30)  elmd,vlmd
   30    format (' Free Energy Perturbation Step for Electrostatics : ',
     &           f7.4,2x,'for VDW : ',f7.4)
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
            if (alter(i)) then
               pchg(i) = pchg(i) * elmd
            end if
         end do
      end if
c
c     set electrostatic parameters for polarizable multipole models
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (alter(k)) then
               do j = 1, 13
                  pole (j,i) = pole(j,i) * elmd
               end do
            end if
         end do
         do i = 1, npolar
            if (alter(i)) then
               polarity(i) = polarity(i) * elmd
            end if
         end do
      end if
      return
      end
