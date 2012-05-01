c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##  Modified by Chuanjie Wu in 2009              ##
c     ##              All Rights Reserved              ##
c     ###################################################
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
      integer liglist(20)
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
         liglist(i) = 0
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
            read (string,*,err=10)  lambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=10)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
            alter(ihyb) = .true.
         else if(keyword(1:7) .eq. 'ELAMDA ') then
            string = record(next:120) 
            read (string,*,err=10) elmd
         else if(keyword(1:7) .eq. 'VLAMDA ') then
            string = record(next:120)
            read (string,*,err=10) vlmd
         else if(keyword(1:7) .eq. 'LIGAND ') then
            string = record(next:120)
            read (string,*,err=20,end=20) (liglist(k),k=1,20)
  20        continue
            k = 1
            dowhile (liglist(k) .ne. 0) 
               if (liglist(k) .gt. 0) then
                  alter(liglist(k)) = .true. 
                  write(*,*) 'Ligand Atom', liglist(k) 
                  k = k + 1
               else 
                  do j = abs(liglist(k)), abs(liglist(k+1))
                     alter(j) = .true.
                     write(*,*) 'Ligand Atom', j
                  end do   
                  k = k + 2
               end if
            end do   
         end if
   10    continue
      end do
      if (elmd .ge. 0.0d0 .and. elmd .lt. 1.0d0) call altelec
      if (elmd .ne. 1.0d0 .and. vlmd .ne. 1.0d0) write(*,30) elmd,vlmd
   30 format(1x,'Free Energy Perturbation step for Electrostatic: ',
     &       f7.4,' for VDW : ',f7.4)
      return
      end

c     ###############################################################
c     ##                                                           ##
c     ##  subroutine altelec -- alter electrostatic parameters     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "altelec" constructs the mutated electrostatic parameters
c     electrostatic mutation parameter "elmd (ELAMDA)"
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
      integer ii,jj,kk
 
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

      if (use_charge) then 
         do i = 1, nion
            if (alter(i)) then
               pchg(i) = pchg(i) * elmd
            end if
         end do
      end if

      return
      end          
  
