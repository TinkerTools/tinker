c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutate  --  set the chimeric atoms and values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutate" sets the list of chimeric atoms and mutational
c     parameters that are used during a free energy calculation
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
      integer maxlist
      parameter (maxlist=50)
      integer i,j,k,ia
      integer it0,it1,next
      integer list(maxlist)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     zero out the number of mutated atoms and atom list
c
      nmut = 0
      do i = 1, n
         mut(i) = .false.
      end do
c
c     search keywords for free energy perturbation options
c
      do j = 1, nkey
         next = 1
         record = keyline(j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            lambda = 0.0d0
            string = record(next:120)
            read (string,*,err=20)  lambda
            vlambda = lambda
            clambda = lambda
            dlambda = lambda
            mlambda = lambda
            plambda = lambda
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  vlambda
         else if (keyword(1:11) .eq. 'CHG-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  clambda
         else if (keyword(1:11) .eq. 'DPL-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  dlambda
         else if (keyword(1:13) .eq. 'MPOLE-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  mlambda
         else if (keyword(1:13) .eq. 'POLAR-LAMBDA ') then
            string = record(next:120)
            read (string,*,err=20)  plambda
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:120)
            read (string,*,err=20)  ia,it0,it1
            nmut = nmut + 1
            imut(nmut) = ia
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
            mut(ia) = .true.
         else if (keyword(1:7) .eq. 'LIGAND ') then
            do i = 1, maxlist
               list(i) = 0
            end do
            string = record(next:120)
            read (string,*,err=10,end=10)  (list(i),i=1,maxlist)
   10       continue
            k = 1
            do while (list(k) .ne. 0)
               ia = list(k)
               if (list(k) .gt. 0) then
                  if (.not. mut(ia)) then
                     nmut = nmut + 1
                     imut(nmut) = ia
                     type0(nmut) = 0
                     type1(nmut) = type(ia)
                     class0(nmut) = 0
                     class1(nmut) = class(ia)
                     mut(ia) = .true.
                  end if
                  k = k + 1
               else
                  do i = abs(list(k)), abs(list(k+1))
                     if (.not. mut(i)) then
                        nmut = nmut + 1
                        imut(nmut) = i
                        type0(nmut) = 0
                        type1(nmut) = type(i)
                        class0(nmut) = 0
                        class1(nmut) = class(i)
                        mut(i) = .true.
                     end if
                  end do
                  k = k + 2
               end if
            end do
         end if
   20    continue
      end do
      return
      end
