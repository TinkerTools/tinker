c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine krepel  --  Pauli repulsion term assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "krepel" assigns the size values, exponential parameter and
c     number of valence electrons for Pauli repulsion interactions
c     and processes any new or changed values for these parameters
c
c
      subroutine krepel
      use sizes
      use atomid
      use atoms
      use inform
      use iounit
      use krepl
      use keys
      use potent
      use repel
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 spr,apr,epr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing Pauli repulsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'REPULSION ') then
            k = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  spr,apr,epr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Pauli Repulsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Size',11x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (k .le. maxclass) then
                  prsiz(k) = spr
                  prdmp(k) = apr
                  prele(k) = -abs(epr)
                  if (.not. silent) then
                     write (iout,30)  k,spr,apr,epr
   30                format (6x,i6,7x,2f15.4,f15.3)
                  end if
               else
                  write (iout,40)
   40             format (/,' KREPEL  --  Too many Pauli Repulsion',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(sizpr))  deallocate (sizpr)
      if (allocated(dmppr))  deallocate (dmppr)
      if (allocated(elepr))  deallocate (elepr)
      allocate (sizpr(n))
      allocate (dmppr(n))
      allocate (elepr(n))
c
c     assign the repulsion size, alpha and valence parameters 
c     
      do i = 1, n
         ic = class(i)
         sizpr(i) = prsiz(ic)
         dmppr(i) = prdmp(ic)
         elepr(i) = prele(ic)
      end do
c
c     process keywords containing atom specific Pauli repulsion
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:10) .eq. 'REPULSION ') then
            ia = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,spr,apr,epr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Pauli Repulsion Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Size',12x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,spr,apr,epr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
               sizpr(ia) = spr
               dmppr(ia) = apr
               elepr(ia) = -abs(epr)
            end if
   70       continue
         end if
      end do
c
c     remove zero and undefined repulsion sites from the list
c
      nrep = 0
      do i = 1, n
         if (sizpr(i) .ne. 0.0d0) then
            nrep = nrep + 1
         end if
      end do
c
c     turn off the Pauli repulsion potential if not used
c
      if (nrep .eq. 0)  use_repuls = .false.
      return
      end
