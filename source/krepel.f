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
      use atomid
      use atoms
      use inform
      use iounit
      use krepl
      use keys
      use mpole
      use potent
      use repel
      use reppot
      use sizes
      implicit none
      integer i,j,k,ii
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
      if (allocated(irep))  deallocate (irep)
      if (allocated(sizpr))  deallocate (sizpr)
      if (allocated(dmppr))  deallocate (dmppr)
      if (allocated(elepr))  deallocate (elepr)
      if (allocated(repole))  deallocate (repole)
      if (allocated(rrepole))  deallocate (rrepole)
      allocate (irep(n))
      allocate (sizpr(n))
      allocate (dmppr(n))
      allocate (elepr(n))
      allocate (repole(maxpole,n))
      allocate (rrepole(maxpole,n))
c
c     assign the repulsion size, alpha and valence parameters 
c
      nrep = n
      do i = 1, n
         sizpr(i) = 0.0d0
         dmppr(i) = 0.0d0
         elepr(i) = 0.0d0
         ic = class(i)
         if (ic .ne. 0) then
            sizpr(i) = prsiz(ic)
            dmppr(i) = prdmp(ic)
            elepr(i) = prele(ic)
         end if
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
c     condense repulsion sites to the list of multipole sites
c
      if (use_repuls) then
         nrep = npole
         do ii = 1, npole
            i = ipole(ii)
            irep(ii) = i
            sizpr(ii) = sizpr(i)
            dmppr(ii) = dmppr(i)
            elepr(ii) = elepr(i)
            do j = 1, maxpole
               repole(j,ii) = pole(j,ii)
            end do
         end do
      end if
c
c     turn off the Pauli repulsion potential if not used
c
      if (nrep .eq. 0)  use_repuls = .false.
      return
      end
