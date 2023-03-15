c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2022 by Moses KJ Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kxrepel  --  exch repulsion term assignment  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kxrepel" assigns the nuclear charge parameter and exponential
c     parameter for exchange repulsion interaction and processes any
c     new or changed values for these parameters
c
c
      subroutine kxrepel
      use atomid
      use atoms
      use inform
      use iounit
      use kxrepl
      use keys
      use mpole
      use potent
      use xrepel
      use repel
      use reppot
      use sizes
      implicit none
      integer i,j,k
      integer iboys
      integer freeunit
      integer ia,ic,next
      real*8 zpr,apr,cpr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     exit if using "Pauli" repulsion
c
      if (use_repel) return
      use_repel = .true.
c
c     process keywords containing exch repulsion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'XREPULSION ') then
            k = 0
            zpr = 0.0d0
            apr = 0.0d0
            cpr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  zpr,apr,cpr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Repulsion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Core',11x,'Damp',
     &                       6x,'P/S Coeff'/)
               end if
               if (k .le. maxclass) then
                  pxrz(k) = zpr
                  pxrdmp(k) = apr
                  pxrcr(k) = cpr
                  if (.not. silent) then
                     write (iout,30)  k,zpr,apr,cpr
   30                format (6x,i6,7x,3f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KXREPEL  --  Too many Exchange Repulsion',
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
      if (allocated(replist))  deallocate (replist)
      if (allocated(zpxr))  deallocate (zpxr)
      if (allocated(dmppxr))  deallocate (dmppxr)
      if (allocated(crpxr))  deallocate (crpxr)
      if (allocated(cpxr))  deallocate (cpxr)
      if (allocated(rcpxr))  deallocate (rcpxr)
      if (allocated(repole))  deallocate (repole)
      if (allocated(rrepole))  deallocate (rrepole)
      allocate (irep(n))
      allocate (replist(n))
      allocate (zpxr(n))
      allocate (dmppxr(n))
      allocate (crpxr(n))
      allocate (cpxr(4,n))
      allocate (rcpxr(4,n))
      allocate (repole(maxpole,n))
      allocate (rrepole(maxpole,n))
c
c     assign the core, alpha, and coefficient ratio parameters 
c
      do i = 1, n
         irep(i) = 0
         replist(i) = 0
         zpxr(i) = 0.0d0
         dmppxr(i) = 0.0d0
         crpxr(i) = 0.0d0
         ic = class(i)
         if (ic .ne. 0) then
            zpxr(i) = pxrz(ic)
            dmppxr(i) = pxrdmp(ic)
            crpxr(i) = pxrcr(ic)
         end if
      end do
c
c     process keywords containing atom specific exchange repulsion
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'XREPULSION ') then
            ia = 0
            zpr = 0.0d0
            apr = 0.0d0
            cpr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,zpr,apr,cpr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Repulsion Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',15x,'Core',11x,'Damp',
     &                       6x,'P/S Coeff'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,zpr,apr,cpr
   60             format (6x,i6,7x,3f15.4)
               end if
               zpxr(ia) = zpr
               dmppxr(ia) = apr
               crpxr(ia) = cpr
            end if
   70       continue
         end if
      end do
c
c     condense repulsion sites to the list of multipole sites
c
      nrep = 0
      if (use_repel) then
         do i = 1, n
            if (zpxr(i) .ne. 0) then
               nrep = nrep + 1
               irep(nrep) = i
               replist(i) = nrep
               do j = 1, maxpole
                  repole(j,i) = pole(j,i)
               end do
            end if
         end do
      end if
c
c     test multipoles at chiral sites and invert if necessary
c
      call chkpole
c
c     turn off the exchange repulsion potential if not used
c
      if (nrep .eq. 0)  use_repel = .false.
c
c     set repulsion type
c
      if (use_repel) reptyp = 'EXCHANGE'
      return
      end
