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
      if (allocated(ixrep))  deallocate (ixrep)
      if (allocated(xreplist))  deallocate (xreplist)
      if (allocated(zpxr))  deallocate (zpxr)
      if (allocated(dmppxr))  deallocate (dmppxr)
      if (allocated(crpxr))  deallocate (crpxr)
      if (allocated(cpxr))  deallocate (cpxr)
      if (allocated(rcpxr))  deallocate (rcpxr)
      if (allocated(xrepole))  deallocate (xrepole)
      allocate (ixrep(n))
      allocate (xreplist(n))
      allocate (zpxr(n))
      allocate (dmppxr(n))
      allocate (crpxr(n))
      allocate (cpxr(4,n))
      allocate (rcpxr(4,n))
      allocate (xrepole(4,n))
c
c     assign the core, alpha, and coefficient ratio parameters 
c
      do i = 1, n
         ixrep(i) = 0
         xreplist(i) = 0
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
      nxrep = 0
      if (use_xrepel) then
         do i = 1, n
            if (zpxr(i) .ne. 0) then
               nxrep = nxrep + 1
               ixrep(nxrep) = i
               xreplist(i) = nxrep
               do j = 1, 4
                  xrepole(j,i) = pole(j,i)
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
      if (nxrep .eq. 0)  use_xrepel = .false.
      return
      end
