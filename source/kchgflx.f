c
c
c     ##########################################################
c     ##  COPYRIGHT (C) 2020 by Chengwen Liu & Jay W. Ponder  ##
c     ##                 All Rights Reserved                  ##
c     ##########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kchgflx  --  charge flux parameter assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kchgflx" assigns a force constant and ideal bond length
c     to each bond in the structure and processes any new or
c     changed parameter values
c
c     literature reference:
c
c     C. Liu, J.-P. Piquemal and P. Ren, "Implementation of Geometry-
c     Dependent Charge Flux into the Polarizable AMOEBA+ Potential",
c     Journal of Physical Chemistry Letters, 11, 419-426 (2020)
c
c
      subroutine kchgflx
      use sizes
      use angbnd
      use atomid
      use atoms
      use bndstr
      use cflux
      use inform
      use iounit
      use kangs
      use kbonds
      use kcflux
      use keys
      use potent
      use usage
      implicit none
      integer i,j,k
      integer ia,ib,ic
      integer ita,itb,itc
      integer na,nb
      integer size,next
      real*8 fc,bd,fj
      real*8 ja1,ja2
      real*8 jb1,jb2
      logical headerb
      logical headera
      character*4 pa,pb,pc
      character*8 blank8,pt,pt2
      character*8 ptab,ptbc
      character*12 blank12,pt3
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing charge flux parameters
c
      blank8 = '        '
      blank12 = '            '
      size = 4
      headerb = .true.
      headera = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'BNDCFLUX ') then
            ia = 0
            ib = 0
            fj = 0.0d0
            string = record(next:240)
            read (string,*,err=10,end=10) ia,ib,fj
   10       continue
            if (headerb .and. .not.silent) then
               headerb = .false.
               write (iout,20)
   20          format (/,' Additional Bond Charge Flux Parameters :',
     &                 //,5x,'Atom Classes',19x,'K(CFB)',/)
            end if
            if (.not. silent) then
               write (iout,30)  ia,ib,fj
   30          format (6x,2i4,13x,f15.6)
            end if
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt2 = pa//pb
            else
               pt2 = pb//pa
            end if
            do j = 1, maxnbcf
               if (kcfb(j).eq.blank8 .or. kcfb(j).eq.pt2) then
                  kcfb(j) = pt2
                  jbnd(j) = fj
                  goto 40
               end if
            end do
   40       continue
         else if (keyword(1:9) .eq. 'ANGCFLUX ') then
            ia = 0
            ib = 0
            ic = 0
            ja1 = 0.0d0
            ja2 = 0.0d0
            jb1 = 0.0d0
            jb2 = 0.0d0
            string = record(next:240)
            read (string,*,err=50,end=50) ia,ib,ic,ja1,ja2,jb1,jb2
   50       continue
            if (headera .and. .not.silent) then
               headera = .false.
               write (iout,60)
   60          format (/,' Additional Angle Charge Flux Parameters :',
     &                 //,5x,'Atom Classes',10x,'K(CFA1)',
     &                    7x,'K(CFA2)',7x,'K(CFB1)',7x,'K(CFB2)',/)
            end if
            if (.not. silent) then
               write (iout,70)  ia,ib,ic,ja1,ja2,jb1,jb2
   70          format (4x,3i4,4x,4f14.6)
            end if
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            call numeral (ic,pc,size)
            if (ia .le. ic) then
               pt3 = pa//pb//pc
            else
               pt3 = pc//pb//pa
            end if
            do j = 1, maxnacf
               if (kcfa(j).eq.blank12 .or. kcfa(j).eq.pt3) then
                  kcfa(j) = pt3
                  jthetal(1,j) = ja1
                  jthetal(2,j) = ja2
                  jbpl(1,j) = jb1
                  jbpl(2,j) = jb2
                  goto 80
               end if
            end do
   80       continue
         end if
      end do
c
c     determine the total number of forcefield parameters
c
      nb = maxnbcf
      do i = maxnbcf, 1, -1
         if (kcfb(i) .eq. blank8)  nb = i - 1
      end do
      na = maxnacf
      do i = maxnacf, 1, -1
         if (kcfa(i) .eq. blank12)  na = i - 1
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(b0))  deallocate (b0)
      if (allocated(jb))  deallocate (jb)
      allocate (b0(nbond))
      allocate (jb(nbond))
c
c     assign bond charge flux parameters for each bond
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         ita = class(ia)
         itb = class(ib)
         size = 4
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         if (ita .le. itb) then
            pt2 = pa//pb
         else
            pt2 = pb//pa
         end if
         b0(i) = 0.0d0
         jb(i) = 0.0d0
         do j = 1, nb
            if (kcfb(j) .eq. pt2) then
               jb(i) = jbnd(j)
            end if
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(theta0))  deallocate (theta0)
      if (allocated(bp0))  deallocate (bp0)
      if (allocated(jbp))  deallocate (jbp)
      if (allocated(jtheta))  deallocate (jtheta)
      allocate (theta0(nangle))
      allocate (bp0(2,nangle))
      allocate (jbp(2,nangle))
      allocate (jtheta(2,nangle))
c
c    assign angle charge flux parameters for each angle
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         ita = class(ia)
         itb = class(ib)
         itc = class(ic)
         call numeral (ita,pa,size)
         call numeral (itb,pb,size)
         call numeral (itc,pc,size)
         if (ita .le. itc) then
            pt3 = pa//pb//pc
         else
            pt3 = pc//pb//pa
         end if
         if (ita .le. itb) then
            ptab = pa//pb
         else
            ptab = pb//pa
         end if
         if (itb .le. itc) then
            ptbc = pb//pc
         else
            ptbc = pc//pb
         end if
         bp0(1,i) = 0.0d0
         bp0(2,i) = 0.0d0
         jbp(1,i) = 0.0d0
         jbp(2,i) = 0.0d0
         jtheta(1,i) = 0.0d0
         jtheta(2,i) = 0.0d0
         do j = 1, na
            if (kcfa(j) .eq. pt3) then
               jtheta(1,i) = jthetal(1,j)
               jtheta(2,i) = jthetal(2,j)
               jbp(1,i) = jbpl(1,j)
               jbp(2,i) = jbpl(2,j)
            end if
            do k = 1, nbond
               if (kb(k) .eq. ptab) then
                  bp0(1,i) = blen(k)
               end if
               if (kb(k) .eq. ptbc) then
                  bp0(2,i) = blen(k)
               end if
            end do
         end do
      end do
c
c     turn off charge flux if bonds and angles are not used
c
      if (nb.eq.0 .and. na.eq.0)  use_chgflx = .false.
      return
      end
