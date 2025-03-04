c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine rings  --  locate and store small rings  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "rings" searches the structure for small rings and stores
c     their constituent atoms, and optionally reduces rings into
c     their component smaller rings
c
c     note by default reducible rings are not removed since they
c     are needed for force field parameter assignment
c
c
      subroutine rings
      use angbnd
      use atoms
      use bitor
      use bndstr
      use couple
      use inform
      use iounit
      use ring
      use tettor
      use tors
      use tritor
      implicit none
      integer i,j,k,imin
      integer ia,ib,ic,id
      integer ie,ig,ih
      integer list1,list2
      integer list3,list4
      integer, allocatable :: list(:)
      logical reduce
c
c
c     zero out the number of small rings in the structure
c
      reduce = .false.
      nring3 = 0
      nring4 = 0
      nring5 = 0
      nring6 = 0
      nring7 = 0
c
c     parse to find bonds, angles, torsions and bitorsions
c
      if (nbond .eq. 0)  call bonds
      if (nangle .eq. 0)  call angles
      if (ntors .eq. 0)  call torsions
      if (nbitor .eq. 0)  call bitors
      if (ntritor .eq. 0)  call tritors
      if (ntettor .eq. 0)  call tettors
c
c     initial count of the total number of each ring size
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (ia .gt. ic) then
            ia = iang(3,i)
            ic = iang(1,i)
         end if
         imin = min(ia,ib,ic)
         if (ia.eq.imin .and. ib.lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ic) then
                  nring3 = nring3 + 1
                  goto 10
               end if
            end do
   10       continue
         end if
      end do
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (ia .gt. id) then
            ia = itors(4,i)
            ib = itors(3,i)
            ic = itors(2,i)
            id = itors(1,i)
         end if
         imin = min(ia,ib,ic,id)
         if (ia.eq.imin .and. ib.lt.id) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. id) then
                  nring4 = nring4 + 1
                  goto 20
               end if
            end do
   20       continue
         end if
      end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         if (ia .gt. ie) then
            ia = ibitor(5,i)
            ib = ibitor(4,i)
            id = ibitor(2,i)
            ie = ibitor(1,i)
         end if
         imin = min(ia,ib,ic,id,ie)
         if (ia.eq.imin .and. ib.le.ie) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ie) then
                  nring5 = nring5 + 1
                  goto 30
               end if
            end do
   30       continue
         end if
      end do
      do i = 1, ntritor
         ia = itritor(1,i)
         ib = itritor(2,i)
         ic = itritor(3,i)
         id = itritor(4,i)
         ie = itritor(5,i)
         ig = itritor(6,i)
         if (ia .gt. ig) then
            ia = itritor(6,i)
            ib = itritor(5,i)
            ic = itritor(4,i)
            id = itritor(3,i)
            ie = itritor(2,i)
            ig = itritor(1,i)
         end if
         imin = min(ia,ib,ic,id,ie,ig)
         if (ia.eq.imin .and. ib.lt.ig) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ig) then
                  nring6 = nring6 + 1
                  goto 40
               end if
            end do
   40       continue
         end if
      end do
      do i = 1, ntettor
         ia = itettor(1,i)
         ib = itettor(2,i)
         ic = itettor(3,i)
         id = itettor(4,i)
         ie = itettor(5,i)
         ig = itettor(6,i)
         ih = itettor(7,i)
         if (ia .gt. ih) then
            ia = itettor(7,i)
            ib = itettor(6,i)
            ic = itettor(5,i)
            ie = itettor(3,i)
            ig = itettor(2,i)
            ih = itettor(1,i)
         end if
         imin = min(ia,ib,ic,id,ie,ig,ih)
         if (ia.eq.imin .and. ib.lt.ih) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ih) then
                  nring7 = nring7 + 1
                  goto 50
               end if
            end do
   50       continue
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(iring3))  deallocate (iring3)
      if (allocated(iring4))  deallocate (iring4)
      if (allocated(iring5))  deallocate (iring5)
      if (allocated(iring6))  deallocate (iring6)
      if (allocated(iring7))  deallocate (iring7)
      allocate (iring3(3,nring3))
      allocate (iring4(4,nring4))
      allocate (iring5(5,nring5))
      allocate (iring6(6,nring6))
      allocate (iring7(7,nring7))
c
c     search for and store all of the 3-membered rings
c
      if (nring3 .ne. 0) then
         nring3 = 0
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (ia .gt. ic) then
               ia = iang(3,i)
               ic = iang(1,i)
            end if
            imin = min(ia,ib,ic)
            if (ia.eq.imin .and. ib.lt.ic) then
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. ic) then
                     nring3 = nring3 + 1
                     iring3(1,nring3) = ia
                     iring3(2,nring3) = ib
                     iring3(3,nring3) = ic
                     goto 60
                  end if
               end do
   60          continue
            end if
         end do
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     search for and store all of the 4-membered rings
c
      if (nring4 .ne. 0) then
         nring4 = 0
         do i = 1, n
            list(i) = 0
         end do
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (ia .gt. id) then
               ia = itors(4,i)
               ib = itors(3,i)
               ic = itors(2,i)
               id = itors(1,i)
            end if
            imin = min(ia,ib,ic,id)
            if (ia.eq.imin .and. ib.lt.id) then
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. id) then
                     nring4 = nring4 + 1
                     iring4(1,nring4) = ia
                     iring4(2,nring4) = ib
                     iring4(3,nring4) = ic
                     iring4(4,nring4) = id
                     if (reduce) then
                        list(ia) = nring4
                        list(ib) = nring4
                        list(ic) = nring4
                        list(id) = nring4
                        do k = 1, nring3
                           list1 = list(iring3(1,k))
                           list2 = list(iring3(2,k))
                           list3 = list(iring3(3,k))
                           if (list1.eq.nring4 .and. list2.eq.nring4
     &                            .and. list3.eq.nring4) then
                              nring4 = nring4 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              goto 70
                           end if
                        end do
                     end if
                     goto 70
                  end if
               end do
   70          continue
            end if
         end do
      end if
c
c     search for and store all of the 5-membered rings
c
      if (nring5 .ne. 0) then
         nring5 = 0
         do i = 1, n
             list(i) = 0
         end do
         do i = 1, nbitor
            ia = ibitor(1,i)
            ib = ibitor(2,i)
            ic = ibitor(3,i)
            id = ibitor(4,i)
            ie = ibitor(5,i)
            if (ia .gt. ie) then
               ia = ibitor(5,i)
               ib = ibitor(4,i)
               id = ibitor(2,i)
               ie = ibitor(1,i)
            end if
            imin = min(ia,ib,ic,id,ie)
            if (ia.eq.imin .and. ib.lt.ie) then
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. ie) then
                     nring5 = nring5 + 1
                     iring5(1,nring5) = ia
                     iring5(2,nring5) = ib
                     iring5(3,nring5) = ic
                     iring5(4,nring5) = id
                     iring5(5,nring5) = ie
                     if (reduce) then
                        list(ia) = nring5
                        list(ib) = nring5
                        list(ic) = nring5
                        list(id) = nring5
                        list(ie) = nring5
                        do k = 1, nring3
                           list1 = list(iring3(1,k))
                           list2 = list(iring3(2,k))
                           list3 = list(iring3(3,k))
                           if (list1.eq.nring5 .and. list2.eq.nring5
     &                            .and. list3.eq.nring5) then
                              nring5 = nring5 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              goto 80
                           end if
                        end do
                     end if
                     goto 80
                  end if
               end do
   80          continue
            end if
         end do
      end if
c
c     search for and store all of the 6-membered rings
c
      if (nring6 .ne. 0) then
         nring6 = 0
         do i = 1, n
            list(i) = 0
         end do
         do i = 1, ntritor
            ia = itritor(1,i)
            ib = itritor(2,i)
            ic = itritor(3,i)
            id = itritor(4,i)
            ie = itritor(5,i)
            ig = itritor(6,i)
            if (ia .gt. ig) then
               ia = itritor(6,i)
               ib = itritor(5,i)
               ic = itritor(4,i)
               id = itritor(3,i)
               ie = itritor(2,i)
               ig = itritor(1,i)
            end if
            imin = min(ia,ib,ic,id,ie,ig)
            if (ia.eq.imin .and. ib.lt.ig) then
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. ig) then
                     nring6 = nring6 + 1
                     iring6(1,nring6) = ia
                     iring6(2,nring6) = ib
                     iring6(3,nring6) = ic
                     iring6(4,nring6) = id
                     iring6(5,nring6) = ie
                     iring6(6,nring6) = ig
                     if (reduce) then
                        list(ia) = nring6
                        list(ib) = nring6
                        list(ic) = nring6
                        list(id) = nring6
                        list(ie) = nring6
                        list(ig) = nring6
                        do k = 1, nring3
                           list1 = list(iring3(1,k))
                           list2 = list(iring3(2,k))
                           list3 = list(iring3(3,k))
                           if (list1.eq.nring6 .and. list2.eq.nring6
     &                            .and. list3.eq.nring6) then
                              nring6 = nring6 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              goto 90
                           end if
                        end do
                        do k = 1, nring4
                           list1 = list(iring4(1,k))
                           list2 = list(iring4(2,k))
                           list3 = list(iring4(3,k))
                           list4 = list(iring4(4,k))
                           if (list1.eq.nring6 .and. list2.eq.nring6
     &                            .and. list3.eq.nring6
     &                            .and. list4.eq.nring6) then
                              nring6 = nring6 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              goto 90
                           end if
                        end do
                     end if
                     goto 90
                  end if
               end do
   90          continue
            end if
         end do
      end if
c
c     search for and store all of the 7-membered rings
c
      if (nring7 .ne. 0) then
         nring7 = 0
         do i = 1, n
            list(i) = 0
         end do
         do i = 1, ntettor
            ia = itettor(1,i)
            ib = itettor(2,i)
            ic = itettor(3,i)
            id = itettor(4,i)
            ie = itettor(5,i)
            ig = itettor(6,i)
            ih = itettor(7,i)
            if (ia .gt. ih) then
               ia = itettor(7,i)
               ib = itettor(6,i)
               ic = itettor(5,i)
               ie = itettor(3,i)
               ig = itettor(2,i)
               ih = itettor(1,i)
            end if
            imin = min(ia,ib,ic,id,ie,ig,ih)
            if (ia.eq.imin .and. ib.lt.ih) then
               do j = 1, n12(ia)
                  if (i12(j,ia) .eq. ih) then
                     nring7 = nring7 + 1
                     iring7(1,nring7) = ia
                     iring7(2,nring7) = ib
                     iring7(3,nring7) = ic
                     iring7(4,nring7) = id
                     iring7(5,nring7) = ie
                     iring7(6,nring7) = ig
                     iring7(7,nring7) = ih
                     if (reduce) then
                        list(ia) = nring7
                        list(ib) = nring7
                        list(ic) = nring7
                        list(id) = nring7
                        list(ie) = nring7
                        list(ig) = nring7
                        list(ih) = nring7
                        do k = 1, nring3
                           list1 = list(iring3(1,k))
                           list2 = list(iring3(2,k))
                           list3 = list(iring3(3,k))
                           if (list1.eq.nring7 .and. list2.eq.nring7
     &                            .and. list3.eq.nring7) then
                              nring7 = nring7 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              list(ih) = 0
                              goto 100
                           end if
                        end do
                        do k = 1, nring4
                           list1 = list(iring4(1,k))
                           list2 = list(iring4(2,k))
                           list3 = list(iring4(3,k))
                           list4 = list(iring4(4,k))
                           if (list1.eq.nring7 .and. list2.eq.nring7
     &                            .and. list3.eq.nring7
     &                            .and. list4.eq.nring7) then
                              nring7 = nring7 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              list(ih) = 0
                              goto 100
                           end if
                        end do
                     end if
                     goto 100
                  end if
               end do
  100          continue
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
c
c     print out lists of the small rings in the structure
c
      if (debug) then
         if (nring3 .gt. 0) then
            write (iout,110)
  110       format (/,' Three-Membered Rings in the Structure :',
     &              //,3x,'Ring',14x,'Atoms in Ring',/)
            do i = 1, nring3
               write (iout,120)  i,(iring3(j,i),j=1,3)
  120          format (i6,7x,3i7)
            end do
         end if
         if (nring4 .gt. 0) then
            write (iout,130)
  130       format (/,' Four-Membered Rings in the Structure :',
     &              //,3x,'Ring',17x,'Atoms in Ring',/)
            do i = 1, nring4
               write (iout,140)  i,(iring4(j,i),j=1,4)
  140          format (i6,7x,4i7)
            end do
         end if
         if (nring5 .gt. 0) then
            write (iout,150)
  150       format (/,' Five-Membered Rings in the Structure :',
     &              //,3x,'Ring',20x,'Atoms in Ring',/)
            do i = 1, nring5
               write (iout,160)  i,(iring5(j,i),j=1,5)
  160          format (i6,7x,5i7)
            end do
         end if
         if (nring6 .gt. 0) then
            write (iout,170)
  170       format (/,' Six-Membered Rings in the Structure :',
     &              //,3x,'Ring',23x,'Atoms in Ring',/)
            do i = 1, nring6
               write (iout,180)  i,(iring6(j,i),j=1,6)
  180          format (i6,7x,6i7)
            end do
         end if
         if (nring7 .gt. 0) then
            write (iout,190)
  190       format (/,' Seven-Membered Rings in the Structure :',
     &              //,3x,'Ring',26x,'Atoms in Ring',/)
            do i = 1, nring7
               write (iout,200)  i,(iring7(j,i),j=1,7)
  200          format (i6,7x,7i7)
            end do
         end if
      end if
      return
      end
