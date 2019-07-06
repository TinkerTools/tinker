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
c     their constituent atoms, and optionally reduces large rings
c     into their component smaller rings
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
      use tors
      implicit none
      integer i,j,k,m
      integer kk,imax
      integer ia,ib,ic,id
      integer ie,ig,ih
      integer list1,list2
      integer list3,list4
      integer maxring
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
c
c     perform dynamic allocation of some global arrays
c
      maxring = 10000
      if (.not. allocated(iring3))  allocate (iring3(3,maxring))
      if (.not. allocated(iring4))  allocate (iring4(4,maxring))
      if (.not. allocated(iring5))  allocate (iring5(5,maxring))
      if (.not. allocated(iring6))  allocate (iring6(6,maxring))
      if (.not. allocated(iring7))  allocate (iring7(7,maxring))
c
c     search for and store all of the 3-membered rings
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (ib.lt.ia .and. ib.lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ic) then
                  nring3 = nring3 + 1
                  if (nring3 .gt. maxring) then
                     write (iout,10)
   10                format (/,' RINGS  --  Too many 3-Membered Rings;',
     &                          ' Increase MAXRING')
                     call fatal
                  end if
                  iring3(1,nring3) = ia
                  iring3(2,nring3) = ib
                  iring3(3,nring3) = ic
                  goto 20
               end if
            end do
   20       continue
         end if
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     search for and store all of the 4-membered rings
c
      do i = 1, n
         list(i) = 0
      end do
      do i = 1, ntors
         ia = itors(1,i)
         ib = itors(2,i)
         ic = itors(3,i)
         id = itors(4,i)
         if (ia.lt.ic .and. id.lt.ib) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. id) then
                  nring4 = nring4 + 1
                  if (nring4 .gt. maxring) then
                     write (iout,30)
   30                format (/,' RINGS  --  Too many 4-Membered Rings;',
     &                          ' Increase MAXRING')
                     call fatal
                  end if
                  iring4(1,nring4) = ia
                  iring4(2,nring4) = ib
                  iring4(3,nring4) = ic
                  iring4(4,nring4) = id
c
c     remove the ring if it is reducible into smaller rings
c
                  if (reduce) then
                     list(ia) = nring4
                     list(ib) = nring4
                     list(ic) = nring4
                     list(id) = nring4
                     do m = 1, nring3
                        list1 = list(iring3(1,m))
                        list2 = list(iring3(2,m))
                        list3 = list(iring3(3,m))
                        if (list1.eq.nring4 .and.
     &                      list2.eq.nring4 .and.
     &                      list3.eq.nring4) then
                           nring4 = nring4 - 1
                           list(ia) = 0
                           list(ib) = 0
                           list(ic) = 0
                           list(id) = 0
                           goto 40
                        end if
                     end do
                  end if
                  goto 40
               end if
            end do
   40       continue
         end if
      end do
c
c     search for and store all of the 5-membered rings
c
      do i = 1, n
         list(i) = 0
      end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         if (ia.lt.id .and. ie.lt.ib .and. min(ia,ie).lt.ic) then
            do j = 1, n12(ia)
               if (i12(j,ia) .eq. ie) then
                  nring5 = nring5 + 1
                  if (nring5 .gt. maxring) then
                     write (iout,50)
   50                format (/,' RINGS  --  Too many 5-Membered Rings;',    
     &                          ' Increase MAXRING')
                     call fatal
                  end if
                  iring5(1,nring5) = ia
                  iring5(2,nring5) = ib
                  iring5(3,nring5) = ic
                  iring5(4,nring5) = id
                  iring5(5,nring5) = ie
c
c     remove the ring if it is reducible into smaller rings
c
                  if (reduce) then
                     list(ia) = nring5
                     list(ib) = nring5
                     list(ic) = nring5
                     list(id) = nring5
                     list(ie) = nring5
                     do m = 1, nring3
                        list1 = list(iring3(1,m))
                        list2 = list(iring3(2,m))
                        list3 = list(iring3(3,m))
                        if (list1.eq.nring5 .and.
     &                      list2.eq.nring5 .and.
     &                      list3.eq.nring5) then
                           nring5 = nring5 - 1
                           list(ia) = 0
                           list(ib) = 0
                           list(ic) = 0
                           list(id) = 0
                           list(ie) = 0
                           goto 60
                        end if
                     end do
                  end if
                  goto 60
               end if
            end do
   60       continue
         end if
      end do
c
c     search for and store all of the 6-membered rings
c
      do i = 1, n
         list(i) = 0
      end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         imax = max(ia,ib,ic,id,ie)
         do j = 1, n12(ia)
            ig = i12(j,ia)
            if (ig .gt. imax) then
               do k = 1, n12(ie)
                  if (i12(k,ie) .eq. ig) then
                     nring6 = nring6 + 1
                     if (nring6 .gt. maxring) then
                        write (iout,70)
   70                   format (/,' RINGS  --  Too many 6-Membered',
     &                             ' Rings; Increase MAXRING')
                        call fatal
                     end if
                     iring6(1,nring6) = ia
                     iring6(2,nring6) = ib
                     iring6(3,nring6) = ic
                     iring6(4,nring6) = id
                     iring6(5,nring6) = ie
                     iring6(6,nring6) = ig
c
c     remove the ring if it is reducible into smaller rings
c
                     if (reduce) then
                        list(ia) = nring6
                        list(ib) = nring6
                        list(ic) = nring6
                        list(id) = nring6
                        list(ie) = nring6
                        list(ig) = nring6
                        do m = 1, nring3
                           list1 = list(iring3(1,m))
                           list2 = list(iring3(2,m))
                           list3 = list(iring3(3,m))
                           if (list1.eq.nring6 .and.
     &                         list2.eq.nring6 .and.
     &                         list3.eq.nring6) then
                              nring6 = nring6 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              goto 80
                           end if
                        end do
                        do m = 1, nring4
                           list1 = list(iring4(1,m))
                           list2 = list(iring4(2,m))
                           list3 = list(iring4(3,m))
                           list4 = list(iring4(4,m))
                           if (list1.eq.nring6 .and.
     &                         list2.eq.nring6 .and.
     &                         list3.eq.nring6 .and.
     &                         list4.eq.nring6) then
                              nring6 = nring6 - 1
                              list(ia) = 0
                              list(ib) = 0
                              list(ic) = 0
                              list(id) = 0
                              list(ie) = 0
                              list(ig) = 0
                              goto 80
                           end if
                        end do
                     end if
   80                continue
                  end if
               end do
            end if
         end do
      end do
c
c     search for and store all of the 7-membered rings
c
      do i = 1, n
         list(i) = 0
      end do
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         imax = max(ia,ib,ic,id,ie)
         do j = 1, n12(ia)
            ih = i12(j,ia)
            do k = 1, n12(ie)
               ig = i12(k,ie)
               if ((ig.gt.imax.and.ih.gt.ie) .or.
     &             (ih.gt.imax.and.ig.gt.ia)) then
                  do kk = 1, n12(ig)
                     if (i12(kk,ig) .eq. ih) then
                        nring7 = nring7 + 1
                        if (nring7 .gt. maxring) then
                           write (iout,90)
   90                      format (/,' RINGS  --  Too many 7-Membered',
     &                                ' Rings; Increase MAXRING')
                           call fatal
                        end if
                        iring7(1,nring7) = ia
                        iring7(2,nring7) = ib
                        iring7(3,nring7) = ic
                        iring7(4,nring7) = id
                        iring7(5,nring7) = ie
                        iring7(6,nring7) = ig
                        iring7(7,nring7) = ih
c
c     remove the ring if it is reducible into smaller rings
c
                        if (reduce) then
                           list(ia) = nring7
                           list(ib) = nring7
                           list(ic) = nring7
                           list(id) = nring7
                           list(ie) = nring7
                           list(ig) = nring7
                           list(ih) = nring7
                           do m = 1, nring3
                              list1 = list(iring3(1,m))
                              list2 = list(iring3(2,m))
                              list3 = list(iring3(3,m))
                              if (list1.eq.nring7 .and.
     &                            list2.eq.nring7 .and.
     &                            list3.eq.nring7) then
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
                           do m = 1, nring4
                              list1 = list(iring4(1,m))
                              list2 = list(iring4(2,m))
                              list3 = list(iring4(3,m))
                              list4 = list(iring4(4,m))
                              if (list1.eq.nring7 .and.
     &                            list2.eq.nring7 .and.
     &                            list3.eq.nring7 .and.
     &                            list4.eq.nring7) then
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
  100                   continue
                     end if
                  end do
               end if
            end do
         end do
      end do
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
