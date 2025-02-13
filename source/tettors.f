c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine tettors  --  locate and store tetratorsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "tettors" finds the total number of tetratorsions as tuples
c     of adjacent torsional angles, and the numbers of the seven
c     atoms defining each tetratorsion
c
c
      subroutine tettors
      use atoms
      use bitor
      use couple
      use tettor
      implicit none
      integer i,j,k
      integer ia,ib,ic,id
      integer ie,ig,ih
c
c
c     initial count of the total number of tetratorsions
c
      ntettor = 0
      do i = 1, nbitor
         ib = ibitor(1,i)
         ic = ibitor(2,i)
         id = ibitor(3,i)
         ie = ibitor(4,i)
         ig = ibitor(5,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id .and.
     &          ia.ne.ie .and. ia.ne.ig) then
               do k = 1, n12(ig)
                  ih = i12(k,ig)
                  if (ih.ne.ie .and. ih.ne.id .and. ih.ne.ic
     &                .and. ih.ne.ib .and. ih.ne.ia) then
                     ntettor = ntettor + 1
                  end if
               end do
            end if
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(itettor))  deallocate (itettor)
      allocate (itettor(7,ntettor))
c
c     store the list of atoms involved in each tetratorsion
c
      ntettor = 0
      do i = 1, nbitor
         ib = ibitor(1,i)
         ic = ibitor(2,i)
         id = ibitor(3,i)
         ie = ibitor(4,i)
         ig = ibitor(5,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id .and.
     &          ia.ne.ie .and. ia.ne.ig) then
               do k = 1, n12(ig)
                  ih = i12(k,ig)
                  if (ih.ne.ie .and. ih.ne.id .and. ih.ne.ic
     &                .and. ih.ne.ib .and. ih.ne.ia) then
                     ntettor = ntettor + 1
                     itettor(1,ntettor) = ia
                     itettor(2,ntettor) = ib
                     itettor(3,ntettor) = ic
                     itettor(4,ntettor) = id
                     itettor(5,ntettor) = ie
                     itettor(6,ntettor) = ig
                     itettor(7,ntettor) = ih
                  end if
               end do
            end if
         end do
      end do
      return
      end
