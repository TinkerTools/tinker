c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine tritors  --  locate and store tritorsions  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "tritors" finds the total number of tritorsions as triples
c     of adjacent torsional angles, and the numbers of the six
c     atoms defining each tritorsion
c
c
      subroutine tritors
      use atoms
      use couple
      use tors
      use tritor
      implicit none
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,ig
c
c
c     initial count of the total number of tritorsions
c
      ntritor = 0
      do i = 1, ntors
         ib = itors(1,i)
         ic = itors(2,i)
         id = itors(3,i)
         ie = itors(4,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id .and. ia.ne.ie) then
               do k = 1, n12(ie)
                  ig = i12(k,ie)
                  if (ig.ne.id .and. ig.ne.ic .and.
     &                ig.ne.ib .and. ig.ne.ia) then
                     ntritor = ntritor + 1
                  end if
               end do
            end if
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(itritor))  deallocate (itritor)
      allocate (itritor(6,ntritor))
c
c     store the list of atoms involved in each tritorsion
c
      ntritor = 0
      do i = 1, ntors
         ib = itors(1,i)
         ic = itors(2,i)
         id = itors(3,i)
         ie = itors(4,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id .and. ia.ne.ie) then
               do k = 1, n12(ie)
                  ig = i12(k,ie)
                  if (ig.ne.id .and. ig.ne.ic .and.
     &                ig.ne.ib .and. ig.ne.ia) then
                     ntritor = ntritor + 1
                     itritor(1,ntritor) = ia
                     itritor(2,ntritor) = ib
                     itritor(3,ntritor) = ic
                     itritor(4,ntritor) = id
                     itritor(5,ntritor) = ie
                     itritor(6,ntritor) = ig
                  end if
               end do
            end if
         end do
      end do
      return
      end
