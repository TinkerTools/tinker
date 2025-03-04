c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2003  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine bitors  --  locate and store bitorsions  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "bitors" finds the total number of bitorsions as pairs
c     of adjacent torsional angles, and the numbers of the five
c     atoms defining each bitorsion
c
c
      subroutine bitors
      use angbnd
      use atoms
      use bitor
      use couple
      implicit none
      integer i,j,k
      integer ia,ib,ic,id,ie
c
c
c     initial count of the total number of bitorsions
c
      nbitor = 0
      do i = 1, nangle
         ib = iang(1,i)
         ic = iang(2,i)
         id = iang(3,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id) then
               do k = 1, n12(id)
                  ie = i12(k,id)
                  if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                     nbitor = nbitor + 1
                  end if
               end do
            end if
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ibitor))  deallocate (ibitor)
      allocate (ibitor(5,nbitor))
c
c     store the list of atoms involved in each bitorsion
c
      nbitor = 0
      do i = 1, nangle
         ib = iang(1,i)
         ic = iang(2,i)
         id = iang(3,i)
         do j = 1, n12(ib)
            ia = i12(j,ib)
            if (ia.ne.ic .and. ia.ne.id) then
               do k = 1, n12(id)
                  ie = i12(k,id)
                  if (ie.ne.ic .and. ie.ne.ib .and. ie.ne.ia) then
                     nbitor = nbitor + 1
                     ibitor(1,nbitor) = ia
                     ibitor(2,nbitor) = ib
                     ibitor(3,nbitor) = ic
                     ibitor(4,nbitor) = id
                     ibitor(5,nbitor) = ie
                  end if
               end do
            end if
         end do
      end do
      return
      end
