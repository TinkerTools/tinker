c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bonds  --  locate and store covalent bonds  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bonds" finds the total number of covalent bonds and
c     stores the atom numbers of the atoms defining each bond
c
c
      subroutine bonds
      use atmlst
      use atoms
      use bndstr
      use couple
      implicit none
      integer i,j,k,m
c
c
c     initial count of the total number of bonds
c
      nbond = 0
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            if (i .lt. k) then
               nbond = nbond + 1
            end if
         end do
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ibnd))  deallocate (ibnd)
      if (allocated(bndlist))  deallocate (bndlist)
      allocate (ibnd(2,nbond))
      allocate (bndlist(maxval,n))
c
c     store the list of atoms involved in each bond
c
      nbond = 0
      do i = 1, n
         do j = 1, n12(i)
            k = i12(j,i)
            if (i .lt. k) then
               nbond = nbond + 1
               ibnd(1,nbond) = i
               ibnd(2,nbond) = k
c
c     store the numbers of the bonds involving each atom
c
               bndlist(j,i) = nbond
               do m = 1, n12(k)
                  if (i .eq. i12(m,k)) then
                     bndlist(m,k) = nbond
                     goto 10
                  end if
               end do
   10          continue
            end if
         end do
      end do
      return
      end
