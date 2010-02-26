c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine molecule  --  assign atoms to molecules  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "molecule" counts the molecules, assigns each atom to
c     its molecule and computes the mass of each molecule
c
c
      subroutine molecule
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'molcul.i'
      integer i,j,k
      integer mi,mj,mk
      integer iattach
      integer list(maxatm)
c
c
c     zero number of molecules and molecule membership list
c
      nmol = 0
      do i = 1, n
         molcule(i) = 0
      end do
c
c     assign each atom to its respective molecule
c
      do i = 1, n
         if (molcule(i) .eq. 0) then
            nmol = nmol + 1
            molcule(i) = nmol
         end if
         mi = molcule(i)
         do iattach = 1, n12(i)
            j = i12(iattach,i)
            mj = molcule(j)
            if (mj .eq. 0) then
               molcule(j) = mi
            else if (mi .lt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mj) then
                     molcule(k) = mi
                  else if (mk .gt. mj) then
                     molcule(k) = mk - 1
                  end if
               end do
            else if (mi .gt. mj) then
               nmol = nmol - 1
               do k = 1, n
                  mk = molcule(k)
                  if (mk .eq. mi) then
                     molcule(k) = mj
                  else if (mk .gt. mi) then
                     molcule(k) = mk - 1
                  end if
               end do
               mi = mj
            end if
         end do
      end do
c
c     pack atoms of each molecule into a contiguous indexed list
c
      do i = 1, n
         list(i) = molcule(i)
      end do
      call sort3 (n,list,kmol)
c
c     find the first and last atom in each molecule
c
      k = 1
      imol(1,1) = 1
      do i = 2, n
         j = list(i)
         if (j .ne. k) then
            imol(2,k) = i - 1
            k = j
            imol(1,k) = i
         end if
      end do
      imol(2,nmol) = n
c
c     sort the list of atoms in each molecule by atom number
c
      do i = 1, nmol
         k = imol(2,i) - imol(1,i) + 1
         call sort (k,kmol(imol(1,i)))
      end do
c
c     if all atomic masses are zero, set them all to unity
c
      do i = 1, n
         if (mass(i) .ne. 0.0d0)  goto 10
      end do
      do i = 1, n
         mass(i) = 1.0d0
      end do
   10 continue
c
c     compute the mass of each molecule and the total mass
c
      totmass = 0.0d0
      do i = 1, nmol
         molmass(i) = 0.0d0
         do k = imol(1,i), imol(2,i)
            molmass(i) = molmass(i) + mass(kmol(k))
         end do
         totmass = totmass + molmass(i)
      end do
      return
      end
