c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine mlgrid  --  maintain monotonic logical grid  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "mlgrid" sets up a monotonic logical grid for all atoms and
c     computes some statistics on physical vs. logical separations
c
c
      subroutine mlgrid
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'iounit.i'
      integer maxgrid
      parameter (maxgrid=20)
      integer i,j,k,nx,ny,nz
      integer iatom,ilist,next
      integer inside(maxgrid),outside(maxgrid)
      integer grid(maxgrid,maxgrid,maxgrid)
      integer list(maxatm),mesh(3,maxatm)
      real*8 elapsed,vdwcut2,dist2
c
c
c     get smallest possible size for the grid
c
      call setime
      nx = int(dble(n)**(1.0d0/3.0d0))
      ny = nx
      nz = nx
      if (nx*ny*nz .lt. n) then
         nx = nx + 1
         if (nx*ny*nz .lt. n) then
            ny = ny + 1
            if (nx*ny*nz .lt. n) then
               nz = nz + 1
            end if
         end if
      end if
c
c     initial assignment of atoms to the list
c
      do i = 1, n
         list(i) = i
      end do
      do i = n+1, nx*ny*nz
         list(i) = 0
      end do
c
c     use repeated sorts to correctly order the list
c
      next = 1
      call gsort (nx*ny*nz,list(next),z)
      do i = 1, nz
         call gsort (nx*ny,list(next),y)
         do j = 1, ny
            call gsort (nx,list(next),x)
            next = next + nx
         end do
      end do
c
c     transfer the list to the monotonic logical grid
c
      ilist = 0
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               ilist = ilist + 1
               iatom = list(ilist)
               grid(i,j,k) = iatom
               mesh(1,iatom) = i
               mesh(2,iatom) = j
               mesh(3,iatom) = k
            end do
         end do
      end do
      call getime (elapsed)
      write (iout,10)  elapsed
   10 format (/,' MLGRID  --  Time to Construct MLG :  ',f12.3,/)
c
c     print out the monotonic logical grid
c
      if (nx*ny*nz .le. 1000) then
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  iatom = grid(i,j,k)
                  if (iatom .ne. 0) then
                     write (iout,20)  i,j,k,iatom,x(iatom),
     &                                y(iatom),z(iatom)
   20                format (3i5,i10,3f10.4)
                  else
                     write (iout,30)  i,j,k,iatom
   30                format (3i5,i10)
                  end if
               end do
            end do
         end do
      end if
c
c     test grid separation vs. distance for all atom pairs
c
      do i = 1, maxgrid
         inside(i) = 0
         outside(i) = 0
      end do
      vdwcut2 = vdwcut * vdwcut
      do i = 1, n-1
         do j = i+1, n
            k = max(abs(mesh(1,i)-mesh(1,j)),abs(mesh(2,i)-mesh(2,j)),
     &                        abs(mesh(3,i)-mesh(3,j)))
            dist2 = ((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
            if (dist2 .le. vdwcut2) then
               inside(k) = inside(k) + 1
            else
               outside(k) = outside(k) + 1
            end if
         end do
      end do
      do i = 1, maxgrid
         write (iout,40)  i,inside(i),outside(i)
   40    format (' MLGRID Separation :',i6,5x,i10,
     &              ' Inside',i10,' Outside')
      end do
c
c     list all neighbors of a central grid vertex
c
      j = grid(nx/2+1,ny/2+1,nz/2+1)
      write (iout,50)  j,(mesh(k,j),k=1,3)
   50 format (/,' Central Atom : ',i6,5x,'Grid :',3i6)
      do i = 1, n
         dist2 = ((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
         if (dist2 .le. vdwcut2) then
            write (iout,60)  i,(mesh(k,i),k=1,3),sqrt(dist2)
   60       format (' Neighbor Atom :',i6,5x,'Grid :',3i6,f12.4)
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine gsort  --  heapsort of monotonic logical grid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gsort" uses the heapsort algorithm to sort part of a
c     monotonic logical grid based upon coordinate values
c
c
      subroutine gsort (n,grid,coord)
      implicit none
      include 'sizes.i'
      integer i,j,k,n,index,mid
      integer grid(*),grids,gridj,gridj1
      real*8 neg,pos,value
      real*8 coord(*),coordj,coordj1
c
c
c     perform the heapsort of the input list
c
      neg = -1000000.0d0
      pos = 1000000.0d0
      mid = n/2 + 1
      k = mid
      index = n
      dowhile (n .gt. 1)
         if (k .gt. 1) then
            k = k - 1
            grids = grid(k)
            if (grids .ne. 0) then
               value = coord(grids)
            else if (k .lt. mid) then
               value = neg
            else
               value = pos
            end if
         else
            grids = grid(index)
            if (grids .ne. 0) then
               value = coord(grids)
            else if (index .lt. mid) then
               value = neg
            else
               value = pos
            end if
            grid(index) = grid(1)
            index = index - 1
            if (index .eq. 1) then
               grid(1) = grids
               return
            end if
         end if
         i = k
         j = k + k
         dowhile (j .le. index)
            if (j .lt. index) then
               gridj = grid(j)
               if (gridj .ne. 0) then
                  coordj = coord(gridj)
               else if (gridj .lt. mid) then
                  coordj = neg
               else
                  coordj = pos
               end if
               gridj1 = grid(j+1)
               if (gridj1 .ne. 0) then
                  coordj1 = coord(gridj1)
               else if (gridj1 .lt. mid) then
                  coordj1 = neg
               else
                  coordj1 = pos
               end if
               if (coordj .lt. coordj1)  j = j + 1
            end if
            gridj = grid(j)
            if (gridj .ne. 0) then
               coordj = coord(gridj)
            else if (gridj .lt. mid) then
               coordj = neg
            else
               coordj = pos
            end if
            if (value .lt. coordj) then
               grid(i) = gridj
               i = j
               j = j + j
            else
               j = index + 1
            end if
         end do
         grid(i) = grids
      end do
      end
