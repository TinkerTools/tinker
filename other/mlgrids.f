c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program mlgrids  --  monotonic logical grid demonstration  ##
c     ##                                                             ##
c     #################################################################
c
c
      program mlgrids
      implicit none
      include 'sizes.i'
      integer n,n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer i,j,k,count,ii,ij,ik
      real*8 x,y,z,trans
      real*8 boxlength,boxcut
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm),
     &                boxlength,boxcut,trans(3)
c
c
c     get coordinates file
c
      call initial
      call getxyz
      n_orig = n
c
c     define boundary conditions --> expand
c
      call getbound
      call setime
      call expand
      call getime
      n_box = n
      write (*,*) ' n_orig, n_box : ', n_orig, n_box
c
c     setup monotonic logical grid: here is n = n_box
c
      call setime
      call setmlg
      call getime
      write (*,*) ' n_box, n : ', n_box,n
c
c     maintain monotonic logical grid
c
      call setime
      do ii = 1, 50
         call swap
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine setmlg  --  setup a monotonic logical grid  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine setmlg
      implicit none
      include 'sizes.i'
      real*8 x,y,z
      real*8 matx(8),maty(8),matz(8)
      integer n, n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer i,j,k,iatom,array(maxatm)
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(3,maxatm)
c
c
c     constants
c
      iatom = 0
      do i = -1, 1, 2
         do j = -1, 1, 2
            do k = -1, 1, 2
               iatom = iatom + 1
               matx(iatom) = float(i) * 100000.0
               maty(iatom) = float(j) * 100000.0
               matz(iatom) = float(k) * 100000.0
            end do
         end do
      end do
c
c     get smallest possible size of the 'cube-grid'
c
      ngrid = 0
      dowhile (ngrid*ngrid*ngrid .lt. n)
         ngrid = ngrid + 1
      end do
      ngrid2 = ngrid * ngrid
      n = ngrid2 * ngrid
c
c     fill corners of cube with "holes"
c
      do i = n_box+1, n
         iatom = mod(i,8) + 1
         x(i) = matx(iatom)
         y(i) = maty(iatom)
         z(i) = matz(iatom)
      end do
c
c     initial assignment of atoms to the grid
c
      do i = 1, n
         array(i) = i
      end do
c
c     remark: eventually it is more efficient to set up array(i) once
c             and keep it for the next iteration. We would save this
c             simple do-loop and the heapsort is probably faster because
c             the old "array" might be still a good guess.
c             In the same way it might be better to keep the "holes" onc
c             they are set up. This would work, as long as the number of
c             'expanded' atoms does not vary. If it decreases, one has t
c             introduce new "holes".
c
c     use sort to achieve initial monotonic logical grid
c
      call heapsort (n,array(1),z)
      do i = 0, ngrid-1
         call heapsort (ngrid2,array(i*ngrid2+1),y)
         do j = 0, ngrid-1
            call heapsort (ngrid,array(i*ngrid2+j*ngrid+1),x)
         end do
      end do
      iatom = 0
      do k = 1, ngrid
         do j = 1, ngrid
            do i = 1, ngrid
               iatom = iatom + 1
               grid(i,j,k) = array(iatom)
               gridpoint(1,array(iatom)) = i
               gridpoint(2,array(iatom)) = j
               gridpoint(3,array(iatom)) = k
            end do
         end do
      end do
c
c     print out monotonic logical grid
c
c      do i = 1, ngrid
c         do j = 1, ngrid
c            do k = 1, ngrid
c               iatom = grid(i,j,k)
c               write (*,20) i,j,k,iatom,x(iatom),y(iatom),z(iatom)
c            end do
c         end do
c      end do
c   20 format (3i5,i10,3f12.4)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine heapsort  --  get ascending order via Heapsort  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine heapsort (size,grid,coord)
      implicit none
      include 'sizes.i'
      integer i,j,l,size,grid(size)
      integer temp,gridj,gridj1
      real*8 coord(maxatm)
      real*8 value,coordj,coordj1
c
c
      l = size/2 + 1
   10 continue
      if (l .gt. 1) then
         l = l - 1
         temp = grid(l)
         value = coord(temp)
      else
         temp = grid(size)
         value = coord(temp)
         grid(size) = grid(1)
         size = size - 1
         if (size .eq. 1) then
            grid(1) = temp
            return
         end if
      end if
      i = l
      j = l + l
      dowhile (j .le. size)
         if (j .lt. size) then
            gridj = grid(j)
            coordj = coord(gridj)
            gridj1 = grid(j+1)
            coordj1 = coord(gridj1)
            if (coordj .lt. coordj1)  j = j + 1
         end if
         gridj = grid(j)
         coordj = coord(gridj)
         if (value .lt. coordj) then
            grid(i) = gridj
            i = j
            j = j + j
         else
            j = size + 1
         end if
      end do
      grid(i) = temp
      goto 10
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine swap  --  check if mlg is violated -      ##
c     ##                       do one iteration of correction  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine swap
      implicit none
      include 'sizes.i'
      real*8 x,y,z
      integer n,n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer i,j,k,temp,ptemp,nswaps
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm)
c
c
      nswaps = 0
      do k = 1, ngrid-1
         do j = 1, ngrid-1
            do i = 1, ngrid-1
               if (x(grid(i,j,k)) .gt. x(grid(i+1,j,k))) then
                  temp = grid(i,j,k)
                  ptemp = gridpoint(temp)
                  grid(i,j,k) = grid(i+1,j,k)
                  gridpoint(grid(i,j,k)) = gridpoint(grid(i+1,j,k))
                  grid(i+1,j,k) = temp
                  gridpoint(grid(i+1,j,k)) = ptemp
                  nswaps = nswaps + 1
               end if
               if (y(grid(i,j,k)) .gt. y(grid(i,j+1,k))) then
                  temp = grid(i,j,k)
                  ptemp = gridpoint(temp)
                  grid(i,j,k) = grid(i,j+1,k)
                  gridpoint(grid(i,j,k)) = gridpoint(grid(i,j+1,k))
                  grid(i,j+1,k) = temp
                  gridpoint(grid(i,j+1,k)) = ptemp
                  nswaps = nswaps + 1
               end if
               if (z(grid(i,j,k)) .gt. z(grid(i,j,k+1))) then
                  temp = grid(i,j,k)
                  ptemp = gridpoint(temp)
                  grid(i,j,k) = grid(i,j,k+1)
                  gridpoint(grid(i,j,k)) = gridpoint(grid(i,j,k+1))
                  grid(i,j,k+1) = temp
                  gridpoint(grid(i,j,k+1)) = ptemp
                  nswaps = nswaps + 1
               end if
            end do
         end do
      end do
c     write (*,*) ' Number of swaps : ', nswaps
      if (nswaps.eq.0) call lib$show_timer
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine pointer  --  atom-number ---> mlgrid  ##
c     ##                                                   ##
c     #######################################################
c
c
      subroutine pointer (number,i,j,k)
      implicit none
      include 'sizes.i'
      integer n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer number,i,j,k, temp
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm)
c
c
c     get indices i,j,k from atom number
c
      temp = gridpoint(number)-1 - mod(gridpoint(number)-1,ngrid2)
      i = temp / ngrid2
      temp = gridpoint(number) - temp -1
      k = mod(temp,ngrid)
      j = (temp - k)/ ngrid
      k = k + 1
      j = j + 1
      i = i + 1
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine getbound  --  setup boundary conditions  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine getbound
      implicit none
      include 'sizes.i'
      integer n,itype
      integer n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer i,j,k, count
      real*8 x,y,z
      real*8 boxlength,boxcore,boxcut,trans
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm),
     &                boxlength,boxcore,trans(3)
c
c
      write (*,10)
   10 format ('$Enter boxlength [f12.5] : ')
      read (*,20) boxlength
   20 format (f12.5)
   30 continue
      write (*,40)
   40 format ('$Enter box_cutoff [f12.5] : ')
      read (*,20) boxcut
c
c     in future one can look-up vdwcut or dipcut
c
      if (boxcut .gt. boxlength) then
         write (*,50)
   50    format (' box_cutoff too large! Try again. ')
         goto 30
      end if
      boxcore = boxlength/2.0d0 - boxcut
c
c     load trans
c
      trans(1) =  boxlength
      trans(2) =  0.0d0
      trans(3) = -boxlength
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine expand -- surround orig. atoms by box-expansion ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine expand
      implicit none
      include 'sizes.i'
      integer n,itype
      integer n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer ii,i,j,k,ix,iy,iz,ixs,iys,izs
      real*8 x,y,z
      real*8 boxlength,boxcore,trans
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm),
     &                boxlength,boxcore,trans(3)
c
c
      do ii = 1, n_orig
         ix = 2
         ixs = 2
         if (x(ii) .gt. boxcore) ixs = 3
         if (x(ii) .lt. -boxcore) ix = 1
         iy = 2
         iys = 2
         if (y(ii) .gt. boxcore) iys = 3
         if (y(ii) .lt. -boxcore) iy = 1
         iz = 2
         izs = 2
         if (z(ii) .gt. boxcore) izs = 3
         if (z(ii) .lt. -boxcore) iz = 1
         do i = ix, ixs
            do j = iy, iys
               do k = iz, izs
                  if (i.ne.2 .or. j.ne.2 .or. k.ne.2) then
                     n = n + 1
                     x(n) = x(ii) + trans(i)
                     y(n) = y(ii) + trans(j)
                     z(n) = z(ii) + trans(k)
                     itype(n) = itype(ii)
                  end if
               end do
            end do
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine getmol  --  set up molecule list  ##
c     ##                                               ##
c     ###################################################
c
c
      subroutine getmol
      implicit none
      include 'sizes.i'
      integer n_orig,nentry
      integer i,j,n_mol,mol
      integer n12,i12
      logical entry
      common /mlgrid/ n_orig
      common /attach/ n12(maxatm),i12(4,maxatm)
      common /molec / n_mol,mol(2,maxatm)
c
c
c     determine which atoms belong to separate molecules and setup
c     the array mol: mol(1,i) = index of cetral atom in mol #i
c                    mol(2,i) = number of atoms in this molecule
c
      n_mol = 1
      mol(1,n_mol) = 1
      nentry = 0
      do i = 2, n_orig
         entry = .true.
         nentry = nentry + 1
         if (n12(i) .gt. 0 ) then
            do j = 1, n12(i)
               if (i12(j,i) .lt. i) entry = .false.
            end do
         end if
         if (entry) then
            mol(2,n_mol) = nentry
            n_mol = n_mol + 1
            mol(1,n_mol) = i
            nentry = 0
         end if
      end do
      mol(2,n_mol) = nentry + 1
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine correct -- correct for mol. which left the box  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine correct
      implicit none
      include 'sizes.i'
      integer n,itype,n_orig,n_box,ngrid,ngrid2,grid,gridpoint
      integer i,j,n_mol,mol
      real*8 x,y,z,boxlength,halfbox
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /mlgrid/ n_orig,n_box,ngrid,ngrid2,
     &                grid(maxgrd,maxgrd,maxgrd),gridpoint(maxatm),
     &                boxlength
      common /molec / n_mol,mol(2,maxatm)
c
c
c     if an atoms moves out of the box it has to be put in at the
c     opposite edge. Since most of the atoms are not independent, we
c     have to move an entire molecule. Therefore a molecule list must
c     have been set up before: mol(1..2,1..maxmol)=number of central ato
c
      halfbox = boxlength / 2.0d0
      do i =  1, n_mol
         if (x(mol(1,i)) .gt. halfbox) then
            do j = 0, mol(2,i)-1
               x(mol(1,i)+j) = x(mol(1,i)+j) - boxlength
            end do
         elseif (x(mol(1,i)) .lt. -halfbox) then
            do j = 0, mol(2,i)-1
               x(mol(1,i)+j) = x(mol(1,i)+j) + boxlength
            end do
         elseif (y(mol(1,i)) .gt. halfbox) then
            do j = 0, mol(2,i)-1
               y(mol(1,i)+j) = y(mol(1,i)+j) - boxlength
            end do
         elseif (y(mol(1,i)) .lt. -halfbox) then
            do j = 0, mol(2,i)-1
               y(mol(1,i)+j) = y(mol(1,i)+j) + boxlength
            end do
         elseif (z(mol(1,i)) .gt. halfbox) then
            do j = 0, mol(2,i)-1
               z(mol(1,i)+j) = z(mol(1,i)+j) - boxlength
            end do
         elseif (z(mol(1,i)) .lt. -halfbox) then
            do j = 0, mol(2,i)-1
               z(mol(1,i)+j) = z(mol(1,i)+j) + boxlength
            end do
         end if
      end do
      return
      end
