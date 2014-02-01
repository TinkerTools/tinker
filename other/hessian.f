c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine hessian  --  atom-by-atom Hessian elements  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "hessian" calls subroutines to calculate the Hessian elements
c     for each atom in turn with respect of Cartesian coordinates
c
c
      subroutine hessian (hess,hinit,hstop,hindex,hdiag)
      implicit none
      include 'sizes.for'
      integer*2 hindex(maxhess)
      integer i,j,k,n,ii,input,iout,itype,nvdw,nv14,n12,i12
      integer nhess,hinit(3,maxatm),hstop(3,maxatm)
      real hess(maxhess)
      real*8 x,y,z,vrad,veps,seps,reduc,hesscut
      real*8 percent,filled,redi
      real*8 xr(maxatm),yr(maxatm),zr(maxatm)
      real*8 hessx,hessy,hessz,hdiag(3,maxatm)
      common /atoms / n,x(maxatm),y(maxatm),z(maxatm),itype(maxatm)
      common /attach/ n12(maxatm),i12(4,maxatm)
      common /hescut/ hesscut
      common /hessn / hessx(3,maxatm),hessy(3,maxatm),hessz(3,maxatm)
      common /iounit/ input,iout
      common /vdw   / nvdw,nv14,vrad(maxtyp),veps(maxtyp),seps(maxtyp),
     &                reduc(maxtyp)
c
c
c     zero out total number of indexed hessian elements
c
      nhess = 0
      do i = 1, n
         do j = 1, 3
            hdiag(j,i) = 0.0d0
         end do
      end do
c
c     calculate "reduced" atomic coordinates for hydrogens
c
      do i = 1, n
         redi = reduc(itype(i))
         if (redi .eq. 0.0d0) then
            xr(i) = x(i)
            yr(i) = y(i)
            zr(i) = z(i)
         else
            ii = i12(1,i)
            xr(i) = redi*(x(i)-x(ii)) + x(ii)
            yr(i) = redi*(y(i)-y(ii)) + y(ii)
            zr(i) = redi*(z(i)-z(ii)) + z(ii)
         end if
      end do
c
c     zero out the hessian elements for the current atom
c
      do i = 1, n
         do k = 1, n
            do j = 1, 3
               hessx(j,k) = 0.0d0
               hessy(j,k) = 0.0d0
               hessz(j,k) = 0.0d0
            end do
         end do
c
c     call the individual second derivative routines
c
         call ebond2 (i)
         call eangle2 (i)
         call etors2 (i)
         call evdw2 (i,xr,yr,zr)
         call echarge2 (i)
         call echgdpl2 (i)
         call edipole2 (i)
         call extra2 (i)
c
c     set the diagonal hessian matrix elements
c
         hdiag(1,i) = hdiag(1,i) + hessx(1,i)
         hdiag(2,i) = hdiag(2,i) + hessy(2,i)
         hdiag(3,i) = hdiag(3,i) + hessz(3,i)
c
c     copy selected off-diagonal hessian elements for current
c     atom into an indexed master list of hessian elements;
c     if below cutoff, add into diagonal to maintain row sum
c
         hinit(1,i) = nhess + 1
         do j = 2, 3
            if (abs(hessx(j,i)) .gt. hesscut) then
               nhess = nhess + 1
               hindex(nhess) = 3*i + j - 3
               hess(nhess) = hessx(j,i)
            else
               hdiag(1,i) = hdiag(1,i) + abs(hessx(j,i))
               hdiag(j,i) = hdiag(j,i) + abs(hessx(j,i))
            end if
         end do
         do k = i+1, n
            do j = 1, 3
               if (abs(hessx(j,k)) .gt. hesscut) then
                  nhess = nhess + 1
                  hindex(nhess) = 3*k + j - 3
                  hess(nhess) = hessx(j,k)
               else
                  hdiag(1,i) = hdiag(1,i) + abs(hessx(j,k))
                  hdiag(j,k) = hdiag(j,k) + abs(hessx(j,k))
               end if
            end do
         end do
         hstop(1,i) = nhess
c
         hinit(2,i) = nhess + 1
         if (abs(hessy(3,i)) .gt. hesscut) then
            nhess = nhess + 1
            hindex(nhess) = 3*i
            hess(nhess) = hessy(3,i)
         else
            hdiag(2,i) = hdiag(2,i) + abs(hessy(3,i))
            hdiag(3,i) = hdiag(3,i) + abs(hessy(3,i))
         end if
         do k = i+1, n
            do j = 1, 3
               if (abs(hessy(j,k)) .gt. hesscut) then
                  nhess = nhess + 1
                  hindex(nhess) = 3*k + j - 3
                  hess(nhess) = hessy(j,k)
               else
                  hdiag(2,i) = hdiag(2,i) + abs(hessy(j,k))
                  hdiag(j,k) = hdiag(j,k) + abs(hessy(j,k))
               end if
            end do
         end do
         hstop(2,i) = nhess
c
         hinit(3,i) = nhess + 1
         do k = i+1, n
            do j = 1, 3
               if (abs(hessz(j,k)) .gt. hesscut) then
                  nhess = nhess + 1
                  hindex(nhess) = 3*k + j - 3
                  hess(nhess) = hessz(j,k)
               else
                  hdiag(3,i) = hdiag(3,i) + abs(hessz(j,k))
                  hdiag(j,k) = hdiag(j,k) + abs(hessz(j,k))
               end if
            end do
         end do
         hstop(3,i) = nhess
c
c     check for storage of too many hessian elements
c
         if (nhess .gt. maxhess) then
            write (iout,10)  nhess,maxhess,hesscut
   10       format (' Current required Hessian storage : ',i12,/
     &              ' Maximum allowed Hessian storage :  ',i12,/
     &              ' Minimum significant Hessian value :',f12.6,/
     &              ' -->  Increase MAXHESS and/or HESSCUT')
         end if
      end do
c
c     print message telling how much storage was finally used
c
      percent = 100.0d0 * dfloat(nhess)/dfloat(3*n*(3*n-1)/2)
      filled = 100.0d0 * dfloat(nhess)/dfloat(maxhess)
      write (iout,20)  nhess,percent,filled
   20 format (' HESSIAN :',i11,' Elements',
     &          f9.2,' % Off-diag H',f8.2,' % Storage')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine column  --  access hessian elements by column  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "column" takes the off-diagonal Hessian elements stored
c     as sparse rows and sets up indices to allow column access
c
c
      subroutine column (nvar,h_init,h_stop,h_index,
     &                   c_init,c_stop,c_index,c_value)
      implicit none
      include 'sizes.for'
      integer maxvar
      parameter (maxvar=3*maxatm)
      integer*2 h_index(maxhess),c_index(maxhess)
      integer i,j,k,m,nvar
      integer h_init(maxvar),h_stop(maxvar)
      integer c_init(maxvar),c_stop(maxvar)
      integer c_value(maxhess)
c
c
c     zero out beginning and end marker for each column
c
      do i = 1, nvar
         c_init(i) = 0
         c_stop(i) = 0
      end do
c
c     count the number of elements in each column
c
      do i = 1, nvar
         do j = h_init(i), h_stop(i)
            k = h_index(j)
            c_stop(k) = c_stop(k) + 1
         end do
      end do
c
c     set each beginning marker just beyond
c     the last element for its column
c
      c_init(1) = c_stop(1) + 1
      do i = 2, nvar
         c_init(i) = c_init(i-1) + c_stop(i)
      end do
c
c     set column index by scanning rows in reverse order
c
      do i = nvar, 1, -1
         do j = h_init(i), h_stop(i)
            k = h_index(j)
            m = c_init(k) - 1
            c_init(k) = m
            c_index(m) = i
            c_value(m) = j
         end do
      end do
c
c     convert 'c_stop' from number of elements
c     in column to end marker for the column
c
      do i = 1, nvar
         c_stop(i) = c_init(i) + c_stop(i) - 1
      end do
      return
      end
