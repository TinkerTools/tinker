c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine trngl  --  triangle inequality via Dijkstra  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "trngl" smooths the upper distance bounds using the triangle
c     inequality implemented via Dijkstra's shortest path algorithm
c
c     literature reference:
c
c     T. F. Havel, I. D. Kuntz and G. M. Crippen, "The Theory and
c     Practice of Distance Geometry", Bulletin of Mathematical
c     Biology, 45, 665-720 (1983); see algorithm 2.1 on page 685
c
c
      subroutine trngl
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'disgeo.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,m,nalter
      real*8 bndmin,bndmax
      real*8 eps,big
      logical s(maxgeo)
c
c
c     find the maximal upper distance bound value
c
      eps = 1.0d-10
      big = 0.0d0
      do i = 1, n-1
         do j = i+1, n
            big = max(big,bnd(i,j))
         end do
      end do
c
c     apply the shortest path algorithm to each atom
c
      nalter = 0
      do i = 1, n
         do k = i+1, n
            s(k) = .true.
         end do
         do j = 1, n-i
            bndmin = big
            do k = i+1, n
               if (s(k) .and. bnd(i,k).le.bndmin) then
                  m = k
                  bndmin = bnd(i,k)
               end if
            end do
            s(m) = .false.
            do k = i+1, n
               if (s(k)) then
                  if (m .lt. k) then
                     bndmax = bnd(i,m) + bnd(m,k)
                  else
                     bndmax = bnd(i,m) + bnd(k,m)
                  end if
                  if (bnd(i,k) .gt. bndmax+eps) then
                     nalter = nalter + 1
                  end if
                  bnd(i,k) = min(bnd(i,k),bndmax)
               end if
            end do
         end do
         do j = i+1, n-1
            do k = j+1, n
               bndmax = bnd(i,k) + bnd(i,j)
               if (bnd(j,k) .gt. bndmax+eps) then
                  nalter = nalter + 1
               end if
               bnd(j,k) = min(bnd(j,k),bndmax)
            end do
         end do
      end do
c
c     write the results of the triangle inequality smoothing
c
      if (verbose) then
         write (iout,10)  nalter
   10    format (' TRNGL  --  Number of Upper Bound Changes :',i12)
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine trinv  --  inverse triangle via Dijkstra  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "trinv" smooths the lower distance bounds using the inverse
c     triangle inequality implemented via Dijkstra's shortest path
c     algorithm; this routine must be preceded by "trngl"
c
c     literature reference:
c
c     T. F. Havel, I. D. Kuntz and G. M. Crippen, "The Theory and
c     Practice of Distance Geometry", Bulletin of Mathematical
c     Biology, 45, 665-720 (1983); see algorithm 2.2 on page 688
c
c
      subroutine trinv
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'disgeo.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,m,nalter
      real*8 bndmin,bndmax,eps
      logical s(maxgeo)
c
c
c     first, apply the shortest path algorithm to each atom
c
      eps = 1.0d-10
      nalter = 0
      do i = 1, n
         do k = i+1, n
            s(k) = .true.
         end do
         do j = 1, n-i
            bndmax = 0.0d0
            do k = i+1, n
               if (s(k) .and. bnd(i,k).ge.bndmax) then
                  m = k
                  bndmax = bnd(i,k)
               end if
            end do
            s(m) = .false.
            do k = i+1, n
               if (s(k)) then
                  if (m .lt. k) then
                     bndmin = bnd(m,i) - bnd(m,k)
                  else
                     bndmin = bnd(m,i) - bnd(k,m)
                  end if
                  if (bnd(k,i) .lt. bndmin-eps) then
                     nalter = nalter + 1
                  end if
                  bnd(k,i) = max(bnd(k,i),bndmin)
               end if
            end do
         end do
         do j = 1, n-1
            do k = j+1, n
               if (j.ne.i .and. k.ne.i) then
                  if (i .lt. j) then
                     bndmin = max(bnd(k,i)-bnd(i,j),bnd(j,i)-bnd(i,k))
                  else if (i .lt. k) then
                     bndmin = max(bnd(k,i)-bnd(j,i),bnd(i,j)-bnd(i,k))
                  else
                     bndmin = max(bnd(i,k)-bnd(j,i),bnd(i,j)-bnd(k,i))
                  end if
                  if (bnd(k,j) .lt. bndmin-eps) then
                     nalter = nalter + 1
                  end if
                  bnd(k,j) = max(bnd(k,j),bndmin)
               end if
            end do
         end do
      end do
c
c     now, repeat the entire procedure in reverse atom order
c
      do i = n, 1, -1
         do k = 1, i-1
            s(k) = .true.
         end do
         do j = 1, i-1
            bndmax = 0.0d0
            do k = 1, i-1
               if (s(k) .and. bnd(k,i).ge.bndmax) then
                  m = k
                  bndmax = bnd(k,i)
               end if
            end do
            s(m) = .false.
            do k = 1, i-1
               if (s(k)) then
                  if (m .lt. k) then
                     bndmin = bnd(i,m) - bnd(m,k)
                  else
                     bndmin = bnd(i,m) - bnd(k,m)
                  end if
                  if (bnd(i,k) .lt. bndmin-eps) then
                     nalter = nalter + 1
                  end if
                  bnd(i,k) = max(bnd(i,k),bndmin)
               end if
            end do
         end do
         do j = 1, n-1
            do k = j+1, n
               if (j.ne.i .and. k.ne.i) then
                  if (i .lt. j) then
                     bndmin = max(bnd(k,i)-bnd(i,j),bnd(j,i)-bnd(i,k))
                  else if (i .lt. k) then
                     bndmin = max(bnd(k,i)-bnd(j,i),bnd(i,j)-bnd(i,k))
                  else
                     bndmin = max(bnd(i,k)-bnd(j,i),bnd(i,j)-bnd(k,i))
                  end if
                  if (bnd(k,j) .lt. bndmin-eps) then
                     nalter = nalter + 1
                  end if
                  bnd(k,j) = max(bnd(k,j),bndmin)
               end if
            end do
         end do
      end do
c
c     write the results of inverse triangle inequality smoothing
c
      if (verbose) then
         write (iout,10)  nalter
   10    format (' TRINV  --  Number of Lower Bound Changes :',i12)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine trngl2  --  triangle inequality via all triples  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "trngl2" smooths the upper distance bounds using the triangle
c     inequality via the Floyd-Warshall shortest path algorithm
c
c
      subroutine trngl2
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'disgeo.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k,nalter
      real*8 bndmax,eps
c
c
c     apply the triple loop structure over all atoms
c
      eps = 1.0d-10
      nalter = 0
      do i = 1, n
         do j = 1, n
            if (j .lt. i) then
               do k = 1, n
                  if (k .lt. j) then
                     bndmax = bnd(j,i) + bnd(k,i)
                     if (bnd(k,j) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(k,j) = min(bnd(k,j),bndmax)
                  else if (k.gt.j .and. k.lt.i) then
                     bndmax = bnd(j,i) + bnd(k,i)
                     if (bnd(j,k) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(j,k) = min(bnd(j,k),bndmax)
                  else if (k .gt. i) then
                     bndmax = bnd(j,i) + bnd(i,k)
                     if (bnd(j,k) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(j,k) = min(bnd(j,k),bndmax)
                  end if
               end do
            else if (i .lt. j) then
               do k = 1, n
                  if (k .lt. i) then
                     bndmax = bnd(i,j) + bnd(k,i)
                     if (bnd(k,j) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(k,j) = min(bnd(k,j),bndmax)
                  else if (k.gt.i .and. k.lt.j) then
                     bndmax = bnd(i,j) + bnd(i,k)
                     if (bnd(k,j) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(k,j) = min(bnd(k,j),bndmax)
                  else if (k .gt. j) then
                     bndmax = bnd(i,j) + bnd(i,k)
                     if (bnd(j,k) .gt. bndmax+eps) then
                        nalter = nalter + 1
                     end if
                     bnd(j,k) = min(bnd(j,k),bndmax)
                  end if
               end do
            end if
         end do
      end do
c
c     write the results of the triangle inequality smoothing
c
      if (verbose) then
         write (iout,10)  nalter
   10    format (' TRNGL2  --  Number of Upper Bound Changes :',i12)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine trinv2  --  inverse triangle via all triples  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "trinv2" iterates over all atom triples, raising the lower
c     bounds whenever these are less than possible by the inverse
c     triangle inequality; this routine must be preceded by "trngl2"
c
c
      subroutine trinv2
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'disgeo.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,k
      integer ialter,nalter
      integer niter,maxiter
      real*8 bndmin,eps
      logical done
c
c
c     initialize iteration counter and termination flag
c
      done = .false.
      maxiter = 100
      niter = 0
      eps = 1.0d-10
      nalter = 0
c
c     iterate in a triple loop structure over all atoms
c
      dowhile (.not. done)
         niter = niter + 1
         ialter = 0
         do i = 1, n-1
            do j = i+1, n
               do k = 1, i-1
                  bndmin = bnd(j,i) - bnd(k,j)
                  if (bnd(i,k) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(i,k) = max(bnd(i,k),bndmin)
                  bndmin = bnd(j,i) - bnd(k,i)
                  if (bnd(j,k) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(j,k) = max(bnd(j,k),bndmin)
               end do
               do k = i+1, j-1
                  bndmin = bnd(j,i) - bnd(k,j)
                  if (bnd(k,i) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(k,i) = max(bnd(k,i),bndmin)
                  bndmin = bnd(j,i) - bnd(i,k)
                  if (bnd(j,k) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(j,k) = max(bnd(j,k),bndmin)
               end do
               do k = j+1, n
                  bndmin = bnd(j,i) - bnd(j,k)
                  if (bnd(k,i) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(k,i) = max(bnd(k,i),bndmin)
                  bndmin = bnd(j,i) - bnd(i,k)
                  if (bnd(k,j) .lt. bndmin-eps) then
                     ialter = ialter + 1
                  end if
                  bnd(k,j) = max(bnd(k,j),bndmin)
               end do
            end do
         end do
         if (ialter.eq.0 .or. niter.eq.maxiter)  done = .true.
         nalter = nalter + ialter
      end do
c
c     write the results of inverse triangle inequality smoothing
c
      if (verbose) then
         write (iout,10)  nalter,niter
   10    format (' TRINV2  --  Number of Lower Bound Changes :',i12,
     &              3x,'at',i3,' Iterations')
      end if
      return
      end
