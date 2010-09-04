c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program vibrot  --  vibrational analysis over torsions  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "vibrot" computes the eigenvalues and eigenvectors of the
c     torsional Hessian matrix
c
c     literature reference:
c
c     M. Levitt, C. Sander and P. S. Stern, "Protein Normal-mode
c     Dynamics: Trypsin Inhibitor, Crambin, Ribonuclease and Lysozyme",
c     Journal of Molecular Biology, 181, 423-447 (1985)
c
c
      program vibrot
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'omega.i'
      integer i,j,ihess
      real*8 hrot(maxrot,maxrot)
      real*8, allocatable :: a(:)
      real*8, allocatable :: b(:)
      real*8, allocatable :: p(:)
      real*8, allocatable :: w(:)
      real*8, allocatable :: ta(:)
      real*8, allocatable :: tb(:)
      real*8, allocatable :: ty(:)
      real*8, allocatable :: eigen(:)
      real*8, allocatable :: matrix(:)
      real*8, allocatable :: vects(:,:)
c
c
c     set up the mechanics calculation and rotatable bonds
c
      call initial
      call getint
      call mechanic
      call initrot
c
c     perform dynamic allocation of some local arrays
c
      allocate (a(nomega))
      allocate (b(nomega))
      allocate (p(nomega))
      allocate (w(nomega))
      allocate (ta(nomega))
      allocate (tb(nomega))
      allocate (ty(nomega))
      allocate (eigen(nomega))
      allocate (matrix(nomega*(nomega+1)/2))
      allocate (vects(nomega,nomega))
c
c     calculate the full torsional Hessian matrix
c
      call hessrot ('FULL',hrot)
c
c     write out the torsional Hessian diagonal
c
      write (iout,10)
   10 format (/,' Diagonal of the Torsional Hessian :',/)
      write (iout,20)  (i,hrot(i,i),i=1,nomega)
   20 format (4(i8,f11.3))
c
c     write out the torsional Hessian elements
c
      if (nomega .le. 30) then
         write (iout,30)
   30    format (/,' Torsional Hessian Matrix Elements :')
         do i = 1, nomega
            write (iout,40)
   40       format ()
            write (iout,50)  (hrot(j,i),j=1,nomega)
   50       format (6f13.4)
         end do
      end if
c
c     place Hessian elements into triangular form
c
      ihess = 0
      do i = 1, nomega
         do j = i, nomega
            ihess = ihess + 1
            matrix(ihess) = hrot(i,j)
         end do
      end do
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nomega,nomega,nomega,matrix,eigen,vects,
     &                     a,b,p,w,ta,tb,ty)
      write (iout,60)
   60 format (/,' Eigenvalues of the Hessian Matrix :',/)
      write (iout,70)  (i,eigen(i),i=1,nomega)
   70 format (4(i8,f11.3))
c
c     perform deallocation of some local arrays
c
      deallocate (a)
      deallocate (b)
      deallocate (p)
      deallocate (w)
      deallocate (ta)
      deallocate (tb)
      deallocate (ty)
      deallocate (eigen)
      deallocate (matrix)
      deallocate (vects)
c
c     perform any final tasks before program exit
c
      call final
      end
