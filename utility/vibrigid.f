c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2002  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program vibrigid  --  rigid body vibrational analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "vibrigid" computes the eigenvalues and eigenvectors of the
c     Hessian matrix over rigid body degrees of freedom
c
c
      program vibrigid
      implicit none
      use group
      use iounit
      use rigid
      use sizes
      integer maxrgd
      parameter (maxrgd=6*maxgrp)
      integer i,j,ihess,nvar
      real*8 a(maxrgd+1)
      real*8 b(maxrgd+1)
      real*8 p(maxrgd+1)
      real*8 w(maxrgd+1)
      real*8 ta(maxrgd+1)
      real*8 tb(maxrgd+1)
      real*8 ty(maxrgd+1)
      real*8 eigen(maxrgd)
      real*8 matrix((maxrgd+1)*maxrgd/2)
      real*8 hrigid(maxrgd,maxrgd)
      real*8 vects(maxrgd,maxrgd)
c
c
c     set up the molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     set up the use of rigid body coordinate system
c
      use_rigid = .true.
      call orient
c
c     compute rigid body Hessian matrix elements
c
      call hessrgd (hrigid)
c
c     place Hessian elements into triangular form
c
      nvar = 6 * ngrp
      ihess = 0
      do i = 1, nvar
         do j = i, nvar
            ihess = ihess + 1
            matrix(ihess) = hrigid(j,i)
         end do
      end do
c
c     perform diagonalization to get Hessian eigenvalues
c
      call diagq (nvar,maxrgd,nvar,matrix,eigen,vects,
     &                    a,b,p,w,ta,tb,ty)
      write (iout,50)
   50 format (/,' Eigenvalues of the Hessian Matrix :',/)
      write (iout,60)  (eigen(i),i=1,nvar)
   60 format (6d13.4)
c
c     perform any final tasks before program exit
c
      call final
      end
