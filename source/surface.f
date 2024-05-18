c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Patrice Koehl & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine surface  --  alpha shapes accessible surface  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "surface" computes the weighted solvent accessible surface
c     area each atom via the inclusion-exclusion method of Herbert
c     Edelsbrunner based on alpha shapes
c
c     note for small or symmetric structures where alpha shapes
c     may fail, switch to the Richmond method
c
c     developed to facilitate calling UnionBall from Tinker by
c     Jay W. Ponder, Washington University, October 2023
c
c     literature references:
c
c     P. Mach and P. Koehl, "Geometric Measures of Large Biomolecules:
c     Surface, Volume, and Pockets", Journal of Computational Chemistry,
c     32, 3023-3038 (2011)
c
c     P. Koehl, A. Akopyan and H. Edelsbrunner, "Computing the Volume,
c     Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
c     Derivatives", Journal of Chemical Information and Modeling, 63,
c     973-985 (2023)
c
c     variables and parameters:
c
c     nsphere    number of spheres/balls in the system
c     coords     coordinates of the center of each sphere
c     radii      radius value for each sphere
c     weight     weight value for each sphere
c     probe      radius value of the probe sphere
c     surf       weighted surface area of union of spheres
c     usurf      unweighted surface area of union of spheres
c     asurf      weighted area contribution of each sphere
c
c
      subroutine surface (rad,weight,probe,surf,asurf)
      use atoms
      implicit none
      integer i,nsphere
      integer nsize,nfudge
      integer nredundant
      integer, allocatable :: redlist(:)
      real*8 surf,usurf,eps
      real*8 probe,alpha
      real*8 rad(*)
      real*8 weight(*)
      real*8 asurf(*)
      real*8, allocatable :: radii(:)
      real*8, allocatable :: asurfx(:)
      real*8, allocatable :: coords(:,:)
      logical dowiggle
      character*6 symmtyp
c
c
c     check coordinates for linearity, planarity and symmetry
c
      symmtyp = 'NONE'
      call chksymm (symmtyp)
      dowiggle = .false.
      if (n.gt.2 .and. symmtyp.eq.'LINEAR')  dowiggle = .true.
      if (n.gt.3 .and. symmtyp.eq.'PLANAR')  dowiggle = .true.
c
c     use Richmond method for small symmetric structures
c
      if (dowiggle) then
         call richmond (n,x,y,z,rad,weight,probe,surf,asurf)
         return
      end if
c
c     perform dynamic allocation of some local arrays
c
      nfudge = 10
      nsize = n + nfudge
      allocate (radii(nsize))
      allocate (asurfx(nsize))
      allocate (coords(3,nsize))
      allocate (redlist(nsize))
c
c     set the coordinates and sphere radii plus probe`
c
      nsphere = n
      do i = 1, n
         coords(1,i) = x(i)
         coords(2,i) = y(i)
         coords(3,i) = z(i)
         radii(i) = 0.0d0
         if (rad(i) .ne. 0.0d0)  radii(i) = rad(i) + probe
      end do
c
c     random coordinate perturbation to avoid numerical issues
c
      if (dowiggle) then
         eps = 0.001d0
         call wiggle (n,coords,eps)
      end if
c
c     transfer coordinates, complete to minimum of four spheres
c     if needed, set Delaunay and alpha complex arrays
c
      call setunion (nsphere,coords,radii)
c
c     compute the weighted Delaunay triangulation
c
      call regular3 (nredundant,redlist)
c
c     compute the alpha complex for fixed value of alpha
c
      alpha = 0.0d0
      call alfcx (alpha,nredundant,redlist)
c
c     if fewer than four balls, set artificial spheres as redundant
c
      call readjust_sphere (nsphere,nredundant,redlist)
c
c     get accessible surface area via the UnionBall method
c
      call ball_surf (weight,surf,usurf,asurfx)
c
c     copy surface area of each sphere into Tinker array
c
      do i = 1, n
         asurf(i) = asurfx(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (radii)
      deallocate (asurfx)
      deallocate (coords)
      deallocate (redlist)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine surface1  --  alpha shapes surface & derivs  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "surface1" computes the weighted solvent accessible surface
c     area of each atom and the first derivatives of the area with
c     respect to Cartesian coordinates via the inclusion-exclusion
c     method of Herbert Edelsbrunner based on alpha shapes
c
c     note for small or symmetric structures where alpha shapes
c     may fail, switch to the Richmond method
c
c     developed to facilitate calling UnionBall from Tinker by
c     Jay W. Ponder, Washington University, October 2023
c
c     literature references:
c
c     P. Mach and P. Koehl, "Geometric Measures of Large Biomolecules:
c     Surface, Volume, and Pockets", Journal of Computational Chemistry,
c     32, 3023-3038 (2011)
c
c     P. Koehl, A. Akopyan and H. Edelsbrunner, "Computing the Volume,
c     Surface Area, Mean, and Gaussian Curvatures of Molecules and Their
c     Derivatives", Journal of Chemical Information and Modeling, 63,
c     973-985 (2023)
c
c     variables and parameters:
c
c     nsphere    number of spheres/balls in the system
c     coords     coordinates of the center of each sphere
c     radii      radius value for each sphere
c     weight     weight value for each sphere
c     probe      radius value of the probe sphere
c     surf       weighted surface area of union of spheres
c     usurf      unweighted surface area of union of spheres
c     asurf      weighted area contribution of each sphere
c     dsurf      derivatives of the weighted surface area over
c                  coordinates of the sphere centers
c
c
      subroutine surface1 (rad,weight,probe,surf,asurf,dsurf)
      use atoms
      implicit none
      integer i,nsphere
      integer nsize,nfudge
      integer nredundant
      integer, allocatable :: redlist(:)
      real*8 surf,usurf,eps
      real*8 probe,alpha
      real*8 rad(*)
      real*8 weight(*)
      real*8 asurf(*)
      real*8 dsurf(3,*)
      real*8, allocatable :: radii(:)
      real*8, allocatable :: asurfx(:)
      real*8, allocatable :: coords(:,:)
      real*8, allocatable :: dsurfx(:,:)
      logical dowiggle
      character*6 symmtyp
c
c
c     check coordinates for linearity, planarity and symmetry
c
      symmtyp = 'NONE'
      call chksymm (symmtyp)
      dowiggle = .false.
      if (n.gt.2 .and. symmtyp.eq.'LINEAR')  dowiggle = .true.
      if (n.gt.3 .and. symmtyp.eq.'PLANAR')  dowiggle = .true.
c
c     use Richmond method for small symmetric structures
c
      if (dowiggle) then
         call richmond1 (n,x,y,z,rad,weight,probe,surf,asurf,dsurf)
         return
      end if
c
c     perform dynamic allocation of some local arrays
c
      nfudge = 10
      nsize = n + nfudge
      allocate (radii(nsize))
      allocate (asurfx(nsize))
      allocate (coords(3,nsize))
      allocate (dsurfx(3,nsize))
      allocate (redlist(nsize))
c
c     set the coordinates and sphere radii plus probe`
c
      nsphere = n
      do i = 1, n
         coords(1,i) = x(i)
         coords(2,i) = y(i)
         coords(3,i) = z(i)
         radii(i) = 0.0d0
         if (rad(i) .ne. 0.0d0)  radii(i) = rad(i) + probe
      end do
c
c     random coordinate perturbation to avoid numerical issues
c
      if (dowiggle) then
         eps = 0.001d0
         call wiggle (n,coords,eps)
      end if
c
c     transfer coordinates, complete to minimum of four spheres
c     if needed, set Delaunay and alpha complex arrays
c
      call setunion (nsphere,coords,radii)
c
c     compute the weighted Delaunay triangulation
c
      call regular3 (nredundant,redlist)
c
c     compute the alpha complex for fixed value of alpha
c
      alpha = 0.0d0
      call alfcx (alpha,nredundant,redlist)
c
c     if fewer than four balls, set artificial spheres as redundant
c
      call readjust_sphere (nsphere,nredundant,redlist)
c
c     get accessible surface area via the UnionBall method
c
      call ball_dsurf (weight,surf,usurf,asurfx,dsurfx)
c
c     copy surface area of each sphere into Tinker array
c
      do i = 1, n
         asurf(i) = asurfx(i)
         dsurf(1,i) = dsurfx(1,i)
         dsurf(2,i) = dsurfx(2,i)
         dsurf(3,i) = dsurfx(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (radii)
      deallocate (asurfx)
      deallocate (coords)
      deallocate (dsurfx)
      deallocate (redlist)
      return
      end
