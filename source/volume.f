c
c
c     ################################################################
c     ##  COPYRIGHT (C) 1990 by Patrice Koehl & Jay William Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine volume  --  alpha shapes excluded volume  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "volume" computes the weighted solvent excluded volume via
c     the inclusion-exclusion method of Herbert Edelsbrunner based
c     on alpha shapes; also finds the accessible surface area
c
c     note for small or symmetric structures where alpha shapes
c     may fail, switch to the Connolly method
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
c     vol        weighted excluded volume of union of spheres
c     usurf      unweighted surface area of union of spheres
c     uvol       unweighted excluded volume of union of spheres
c     asurf      weighted area contribution of each sphere
c     avol       weighted volume contribution of each sphere
c
c
      subroutine volume (surf,vol,asurf,avol,rad,weight,probe)
      use atoms
      implicit none
      integer i,nsphere
      integer nsize,nfudge
      integer nredundant
      integer, allocatable :: listredundant(:)
      real*8 surf,usurf
      real*8 vol,uvol
      real*8 reentrant
      real*8 probe,alpha
      real*8 rad(*)
      real*8 weight(*)
      real*8 asurf(*)
      real*8 avol(*)
      real*8, allocatable :: radii(:)
      real*8, allocatable :: asurfx(:)
      real*8, allocatable :: avolx(:)
      real*8, allocatable :: coords(:,:)
      character*6 symmtyp
c
c
c     use Connolly method for small symmetric structures
c
      call chksymm (symmtyp)
      if (n.le.3 .or. (n.le.50.and.symmtyp.ne.'NONE')) then
         reentrant = 0.0d0
         call connolly (vol,surf,rad,reentrant,probe)
         return
      end if
c
c     perform dynamic allocation of some local arrays
c
      nfudge = 10
      nsize = n + nfudge
      allocate (radii(nsize))
      allocate (asurfx(nsize))
      allocate (avolx(nsize))
      allocate (coords(3,nsize))
      allocate (listredundant(nsize))
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
c     transfer coordinates, complete to minimum of four spheres
c     if needed, set Delaunay and alpha complex arrays
c
      call setunion (nsphere,coords,radii)
c
c     compute the weighted Delaunay triangulation
c
      call regular3 (nredundant,listredundant)
c
c     compute the alpha complex for fixed value of alpha
c
      alpha = 0.0d0
      call alfcx (alpha,nredundant,listredundant)
c
c     adjust number of vertices, if artificially set to four
c
      call readjust_nsphere (nsphere,nredundant,listredundant)
c
c     get the accessible surface area and excluded volume
c
      call ball_vol (weight,surf,vol,usurf,uvol,asurfx,avolx)
c
c     copy area and volume of each sphere into Tinker array
c
      do i = 1, n
         asurf(i) = asurfx(i)
         avol(i) = avolx(i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (radii)
      deallocate (asurfx)
      deallocate (avolx)
      deallocate (coords)
      deallocate (listredundant)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine volume1  --  alpha shapes volume & derivs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "volume1" computes the weighted solvent excluded volume
c     and the first derivatives of the volume with respect to
c     Cartesian coordinates via the inclusion-exclusion method
c     of Herbert Edelsbrunner based on alpha shapes; also finds
c     the accessible surface area and first derivatives
c
c     note for small or symmetric structures where alpha shapes
c     may fail, swith to Richmond, Connolly and Kundrot methods
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
c     vol        weighted excluded volume of union of spheres
c     usurf      unweighted surface area of union of spheres
c     uvol       unweighted excluded volume of union of spheres
c     asurf      weighted area contribution of each sphere
c     avol       weighted volume contribution of each sphere
c     dsurf      derivatives of the weighted surface area over
c                  coordinates of the sphere centers
c     dvol       derivatives of the weighted volume over
c                  coordinates of the sphere centers
c
c
      subroutine volume1 (surf,vol,asurf,avol,dsurf,dvol,
     &                       rad,weight,probe)
      use atoms
      implicit none
      integer i,nsphere
      integer nsize,nfudge
      integer nredundant
      integer, allocatable :: listredundant(:)
      real*8 surf,usurf
      real*8 vol,uvol
      real*8 reentrant
      real*8 probe,alpha
      real*8 rad(*)
      real*8 weight(*)
      real*8 asurf(*)
      real*8 avol(*)
      real*8 dsurf(3,*)
      real*8 dvol(3,*)
      real*8, allocatable :: radii(:)
      real*8, allocatable :: asurfx(:)
      real*8, allocatable :: avolx(:)
      real*8, allocatable :: coords(:,:)
      real*8, allocatable :: dsurfx(:,:)
      real*8, allocatable :: dvolx(:,:)
      character*6 symmtyp
c
c
c     use other methods for small symmetric structures
c
      call chksymm (symmtyp)
      if (n.le.3 .or. (n.le.50.and.symmtyp.ne.'NONE')) then
         reentrant = 0.0d0
         call richmond1 (surf,asurf,dsurf,rad,weight,probe)
         call connolly (vol,surf,rad,reentrant,probe)
         call kundrot1 (rad,probe,dvol)
         return
      end if
c
c     perform dynamic allocation of some local arrays
c
      nfudge = 10
      nsize = n + nfudge
      allocate (radii(nsize))
      allocate (asurfx(nsize))
      allocate (avolx(nsize))
      allocate (coords(3,nsize))
      allocate (dsurfx(3,nsize))
      allocate (dvolx(3,nsize))
      allocate (listredundant(nsize))
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
c     transfer coordinates, complete to minimum of four spheres
c     if needed, set Delaunay and alpha complex arrays
c
      call setunion (nsphere,coords,radii)
c
c     compute the weighted Delaunay triangulation
c
      call regular3 (nredundant,listredundant)
c
c     compute the alpha complex for fixed value of alpha
c
      alpha = 0.0d0
      call alfcx (alpha,nredundant,listredundant)
c
c     adjust number of vertices, if artificially set to four
c
      call readjust_nsphere (nsphere,nredundant,listredundant)
c
c     get the accessible surface area and excluded volume
c
      call ball_dvol (weight,surf,vol,usurf,uvol,asurfx,avolx,
     &                   dsurfx,dvolx)
c
c     copy area and volume of each sphere into Tinker array
c
      do i = 1, n
         asurf(i) = asurfx(i)
         avol(i) = avolx(i)
         dsurf(1,i) = dsurfx(1,i)
         dsurf(2,i) = dsurfx(2,i)
         dsurf(3,i) = dsurfx(3,i)
         dvol(1,i) = dvolx(1,i)
         dvol(2,i) = dvolx(2,i)
         dvol(3,i) = dvolx(3,i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (radii)
      deallocate (asurfx)
      deallocate (avolx)
      deallocate (coords)
      deallocate (dsurfx)
      deallocate (dvolx)
      deallocate (listredundant)
      return
      end
