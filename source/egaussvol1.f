c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egaussvol1  --  gaussvol volume and gradient  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "gaussvol" uses the tree algorithms from Emilio Gallicchio's
c     group to compute the molecular surface area and volume of a
c     collection of spherical atoms described as gaussians; also
c     computes volume and surface area gradients
c
c     literature references:
c
c     J.A. Grant and B.T. Pickup, "A Gaussian Description of Molecular
c     Shape", Journal of Physical Chemistry, 99, 3503-3510 (1995).
c
c     E. Gallicchio, and R.M. Levy, "AGBNP, An Analytic Implicit Solvent
c     Model Suitable for Molecular Dynamics Simulations and
c     High-Resolution Modeling", Journal of Computational Chememistry,
c     25, 479-499 (2004).
c
c     B. Zhang, D. Kilburg, P. Eastman, V.S. Pande, E. Gallicchio,
c     "Efficient Gaussian Density Formulation of Volume and Surface
c     Areas of Macromolecules on Graphical Processing Units", Journal 
c     of Computational Chemistry, 38, 740-752 (2017).
c
c     TODO:
c
c     1. Neighborlist implementation
c     2. Cutoff implementation
c     3. OpenMP parallelization
c
c
      subroutine egaussvol1 (volume, area, dvol, dsurf)
      use atoms
      use gaussvolconst
      use gssvol
      implicit none
      real*8 area
      real*8 energy
      real*8 volume,volume2
      real*8 dvol(3,*)
      real*8 dsurf(3,*)
      integer i
c
c
c     initialize gaussvol
c
      call gvol%GaussVol_init(n)
c
c     compute gaussvol volume
c
      call gvol%setRadii(gvradius)
      call gvol%setVolumes(gvvol)
      call gvol%compute_tree()
c      call gvol%GaussVol_print_tree()
      call gvol%compute_volume(volume, energy, gvdr, gvdv,
     &          gvfree_volume, gvself_volume)
c
c     copy volume gradient
c
      do i = 1, n
         dvol(1,i) = gvdr(1,i)
         dvol(2,i) = gvdr(2,i)
         dvol(3,i) = gvdr(3,i)
         dsurf(1,i) = gvdr(1,i)
         dsurf(2,i) = gvdr(2,i)
         dsurf(3,i) = gvdr(3,i)
      end do
c
c     compute volume with offset
c
      call gvol%setRadii(gvradius2)
      call gvol%setVolumes(gvvol2)
      call gvol%rescan_tree_volumes()
c       call gvol%GaussVol_print_tree()
      call gvol%compute_volume(volume2, energy, gvdr, gvdv,
     &          gvfree_volume, gvself_volume)
c
c     compute surface area gradient via finite difference
c
      area = (volume2 - volume)/rad_offset
      do i = 1, n
         dsurf(1,i) = (gvdr(1,i) - dsurf(1,i)) / rad_offset
         dsurf(2,i) = (gvdr(2,i) - dsurf(2,i)) / rad_offset
         dsurf(3,i) = (gvdr(3,i) - dsurf(3,i)) / rad_offset
      end do
c      write(*,*) "GaussVol Volume:  ", volume
c      write(*,*) "GaussVol Area:    ", area
c
c     perform deallocation of gaussvol objects
c
      call gvol%destroy()
      return
      end
