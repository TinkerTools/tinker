c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine egaussvol3  --  gaussvol volume and analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "gaussvol" uses the tree algorithms from Emilio Gallicchio's
c     group to compute the molecular surface area and volume of a
c     collection of spherical atoms described as gaussians; also
c     partitions the volume and surface area among the atoms
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
      subroutine egaussvol3 (volume, area, self_volume, self_area)
      use atoms
      use gaussvolconst
      use gssvol
      use inform
      use iounit
      implicit none
      real*8 area
      real*8 energy
      real*8 volume,volume2
      real*8 self_volume(*)
      real*8 self_area(*)
      integer i,ia,k
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
c     partition self volume to each atom
c
      do i = 1, n
         self_volume(i) = gvself_volume(i)
         self_area(i) = gvself_volume(i)
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
c     compute surface area via finite difference
c
      area = (volume2 - volume)/rad_offset
      do i = 1, n
         self_area(i) = (gvself_volume(i) - self_area(i)) / rad_offset
      end do
c      write(*,*) "GaussVol Volume:  ", volume
c      write(*,*) "GaussVol Area:    ", area
c
c     print out the decomposition of the volume and area
c
      if (debug) then
         write (iout,10)
  10     format (/,' Self Volume for Individual Atoms :',/)
         k = 1
         do while (k .le. n)
            write (iout,20)  (ia,self_volume(ia),ia=k,min(k+4,n))
  20        format (1x,5(i7,f8.3))
            k = k + 5
         end do
         write (iout,30)
  30     format (//,' Self Area for Individual Atoms :',/)
         k = 1
         do while (k .le. n)
            write (iout,40)  (ia,self_area(ia),ia=k,min(k+4,n))
  40        format (1x,5(i7,f8.3))
            k = k + 5
         end do
      end if
c
c     perform deallocation of gaussvol objects
c
      call gvol%destroy()
      return
      end
