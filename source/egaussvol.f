c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine egaussvol  --  gaussvol surface area & volume  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "gaussvol" uses the tree algorithms from Emilio Gallicchio's
c     group to compute the molecular surface area and volume of a
c     collection of spherical atoms described as gaussians
c
c     literature references:
c
c     first GaussVol paper
c
c     second GaussVol paper
c
c
      subroutine egaussvol (volume,area)
      use atoms
      use gaussvolconst
      use gssvol
      implicit none
      real*8 area
      real*8 energy
      real*8 volume,volume2
c
c
c     initialize gaussvol
c
      call gvol%GaussVol_init(n)
c
c     compute gaussvol volume and area
c
      call gvol%setRadii(gvradius)	
      call gvol%setVolumes(gvvol)	
      call gvol%compute_tree()	
c       call gvol%GaussVol_print_tree()	
      call gvol%compute_volume(volume, energy, gvdr, gvdv,	
     &          gvfree_volume, gvself_volume)	
      call gvol%setRadii(gvradius2)	
      call gvol%setVolumes(gvvol2)	
      call gvol%rescan_tree_volumes()	
c       call gvol%GaussVol_print_tree()	
      call gvol%compute_volume(volume2, energy, gvdr, gvdv,	
     &          gvfree_volume, gvself_volume)	
      area = (volume2 - volume)/rad_offset
      write(*,*) "GaussVol Volume:  ", volume
      write(*,*) "GaussVol Volume2:  ", volume2
      write(*,*) "GaussVol Area:    ", area
c
c     perform deallocation of gaussvol objects
c
      call gvol%destroy()
      return
      end
