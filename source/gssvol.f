c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ########################################################
c     ##                                                    ##
c     ##  module gssvol  --  GaussVol in current structure  ##
c     ##                                                    ##
c     ########################################################
c
c
c     gvol             gaussvol object
c     gvradius         vdw radius at each site
c     gvradius2        vdw radius + offset at each site
c     gvvol            vdw volume at each site
c     gvvol2           vdw volume offset at each site
c     gvdv             derivative wrt to volume at each site
c     gvfree_volume    free volume at each site
c     gvself_volume    self volume at each site
c     gvdr             derivative wrt to position at each site
c
c
      module gssvol
      use gaussvolmodule
      implicit none
      type(GaussVol) gvol
      real*8, allocatable :: gvradius(:)
      real*8, allocatable :: gvradius2(:)
      real*8, allocatable :: gvvol(:)
      real*8, allocatable :: gvvol2(:)
      real*8, allocatable :: gvdv(:)
      real*8, allocatable :: gvfree_volume(:)
      real*8, allocatable :: gvself_volume(:)
      real*8, allocatable :: gvdr(:,:)
      save
      end
