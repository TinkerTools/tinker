c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module gaussvolmodule  --  module to compute GaussVol  ##
c     ##                                                         ##
c     #############################################################
c
c
c     tree         GOverlap_Tree object
c     natoms       number of atoms in the system
c     ishydrogen   d
c     radii        atomic radii
c     volumes      atomic volume
c     gammas       atomic gamma
c
c
      module gaussvolmodule
      use goverlapmodule
      use goverlaptreemodule
      implicit none
c
c
c     A class that implements the Gaussian description of an object
c     (molecule) made of overlapping spheres
c
      type GaussVol
      private
      type(GOverlap_Tree) tree
      integer natoms
      integer, allocatable :: ishydrogen(:)
      real*8, allocatable :: radii(:)
      real*8, allocatable :: volumes(:)
      real*8, allocatable :: gammas(:)
      contains
c
c     Creates/Initializes a GaussVol instance
c
      procedure, pass(this) :: GaussVol_init

c
c     Destructor of GaussVol object
c
      procedure, pass(this) :: destroy
c
c     Public methods of GaussVol object
c
      procedure, pass(this) :: setRadii
      procedure, pass(this) :: setVolumes
      procedure, pass(this) :: setGammas
      procedure, pass(this) :: compute_tree
      procedure, pass(this) :: compute_volume
      procedure, pass(this) :: rescan_tree_volumes
      procedure, pass(this) :: rescan_tree_gammas
      procedure, pass(this) :: getstat
      procedure, pass(this) :: GaussVol_print_tree
      end type GaussVol
      contains
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine GaussVol_init  --  initialize GaussVol object  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "GaussVol_init" initializes attributes of a GaussVol class object
c
c
      subroutine GaussVol_init(this, natoms, ishydrogen)
      class(GaussVol) this
      integer natoms, i
      integer ishydrogen(*)
c
c
c     set number of atoms and root of GOverlap_Tree class object
c
      this%natoms = natoms
      call this%tree%create_GOverlap_Tree(natoms)
c
c     perform dynamic allocation of GaussVol arrays
c
      allocate(this%radii(natoms))
      allocate(this%volumes(natoms))
      allocate(this%gammas(natoms))
      allocate(this%ishydrogen(natoms))
c
c     initialize GaussVol radii, volumes, gammas, and ishydrogen
c
      do i = 1, natoms
         this%radii(i) = 1.0d0
         this%volumes(i) = 0.0d0
         this%gammas(i) = 1.0d0
         this%ishydrogen(i) = ishydrogen(i)
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine destroy  --  deallocate GaussVol object  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "destroy" deallocates GaussVol object
c
c
      subroutine destroy(this)
      class(GaussVol) this
      call this%tree%overlaps%clear() 
      deallocate(this%radii)
      deallocate(this%volumes)
      deallocate(this%gammas)
      deallocate(this%ishydrogen)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine setRadii  --  set GaussVol radii  ##
c     ##                                               ##
c     ###################################################
c
c
c     "setRadii" sets GaussVol radii
c
c
      subroutine setRadii(this, radii)
      class(GaussVol) this
      integer i
      real*8 radii(*)
c
c
c     copy radii into GaussVol%radii
c
      do i = 1, this%natoms
         this%radii(i) = radii(i)
      end do
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine setVolumes  --  set GaussVol volumes  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "setVolumes" sets GaussVol volumes
c
c
      subroutine setVolumes(this, volumes)
      class(GaussVol) this
      integer i
      real*8 volumes(*)
c
c
c     copy volumes into GaussVol%volumes
c
      do i = 1, this%natoms
         this%volumes(i) = volumes(i)
      end do
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine setGammas  --  set GaussVol gammas  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "setGammas" sets GaussVol gammas
c
c
      subroutine setGammas(this, gammas)
      class(GaussVol) this
      integer i
      real*8 gammas(*)
c
c
c     copy gammas into GaussVol%gammas
c
      do i = 1, this%natoms
         this%gammas(i) = gammas(i)
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine compute_tree  --  constructs GOverlap_Tree  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "compute_tree" constructs the GOverlap_Tree object
c
c
      subroutine compute_tree(this, positions)
      class(GaussVol) this
      integer ret
      real*8 positions(3,*)
c
c
c     call recursive overlap tree method in GOverlap_Tree object
c
      ret = this%tree%compute_overlap_tree_r(positions,this%radii,
     &                 this%volumes, this%gammas, this%ishydrogen)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine compute_volume  --  compute GaussVol volume  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "compute_volume" computes the volume, energy, and forces of
c     the GaussVol object
c
c
      subroutine compute_volume(this, positions, volume, energy, 
     &                  force, gradV, free_volume, self_volume)
      class(GaussVol) this
      integer ret
      real*8 volume, energy
      real*8 free_volume(*)
      real*8 gradV(*)
      real*8 self_volume(*)
      real*8 force(3,*)
      real*8 positions(3,*)
c
c
c     call recursive volume method in GOverlap_Tree object
c
      ret = this%tree%compute_volume2_r(positions, volume, energy,
     &                    force, gradV, free_volume, self_volume)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine rescan_tree_volumes  --  rescans tree volumes  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "rescan_tree_volumes" rescans the tree after resetting volumes
c
c
      subroutine rescan_tree_volumes(this, positions)
      class(GaussVol) this
      real*8 positions(3,*)
      integer ret
c
c
c     call recursive rescan tree volume method in GOverlap_Tree object
c
      ret = this%tree%rescan_tree_v(positions, this%radii,
     &               this%volumes, this%gammas, this%ishydrogen)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine rescan_tree_gammas  --  rescans tree gammas  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "rescan_tree_gammas" rescans the tree after resetting gammas
c
c
      subroutine rescan_tree_gammas(this)
      class(GaussVol) this
      integer ret
c
c
c     call recursive rescan tree gamma method in GOverlap_Tree object
c
      ret = this%tree%rescan_tree_g(this%gammas)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine getstat  --  returns the number of overlaps  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "getstat" returns the number of overlaps for each atom
c
c
      subroutine getstat(this, nov)
      class(GaussVol) this
      integer atom
      integer nov(*)
      do atom = 1, this%natoms
         nov(atom) = 0
         nov(atom) = nchildren_under_slot_r(this%tree, atom)
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine GaussVol_print_tree  --  print GaussVol tree  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "GaussVol_print_tree" prints GaussVol tree
c
c
      subroutine GaussVol_print_tree(this)
      class(GaussVol) this
      call this%tree%print_tree()
      return
      end
      end module gaussvolmodule
