      module gaussvolmodule
        use gaussianvcamodule
        use goverlapmodule
        use goverlaptreemodule
        implicit none

        ! A class that implements the Gaussian description of an object
        ! (molecule) made of overlapping spheres
        type :: GaussVol
          private
          integer :: natoms
          type(GOverlap_Tree) :: tree
          real*8, dimension(:), allocatable :: radii
          real*8, dimension(:), allocatable :: volumes
          real*8, dimension(:), allocatable :: gammas
          integer, dimension(:), allocatable :: ishydrogen
        contains
          ! Creates/Initializes a GaussVol instance
          procedure, pass(this) :: init => GaussVol_init

          ! Destructor
          procedure, pass(this) :: destroy

          ! Public methods
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

        ! GaussVol initialization
         subroutine GaussVol_init(this, natoms, ishydrogen)
           class(GaussVol), intent(inout) :: this
           integer, intent(in) :: natoms
           integer, dimension(:), intent(in) :: ishydrogen
 
           this%natoms = natoms
           call this%tree%create_GOverlap_Tree(natoms)
   
           allocate(this%radii(natoms))
           this%radii = 1.0d0
   
           allocate(this%volumes(natoms))
           this%volumes = 0.0d0
   
           allocate(this%gammas(natoms))
           this%gammas = 1.0d0
   
           allocate(this%ishydrogen(natoms))
           this%ishydrogen = ishydrogen
         end subroutine GaussVol_init

        ! GaussVol destructor
        subroutine destroy(this)
          class(GaussVol), intent(inout) :: this
          call this%tree%overlaps%clear() 
          deallocate(this%radii)
          deallocate(this%volumes)
          deallocate(this%gammas)
          deallocate(this%ishydrogen)
        end subroutine destroy

        ! Sets the radii of the GaussVol
        subroutine setRadii(this, radii)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:), intent(in) :: radii

          if (this%natoms == size(radii)) then
            this%radii = radii
          else
            error stop "setRadii: number of atoms does not match"
          end if
        end subroutine setRadii

        ! Sets the volumes of the GaussVol
        subroutine setVolumes(this, volumes)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:), intent(in) :: volumes

          if (this%natoms == size(volumes)) then
            this%volumes = volumes
          else
            error stop "setVolumes: number of atoms does not match"
          end if
        end subroutine setVolumes

        ! Sets the gammas of the GaussVol
        subroutine setGammas(this, gammas)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:), intent(in) :: gammas
  
          if (this%natoms == size(gammas)) then
            this%gammas = gammas
          else
            error stop "setGammas: number of atoms does not match"
          end if
        end subroutine setGammas

        ! Constructs the tree
        subroutine compute_tree(this, positions)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:,:), intent(in) :: positions
          integer :: ret

          ret = this%tree%compute_overlap_tree_r(positions, this%radii,
     &       this%volumes, this%gammas, this%ishydrogen)
        end subroutine compute_tree

      ! Computes the volume, energy, and forces of the GaussVol
        subroutine compute_volume(this, positions, volume, energy, 
     &     force, gradV, free_volume, self_volume)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:,:), intent(in) :: positions
          real*8, intent(out) :: volume, energy
          real*8, dimension(:,:), intent(out) :: force
          real*8, dimension(:), intent(out) :: gradV, free_volume
          real*8, dimension(:), intent(out) :: self_volume
          integer :: ret

          ret = this%tree%compute_volume2_r(positions, volume, energy,
     &      force, gradV, free_volume, self_volume)
        end subroutine compute_volume

        ! Rescans the tree after resetting volumes
        subroutine rescan_tree_volumes(this, positions)
          class(GaussVol), intent(inout) :: this
          real*8, dimension(:,:), intent(in) :: positions
          integer :: ret

          ret = this%tree%rescan_tree_v(positions, this%radii,
     &     this%volumes, this%gammas, this%ishydrogen)
        end subroutine rescan_tree_volumes

        ! Rescans the tree after resetting gammas
        subroutine rescan_tree_gammas(this)
          class(GaussVol), intent(inout) :: this
          integer :: ret

          ret = this%tree%rescan_tree_g(this%gammas)
        end subroutine rescan_tree_gammas

        ! Returns the number of overlaps for each atom
        subroutine getstat(this, nov)
          class(GaussVol), intent(inout) :: this
          integer, intent(inout), allocatable :: nov(:)
          integer :: atom, slot

          allocate(nov(this%natoms))
          nov = 0

          do atom = 0, this%natoms - 1
            slot = atom + 1
            nov(atom) = nchildren_under_slot_r(this%tree, slot)
          end do

        end subroutine getstat

        ! Prints the tree
        subroutine GaussVol_print_tree(this)
          class(GaussVol), intent(inout) :: this

          call this%tree%print_tree()
        end subroutine GaussVol_print_tree

      end module GaussVolModule
