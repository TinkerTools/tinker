      module gaussvolmodule
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
            procedure, pass(this) :: GaussVol_init

c            ! Destructor
c            procedure, pass(this) :: destroy
c   
            ! Public methods
            procedure, pass(this) :: setRadii
            procedure, pass(this) :: setVolumes
c            procedure, pass(this) :: setGammas
            procedure, pass(this) :: compute_tree
            procedure, pass(this) :: compute_volume
c            procedure, pass(this) :: rescan_tree_volumes
c            procedure, pass(this) :: rescan_tree_gammas
c            procedure, pass(this) :: getstat
c            procedure, pass(this) :: GaussVol_print_tree
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
c
c        ! GaussVol destructor
c        subroutine destroy(this)
c          class(GaussVol), intent(inout) :: this
c          call this%tree%overlaps%clear() 
c          deallocate(this%radii)
c          deallocate(this%volumes)
c          deallocate(this%gammas)
c          deallocate(this%ishydrogen)
c        end subroutine destroy
c
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
c
c        ! Sets the gammas of the GaussVol
c        subroutine setGammas(this, gammas)
c          class(GaussVol), intent(inout) :: this
c          real*8, dimension(:), intent(in) :: gammas
c  
c          if (this%natoms == size(gammas)) then
c            this%gammas = gammas
c          else
c            error stop "setGammas: number of atoms does not match"
c          end if
c        end subroutine setGammas
c
         ! Constructs the tree
         subroutine compute_tree(this, positions)
            class(GaussVol), intent(inout) :: this
            real*8, dimension(:,:), intent(in) :: positions
            integer :: ret
            ret = this%tree%compute_overlap_tree_r(positions,this%radii,
     &                 this%volumes, this%gammas, this%ishydrogen)
         end subroutine compute_tree

         ! Computes the volume, energy, and forces of the GaussVol
         subroutine compute_volume(this, positions, volume, energy, 
     &          force, gradV, free_volume, self_volume)
            class(GaussVol), intent(inout) :: this
            real*8, dimension(:,:), intent(in) :: positions
            real*8, intent(out) :: volume, energy
            real*8, dimension(:,:), intent(out) :: force
            real*8, dimension(:), intent(out) :: gradV, free_volume
            real*8, dimension(:), intent(out) :: self_volume
            integer :: ret

            ret = this%tree%compute_volume2_r(positions, volume, energy,
     &              force, gradV, free_volume, self_volume)
         end subroutine compute_volume
c
c        ! Rescans the tree after resetting volumes
c        subroutine rescan_tree_volumes(this, positions)
c          class(GaussVol), intent(inout) :: this
c          real*8, dimension(:,:), intent(in) :: positions
c          integer :: ret
c
c          ret = this%tree%rescan_tree_v(positions, this%radii,
c     &     this%volumes, this%gammas, this%ishydrogen)
c        end subroutine rescan_tree_volumes
c
c        ! Rescans the tree after resetting gammas
c        subroutine rescan_tree_gammas(this)
c          class(GaussVol), intent(inout) :: this
c          integer :: ret
c
c          ret = this%tree%rescan_tree_g(this%gammas)
c        end subroutine rescan_tree_gammas
c
c        ! Returns the number of overlaps for each atom
c        subroutine getstat(this, nov)
c          class(GaussVol), intent(inout) :: this
c          integer, intent(inout), allocatable :: nov(:)
c          integer :: atom, slot
c
c          allocate(nov(this%natoms))
c          nov = 0
c
c          do atom = 0, this%natoms - 1
c            slot = atom + 1
c            nov(atom) = nchildren_under_slot_r(this%tree, slot)
c          end do
c
c        end subroutine getstat
c
c        ! Prints the tree
c        subroutine GaussVol_print_tree(this)
c          class(GaussVol), intent(inout) :: this
c
c          call this%tree%print_tree()
c        end subroutine GaussVol_print_tree
c
      end module gaussvolmodule
