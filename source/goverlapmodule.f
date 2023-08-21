      module goverlapmodule
        use gaussianvcamodule
        implicit none

        type :: GOverlap
          integer :: level               ! level (0=root, 1=atoms, 2=2-body, 3=3-body, etc.)
          type(GaussianVca) :: g         ! Gaussian representing overlap 
          real*8 :: volume               ! volume of overlap (also stores Psi1..i in GPU version)
          real*8 :: dvv1                 ! derivative wrt volume of first atom (also stores F1..i in GPU version)
          real*8, dimension(3) :: dv1    ! derivative wrt position of first atom (also stores P1..i in GPU version)
          real*8 :: gamma1i              ! sum gammai for this overlap
          real*8 :: self_volume          ! self volume accumulator (also stores Psi1..i in GPU version)
          real*8 :: sfp                  ! switching function derivatives
          integer :: atom                ! the atomic index of the last atom of the overlap list (i, j, k, ..., atom) ! = (Parent, atom)
          integer :: parent_index        ! index in tree list of parent overlap
          integer :: children_startindex ! start index in tree array of children
          integer :: children_count      ! number of children contains
          contains
            procedure, pass(this) :: create_GOverlap
            procedure, pass(this) :: print_overlap 
        end type GOverlap

        contains

        subroutine print_overlap(this)
           class(GOverlap), intent(in) :: this
           ! Print the overlap details
  
           write(*, "('           ', I6, I7, I7, I8, I8, 11F11.6)")
     &       this%level, this%atom, this%parent_index,
     &       this%children_startindex, this%children_count,
     &       this%volume, this%gamma1i, this%g%a,
     &       this%g%v, this%g%c(1), this%g%c(2), this%g%c(3), 
     &       this%dv1(1), this%dv1(2), this%dv1(3), this%sfp
        end subroutine print_overlap

        subroutine create_GOverlap(this, level, g, volume, dvv1, dv1,
     &    gamma1i, self_volume, sfp, atom, parent_index,
     &    children_startindex, children_count)
          class(GOverlap), intent(inout) :: this
          integer, intent(in) :: level
          type(GaussianVca), intent(in) :: g
          real*8, intent(in) :: volume, dvv1, gamma1i
          real*8, intent(in) :: self_volume, sfp
          real*8, dimension(3), intent(in) :: dv1
          integer, intent(in) :: atom, parent_index
          integer, intent(in) :: children_startindex, children_count
          ! Initialize the overlap properties
          this%level = level
          this%g = g
          this%volume = volume
          this%dvv1 = dvv1
          this%dv1 = dv1
          this%gamma1i = gamma1i
          this%self_volume = self_volume
          this%sfp = sfp
          this%atom = atom
          this%parent_index = parent_index
          this%children_startindex = children_startindex
          this%children_count = children_count
        end subroutine create_GOverlap

      end module goverlapmodule
