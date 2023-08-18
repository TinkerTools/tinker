      module goverlaptreemodule
        use gaussvolconst
        use gaussianvcamodule
        use goverlapmodule
        use vectormodule
        implicit none

        type :: GOverlap_Tree
          integer :: natoms
          type(Vector) :: overlaps
          contains
            procedure, pass(this) :: create_GOverlap_Tree 
            procedure, pass(this) :: init_overlap_tree 
            procedure, pass(this) :: add_children 
            procedure, pass(this) :: compute_children 
            procedure, pass(this) :: compute_andadd_children_r 
            procedure, pass(this) :: compute_overlap_tree_r 
            procedure, pass(this) :: compute_volume_underslot2_r 
            procedure, pass(this) :: compute_volume2_r
            procedure, pass(this) :: rescan_r 
            procedure, pass(this) :: rescan_tree_v 
            procedure, pass(this) :: rescan_gamma_r 
            procedure, pass(this) :: rescan_tree_g 
            procedure, pass(this) :: print_tree 
            procedure, pass(this) :: print_tree_r 
            procedure, pass(this) :: nchildren_under_slot_r 
        end type GOverlap_Tree

        contains

          subroutine create_GOverlap_Tree(this, natoms)
            class(GOverlap_Tree), intent(inout) :: this
            integer, intent(in) :: natoms

            this%natoms = natoms
            call this%overlaps%create_vector()  ! Use create_vector subroutine from VectorModule

            ! The root is at index 0, atoms are at 1..natoms
            call this%overlaps%push_back(GOverlap(level=0,
     &        g=GaussianVca(v=0.0d0,a=0.0d0,c=(/ 0.0d0,0.0d0,0.0d0 /)),
     &        volume=0.0d0, dvv1=0.0d0, dv1=(/ 0.0d0, 0.0d0, 0.0d0 /), 
     &        gamma1i=0.0d0, self_volume=0.0d0, sfp=0.0d0, atom=0, 
     &        parent_index=0, children_startindex=0, children_count=0))
        end subroutine create_GOverlap_Tree

        subroutine init_overlap_tree(this, pos, radius, volume, gama,
     &          ishydrogen)
          class(GOverlap_Tree), intent(inout) :: this
          real*8, dimension(:,:), intent(in) :: pos
          real*8, dimension(:), intent(in) :: radius, volume, gama
          integer, dimension(:), intent(in) :: ishydrogen
          type(GOverlap), pointer :: overlap
          integer :: iat
          real*8 :: a, vol

          ! Reset tree
          call this%overlaps%clear()  ! Clear the overlaps vector

          ! Slot 0 contains the master tree information, children = all of the atoms
          allocate(overlap)
          overlap%level = 0
          overlap%volume = 0.0d0
          overlap%dv1 = [0.0d0, 0.0d0, 0.0d0]
          overlap%dvv1 = 0.0d0
          overlap%self_volume = 0.0d0
          overlap%sfp = 1.0d0
          overlap%gamma1i = 0.0d0
          overlap%parent_index = -1
          overlap%atom = -1
          overlap%children_startindex = 1
          overlap%children_count = this%natoms

          ! Add the root overlap to the vector.
          call this%overlaps%push_back(overlap) 

          ! List of atoms start at slot 1.
          do iat = 1, this%natoms
            allocate(overlap)
            overlap%level = 1
            a = KFC / (radius(iat) * radius(iat))
            if (ishydrogen(iat) > 0) then
              a = 0.0d0
              vol = 0.0d0
            else
              vol = volume(iat)
            end if
            overlap%g%v = vol
            overlap%g%a = a
            overlap%g%c = pos(:, iat)
            overlap%volume = vol
            overlap%dv1 = [0.0d0, 0.0d0, 0.0d0]
            overlap%dvv1 = 1.0d0
            overlap%self_volume = 0.0d0
            overlap%sfp = 1.0d0
            overlap%gamma1i = gama(iat)
            overlap%parent_index = 0
            overlap%atom = iat
            overlap%children_startindex = -1
            overlap%children_count = -1
            ! Add the atom to the vector
            call this%overlaps%push_back(overlap)
          end do

        end subroutine init_overlap_tree
  
        ! Adds to the tree the children of overlap identified by
        ! "parent_index" in the tree.
        function add_children(this, parent_index, children_overlaps)
     &        result(start_index)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: parent_index
          type(Vector), intent(inout) :: children_overlaps
          integer :: i, ip, slot, start_index, noverlaps
          integer :: parent_level
          ! Retrieves address of root overlap
          type(GOverlap) :: parent, child_overlap

          ! Adds children starting at the last slot
          start_index = this%overlaps%vector_size()
          noverlaps = children_overlaps%vector_size()
          parent = this%overlaps%get_element(parent_index)

          ! Registers list of children
          parent%children_startindex = start_index
          parent%children_count = noverlaps
          ! Write back the changes,
          call this%overlaps%set_element(parent_index, parent)

          ! Sort neighbors by overlap volume
          call sort(children_overlaps)

          parent_level = parent%level

          ! Now copies the children overlaps from temp buffer
          do ip = 0, noverlaps - 1
            child_overlap = children_overlaps%get_element(ip)
            child_overlap%level = parent_level + 1
            ! Connect overlap to parent
            child_overlap%parent_index = parent_index
            ! Reset its children indexes
            child_overlap%children_startindex = -1
            child_overlap%children_count = -1
            ! Add to tree
            call this%overlaps%push_back(child_overlap)
          end do
        end function add_children

        ! Scans the siblings of overlap identified by "root_index" to create children overlaps, 
        ! returns them into the "children_overlaps" buffer: (root) + (atom) -> (root, atom).
        function compute_children(this, root_index, 
     &                            children_overlaps) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: root_index
          type(Vector), intent(inout) :: children_overlaps
          type(Vector) :: overlaps
          type(GOverlap) :: root, parent, sibling, ov2
          type(GOverlap), pointer :: ov
          type(GaussianVca) :: g1, g2, g12
          integer :: parent_index, sibling_start, sibling_count
          integer :: j, atom2, ret
          real*8 :: gvol, dVdr, dVdV, sfp

          ! Reset output buffer
          call children_overlaps%clear()

          ! Retrieves overlap
          overlaps = this%overlaps
          root = overlaps%get_element(root_index)

          ! Retrieves parent overlap
          parent_index = root%parent_index
          if (parent_index < 0) then
            write(*,*) "Error -- invalid parent index"
            ret = 1 ! master root? can't compute children
            return
          end if
          if (root%level >= MAX_ORDER) then
            write(*,*) "Reached max level", root%level
            ret = 1
            return
          end if
          parent = overlaps%get_element(parent_index)

          ! Retrieves start index and count of siblings
          sibling_start = parent%children_startindex
          sibling_count = parent%children_count
          if (sibling_start < 0 .or. sibling_count < 0) then
           ! Parent is not initialized?
            write(*,*) "Error -- parent is not initialized"
            ret = -1
            return
          end if
          ! Check that this overlap is the child of the registered parent
          if (root_index < sibling_start .or. 
     &        root_index > sibling_start + sibling_count - 1) then
            write(*,*) "Error -- overlap parent is not registered"
            ret = -1 
            return
          end if

          ! Loop over "younger" siblings (i < j loop) to compute new overlaps
          do j = root_index + 1, sibling_start + sibling_count - 1
            sibling = overlaps%get_element(j)
            ! Atomic gaussian of the last atom of the sibling
            atom2 = sibling%atom
            g1 = root%g
            ov2 = overlaps%get_element(atom2)
            g2 = ov2%g
            gvol = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp)

            ! Create child if overlap volume is not zero
            if (gvol > MIN_GVOL) then
              allocate(ov)
              ov%g = g12
              ov%volume = gvol
              ov%self_volume = 0
              ov%atom = atom2
              ! dv1 is the gradient of V(123..)n with respect to the position of 1
              ov%dv1 = (g2%c - g1%c) * (-dVdr)
              ! dvv1 is the derivative of V(123...)n with respect to V(123...)
              ov%dvv1 = dVdV
              ov%sfp = sfp
              ov%gamma1i = root%gamma1i + ov2%gamma1i
              call children_overlaps%push_back(ov)
            end if
          end do

          ret = 1
        end function compute_children

        recursive function compute_andadd_children_r(this, root) 
     &      result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: root
          integer :: ret
          type(Vector) :: children_overlaps
          type(GOverlap) :: child, parent
          integer :: noverlaps, start_slot, ichild

          ! Compute children overlaps
          call children_overlaps%create_vector()
          ret = this%compute_children(root, children_overlaps)

          ! Get the number of overlaps
          noverlaps = children_overlaps%vector_size()

          if (noverlaps > 0) then
            ! Add children overlaps and get the starting slot
            start_slot = this%add_children(root, children_overlaps)

            ! Recursively compute and add children for each child overlap
            do ichild = start_slot, start_slot + noverlaps - 1
              ret = this%compute_andadd_children_r(ichild)
            end do
          end if

          ret = 1
        end function compute_andadd_children_r

        function compute_overlap_tree_r(this, pos, radius, volume, gama, 
     &      ishydrogen) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          real*8, dimension(:,:), intent(in) :: pos
          real*8, dimension(:), intent(in) :: radius, volume, gama
          integer, dimension(:), intent(in) :: ishydrogen
          integer :: ret, slot

          ! Initialize overlap tree
          call this%init_overlap_tree(pos, radius, volume, gama, 
     &      ishydrogen)

          ! Compute and add children overlaps recursively for each atom
          do slot = 1, this%natoms
            ret = this%compute_andadd_children_r(slot)
          end do

          ret = 1
        end function compute_overlap_tree_r

        ! Compute volumes, energy of this volume and calls itself to get the volumes of the children
        recursive function compute_volume_underslot2_r(this, slot, 
     &      psi1i, f1i, p1i, psip1i, fp1i, pp1i, energy1i, fenergy1i, 
     &      penergy1i, dr, dv, free_volume, self_volume) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: slot
          real*8, intent(out) :: psi1i, f1i, psip1i, fp1i
          real*8, intent(out) :: energy1i, fenergy1i
          real*8, dimension(3), intent(out) :: p1i, pp1i
          real*8, dimension(3), intent(out) :: penergy1i
          real*8, dimension(:,:), intent(inout) :: dr
          real*8, dimension(:), intent(inout) :: dv, free_volume
          real*8, dimension(:), intent(inout) :: self_volume
          integer :: ret

          type(GOverlap) :: ov, ov2
          real*8 :: cf, volcoeff, volcoeffp
          integer :: atom, i, j, sloti
          real*8 :: ai, a1i, a1, c1, c1p, c2
          real*8 :: psi1it, f1it, psip1it, fp1it
          real*8 :: energy1it, fenergy1it
          real*8, dimension(3) :: p1it, pp1it, penergy1it

          ov = this%overlaps%get_element(slot)
          cf = 1.0d0
          if (mod(ov%level, 2) == 0) then
            cf = -1.0d0
          end if
          volcoeff = 0.0d0
          if (ov%level > 0) then
            volcoeff = cf
          end if
          volcoeffp = 0.0d0
          if (ov%level > 0) then
            volcoeffp = volcoeff / real(ov%level, kind=8)
          end if

          atom = ov%atom
          ov2 = this%overlaps%get_element(atom)
          ai = ov2%g%a
          a1i = ov%g%a
          a1 = a1i - ai

          psi1i = volcoeff * ov%volume
          f1i = volcoeff * ov%sfp
          p1i = [0.0d0, 0.0d0, 0.0d0]

          psip1i = volcoeffp * ov%volume
          fp1i = volcoeffp * ov%sfp
          pp1i = [0.0d0, 0.0d0, 0.0d0]

          energy1i = volcoeffp * ov%gamma1i * ov%volume
          fenergy1i = volcoeffp * ov%sfp * ov%gamma1i
          penergy1i = [0.0d0, 0.0d0, 0.0d0]

          if (ov%children_startindex >= 0) then
            do sloti = ov%children_startindex, ov%children_startindex 
     &          + ov%children_count - 1
              ret = compute_volume_underslot2_r(this, sloti, psi1it, 
     &          f1it, p1it, psip1it, fp1it, pp1it, energy1it, 
     &          fenergy1it, penergy1it, dr, dv, free_volume, 
     &          self_volume)

              psi1i = psi1i + psi1it
              f1i = f1i + f1it
              p1i = p1i + p1it

              psip1i = psip1i + psip1it
              fp1i = fp1i + fp1it
              pp1i = pp1i + pp1it

              energy1i = energy1i + energy1it
              fenergy1i = fenergy1i + fenergy1it
              penergy1i = penergy1i + penergy1it
            end do
          end if

          if (ov%level > 0) then
            ! Contributions to free and self volume of last atom
            free_volume(atom) = free_volume(atom) + psi1i
            self_volume(atom) = self_volume(atom) + psip1i

            ! Contributions to energy gradients
            c2 = ai / a1i
            dr(atom, :) = dr(atom, :) 
     &                  + (-ov%dv1) * fenergy1i + penergy1i * c2
            dv(atom) = dv(atom) + ov%g%v * fenergy1i ! Will be divided by Vatom later

            ! Update subtree P1..i's for parent
            c2 = a1 / a1i
            p1i = (ov%dv1) * f1i + p1i * c2
            pp1i = (ov%dv1) * fp1i + pp1i * c2
            penergy1i = (ov%dv1) * fenergy1i + penergy1i * c2

            ! Update subtree F1..i's for parent
            f1i = ov%dvv1 * f1i
            fp1i = ov%dvv1 * fp1i
            fenergy1i = ov%dvv1 * fenergy1i
          end if

          ret = 1
        end function compute_volume_underslot2_r

        ! Traverses tree and computes volumes, etc.
        function compute_volume2_r(this, pos, volume, energy, dr, dv,
     &      free_volume, self_volume) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          real*8, intent(in) :: pos(:,:)
          real*8, intent(out) :: volume, energy
          real*8, dimension(:,:), intent(inout) :: dr
          real*8, dimension(:), intent(inout) :: dv, free_volume
          real*8, dimension(:), intent(inout) :: self_volume
          integer :: ret, slot, i, j
          real*8 :: psi1i, f1i, psip1i, fp1i, energy1i, fenergy1i
          real*8, dimension(3) :: p1i, pp1i, penergy1i
          real*8, dimension(3) :: zero3

          slot = 0

          ! Reset volumes, gradients
          zero3 = [0.0d0, 0.0d0, 0.0d0]
          dr = reshape(zero3, shape(dr))
          dv = 0.0d0
          free_volume = 0.0d0
          self_volume = 0.0d0

          ret = compute_volume_underslot2_r(this, slot, psi1i, f1i, p1i,
     &      psip1i, fp1i, pp1i, energy1i, fenergy1i, penergy1i, dr, dv, 
     &      free_volume, self_volume)

          volume = psi1i
          energy = energy1i
          ret = 1
        end function compute_volume2_r

        ! Rescan the sub-tree to recompute the volumes, does not modify the tree.
        recursive function rescan_r(this, slot) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: slot
          integer :: parent_index, sibling_start, sibling_count
          integer :: atom, slot_child, ret
          type(GOverlap) :: ov, ov2, parent
          type(GaussianVca) :: g1, g2, g12
          real*8 :: dVdr, dVdV, dVdalpha, d2Vdalphadr, d2VdVdr
          real*8 :: sfp, gvol

          ov = this%overlaps%get_element(slot)

          parent_index = ov%parent_index
          if (parent_index > 0) then
            atom = ov%atom
            parent = this%overlaps%get_element(parent_index)
            g1 = parent%g
            ov2 = this%overlaps%get_element(atom)
            g2 = ov2%g
            gvol = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp)
            ov%g = g12
            ov%volume = gvol
            ov%dv1 = (g2%c - g1%c) * (-dVdr)
            ov%dvv1 = dVdV
            ov%sfp = sfp
            ov%gamma1i = parent%gamma1i + ov2%gamma1i
            ! Write back the changes.
            call this%overlaps%set_element(slot, ov)
          end if

          do slot_child = ov%children_startindex, 
     &       ov%children_startindex + ov%children_count - 1
             ret = rescan_r(this, slot_child)
          end do

          ret = 1
        end function rescan_r

        ! Rescan the tree to recompute the volumes, does not modify the tree.
        recursive function rescan_tree_v(this, pos, radius, volume, 
     &      gama, ishydrogen) result(ret)
          class(GOverlap_Tree), intent(inout) :: this
          real*8, intent(in) :: pos(:,:)
          real*8, intent(in) :: radius(:)
          real*8, intent(in) :: volume(:)
          real*8, intent(in) :: gama(:)
          integer, intent(in) :: ishydrogen(:)
          integer :: iat
          type(GOverlap) :: ov
          real*8 :: a, vol, ret

          ov = this%overlaps%get_element(0)
          ov%level = 0
          ov%volume = 0
          ov%dv1 = [0.0d0, 0.0d0, 0.0d0]
          ov%dvv1 = 0.0d0
          ov%self_volume = 0
          ov%sfp = 1.0d0
          ov%gamma1i = 0.0d0
          ! Write back updates.
          call this%overlaps%set_element(0, ov)

          do iat = 1, this%natoms
            ov = this%overlaps%get_element(iat)
            a = KFC / (radius(iat) * radius(iat))
            vol = merge(0.0d0, volume(iat), ishydrogen(iat) > 0)
            ov%level = 1
            ov%g%v = vol
            ov%g%a = a
            ov%g%c = pos(:, iat)
            ov%volume = vol
            ov%dv1 = [0.0d0, 0.0d0, 0.0d0]
            ov%dvv1 = 1.0d0
            ov%self_volume = 0.0d0
            ov%sfp = 1.0d0
            ov%gamma1i = gama(iat)
            ! Write back updates.
            call this%overlaps%set_element(iat, ov)
          end do

          ret = this%rescan_r(0)

          ret = 1
        end function rescan_tree_v

        ! Rescan the sub-tree to recompute the gammas, does not modify
        ! the volumes nor the tree.
        recursive function rescan_gamma_r(this, slot) result(res)
          class(GOverlap_Tree), intent(inout) :: this
          integer, intent(in) :: slot
          integer :: parent_index, sibling_start, sibling_count
          logical :: res
          type(GOverlap) :: ov, ov2, parent
          integer :: atom, slot_child

          ! This overlap
          ov = this%overlaps%get_element(slot)

          ! Recompute its own overlap by merging parent and last atom
          parent_index = ov%parent_index
          if (parent_index > 0) then
            atom = ov%atom
            parent = this%overlaps%get_element(parent_index)
            ov2 = this%overlaps%get_element(atom)
            ov%gamma1i = parent%gamma1i + ov2%gamma1i
            ! Write back updates.
            call this%overlaps%set_element(slot, ov)
          endif

          ! Calls itself recursively on the children
          do slot_child = ov%children_startindex, 
     &        ov%children_startindex + ov%children_count-1
            res = rescan_gamma_r(this, slot_child)
          end do

          res = .true.
        end function rescan_gamma_r

        ! Rescan the tree to recompute the gammas only, does not modify
        ! volumes and the tree.
        function rescan_tree_g(this, gama) result(res)
          class(GOverlap_Tree), intent(inout) :: this
          real*8, intent(in) :: gama(:)
          logical :: res
          integer :: iat
          type(GOverlap) :: ov

          ov = this%overlaps%get_element(0)
          ov%gamma1i = 0.0d0
          call this%overlaps%set_element(0, ov)

          do iat = 1, this%natoms
            ov = this%overlaps%get_element(iat)
            ov%gamma1i = gama(iat)
            call this%overlaps%set_element(iat, ov)
          end do

          res = this%rescan_gamma_r(0)

          res = .true.
        end function rescan_tree_g

        subroutine print_tree(this)
          class(GOverlap_Tree), intent(in) :: this
          integer :: i

          write(*, '(A)', advance='no') "slot level LastAtom parent 
     &      ChStart ChCount SelfV V gamma a x y z dedx dedy dedz sfp"
          write(*, '(A)', advance='no') new_line('A')
  
          call print_tree_r(this, 0)
!          do i = 1, this%natoms
!            call print_tree_r(this, i)
!          end do

        end subroutine print_tree

        recursive subroutine print_tree_r(this, slot)
          class(GOverlap_Tree), intent(in) :: this
          integer, intent(in) :: slot
          type(GOverlap) :: ov
          integer :: i

          ov = this%overlaps%get_element(slot)

          write(*, '(A,I6,A)', advance='no') "tg: ", slot, " "
          call ov%print_overlap()

          do i = ov%children_startindex, 
     &           ov%children_startindex + ov%children_count - 1
            call print_tree_r(this, i)
          end do
        end subroutine print_tree_r

        recursive function nchildren_under_slot_r(this, slot) result(n)
          class(GOverlap_Tree), intent(in) :: this
          integer, intent(in) :: slot
          integer :: i, n
          type(GOverlap) :: ov

          n = 0
          ov = this%overlaps%get_element(slot)
          if (ov%children_count > 0) then
            n = n + ov%children_count
            ! now calls itself on the children
            do i = 0, ov%children_count - 1
              n = n + this%nchildren_under_slot_r( 
     &                  ov%children_startindex + i)
            end do
          end if

          return
        end function nchildren_under_slot_r

      end module goverlaptreemodule
