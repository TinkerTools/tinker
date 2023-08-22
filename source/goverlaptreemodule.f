      ! Overlap volume switching function + 1st derivative 
      subroutine pol_switchfunc(gvol, volmina, volminb, sp, s)
         real*8 gvol, volmina, volminb, sp, s
         real*8 swf, swfp, swd, swu, swu2, swu3

         swf = 0.0d0
         swfp = 1.0d0

         if (gvol > volminb) then
            swf = 1.0d0
            swfp = 0.0d0
         else if (gvol < volmina) then
            swf = 0.0d0
            swfp = 0.0d0
         end if

         swd = 1.0d0 / (volminb - volmina)
         swu = (gvol - volmina) * swd
         swu2 = swu * swu
         swu3 = swu * swu2
         s = swf + swfp * swu3 * (10.0d0 - 15.0d0 * swu + 6.0d0 * swu2)
         sp = swfp * swd * 30.0d0 * swu2 * (1.0d0 - 2.0d0 * swu + swu2)

         ! turn off switching function
         !sp = 0.0d0
         !s = 1.0d0
      return
      end

      ! overlap between two Gaussians represented by a (V,c,a) triplet
      !   V: volume of Gaussian
      !   c: position of Gaussian
      !   a: exponential coefficient
      !
      !   g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
      !
      !   this version is based on V=V(V1,V2,r1,r2,alpha)
      !   alpha = (a1 + a2)/(a1 a2)
      !
      !   dVdr is (1/r)*(dV12/dr)
      !   dVdV is dV12/dV1 
      !   dVdalpha is dV12/dalpha
      !   d2Vdalphadr is (1/r)*d^2V12/dalpha dr
      !   d2VdVdr is (1/r) d^2V12/dV1 dr
      subroutine ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp, sgvol)
         use gaussianvcamodule
         use gaussvolconst
         use math
         implicit none
         type(GaussianVca) :: g1, g2, g12
         real*8 dVdr, dVdV, sfp, sgvol
         real*8 d2, deltai, gvol, p12, a12
         real*8 s, sp, df, dgvol, dgvolv, ef, dgvola2, dgvola1
         real*8 dgalpha, dgalpha2, dgvolvdr
         real*8, dimension(3) :: c1, c2, dist

         c1 = g1%c
         c2 = g2%c
         dist = c2 - c1
         d2 = dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3)

         a12 = g1%a + g2%a
         deltai = 1.0d0 / a12
         df = (g1%a) * (g2%a) * deltai

         ef = exp(-df * d2)
         gvol = ((g1%v * g2%v) / (PI / df)**1.5) * ef
         dgvol = -2.0d0 * df * gvol
         dgvolv = 0.0d0
         if (g1%v > 0.0d0) then
            dgvolv = gvol / g1%v
         end if

         g12%c = ((c1 * g1%a) + (c2 * g2%a)) * deltai
         g12%a = a12
         g12%v = gvol

         call pol_switchfunc(gvol, volmina, volminb, sp, s)
         sfp = sp * gvol + s
         dVdr = dgvol
         dVdV = dgvolv

         sgvol = s * gvol
      return
      end

      module goverlaptreemodule
         use gaussianvcamodule
         use gaussvolconst
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
c            procedure, pass(this) :: rescan_r 
c            procedure, pass(this) :: rescan_tree_v 
c            procedure, pass(this) :: rescan_gamma_r 
c            procedure, pass(this) :: rescan_tree_g 
c            procedure, pass(this) :: print_tree 
c            procedure, pass(this) :: print_tree_r 
c            procedure, pass(this) :: nchildren_under_slot_r 
         end type GOverlap_Tree

         contains

         subroutine create_GOverlap_Tree(this, natoms)
            class(GOverlap_Tree), intent(inout) :: this
            integer, intent(in) :: natoms

            this%natoms = natoms
            call this%overlaps%create_vector()  ! Use create_vector subroutine from vectormodule
         end subroutine create_GOverlap_Tree

         subroutine init_overlap_tree(this, pos, radius, volume, gamma,
     &          ishydrogen)
            class(GOverlap_Tree), intent(inout) :: this
            real*8, dimension(:,:), intent(in) :: pos
            real*8, dimension(:), intent(in) :: radius, volume, gamma
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
               a = kfc / (radius(iat) * radius(iat))
               vol = merge(0.0d0, volume(iat), ishydrogen(iat) > 0)
               overlap%g%v = vol
               overlap%g%a = a
               overlap%g%c = pos(:, iat)
               overlap%volume = vol
               overlap%dv1 = [0.0d0, 0.0d0, 0.0d0]
               overlap%dvv1 = 1.0d0
               overlap%self_volume = 0.0d0
               overlap%sfp = 1.0d0
               overlap%gamma1i = gamma(iat)
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
     &                                  result(start_index)
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
     &                                   children_overlaps) result(ret)
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
            if (root%level >= max_order) then
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
            if (root_index < sibling_start .and.
     &            root_index > sibling_start + sibling_count - 1) then
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
               call ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp, gvol)
               print*, "sfp, gvol", sfp, gvol
               print*, "g1v g1a g1c",g1%v,g1%a,g1%c(1),g1%c(2),g1%c(3)
               print*, "g2v g2a g2c",g2%v,g2%a,g2%c(1),g2%c(2),g2%c(3)
               print*, "g12v g12a g12c",g12%v,g12%a,g12%c(1),g12%c(2),
     &                                                       g12%c(3)
               print*, ""

               ! Create child if overlap volume is not zero
               if (gvol > min_gvol) then
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
     &                                          result(ret)
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

         function compute_overlap_tree_r(this, pos, radius, volume,
     &                    gamma, ishydrogen) result(ret)
            class(GOverlap_Tree), intent(inout) :: this
            real*8, dimension(:,:), intent(in) :: pos
            real*8, dimension(:), intent(in) :: radius, volume, gamma
            integer, dimension(:), intent(in) :: ishydrogen
            integer :: ret, slot
            ! Initialize overlap tree
            call this%init_overlap_tree(pos, radius, volume, gamma, 
     &                              ishydrogen)

            ! Compute and add children overlaps recursively for each atom
            do slot = 1, this%natoms
               ret = this%compute_andadd_children_r(slot)
            end do

            ret = 1
         end function compute_overlap_tree_r

         ! Compute volumes, energy of this volume and calls itself to get the volumes of the children
         recursive function compute_volume_underslot2_r(this, slot, 
     &        psi1i, f1i, p1i, psip1i, fp1i, pp1i, energy1i, fenergy1i, 
     &        penergy1i, dr, dv, free_volume, self_volume) result(ret)
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
     &                 + ov%children_count - 1
                  ret = compute_volume_underslot2_r(this, sloti, psi1it, 
     &                    f1it, p1it, psip1it, fp1it, pp1it, energy1it, 
     &                    fenergy1it, penergy1it, dr, dv, free_volume, 
     &                    self_volume)

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
     &                       + (-ov%dv1) * fenergy1i + penergy1i * c2
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
     &                     free_volume, self_volume) result(ret)
            class(GOverlap_Tree), intent(inout) :: this
            real*8, intent(in) :: pos(:,:)
            real*8, intent(out) :: volume, energy
            real*8, dimension(:,:), intent(inout) :: dr
            real*8, dimension(:), intent(inout) :: dv, free_volume
            real*8, dimension(:), intent(inout) :: self_volume
            integer :: ret, slot, i, j
            real*8 :: psi1i, f1i, psip1i, fp1i, energy1i, fenergy1i
            real*8, dimension(3) :: p1i, pp1i, penergy1i

            slot = 0

            ! Reset volumes, gradients
            do i = 1, this%natoms
               dr(1,i) = 0.0d0
               dr(2,i) = 0.0d0
               dr(3,i) = 0.0d0
               dv(i) = 0.0d0
               free_volume(i) = 0.0d0
               self_volume(i) = 0.0d0
            end do
          
            ret = compute_volume_underslot2_r(this, slot, psi1i, f1i,
     &                    p1i, psip1i, fp1i, pp1i, energy1i, fenergy1i,
     &                    penergy1i, dr, dv, free_volume, self_volume)

            volume = psi1i
            energy = energy1i
            ret = 1
         end function compute_volume2_r

c         ! Rescan the sub-tree to recompute the volumes, does not modify the tree.
c         recursive function rescan_r(this, slot) result(ret)
c           class(GOverlap_Tree), intent(inout) :: this
c           integer, intent(in) :: slot
c           integer :: parent_index, sibling_start, sibling_count
c           integer :: atom, slot_child, ret
c           type(GOverlap) :: ov, ov2, parent
c           type(GaussianVca) :: g1, g2, g12
c           real*8 :: dVdr, dVdV, dVdalpha, d2Vdalphadr, d2VdVdr
c           real*8 :: sfp, gvol

c           ov = this%overlaps%get_element(slot)

c           parent_index = ov%parent_index
c           if (parent_index > 0) then
c             atom = ov%atom
c             parent = this%overlaps%get_element(parent_index)
c             g1 = parent%g
c             ov2 = this%overlaps%get_element(atom)
c             g2 = ov2%g
c             gvol = ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp)
c             ov%g = g12
c             ov%volume = gvol
c             ov%dv1 = (g2%c - g1%c) * (-dVdr)
c             ov%dvv1 = dVdV
c             ov%sfp = sfp
c             ov%gamma1i = parent%gamma1i + ov2%gamma1i
c             ! Write back the changes.
c             call this%overlaps%set_element(slot, ov)
c           end if

c           do slot_child = ov%children_startindex, 
c      &       ov%children_startindex + ov%children_count - 1
c              ret = rescan_r(this, slot_child)
c           end do

c           ret = 1
c         end function rescan_r

c         ! Rescan the tree to recompute the volumes, does not modify the tree.
c         recursive function rescan_tree_v(this, pos, radius, volume, 
c      &      gamma, ishydrogen) result(ret)
c           class(GOverlap_Tree), intent(inout) :: this
c           real*8, intent(in) :: pos(:,:)
c           real*8, intent(in) :: radius(:)
c           real*8, intent(in) :: volume(:)
c           real*8, intent(in) :: gamma(:)
c           integer, intent(in) :: ishydrogen(:)
c           integer :: iat
c           type(GOverlap) :: ov
c           real*8 :: a, vol, ret

c           ov = this%overlaps%get_element(0)
c           ov%level = 0
c           ov%volume = 0
c           ov%dv1 = [0.0d0, 0.0d0, 0.0d0]
c           ov%dvv1 = 0.0d0
c           ov%self_volume = 0
c           ov%sfp = 1.0d0
c           ov%gamma1i = 0.0d0
c           ! Write back updates.
c           call this%overlaps%set_element(0, ov)

c           do iat = 1, this%natoms
c             ov = this%overlaps%get_element(iat)
c             a = kfc / (radius(iat) * radius(iat))
c             vol = merge(0.0d0, volume(iat), ishydrogen(iat) > 0)
c             ov%level = 1
c             ov%g%v = vol
c             ov%g%a = a
c             ov%g%c = pos(:, iat)
c             ov%volume = vol
c             ov%dv1 = [0.0d0, 0.0d0, 0.0d0]
c             ov%dvv1 = 1.0d0
c             ov%self_volume = 0.0d0
c             ov%sfp = 1.0d0
c             ov%gamma1i = gamma(iat)
c             ! Write back updates.
c             call this%overlaps%set_element(iat, ov)
c           end do

c           ret = this%rescan_r(0)

c           ret = 1
c         end function rescan_tree_v

c         ! Rescan the sub-tree to recompute the gammas, does not modify
c         ! the volumes nor the tree.
c         recursive function rescan_gamma_r(this, slot) result(res)
c           class(GOverlap_Tree), intent(inout) :: this
c           integer, intent(in) :: slot
c           integer :: parent_index, sibling_start, sibling_count
c           logical :: res
c           type(GOverlap) :: ov, ov2, parent
c           integer :: atom, slot_child

c           ! This overlap
c           ov = this%overlaps%get_element(slot)

c           ! Recompute its own overlap by merging parent and last atom
c           parent_index = ov%parent_index
c           if (parent_index > 0) then
c             atom = ov%atom
c             parent = this%overlaps%get_element(parent_index)
c             ov2 = this%overlaps%get_element(atom)
c             ov%gamma1i = parent%gamma1i + ov2%gamma1i
c             ! Write back updates.
c             call this%overlaps%set_element(slot, ov)
c           endif

c           ! Calls itself recursively on the children
c           do slot_child = ov%children_startindex, 
c      &        ov%children_startindex + ov%children_count-1
c             res = rescan_gamma_r(this, slot_child)
c           end do

c           res = .true.
c         end function rescan_gamma_r

c         ! Rescan the tree to recompute the gammas only, does not modify
c         ! volumes and the tree.
c         function rescan_tree_g(this, gamma) result(res)
c           class(GOverlap_Tree), intent(inout) :: this
c           real*8, intent(in) :: gamma(:)
c           logical :: res
c           integer :: iat
c           type(GOverlap) :: ov

c           ov = this%overlaps%get_element(0)
c           ov%gamma1i = 0.0d0
c           call this%overlaps%set_element(0, ov)

c           do iat = 1, this%natoms
c             ov = this%overlaps%get_element(iat)
c             ov%gamma1i = gamma(iat)
c             call this%overlaps%set_element(iat, ov)
c           end do

c           res = this%rescan_gamma_r(0)

c           res = .true.
c         end function rescan_tree_g

c         subroutine print_tree(this)
c           class(GOverlap_Tree), intent(in) :: this
c           integer :: i
c           write(*, '(A)', advance='no') "      slot  level   Atom parent
c      & ChStart ChCount          V      gamma          a          v
c      &    x          y          z       dedx       dedy       dedz
c      &  sfp"
c           write(*, '(A)', advance='no') new_line('A')
  
c           call print_tree_r(this, 0)
c !          do i = 1, this%natoms
c !            call print_tree_r(this, i)
c !          end do

c         end subroutine print_tree

c         recursive subroutine print_tree_r(this, slot)
c           class(GOverlap_Tree), intent(in) :: this
c           integer, intent(in) :: slot
c           type(GOverlap) :: ov
c           integer :: i

c           ov = this%overlaps%get_element(slot)

c           write(*, '(A,I6,A)', advance='yes') "tg: ", slot, " "
c           call ov%print_overlap()

c           do i = ov%children_startindex, 
c      &           ov%children_startindex + ov%children_count - 1
c             call print_tree_r(this, i)
c           end do
c         end subroutine print_tree_r

c         recursive function nchildren_under_slot_r(this, slot) result(n)
c           class(GOverlap_Tree), intent(in) :: this
c           integer, intent(in) :: slot
c           integer :: i, n
c           type(GOverlap) :: ov

c           n = 0
c           ov = this%overlaps%get_element(slot)
c           if (ov%children_count > 0) then
c             n = n + ov%children_count
c             ! now calls itself on the children
c             do i = 0, ov%children_count - 1
c               n = n + this%nchildren_under_slot_r( 
c      &                  ov%children_startindex + i)
c             end do
c           end if

c           return
c         end function nchildren_under_slot_r

      end module goverlaptreemodule
