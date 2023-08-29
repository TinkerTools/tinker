c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2023 by MKJ Chung, MJ Schnieders & JW Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  module goverlaptreemodule  --  GaussVol tree class  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     overlaps   vector that stores gaussian overlaps
c     natoms     number of atoms in the system
c
c
      module goverlaptreemodule
      use gaussianvcamodule
      use gaussvolconst
      use vectormodule
      implicit none
c
c
c     tree data structure for storing and computing GaussVol
c
      type GOverlap_Tree
      type(Vector) overlaps
      integer natoms
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
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine create_GOverlap_Tree  --  create overlap tree  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "create_GOverlap_Tree" creates the vectors that stores
c     the gaussian overlap trees
c
c
      subroutine create_GOverlap_Tree(this, natoms)
      class(GOverlap_Tree) this
      integer natoms
c
c
c     create vector that stores the gaussian overlap tree
c
      this%natoms = natoms
      call this%overlaps%create_vector()
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine init_overlap_tree  --  initialize overlap tree  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "init_overlap_tree" initializes the overlap tree by appending
c     atomic gaussian volumes to the tree
c
c
      subroutine init_overlap_tree(this, x, y, z, radius, volume, gamma,
     &                             use)
      class(GOverlap_Tree) this
      type(GOverlap), pointer :: overlap
      integer :: iat
      real*8 :: a, vol
      real*8 gamma(*)
      real*8 radius(*)
      real*8 volume(*)
      real*8 x(*)
      real*8 y(*)
      real*8 z(*)
      logical use(*)
c
c
c     reset tree
c
      call this%overlaps%clear()
c
c     slot 0 contains the master tree information
c     children = all of the atoms
c
      allocate(overlap)
      overlap%level = 0
      overlap%volume = 0.0d0
      overlap%dv1(1) = 0.0d0
      overlap%dv1(2) = 0.0d0
      overlap%dv1(3) = 0.0d0
      overlap%dvv1 = 0.0d0
      overlap%self_volume = 0.0d0
      overlap%sfp = 1.0d0
      overlap%gamma1i = 0.0d0
      overlap%parent_index = -1
      overlap%atom = -1
      overlap%children_startindex = 1
      overlap%children_count = this%natoms
c
c     add the root overlap to the vector
c
      call this%overlaps%push_back(overlap)
      deallocate(overlap)
c
c     add all atoms to the tree, where list of atoms start at slot 1
c
      do iat = 1, this%natoms
         allocate(overlap)
         overlap%level = 1
         a = kfc / (radius(iat) * radius(iat))
         vol = volume(iat)
         if (.not. use(iat)) vol = 0.0d0
         overlap%g%v = vol
         overlap%g%a = a
         overlap%g%c(1) = x(iat)
         overlap%g%c(2) = y(iat)
         overlap%g%c(3) = z(iat)
         overlap%volume = vol
         overlap%dv1(1) = 0.0d0
         overlap%dv1(2) = 0.0d0
         overlap%dv1(3) = 0.0d0
         overlap%dvv1 = 1.0d0
         overlap%self_volume = 0.0d0
         overlap%sfp = 1.0d0
         overlap%gamma1i = gamma(iat)
         overlap%parent_index = 0
         overlap%atom = iat
         overlap%children_startindex = -1
         overlap%children_count = -1
c
c     add the atom to the vector
c
         call this%overlaps%push_back(overlap)
         deallocate(overlap)
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function add_children  --  add child gaussian overlap  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "add_children" adds to the tree the children of overlap
c     identified by "parent_index" in the tree.
c
c
      function add_children(this, parent_index, children_overlaps)
     &                                  result(start_index)
      class(GOverlap_Tree) this
      type(GOverlap) parent, child_overlap
      type(Vector) children_overlaps
      integer i, ip, noverlaps
      integer parent_index, parent_level
      integer slot, start_index
c
c
c     adds children starting at the last slot
c
      start_index = this%overlaps%vector_size()
      noverlaps = children_overlaps%vector_size()
c
c     retrieves address of root overlap
c
      parent = this%overlaps%get_element(parent_index)
c
c     registers list of children
c
      parent%children_startindex = start_index
      parent%children_count = noverlaps
c
c     write back the changes
c
      call this%overlaps%set_element(parent_index, parent)
c
c     sort neighbors by overlap volume
c
      call sort(children_overlaps)
      parent_level = parent%level
c
c     now copies the children overlaps from temp buffer
c
      do ip = 0, noverlaps - 1
         child_overlap = children_overlaps%get_element(ip)
         child_overlap%level = parent_level + 1
c
c     connect overlap to parent
c
         child_overlap%parent_index = parent_index
c
c     reset its children indexes
c
         child_overlap%children_startindex = -1
         child_overlap%children_count = -1
c
c     add to tree
c
         call this%overlaps%push_back(child_overlap)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function compute_children  --  computes children overlap  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "compute_children" scans the siblings of overlap identified by
c     "root_index" to create children overlaps, returns them into the
c     "children_overlaps" buffer: (root) + (atom) -> (root, atom).
c
c
      function compute_children(this, root_index, children_overlaps)
     &                                    result(ret)
      class(GOverlap_Tree) this
      type(GaussianVca) g1, g2, g12
      type(GOverlap) root, parent, sibling, ov2
      type(GOverlap), pointer :: ov
      type(Vector) children_overlaps
      type(Vector) overlaps
      integer root_index
      integer parent_index, sibling_start, sibling_count
      integer j, atom2, ret
      real*8 gvol, dVdr, dVdV, sfp
c
c
c     reset output buffer
c
      call children_overlaps%clear()
c
c     retrieves overlap
c
      overlaps = this%overlaps
      root = overlaps%get_element(root_index)
c
c     retrieves parent overlap
c
      parent_index = root%parent_index
c
c     master root, can't compute children
c
      if (parent_index < 0) then
         write(*,*) "Error -- invalid parent index"
         ret = 1 ! 
         return
      end if
c
c     max_order reached, don't compute further
c
      if (root%level >= max_order) then
         write(*,*) "Reached max level", root%level
         ret = 1
         return
      end if
c
c     retrieves start index and count of siblings
c
      parent = overlaps%get_element(parent_index)
      sibling_start = parent%children_startindex
      sibling_count = parent%children_count
c
c     parent is not initialized
c
      if (sibling_start < 0 .or. sibling_count < 0) then
         write(*,*) "Error -- parent is not initialized"
         ret = -1
         return
      end if
c
c     check that this overlap is the child of the registered parent
c
      if (root_index < sibling_start .and. root_index >
     &             sibling_start + sibling_count - 1) then
         write(*,*) "Error -- overlap parent is not registered"
         ret = -1 
         return
      end if
c
c     loop over "younger" siblings (i < j loop) to compute new overlaps
c
      do j = root_index + 1, sibling_start + sibling_count - 1
         sibling = overlaps%get_element(j)
c
c     atomic gaussian of the last atom of the sibling
c
         atom2 = sibling%atom
         g1 = root%g
         ov2 = overlaps%get_element(atom2)
         g2 = ov2%g
         call ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp, gvol)
c
c     create child if overlap volume is not zero
c
         if (gvol > min_gvol) then
            allocate(ov)
            ov%g = g12
            ov%volume = gvol
            ov%self_volume = 0
            ov%atom = atom2
c
c     dv1 is the gradient of V(123..)n with respect to the position of 1
c
            ov%dv1 = (g2%c - g1%c) * (-dVdr)
c
c     dvv1 is the derivative of V(123...)n with respect to V(123...)
c
            ov%dvv1 = dVdV
            ov%sfp = sfp
            ov%gamma1i = root%gamma1i + ov2%gamma1i
            call children_overlaps%push_back(ov)
            deallocate(ov)
         end if
      end do
      ret = 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function compute_andadd_children_r  --  recursive children  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "compute_andadd_children_r" recursively adds children if
c     noverlap > 0
c
c
      recursive function compute_andadd_children_r(this, root)
     &                                          result(ret)
      class(GOverlap_Tree) this
      type(GOverlap) child, parent
      type(Vector) children_overlaps
      integer ret
      integer root
      integer noverlaps, start_slot, ichild
c
c
c     compute children overlaps
c
      call children_overlaps%create_vector()
      ret = this%compute_children(root, children_overlaps)
c
c     get the number of overlaps
c
      noverlaps = children_overlaps%vector_size()
c
c     add children overlaps and get the starting slot
c
      if (noverlaps > 0) then
         start_slot = this%add_children(root, children_overlaps)
c
c     recursively compute and add children for each child overlap
c
         do ichild = start_slot, start_slot + noverlaps - 1
            ret = this%compute_andadd_children_r(ichild)
         end do
      end if
      ret = 1
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function compute_overlap_tree_r  --  recursive overlap  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "compute_overlap_tree_r" generates an overlap tree structure
c     object by recursively computing atomic sphere overlaps
c
c
      function compute_overlap_tree_r(this, x, y, z, radius, volume,
     &                                gamma, use) result(ret)
      class(GOverlap_Tree) this
      integer ret, slot
      real*8 x(*)
      real*8 y(*)
      real*8 z(*)
      real*8 radius(*)
      real*8 volume(*)
      real*8 gamma(*)
      logical use(*)
c
c
c     initialize overlap tree
c
      call this%init_overlap_tree(x, y, z, radius,volume,gamma,use)
c
c     compute and add children overlaps recursively for each atom
c
      do slot = 1, this%natoms
         ret = this%compute_andadd_children_r(slot)
      end do
      ret = 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function compute_volume_underslot2_r  --  recursive volume  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "compute_volume_underslot2_r" compute volumes, energy of this
c     volume and calls itself to get the volumes of the children
c
c
      recursive function compute_volume_underslot2_r(this, slot, 
     &        psi1i, f1i, p1i, psip1i, fp1i, pp1i, energy1i, fenergy1i, 
     &        penergy1i, dr, dv, free_volume, self_volume) result(ret)
      class(GOverlap_Tree) this
      type(GOverlap) ov, ov2
      integer ret
      integer slot
      integer atom, i, j, sloti
      real*8 psi1i, f1i, psip1i, fp1i
      real*8 energy1i, fenergy1i
      real*8 cf, volcoeff, volcoeffp
      real*8 ai, a1i, a1, c1, c1p, c2
      real*8 psi1it, f1it, psip1it, fp1it
      real*8 energy1it, fenergy1it
      real*8 p1i(3), pp1i(3)
      real*8 penergy1i(3)
      real*8 p1it(3), pp1it(3), penergy1it(3)
      real*8 dr(3,*)
      real*8 dv(*)
      real*8 free_volume(*)
      real*8 self_volume(*)
c
c
c     determine if level is even or odd
c
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
c
c     evaluate intermediate terms for volume calculations
c
      atom = ov%atom
      ov2 = this%overlaps%get_element(atom)
      ai = ov2%g%a
      a1i = ov%g%a
      a1 = a1i - ai
      psi1i = volcoeff * ov%volume
      f1i = volcoeff * ov%sfp
      p1i(1) = 0.0d0
      p1i(2) = 0.0d0
      p1i(3) = 0.0d0
      psip1i = volcoeffp * ov%volume
      fp1i = volcoeffp * ov%sfp
      pp1i(1) = 0.0d0
      pp1i(2) = 0.0d0
      pp1i(3) = 0.0d0
      energy1i = volcoeffp * ov%gamma1i * ov%volume
      fenergy1i = volcoeffp * ov%sfp * ov%gamma1i
      penergy1i(1) = 0.0d0
      penergy1i(2) = 0.0d0
      penergy1i(3) = 0.0d0
c
c     add up volume of children
c
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
c
c     contributions to free and self volume of last atom
c
      if (ov%level > 0) then
         free_volume(atom) = free_volume(atom) + psi1i
         self_volume(atom) = self_volume(atom) + psip1i
c
c     contributions to energy gradients
c
         c2 = ai / a1i
         dr(1,atom) = dr(1,atom) + (-ov%dv1(1)) * fenergy1i
     &                     + penergy1i(1) * c2
         dr(2,atom) = dr(2,atom) + (-ov%dv1(2)) * fenergy1i
     &                     + penergy1i(2) * c2
         dr(3,atom) = dr(3,atom) + (-ov%dv1(3)) * fenergy1i
     &                     + penergy1i(3) * c2
         dv(atom) = dv(atom) + ov%g%v * fenergy1i
c
c     update subtree P1..i's for parent
c
         c2 = a1 / a1i
         p1i = (ov%dv1) * f1i + p1i * c2
         pp1i = (ov%dv1) * fp1i + pp1i * c2
         penergy1i = (ov%dv1) * fenergy1i + penergy1i * c2
c
c     update subtree F1..i's for parent
c
         f1i = ov%dvv1 * f1i
         fp1i = ov%dvv1 * fp1i
         fenergy1i = ov%dvv1 * fenergy1i
      end if
      ret = 1
      return
      end 
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function compute_volume2_r  --  traverse and compute  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "compute_volume2_r" traverses tree and computes volumes
c
c
      function compute_volume2_r(this, x, y, z, volume, energy, dr, dv,
     &                     free_volume, self_volume) result(ret)
      class(GOverlap_Tree) this
      integer ret, slot, i, j
      real*8 volume, energy
      real*8 psi1i, f1i, psip1i, fp1i, energy1i, fenergy1i
      real*8 p1i(3), pp1i(3), penergy1i(3)
      real*8 dv(*)
      real*8 free_volume(*)
      real*8 self_volume(*)
      real*8 x(*)
      real*8 y(*)
      real*8 z(*)
      real*8 dr(3,*)
c
c
c     reset volumes, gradients
c
      do i = 1, this%natoms
         dr(1,i) = 0.0d0
         dr(2,i) = 0.0d0
         dr(3,i) = 0.0d0
         dv(i) = 0.0d0
         free_volume(i) = 0.0d0
         self_volume(i) = 0.0d0
      end do
c
c     traverse tree and compute volume
c
      slot = 0
      ret = compute_volume_underslot2_r(this, slot, psi1i, f1i,
     &                    p1i, psip1i, fp1i, pp1i, energy1i, fenergy1i,
     &                    penergy1i, dr, dv, free_volume, self_volume)

      volume = psi1i
      energy = energy1i
      ret = 1
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  function rescan_r  --  rescan sub-tree and compute volume  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rescan_r" rescan the sub-tree to recompute the volumes,
c     does not modify the tree
c
c
      recursive function rescan_r(this, slot) result(ret)
      class(GOverlap_Tree) this
      type(GaussianVca) g1, g2, g12
      type(GOverlap) ov, ov2, parent
      integer slot
      integer parent_index, sibling_start, sibling_count
      integer atom, slot_child, ret
      real*8 dVdr, dVdV, dVdalpha, d2Vdalphadr, d2VdVdr
      real*8 sfp, gvol
c
c
c     choose this overlap
c
      ov = this%overlaps%get_element(slot)
c
c     recompute its own overlap by merging parent and last atom
c
      parent_index = ov%parent_index
      if (parent_index > 0) then
         atom = ov%atom
         parent = this%overlaps%get_element(parent_index)
         g1 = parent%g
         ov2 = this%overlaps%get_element(atom)
         g2 = ov2%g
         call ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp, gvol)
         ov%g = g12
         ov%volume = gvol
         ov%dv1 = (g2%c - g1%c) * (-dVdr)
         ov%dvv1 = dVdV
         ov%sfp = sfp
         ov%gamma1i = parent%gamma1i + ov2%gamma1i
c
c     write back the changes
c
         call this%overlaps%set_element(slot, ov)
      end if
c
c     calls itself recursively on the children
c
      do slot_child = ov%children_startindex, 
     &               ov%children_startindex + ov%children_count - 1
         ret = rescan_r(this, slot_child)
      end do
      ret = 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function rescan_tree_v  --  rescan tree and compute volume  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rescan_tree_v" rescan the tree to recompute the volumes,
c     does not modify the tree
c
c
      recursive function rescan_tree_v(this, x, y, z, radius, volume, 
     &                                gamma, use) result(ret)
      class(GOverlap_Tree) this
      type(GOverlap) :: ov
      integer iat
      real*8 a, vol, ret
      real*8 radius(*)
      real*8 volume(*)
      real*8 gamma(*)
      real*8 x(*)
      real*8 y(*)
      real*8 z(*)
      logical use(*)
c
c
c     initialize root tree
c
      ov = this%overlaps%get_element(0)
      ov%level = 0
      ov%volume = 0
      ov%dv1(1) = 0.0d0
      ov%dv1(2) = 0.0d0
      ov%dv1(3) = 0.0d0
      ov%dvv1 = 0.0d0
      ov%self_volume = 0
      ov%sfp = 1.0d0
      ov%gamma1i = 0.0d0
c
c     write back updates
c
      call this%overlaps%set_element(0, ov)
c
c     initialize first level of tree
c
      do iat = 1, this%natoms
         ov = this%overlaps%get_element(iat)
         a = kfc / (radius(iat) * radius(iat))
         vol = volume(iat)
         if (.not. use(iat)) vol = 0.0d0
         ov%level = 1
         ov%g%v = vol
         ov%g%a = a
         ov%g%c(1) = x(iat)
         ov%g%c(2) = y(iat)
         ov%g%c(3) = z(iat)
         ov%volume = vol
         ov%dv1(1) = 0.0d0
         ov%dv1(2) = 0.0d0
         ov%dv1(3) = 0.0d0
         ov%dvv1 = 1.0d0
         ov%self_volume = 0.0d0
         ov%sfp = 1.0d0
         ov%gamma1i = gamma(iat)
c
c     write back updates
c
         call this%overlaps%set_element(iat, ov)
      end do
      ret = this%rescan_r(0)
      ret = 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function rescan_gamma_r  --  rescan sub-tree compute gamma  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rescan_gamma_r" rescan the sub-tree to recompute the gammas,
c     does not modify volumes nor the tree
c
c
      recursive function rescan_gamma_r(this, slot) result(ret)
      class(GOverlap_Tree) this
      type(GOverlap) ov, ov2, parent
      integer slot
      integer parent_index, sibling_start, sibling_count
      integer ret
      integer atom, slot_child
c
c
c     choose this overlap
c
      ov = this%overlaps%get_element(slot)
c
c     recompute its own overlap by merging parent and last atom
c
      parent_index = ov%parent_index
      if (parent_index > 0) then
         atom = ov%atom
         parent = this%overlaps%get_element(parent_index)
         ov2 = this%overlaps%get_element(atom)
         ov%gamma1i = parent%gamma1i + ov2%gamma1i
c
c     write back the changes
c
         call this%overlaps%set_element(slot, ov)
      endif
c
c     calls itself recursively on the children
c
      do slot_child = ov%children_startindex, 
     &        ov%children_startindex + ov%children_count-1
         ret = rescan_gamma_r(this, slot_child)
      end do

      ret = 1
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function rescan_tree_g  --  rescan tree and compute gammas  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rescan_tree_g" rescan the tree to recompute the gammas,
c     does not modify volumes nor the tree
c
c
      function rescan_tree_g(this, gamma) result(ret)
      class(GOverlap_Tree) this
      type(GOverlap) :: ov
      integer ret, iat
      real*8 gamma(*)
c
c
c     initialize root of tree
c
      ov = this%overlaps%get_element(0)
      ov%gamma1i = 0.0d0
      call this%overlaps%set_element(0, ov)
c
c     initialize first level of tree
c
      do iat = 1, this%natoms
         ov = this%overlaps%get_element(iat)
         ov%gamma1i = gamma(iat)
         call this%overlaps%set_element(iat, ov)
      end do
      ret = this%rescan_gamma_r(0)
      ret = 1
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine print_tree  --  print GaussVol tree  ##
c     ##                                                  ##
c     ######################################################
c
c
c     "print_tree" prints the GaussVol tree
c
c
      subroutine print_tree(this)
      class(GOverlap_Tree) this
      integer i
c
c
c     print header
c
      write(*, '(A)', advance='no') "      slot  level   Atom parent"
      write(*, '(A)', advance='no') " ChStart ChCount          V"
      write(*, '(A)', advance='no') "      gamma          a          v"
      write(*, '(A)', advance='no') "          x          y          z"
      write(*, '(A)', advance='no') "       dedx       dedy       dedz"
      write(*, '(A)', advance='no') "        sfp"
c      write(*, '(A)', advance='no') "  level          V"
      write(*, '(A)', advance='no') new_line('A')
c
c     recursive call to print tree
c
      do i = 1, this%natoms
         call print_tree_r(this, i)
      end do
c      call print_tree_r(this, 0)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine print_tree_r  --  recursive print GaussVol tree  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "print_tree_r" recursively prints the GaussVol tree
c
c
      recursive subroutine print_tree_r(this, slot)
      class(GOverlap_Tree) this
      type(GOverlap) ov
      integer slot
      integer i
      ov = this%overlaps%get_element(slot)

      write(*, '(A,I6,A)', advance='yes') "tg: ", slot, " "
      call ov%print_overlap()
      do i = ov%children_startindex, 
     &           ov%children_startindex + ov%children_count - 1
         call print_tree_r(this, i)
      end do
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine nchildren_under_slot_r  --  recursive nchildren  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "nchildren_under_slot_r" recursively calculates total number of
c     children elements under tree node
c
c
      recursive function nchildren_under_slot_r(this, slot) result(n)
      class(GOverlap_Tree) this
      type(GOverlap) ov
      integer slot
      integer i, n
      n = 0
      ov = this%overlaps%get_element(slot)
      if (ov%children_count > 0) then
         n = n + ov%children_count
         do i = 0, ov%children_count - 1
            n = n + this%nchildren_under_slot_r( 
     &                         ov%children_startindex + i)
         end do
      end if
      return
      end
      end module goverlaptreemodule
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine pol_switchfunc  --  overlap switching function  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "pol_switchfunc" overlap volume switching function
c     + 1st derivative 
c
c
      subroutine pol_switchfunc(gvol, volmina, volminb, sp, s)
      implicit none
      real*8 gvol, volmina, volminb, sp, s
      real*8 swf, swfp, swd, swu, swu2, swu3
c
c
c     set swf and swfp
c
      swf = 0.0d0
      swfp = 1.0d0
      if (gvol > volminb) then
         swf = 1.0d0
         swfp = 0.0d0
      else if (gvol < volmina) then
         swf = 0.0d0
         swfp = 0.0d0
      end if
c
c     compute switching function
c
      swd = 1.0d0 / (volminb - volmina)
      swu = (gvol - volmina) * swd
      swu2 = swu * swu
      swu3 = swu * swu2
      s = swf + swfp * swu3 * (10.0d0 - 15.0d0 * swu + 6.0d0 * swu2)
      sp = swfp * swd * 30.0d0 * swu2 * (1.0d0 - 2.0d0 * swu + swu2)
c
c     turn off switching function
c
c      sp = 0.0d0
c      s = 1.0d0
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ogauss_alpha  --  overlap of two gaussians  ##
c     ##                                                         ##
c     #############################################################
c
c
c     ogauss_alpha computes the overlap between two Gaussians
c     represented by a (V,c,a) triplet:
c        V: volume of Gaussian
c        c: position of Gaussian
c        a: exponential coefficient
c        g(x) = V (a/pi)^(3/2) exp(-a(x-c)^2)
c     this version is based on V=V(V1,V2,r1,r2,alpha)
c        alpha = (a1 + a2)/(a1 a2)
c        dVdr is (1/r)*(dV12/dr)
c        dVdV is dV12/dV1 
c        dVdalpha is dV12/dalpha
c        d2Vdalphadr is (1/r)*d^2V12/dalpha dr
c        d2VdVdr is (1/r) d^2V12/dV1 dr
c
c
      subroutine ogauss_alpha(g1, g2, g12, dVdr, dVdV, sfp, sgvol)
      use gaussianvcamodule
      use gaussvolconst
      use math
      implicit none
      type(GaussianVca) g1, g2, g12
      real*8 dVdr, dVdV, sfp, sgvol
      real*8 d2, deltai, gvol, p12, a12
      real*8 s, sp, df, dgvol, dgvolv, ef, dgvola2, dgvola1
      real*8 dgalpha, dgalpha2, dgvolvdr
      real*8 c1(3), c2(3), dist(3)
c
c
c     get center two of gaussians
c
      c1 = g1%c
      c2 = g2%c
      dist = c2 - c1
      d2 = dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3)
c
c     compute gaussian and derivative
c
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
c
c     set new gaussian formed from two gaussians
c
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
