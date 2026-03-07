c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine findpro  --  search for amino acid residues  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "findpro" locates and stores the atoms in amino acid residues
c     based on atomic element and connectivity information
c
c     note this routine assumes an all-atom model with hydrogen
c     atoms explicitly represented
c
c
      subroutine findpro
      use atomid
      use atoms
      use bitor
      use couple
      implicit none
      integer i,j,nhyd
      integer ia,ib,ic
      integer id,ie,ij
      integer nia,nib,nic
      integer nid,nie,nij
      integer aia,aib,aic
      integer aid,aie,aij
      integer icb1,icb2
      integer icg1,icg2
      integer icd1,icd2
      integer ice1,ice2
      integer icz1,icz2
      integer ich,ind
      integer ine,inz
      integer inh1,inh2
      integer iog,ioh
      integer iod1,iod2
      integer ioe1,ioe2
      integer isg,isd
      integer nha,nhb
      integer nhg,nhd
      logical proceed
      character*3 label
c
c
c     initialize the tier name associated with each atom
c
      do i = 1, n
         tier(i) = '   '
      end do
c
c     cycle over bitorsions to search for amino acids
c
      do i = 1, nbitor
         ia = ibitor(1,i)
         ib = ibitor(2,i)
         ic = ibitor(3,i)
         id = ibitor(4,i)
         ie = ibitor(5,i)
         nia = n12(ia)
         nib = n12(ib)
         nic = n12(ic)
         nid = n12(id)
         nie = n12(ie)
         aia = atomic(ia)
         aib = atomic(ib)
         aic = atomic(ic)
         aid = atomic(id)
         aie = atomic(ie)
c
c     check for five-atom stretch of polypeptide backbone
c
         proceed = .false.
         if (nia.eq.3 .and. nib.eq.3 .and. nic.eq.4
     &          .and. nid.eq.3 .and. nie.eq.3) then
            if ((aia.eq.6 .and. aib.eq.7 .and. aic.eq.6
     &              .and. aid.eq.6 .and. aie.eq.7) .or.
     &          (aia.eq.7 .and. aib.eq.6 .and. aic.eq.6
     &              .and. aid.eq.7 .and. aie.eq.6)) then
               if (.not.proceed .and. aia.eq.6) then
                  do j = 1, n12(ia)
                     ij = i12(j,ia)
                     nij = n12(ij)
                     aij = atomic(ij)
                     if (nij.eq.1 .and. aij.eq.8)  proceed = .true.
                  end do 
               end if
               if (.not.proceed .and. aib.eq.6) then
                  do j = 1, n12(ib)
                     ij = i12(j,ib)
                     nij = n12(ij)
                     aij = atomic(ij)
                     if (nij.eq.1 .and. aij.eq.8)  proceed = .true.
                  end do 
               end if
               if (.not.proceed .and. aid.eq.6) then
                  do j = 1, n12(id)
                     ij = i12(j,id)
                     nij = n12(ij)
                     aij = atomic(ij)
                     if (nij.eq.1 .and. aij.eq.8)  proceed = .true.
                  end do 
               end if
               if (.not.proceed .and. aie.eq.6) then
                  do j = 1, n12(ie)
                     ij = i12(j,ie)
                     nij = n12(ij)
                     aij = atomic(ij)
                     if (nij.eq.1 .and. aij.eq.8)  proceed = .true.
                  end do 
               end if
            end if
         end if
c
c     check for N-terminal and C-terminal residues
c
         if (.not. proceed) then
            if (nia.eq.1 .and. nib.eq.4 .and. nic.eq.4
     &             .and. nid.eq.3 .and. nie.eq.3 .and.
     &          aia.eq.1 .and. aib.eq.7 .and. aic.eq.6
     &             .and. aid.eq.6 .and. aie.eq.7)  proceed = .true.
         end if
         if (.not. proceed) then
            if (nia.eq.3 .and. nib.eq.3 .and. nic.eq.4
     &             .and. nid.eq.4 .and. nie.eq.1 .and.
     &          aia.eq.7 .and. aib.eq.6 .and. aic.eq.6
     &             .and. aid.eq.7 .and. aie.eq.1)  proceed = .true.
         end if
         if (.not. proceed) then
            if (nia.eq.3 .and. nib.eq.3 .and. nic.eq.4
     &             .and. nid.eq.3 .and. nie.eq.1 .and.
     &          aia.eq.6 .and. aib.eq.7 .and. aic.eq.6
     &              .and. aid.eq.6 .and. aie.eq.8)  proceed = .true.
         end if
         if (.not. proceed) then
            if (nia.eq.1 .and. nib.eq.3 .and. nic.eq.4
     &             .and. nid.eq.3 .and. nie.eq.3 .and.
     &          aia.eq.8 .and. aib.eq.6 .and. aic.eq.6
     &             .and. aid.eq.7 .and. aie.eq.6)  proceed = .true.
         end if
c
c     check for N-terminal ACE chain capping group
c
         if (.not. proceed) then
            if (nia.eq.4 .and. nib.eq.3 .and. nic.eq.3
     &             .and. nid.eq.4 .and. nie.eq.3 .and.
     &          aia.eq.6 .and. aib.eq.6 .and. aic.eq.7
     &             .and. aid.eq.6 .and. aie.eq.6) then
               nhyd = 0
               do j = 1, n12(ia)
                  ij = i12(j,ia)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  tier(ia) = 'ACE'
                  do j = 1, n12(ia)
                     ij = i12(j,ia)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'ACE'
                  end do
                  tier(ib) = 'ACE'
                  do j = 1, n12(ib)
                     ij = i12(j,ib)
                     aij = atomic(ij)
                     if (aij .eq. 8)  tier(ij) = 'ACE'
                  end do
               end if
            end if
         end if
         if (.not. proceed) then
            if (nia.eq.3 .and. nib.eq.4 .and. nic.eq.3
     &             .and. nid.eq.3 .and. nie.eq.4 .and.
     &          aia.eq.6 .and. aib.eq.6 .and. aic.eq.7
     &             .and. aid.eq.6 .and. aie.eq.6) then
               nhyd = 0
               do j = 1, n12(ie)
                  ij = i12(j,ie)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  tier(ia) = 'ACE'
                  do j = 1, n12(ie)
                     ij = i12(j,ie)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'ACE'
                  end do
                  tier(ib) = 'ACE'
                  do j = 1, n12(id)
                     ij = i12(j,id)
                     aij = atomic(ij)
                     if (aij .eq. 8)  tier(ij) = 'ACE'
                  end do
               end if
            end if
         end if
c
c     check for C-terminal NME chain capping group
c
         if (.not. proceed) then
            if (nia.eq.3 .and. nib.eq.4 .and. nic.eq.3
     &             .and. nid.eq.3 .and. nie.eq.4 .and.
     &          aia.eq.7 .and. aib.eq.6 .and. aic.eq.6
     &             .and. aid.eq.7 .and. aie.eq.6) then
               nhyd = 0
               do j = 1, n12(ie)
                  ij = i12(j,ie)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  tier(id) = 'NME'
                  do j = 1, n12(id)
                     ij = i12(j,id)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'NME'
                  end do
                  tier(ie) = 'NME'
                  do j = 1, n12(ie)
                     ij = i12(j,ie)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'NME'
                  end do
               end if
            end if
         end if
         if (.not. proceed) then
            if (nia.eq.4 .and. nib.eq.3 .and. nic.eq.3
     &             .and. nid.eq.4 .and. nie.eq.3 .and.
     &          aia.eq.6 .and. aib.eq.7 .and. aic.eq.6
     &             .and. aid.eq.6 .and. aie.eq.7) then
               nhyd = 0
               do j = 1, n12(ia)
                  ij = i12(j,ia)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  tier(ia) = 'NME'
                  do j = 1, n12(ia)
                     ij = i12(j,ia)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'NME'
                  end do
                  tier(ib) = 'NME'
                  do j = 1, n12(ib)
                     ij = i12(j,ib)
                     aij = atomic(ij)
                     if (aij .eq. 1)  tier(ij) = 'NME'
                  end do
               end if
            end if
         end if
         if (.not. proceed)  goto 10
c
c     zero out residue name and possible side chain atoms
c
         label = '   '
         icb1 = 0
         icb2 = 0
         icg1 = 0
         icg2 = 0
         iog = 0
         isg = 0
         icd1 = 0
         icd2 = 0
         ind = 0
         iod1 = 0
         iod2 = 0
         isd = 0
         ice1 = 0
         ice2 = 0
         ine = 0
         ioe1 = 0
         ioe2 = 0
         icz1 = 0
         icz2 = 0
         inz = 0
         ich = 0
         inh1 = 0
         inh2 = 0
         ioh = 0
c
c     inspect the beta position of amino acid residue
c
         if (label .eq. '   ') then
            nha = 0
            do j = 1, n12(ic)
               ij = i12(j,ic)
               aij = atomic(ij)
               if (aij .eq. 1)  nha = nha + 1
               if (ij.ne.ib .and. ij.ne.id .and. aij.eq.6) then
                  if (icb1 .ne. 0)  icb2 = ij
                  if (icb1 .eq. 0)  icb1 = ij
               end if
            end do
            if (nha .eq. 2)  label = 'GLY'
            if (icb2 .ne. 0)  label = 'AIB'
         end if
c
c     inspect the gamma position of amino acid residue
c
         if (label.eq.'   ' .and. icb1.ne.0) then
            nhb = 0
            do j = 1, n12(icb1)
               ij = i12(j,icb1)
               aij = atomic(ij)
               if (aij .eq. 1)  nhb = nhb + 1
               if (ij.ne.ic .and. aij.eq.6) then
                  if (icg1 .ne. 0)  icg2 = ij
                  if (icg1 .eq. 0)  icg1 = ij
               end if
               if (aij .eq. 8)  iog = ij
               if (aij .eq. 16)  isg = ij
            end do
            if (nhb .eq. 3)  label = 'ALA'
            if (iog.ne.0 .and. icg1.eq.0)  label = 'SER'
            if (iog.ne.0 .and. icg1.ne.0)  label = 'THR'
            if (isg .ne. 0) then
               if (n12(isg) .eq. 1) then
                  label = 'CYD'
               else
                  do j = 1, n12(isg)
                     ij = i12(j,isg)
                     aij = atomic(ij)
                     if (aij .eq. 1)  label = 'CYS'
                     if (aij .eq. 16)  label = 'CYX'
                  end do
               end if
            end if
            if (min(icg1,icg2) .ne. 0) then
               nhg = 0
               do j = 1, n12(icg1)
                  ij = i12(j,icg1)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhg = nhg + 1
                  if (ij.ne.icb1 .and. aij.eq.6)  icd1 = ij
               end do
               do j = 1, n12(icg2)
                  ij = i12(j,icg2)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhg = nhg + 1
                  if (ij.ne.icb1 .and. aij.eq.6)  icd1 = ij
               end do
               if (nhg .eq. 5)  label = 'ILE'
               if (nhg .eq. 6)  label = 'VAL'
            end if
         end if
c
c     inspect the delta position of amino acid residue
c
         if (label.eq.'   ' .and. icg1.ne.0) then
            do j = 1, n12(icg1)
               ij = i12(j,icg1)
               aij = atomic(ij)
               if (ij.ne.icb1 .and. aij.eq.6) then
                  if (icd1 .ne. 0)  icd2 = ij
                  if (icd1 .eq. 0)  icd1 = ij
               end if
               if (aij .eq. 7)  ind = ij
               if (aij .eq. 8) then
                  if (iod1 .ne. 0)  iod2 = ij
                  if (iod1 .eq. 0)  iod1 = ij
               end if
               if (aij .eq. 16)  isd = ij
            end do
            if (iod1.ne.0 .and. ind.ne.0)  label = 'ASN'
            if (iod2 .ne. 0) then
               if (n12(iod1).eq.1 .and. n12(iod2).eq.1) then
                  label = 'ASP'
               else
                  label = 'ASH'
               end if
            end if
            if (isd .ne. 0) then
               label = 'MET'
               do j = 1, n12(isd)
                  ij = i12(j,isd)
                  aij = atomic(ij)
                  if (ij.ne.icg1 .and. aij.eq.6)  ice1 = ij
               end do
            end if
            if (min(icd1,icd2) .ne. 0) then
               nhd = 0
               do j = 1, n12(icd1)
                  ij = i12(j,icd1)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhd = nhd + 1
               end do
               do j = 1, n12(icd2)
                  ij = i12(j,icd2)
                  aij = atomic(ij)
                  if (aij .eq. 1)  nhd = nhd + 1
               end do
               if (nhd .eq. 6)  label = 'LEU'
            end if
            if (icd1 .ne. 0) then
               do j = 1, n12(icd1)
                  ij = i12(j,icd1)
                  if (ij.eq.ib .or. ij.eq.id)  label = 'PRO'
               end do
            end if
         end if
c
c     inspect the epsilon position of amino acid residue
c
         if (label.eq.'   ' .and. icd1.ne.0) then
            do j = 1, n12(icd1)
               ij = i12(j,icd1)
               aij = atomic(ij)
               if (ij.ne.icg1 .and. aij.eq.6) then
                  if (ice1 .ne. 0)  ice2 = ij
                  if (ice1 .eq. 0)  ice1 = ij
               end if
               if (aij .eq. 7)  ine = ij
               if (aij .eq. 8) then
                  if (ioe1 .ne. 0)  ioe2 = ij
                  if (ioe1 .eq. 0)  ioe1 = ij
               end if
            end do
            if (ine .ne. 0) then
               if (n12(ine) .eq. 4)  label = 'ORN'
            end if
            if (ioe1.ne.0 .and. ine.ne.0)  label = 'GLN'
            if (ioe2 .ne. 0) then
               if (n12(ioe1).eq.1 .and. n12(ioe2).eq.1) then
                  label = 'GLU'
               else
                  label = 'GLH'
               end if
            end if
         end if
         if (label.eq.'   ' .and. icd2.ne.0) then
            do j = 1, n12(icd2)
               ij = i12(j,icd2)
               aij = atomic(ij)
               if (ij.ne.icg1 .and. aij.eq.6) then
                  if (ice1 .ne. 0)  ice2 = ij
                  if (ice1 .eq. 0)  ice1 = ij
               end if
               if (aij .eq. 7)  ine = ij
            end do
         end if
         if (label.eq.'   ' .and. ind.ne.0) then
            do j = 1, n12(ind)
               ij = i12(j,ind)
               aij = atomic(ij)
               if (ij.ne.icg1 .and. aij.eq.6) then
                  if (ice1 .ne. 0)  ice2 = ij
                  if (ice1 .eq. 0)  ice1 = ij
               end if
            end do
         end if
         if (label.eq.'   ' .and. min(ind,ine).ne.0) then
            do j = 1, n12(ine)
               ij = i12(j,ine)
               if (ij .eq. ice1) then
                  label = 'HIS'
                  if (n12(ind) .eq. 2)  label = 'HIE'
                  if (n12(ine) .eq. 2)  label = 'HID'
               end if
            end do
         end if
         if (label.eq.'   ' .and. ine.ne.0) then
            do j = 1, n12(ine)
               ij = i12(j,ine)
               if (ij .eq. ice1)  label = 'TRP'
               if (ij .eq. ice2)  label = 'TRP'
            end do
            if (ice1 .ne.  0) then
               do j = 1, n12(ice1)
                  ij = i12(j,ice1)
                  aij = atomic(ij)
                  if (ij.ne.icd1 .and. ij.ne.icd2 .and. aij.eq.6) then
                     if (icz1 .ne. 0)  icz2 = ij
                     if (icz1 .eq. 0)  icz1 = ij
                  end if
               end do
            end if
            if (ice2 .ne.  0) then
               do j = 1, n12(ice2)
                  ij = i12(j,ice2)
                  aij = atomic(ij)
                  if (ij.ne.icd1 .and. ij.ne.icd2 .and. aij.eq.6) then
                     if (icz1 .ne. 0)  icz2 = ij
                     if (icz1 .eq. 0)  icz1 = ij
                  end if
               end do
            end if
            if (icz1 .ne.  0) then
               do j = 1, n12(icz1)
                  ij = i12(j,icz1)
                  aij = atomic(ij)
                  if (ij.ne.ice1 .and. ij.ne.ice2 .and. aij.eq.6) then
                     ich = ij
                  end if
               end do
            end if
         end if
c
c     inspect the zeta position of amino acid residue
c
         if (label.eq.'   ' .and. ice1.ne.0) then
            do j = 1, n12(ice1)
               ij = i12(j,ice1)
               aij = atomic(ij)
               if (ij.ne.icd1 .and. aij.eq.6) then
                  if (icz1 .ne. 0)  icz2 = ij
                  if (icz1 .eq. 0)  icz1 = ij
               end if
               if (aij .eq. 7)  inz = ij
            end do
            if (inz.ne.0 .and. n12(ice1).eq.4) then
               if (n12(inz) .eq. 3)  label = 'LYD'
               if (n12(inz) .eq. 4)  label = 'LYS'
            end if
         end if
         if (label.eq.'   ' .and. ice2.ne.0) then
            do j = 1, n12(ice2)
               ij = i12(j,ice2)
               aij = atomic(ij)
               if (ij.ne.icd2 .and. aij.eq.6) then
                  if (icz1 .ne. 0)  icz2 = ij
                  if (icz1 .eq. 0)  icz1 = ij
               end if
            end do
            if (icz1 .eq. icz2) then
               icz2 = 0
               label = 'PHE'
               do j = 1, n12(icz1)
                  ij = i12(j,icz1)
                  aij = atomic(ij)
                  if (aij .eq. 8) then
                     ioh = ij
                     if (n12(ij) .eq. 1)  label = 'TYD'
                     if (n12(ij) .eq. 2)  label = 'TYR'
                  end if
               end do
            end if
         end if
         if (label.eq.'   ' .and. ine.ne.0) then
            do j = 1, n12(ine)
               ij = i12(j,ine)
               aij = atomic(ij)
               if (ij.ne.icd1 .and. aij.eq.6) then
                  if (icz1 .ne. 0)  icz2 = ij
                  if (icz1 .eq. 0)  icz1 = ij
               end if
            end do
            if (icz1 .ne. 0) then
               label = 'ARG'
               do j = 1, n12(icz1)
                  ij = i12(j,icz1)
                  aij = atomic(ij)
                  if (aij .ne. 7) then
                     label = '   '
                  else if (ij .ne. ine) then
                     if (inh1 .ne. 0)  inh2 = ij
                     if (inh1 .eq. 0)  inh1 = ij
                  end if
               end do
            end if
         end if
c
c     propagate the tier name to all atoms of the residue
c
         if (label .ne. '   ') then
            tier(ic) = label
            do j = 1, n12(ic)
               ij = i12(j,ic)
               tier(ij) = label
            end do
            tier(ib) = label
            do j = 1, n12(ib)
               ij = i12(j,ib)
               aij = atomic(ij)
               if (aij .eq. 1)  tier(ij) = label
               if (aij .eq. 8)  tier(ij) = label
            end do
            tier(id) = label
            do j = 1, n12(id)
               ij = i12(j,id)
               aij = atomic(ij)
               if (aij .eq. 1)  tier(ij) = label
               if (aij .eq. 8)  tier(ij) = label
            end do
            if (icb1 .ne. 0) then
               tier(icb1) = label
               do j = 1, n12(icb1)
                  ij = i12(j,icb1)
                  tier(ij) = label
               end do
            end if
            if (icb2 .ne. 0) then
               tier(icb2) = label
               do j = 1, n12(icb2)
                  ij = i12(j,icb2)
                  tier(ij) = label
               end do
            end if
            if (icg1 .ne. 0) then
               tier(icg1) = label
               do j = 1, n12(icg1)
                  ij = i12(j,icg1)
                  tier(ij) = label
               end do
            end if
            if (icg2 .ne. 0) then
               tier(icg2) = label
               do j = 1, n12(icg2)
                  ij = i12(j,icg2)
                  tier(ij) = label
               end do
            end if
            if (icd1 .ne. 0) then
               tier(icd1) = label
               do j = 1, n12(icd1)
                  ij = i12(j,icd1)
                  tier(ij) = label
               end do
            end if
            if (icd2 .ne. 0) then
               tier(icd2) = label
               do j = 1, n12(icd2)
                  ij = i12(j,icd2)
                  tier(ij) = label
               end do
            end if
            if (ice1 .ne. 0) then
               tier(ice1) = label
               do j = 1, n12(ice1)
                  ij = i12(j,ice1)
                  tier(ij) = label
               end do
            end if
            if (ice2 .ne. 0) then
               tier(ice2) = label
               do j = 1, n12(ice2)
                  ij = i12(j,ice2)
                  tier(ij) = label
               end do
            end if
            if (icz1 .ne. 0) then
               tier(icz1) = label
               do j = 1, n12(icz1)
                  ij = i12(j,icz1)
                  tier(ij) = label
               end do
            end if
            if (icz2 .ne. 0) then
               tier(icz2) = label
               do j = 1, n12(icz2)
                  ij = i12(j,icz2)
                  tier(ij) = label
               end do
            end if
            if (ich .ne. 0) then
               tier(ich) = label
               do j = 1, n12(ich)
                  ij = i12(j,ich)
                  tier(ij) = label
               end do
            end if
            if (ind .ne. 0) then
               tier(ind) = label
               do j = 1, n12(ind)
                  ij = i12(j,ind)
                  tier(ij) = label
               end do
            end if
            if (ine .ne. 0) then
               tier(ine) = label
               do j = 1, n12(ine)
                  ij = i12(j,ine)
                  tier(ij) = label
               end do
            end if
            if (inz .ne. 0) then
               tier(inz) = label
               do j = 1, n12(inz)
                  ij = i12(j,inz)
                  tier(ij) = label
               end do
            end if
            if (inh1 .ne. 0) then
               tier(inh1) = label
               do j = 1, n12(inh1)
                  ij = i12(j,inh1)
                  tier(ij) = label
               end do
            end if
            if (inh2 .ne. 0) then
               tier(inh2) = label
               do j = 1, n12(inh2)
                  ij = i12(j,inh2)
                  tier(ij) = label
               end do
            end if
            if (iog .ne. 0) then
               tier(iog) = label
               do j = 1, n12(iog)
                  ij = i12(j,iog)
                  tier(ij) = label
               end do
            end if
            if (iod1 .ne. 0) then
               tier(iod1) = label
               do j = 1, n12(iod1)
                  ij = i12(j,iod1)
                  tier(ij) = label
               end do
            end if
            if (iod2 .ne. 0) then
               tier(iod2) = label
               do j = 1, n12(iod2)
                  ij = i12(j,iod2)
                  tier(ij) = label
               end do
            end if
            if (ioe1 .ne. 0) then
               tier(ioe1) = label
               do j = 1, n12(ioe1)
                  ij = i12(j,ioe1)
                  tier(ij) = label
               end do
            end if
            if (ioe2 .ne. 0) then
               tier(ioe2) = label
               do j = 1, n12(ioe2)
                  ij = i12(j,ioe2)
                  tier(ij) = label
               end do
            end if
            if (ioh .ne. 0) then
               tier(ioh) = label
               do j = 1, n12(ioh)
                  ij = i12(j,ioh)
                  tier(ij) = label
               end do
            end if
            if (isg .ne. 0) then
               tier(isg) = label
               do j = 1, n12(isg)
                  ij = i12(j,isg)
                  tier(ij) = label
               end do
            end if
            if (isd .ne. 0) then
               tier(isd) = label
               do j = 1, n12(isd)
                  ij = i12(j,isd)
                  tier(ij) = label
               end do
            end if
         end if
   10    continue
      end do
      return
      end
