c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine findnuc  --  search for RNA & DNA nucleotides  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "findnuc" locates and stores the atoms in nucleotide units
c     based on atomic element and connectivity information
c
c     note this routine assumes an all-atom model with hydrogen 
c     atoms explicitly represented
c
c
      subroutine findnuc
      use atomid
      use atoms
      use bitor
      use couple
      implicit none
      integer i,j,k
      integer nhyd,nphos
      integer ia,ib,ic
      integer id,ie,ij,ik
      integer nia,nib,nic
      integer nid,nie,nij
      integer aia,aib,aic
      integer aid,aie
      integer aij,aik
      integer ic1s,ic2s,ic3s
      integer ic4s,ic5s,io2s
      integer io3s,io4s,io5s
      integer ipo,iop1
      integer iop2,iop3
      integer in1,in3,in4,in5
      integer in6,in7,in8,in9
      integer ic2,ic3,ic4,ic5
      integer ic6,ic7,ic9
      integer icm,io4,io9
      logical proceed
      logical deoxy
      character*3 label
c
c
c     cycle over bitorsions to search for nucleotides
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
c     check for five-atom stretch from 5' to 3' sugar oxygens
c
         proceed = .true.
         if (aia.ne.8 .or. aib.ne.6 .or. aic.ne.6
     &          .or. aid.ne.6 .or. aie.ne.8)  proceed = .false.
         if (proceed) then
            if (nia.ne.2 .or. nib.ne.4 .or. nic.ne.4
     &             .or. nid.ne.4 .or. nie.ne.2)  proceed = .false.
         end if
c
c     test for phosphorus atoms including uncapped termini
c
         if (proceed) then
            proceed = .false.
            nhyd = 0
            nphos = 0
            do j = 1, n12(ia)
               ij = i12(j,ia)
               aij = atomic(ij)
               if (aij .eq. 1)  nhyd = nhyd + 1
               if (aij .eq. 15)  nphos = nphos + 1
            end do 
            do j = 1, n12(ie)
               ij = i12(j,ie)
               aij = atomic(ij)
               if (aij .eq. 1)  nhyd = nhyd + 1
               if (aij .eq. 15)  nphos = nphos + 1
            end do
            if (nphos.eq.1 .and. nhyd.eq.1)  proceed = .true.
            if (nphos.eq.2 .and. nhyd.eq.0)  proceed = .true.
         end if
         if (.not. proceed)  goto 10
c
c     zero out residue name and possible nucleotide atoms
c
         label = '   '
         ic1s = 0
         ic2s = 0
         ic3s = 0
         ic4s = 0
         ic5s = 0
         io2s = 0
         io3s = 0
         io4s = 0
         io5s = 0
         ipo = 0
         iop1 = 0
         iop2 = 0
         iop3 = 0
         ic2 = 0
         ic3 = 0
         ic4 = 0
         ic5 = 0
         ic6 = 0
         ic7 = 0
         ic9 = 0
         icm = 0
         in1 = 0
         in3 = 0
         in4 = 0
         in5 = 0
         in6 = 0
         in7 = 0
         in8 = 0
         in9 = 0
         io4 = 0
         io9 = 0
c
c     locate sugar ring and assign corresponding atom names
c
         nhyd = 0
         do j = 1, n12(ib)
            ij = i12(j,ib)
            aij = atomic(ij)
            if (aij .eq. 1)  nhyd = nhyd + 1
         end do
         if (nhyd .eq. 2) then
            io5s = ia
            ic5s = ib
            ic4s = ic
            ic3s = id
            io3s = ie
         else
            io5s = ie
            ic5s = id
            ic4s = ic
            ic3s = ib
            io3s = ia
         end if
         do j = 1, n12(ic4s)
            ij = i12(j,ic4s)
            aij = atomic(ij)
            if (aij .eq. 8)  io4s = ij
         end do
         do j = 1, n12(ic3s)
            ij = i12(j,ic3s)
            aij = atomic(ij)
            if (ij.ne.ic4s .and. aij.eq.6)  ic2s = ij
         end do
         deoxy = .true.
         if (ic2s .ne. 0) then
            do j = 1, n12(ic2s)
               ij = i12(j,ic2s)
               aij = atomic(ij)
               if (aij .eq. 8) then
                  io2s = ij
                  deoxy = .false.
               end if
            end do
         end if
         if (io4s .ne. 0) then
            do j = 1, n12(io4s)
               ij = i12(j,io4s)
               aij = atomic(ij)
               if (ij.ne.ic4s .and. aij.eq.6)  ic1s = ij
            end do
         end if
c
c     find phosphate group attached at 5' sugar oxygen
c
         do j = 1, n12(io5s)
            ij = i12(j,io5s)
            aij = atomic(ij)
            if (aij .eq. 15) then
               ipo = ij
            end if
         end do
         if (ipo .ne. 0) then
            do j = 1, n12(ipo)
               ij = i12(j,ipo)
               nij = n12(ij)
               aij = atomic(ij)
               if (nij.eq.1 .and. aij.eq.8) then
                  if (iop2 .ne. 0)  iop3 = ij
                  if (iop1 .ne. 0)  iop2 = ij
                  if (iop1 .eq. 0)  iop1 = ij
               end if
            end do
         end if
c
c     decide if the nucleobase is purine or pyrimidine
c
         if (ic1s .ne. 0) then
            do j = 1, n12(ic1s)
               ij = i12(j,ic1s)
               aij = atomic(ij)
               if (aij .eq. 7)  in1 = ij
            end do
         end if
         if (in1 .ne. 0) then 
            do j = 1, n12(in1)
               ij = i12(j,in1)
               nij = n12(ij)
               aij = atomic(ij)
               if (ij.ne.ic1s .and. aij.eq.6) then
                  do k = 1, n12(ij)
                     ik = i12(k,ij)
                     aik = atomic(ik)
                     if (aik .eq. 1)  ic2 = ij
                     if (aik .eq. 8) then
                        label= 'PYR'
                        ic6 = ij
                     end if
                  end do
                  if (ic6 .eq. 0) then
                     label = 'PUR'
                     ic5 = ij
                  end if
               end if
            end do
         end if
c
c     identify additional atoms of a purine base    
c
         if (label .eq. 'PUR') then
            if (ic2 .ne. 0) then
               do j = 1, n12(ic2)
                  ij = i12(j,ic2)
                  aij = atomic(ij)
                  if (ij.ne.in1 .and. aij.eq.7)  in3 = ij
               end do
            end if
            if (in3 .ne. 0) then
               do j = 1, n12(in3)
                  ij = i12(j,in3)
                  aij = atomic(ij)
                  if (ij.ne.ic2 .and. aij.eq.6)  ic4 = ij
               end do
            end if
            if (ic5 .ne. 0) then
               do j = 1, n12(ic5)
                  ij = i12(j,ic5)
                  aij = atomic(ij)
                  if (ij.ne.in1 .and. aij.eq.7)  in6 = ij
               end do
            end if
            if (in6 .ne. 0) then
               do j = 1, n12(in6)
                  ij = i12(j,in6)
                  aij = atomic(ij)
                  if (ij.ne.ic5 .and. aij.eq.6)  ic7 = ij
               end do
            end if
            if (ic7 .ne. 0) then
               do j = 1, n12(ic7)
                  ij = i12(j,ic7)
                  aij = atomic(ij)
                  if (ij.ne.in6 .and. aij.eq.7) then
                     nhyd = 0
                     do k = 1, n12(ij)
                        ik = i12(k,ij)
                        aik = atomic(ik)
                        if (aik .eq. 1)  nhyd = nhyd + 1
                     end do
                     if (nhyd .le. 1)  in8 = ij
                     if (nhyd .eq. 2)  in7 = ij
                  end if
               end do
            end if
            if (in8 .ne. 0) then
               do j = 1, n12(in8)
                  ij = i12(j,in8)
                  aij = atomic(ij)
                  if (ij.ne.ic7 .and. aij.eq.6)  ic9 = ij
               end do
            end if
            if (ic9 .ne. 0) then
               do j = 1, n12(ic9)
                  ij = i12(j,ic9)
                  aij = atomic(ij)
                  if (ij.ne.in8 .and. aij.eq.7)  in9 = ij
                  if (aij .eq. 8)  io9 = ij
               end do
            end if
            if (io9 .ne. 0) then
               label = '  G'
               if (deoxy)  label = ' DG'
            end if
            if (in9 .ne. 0) then
               label = '  A'
               if (deoxy)  label = ' DA'
            end if
         end if
c
c     identify additional atoms of a pyrimidine base    
c
         if (label .eq. 'PYR') then
            if (ic6 .ne. 0) then
               do j = 1, n12(ic6)
                  ij = i12(j,ic6)
                  aij = atomic(ij)
                  if (ij.ne.in1 .and. aij.eq.7)  in5 = ij
               end do
            end if
            if (in5 .ne. 0) then
               do j = 1, n12(in5)
                  ij = i12(j,in5)
                  aij = atomic(ij)
                  if (ij.ne.ic6 .and. aij.eq.6)  ic4 = ij
               end do
            end if
            if (ic4 .ne. 0) then
               do j = 1, n12(ic4)
                  ij = i12(j,ic4)
                  aij = atomic(ij)
                  if (aij .eq. 6)  ic3 = ij
                  if (ij.ne.in5 .and. aij.eq.7)  in4 = ij
                  if (aij .eq. 8)  io4 = ij
               end do
            end if
            label = '  U'
            if (in5 .ne. 0) then
               if (n12(in5) .eq. 2) then
                  label = '  C'
                  if (deoxy)  label = ' DC'
               end if
            end if
            if (ic3 .ne. 0) then
               do j = 1, n12(ic3)
                  ij = i12(j,ic3)
                  aij = atomic(ij)
                  if (ij.ne.ic2 .and. ij.ne.ic4 .and. aij.eq.6) then
                     label = '  T'
                     icm = ij
                  end if
               end do
            end if
         end if
c
c     propagate the tier name to all atoms of the nucleotide
c
         if (ipo .ne. 0)  tier(ipo) = label
         if (iop1 .ne. 0) then
            tier(iop1) = label
            do j = 1, n12(iop1)
               ij = i12(j,iop1)
               tier(ij) = label
            end do
         end if
         if (iop2 .ne. 0) then
            tier(iop2) = label
            do j = 1, n12(iop2)
               ij = i12(j,iop2)
               tier(ij) = label
            end do
         end if
         if (iop3 .ne. 0) then
            tier(iop3) = label
            do j = 1, n12(iop3)
               ij = i12(j,iop3)
               tier(ij) = label
            end do
         end if
         tier(io5s) = label
         do j = 1, n12(io5s)
            ij = i12(j,io5s)
            tier(ij) = label
         end do
         tier(ic5s) = label
         do j = 1, n12(ic5s)
            ij = i12(j,ic5s)
            tier(ij) = label
         end do
         tier(ic4s) = label
         do j = 1, n12(ic4s)
            ij = i12(j,ic4s)
            tier(ij) = label
         end do
         if (io4s .ne. 0) then
            tier(io4s) = label
            do j = 1, n12(io4s)
               ij = i12(j,io4s)
               tier(ij) = label
            end do
         end if
         tier(ic3s) = label
         do j = 1, n12(ic3s)
            ij = i12(j,ic3s)
            tier(ij) = label
         end do
         tier(io3s) = label
         do j = 1, n12(io3s)
            ij = i12(j,io3s)
            tier(ij) = label
         end do
         if (ic2s .ne. 0) then
            tier(ic2s) = label
            do j = 1, n12(ic2s)
               ij = i12(j,ic2s)
               tier(ij) = label
            end do
         end if
         if (io2s .ne. 0) then
            tier(io2s) = label
            do j = 1, n12(io2s)
               ij = i12(j,io2s)
               tier(ij) = label
            end do
         end if
         if (ic1s .ne. 0) then
            tier(ic1s) = label
            do j = 1, n12(ic1s)
               ij = i12(j,ic1s)
               tier(ij) = label
            end do
         end if
         if (in1 .ne. 0) then
            tier(in1) = label
            do j = 1, n12(in1)
               ij = i12(j,in1)
               tier(ij) = label
            end do
         end if
         if (ic2 .ne. 0) then
            tier(ic2) = label
            do j = 1, n12(ic2)
               ij = i12(j,ic2)
               tier(ij) = label
            end do
         end if
         if (ic3 .ne. 0) then
            tier(ic3) = label
            do j = 1, n12(ic3)
               ij = i12(j,ic3)
               tier(ij) = label
            end do
         end if
         if (in3 .ne. 0) then
            tier(in3) = label
            do j = 1, n12(in3)
               ij = i12(j,in3)
               tier(ij) = label
            end do
         end if
         if (ic4 .ne. 0) then
            tier(ic4) = label
            do j = 1, n12(ic4)
               ij = i12(j,ic4)
               tier(ij) = label
            end do
         end if
         if (in4 .ne. 0) then
            tier(in4) = label
            do j = 1, n12(in4)
               ij = i12(j,in4)
               tier(ij) = label
            end do
         end if
         if (ic5 .ne. 0) then
            tier(ic5) = label
            do j = 1, n12(ic5)
               ij = i12(j,ic5)
               tier(ij) = label
            end do
         end if
         if (in5 .ne. 0) then
            tier(in5) = label
            do j = 1, n12(in5)
               ij = i12(j,in5)
               tier(ij) = label
            end do
         end if
         if (ic6 .ne. 0) then
            tier(ic6) = label
            do j = 1, n12(ic6)
               ij = i12(j,ic6)
               tier(ij) = label
            end do
         end if
         if (in6 .ne. 0) then
            tier(in6) = label
            do j = 1, n12(in6)
               ij = i12(j,in6)
               tier(ij) = label
            end do
         end if
         if (ic7 .ne. 0) then
            tier(ic7) = label
            do j = 1, n12(ic7)
               ij = i12(j,ic7)
               tier(ij) = label
            end do
         end if
         if (in7 .ne. 0) then
            tier(in7) = label
            do j = 1, n12(in7)
               ij = i12(j,in7)
               tier(ij) = label
            end do
         end if
         if (in8 .ne. 0) then
            tier(in8) = label
            do j = 1, n12(in8)
               ij = i12(j,in8)
               tier(ij) = label
            end do
         end if
         if (ic9 .ne. 0) then
            tier(ic9) = label
            do j = 1, n12(ic9)
               ij = i12(j,ic9)
               tier(ij) = label
            end do
         end if
         if (in9 .ne. 0) then
            tier(in9) = label
            do j = 1, n12(in9)
               ij = i12(j,in9)
               tier(ij) = label
            end do
         end if
         if (icm .ne. 0) then
            tier(icm) = label
            do j = 1, n12(icm)
               ij = i12(j,icm)
               tier(ij) = label
            end do
         end if
   10    continue
      end do
      return
      end
