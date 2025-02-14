c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine findseq  --  get protein & nucleotide sequence  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "findseq" locates and stores biopolymer sequences for proteins
c     and nucleic acids from connectivity and residue information
c
c
      subroutine findseq
      use atomid
      use atoms
      use bitor
      use couple
      use inform
      use iounit
      use sequen
      use tettor
      use tritor
      implicit none
      integer i,j,ij
      integer ichar,code
      integer nlist,nhyd
      integer ia,ib,ic,id
      integer ie,ig,ih
      integer nia,nib,nic
      integer nid,nie,nig
      integer aia,aib,aic
      integer aid,aie
      integer aig,aih
      integer, allocatable :: flink(:)
      integer, allocatable :: blink(:)
      integer, allocatable :: clink(:)
      integer, allocatable :: list(:,:)
      logical proceed
      character*1 char
c
c
c     perform dynamic allocation of some local arrays
c
      if (allocated(flink))  deallocate (flink)
      if (allocated(blink))  deallocate (blink)
      if (allocated(clink))  deallocate (clink)
      if (allocated(list))  deallocate (list)
      allocate (flink(n))
      allocate (blink(n))
      allocate (clink(n))
      allocate (list(2,n))
c
c     zero out the sequence length and link arrays
c
      nlist = 0
      do i = 1, n
         flink(i) = 0
         blink(i) = 0
         clink(i) = 0
      end do
c
c     search for alpha carbons along polypeptide backbone
c
      do i = 1, ntritor
         ia = itritor(1,i)
         ib = itritor(2,i)
         ic = itritor(3,i)
         id = itritor(4,i)
         ie = itritor(5,i)
         ig = itritor(6,i)
         nia = n12(ia)
         nib = n12(ib)
         nic = n12(ic)
         nid = n12(id)
         nie = n12(ie)
         nig = n12(ig)
         aia = atomic(ia)
         aib = atomic(ib)
         aic = atomic(ic)
         aid = atomic(id)
         aie = atomic(ie)
         aig = atomic(ig)
         if (aia.eq.7 .and. aib.eq.6 .and. aic.eq.6 .and.
     &       aid.eq.7 .and. aie.eq.6 .and. aig.eq.6) then
            if (nia.ge.3 .and. nib.eq.4 .and. nic.eq.3 .and.
     &          nid.eq.3 .and. nie.eq.4 .and. nig.eq.3) then
               nlist = nlist + 1
               list(1,nlist) = ib
               list(2,nlist) = ie
               flink(ib) = ie
               blink(ie) = ib
               clink(ib) = ib
               clink(ie) = ie
            end if
         end if
         if (aia.eq.6 .and. aib.eq.6 .and. aic.eq.7 .and.
     &       aid.eq.6 .and. aie.eq.6 .and. aig.eq.7) then
            if (nia.eq.3 .and. nib.eq.4 .and. nic.eq.3 .and.
     &          nid.eq.3 .and. nie.eq.4 .and. nig.ge.3) then
               nlist = nlist + 1
               list(1,nlist) = ie
               list(2,nlist) = ib
               flink(ie) = ib
               blink(ib) = ie
               clink(ib) = ib
               clink(ie) = ie
            end if
         end if
      end do
c
c     search for N-terminal ACE capping group at end of chain
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
         if (aia.eq.6 .and. aib.eq.6 .and.aic.eq.7
     &       .and. aid.eq.6 .and. aie.eq.6) then
            if (nia.eq.4 .and. nib.eq.3 .and.nic.eq.3
     &          .and. nid.eq.4 .and. nie.eq.3) then
               nhyd = 0
               do j = 1, nia
                  ij = i12(j,ia)
                  if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  nlist = nlist + 1
                  list(1,nlist) = ia
                  list(2,nlist) = id
                  flink(ia) = id
                  blink(id) = ia
                  clink(ia) = ia
                  clink(id) = id
               end if
            end if
         end if
         if (aia.eq.6 .and. aib.eq.6 .and.aic.eq.7
     &       .and. aid.eq.6 .and. aie.eq.6) then
            if (nia.eq.3 .and. nib.eq.4 .and.nic.eq.3
     &          .and. nid.eq.3 .and. nie.eq.4) then
               nhyd = 0
               do j = 1, nie
                  ij = i12(j,ie)
                  if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  nlist = nlist + 1
                  list(1,nlist) = ie
                  list(2,nlist) = ib
                  flink(ie) = ib
                  blink(ib) = ie
                  clink(ie) = ie
                  clink(ib) = ib
               end if
            end if
         end if
c
c     search for C-terminal NME capping group at end of chain
c
         if (aia.eq.7 .and. aib.eq.6 .and.aic.eq.6
     &       .and. aid.eq.7 .and. aie.eq.6) then
            if (nia.eq.3 .and. nib.eq.4 .and.nic.eq.3
     &          .and. nid.eq.3 .and. nie.eq.4) then
               nhyd = 0
               do j = 1, nie
                  ij = i12(j,ie)
                  if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  nlist = nlist + 1
                  list(1,nlist) = ib
                  list(2,nlist) = ie
                  flink(ib) = ie
                  blink(ie) = ib
                  clink(ib) = ib
                  clink(ie) = ie
               end if
            end if
         end if
         if (aia.eq.6 .and. aib.eq.7 .and.aic.eq.6
     &       .and. aid.eq.6 .and. aie.eq.7) then
            if (nia.eq.4 .and. nib.eq.3 .and.nic.eq.3
     &          .and. nid.eq.4 .and. nie.eq.3) then
               nhyd = 0
               do j = 1, nia
                  ij = i12(j,ia)
                  if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
               end do
               if (nhyd .eq. 3) then
                  nlist = nlist + 1
                  list(1,nlist) = id
                  list(2,nlist) = ia
                  flink(id) = ia
                  blink(ia) = id
                  clink(id) = id
                  clink(ia) = ia
               end if
            end if
         end if
      end do
c
c     search for adjacent phosphates at ends of each nucleotide
c
      do i = 1, ntettor
         ia = itettor(1,i)
         ib = itettor(2,i)
         ic = itettor(3,i)
         id = itettor(4,i)
         ie = itettor(5,i)
         ig = itettor(6,i)
         ih = itettor(7,i)
         aia = atomic(ia)
         aib = atomic(ib)
         aid = atomic(id)
         aig = atomic(ig)
         aih = atomic(ih)
         proceed = .false.
         if (aia.eq.15 .and. aih.eq.15)  proceed = .true.
         if (aia.eq.1 .and. aib.eq.8 .and. aih.eq.15)  proceed = .true.
         if (aia.eq.15 .and. aig.eq.8 .and. aih.eq.1)  proceed = .true.
         if (proceed) then
            nhyd = 0
            do j = 1, n12(ic)
               ij = i12(j,ic)
               if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
            end do
            if (nhyd .eq. 2) then
               nlist = nlist + 1
               list(1,nlist) = ia
               list(2,nlist) = ih
               flink(ia) = ih
               blink(ih) = ia
               clink(ia) = id
            end if
            nhyd = 0
            do j = 1, n12(ie)
               ij = i12(j,ie)
               if (atomic(ij) .eq. 1)  nhyd = nhyd + 1
            end do
            if (nhyd .eq. 2) then
               nlist = nlist + 1
               list(1,nlist) = ih
               list(2,nlist) = ia
               flink(ih) = ia
               blink(ia) = ih
               clink(ih) = id
            end if
         end if
      end do
c
c     place the amino acids or nucleotides in sequence order
c
      nseq = 0
      nchain = 0
      code = ichar('A')
      do i = 1, nlist
         ia = list(1,i)
         ib = blink(ia)
         if (ib .eq. 0) then
            aia = atomic(ia)
            nseq = nseq + 1
            seq(nseq) = tier(ia)
            seqatm(nseq) = clink(ia)
            nchain = nchain + 1
            ichain(1,nchain) = nseq
            chnnam(nchain) = char(code)
            chntyp(nchain) = 'PROTEIN'
            if (aia .eq. 15) then
               chntyp(nchain) = 'NUCLEIC'
            end if
            proceed = .true.
            dowhile (proceed)
               ib = flink(ia)
               ic = clink(ib)
               if (ib.eq.0 .or. ic.eq.0) then
                  proceed = .false.
                  ichain(2,nchain) = nseq
                  code = code + 1
               else
                  aib = atomic(ib)
                  nseq = nseq + 1
                  seq(nseq) = tier(ib)
                  seqatm(nseq) = ic
                  if (aib .eq. 15) then
                     chntyp(nchain) = 'NUCLEIC'
                  end if
                  ia = ib
               end if
            end do
         end if
      end do
      if (nchain .eq. 1)  chnnam(1) = ' '
c
c     print the biopolymer sequence and chain information
c
      if (debug .and. nseq.ne.0) then
         write (iout,10)
   10    format (/,' Biopolymer Sequence Residues :'
     &           //,' Residue',8x,'Name',6x,'Anchor Atom',/)
         do i = 1, nseq
            write (iout,20)  i,seq(i),seqatm(i)
   20       format (i6,11x,a3,6x,i8)
         end do
         write (iout,30)
   30    format (/,' Biopolymer Sequence Chains :'
     &           //,3x,'Chain',8x,'Name',5x,'Residue Range',6x,'Type',/)
         do i = 1, nchain
            if (chnnam(i) .eq. ' ') then
               write (iout,40)  i,(ichain(j,i),j=1,2),chntyp(i)
   40          format (i6,12x,'-',5x,2i6,8x,a7)
            else
               write (iout,50)  i,chnnam(i),(ichain(j,i),j=1,2),
     &                          chntyp(i)
   50          format (i6,12x,a1,5x,2i6,8x,a7)
            end if
         end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (flink)
      deallocate (blink)
      deallocate (clink)
      return
      end
