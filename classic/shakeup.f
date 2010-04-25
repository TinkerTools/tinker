c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine shakeup  --  setup of rattle constraints  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "shakeup" initializes any holonomic constraints for use
c     with the rattle algorithm during molecular dynamics
c
c
      subroutine shakeup
      implicit none
      include 'sizes.i'
      include 'angle.i'
      include 'atmlst.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'bound.i'
      include 'couple.i'
      include 'keys.i'
      include 'math.i'
      include 'molcul.i'
      include 'shake.i'
      include 'usage.i'
      integer i,j,k,m,nh
      integer ia,ib,ic
      integer ja,jb
      integer ilist,next
      real*8 rab,rbc,rac
      real*8 cosine
      logical done
      character*9 rattyp
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     zero out the number of distance and spatial constraints
c
      nrat = 0
      nratx = 0
      use_rattle = .true.
c
c     process keywords containing generic constraint options
c
      do k = 1, nkey
         next = 1
         record = keyline(k)
         call upcase (record)
         call gettext (record,keyword,next)
c
c     get the distance constraint types for the rattle method
c
         if (keyword(1:7) .eq. 'RATTLE ') then
            call getword (record,rattyp,next)
c
c     constrain all bond lengths at their ideal values
c
            if (rattyp(1:5) .eq. 'BONDS') then
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  if (use(ia) .or. use(ib)) then
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = bl(i)
                  end if
               end do
c
c     constrain bonds and independent angles at ideal values
c
            else if (rattyp(1:6) .eq. 'ANGLES') then
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  if (use(ia) .or. use(ib)) then
                     nrat = nrat + 1
                     irat(1,nrat) = ia
                     irat(2,nrat) = ib
                     krat(nrat) = bl(i)
                  end if
               end do
               do i = 1, n
                  if (n12(i) .gt. 1) then
                     do j = 1, 2*n12(i)-3
                        ilist = anglist(j,i)
                        ia = iang(1,ilist)
                        ib = iang(2,ilist)
                        ic = iang(3,ilist)
                        if (use(ia) .or. use(ib) .or. use(ic)) then
                           do m = 1, n12(ib)
                              if (i12(m,ib) .eq. ia) then
                                 rab = bl(bndlist(m,ib))
                              else if (i12(m,ib) .eq. ic) then
                                 rbc = bl(bndlist(m,ib))
                              end if
                           end do
                           cosine = cos(anat(ilist)/radian)
                           rac = sqrt(rab*rab+rbc*rbc
     &                                   -2.0d0*rab*rbc*cosine)
                           nrat = nrat + 1
                           irat(1,nrat) = ia
                           irat(2,nrat) = ic
                           krat(nrat) = rac
                           call chkangle (ia,ib,ic)
                        end if
                     end do
                  end if
               end do
c
c     fix bond length in diatomics to give a rigid molecule
c
            else if (rattyp(1:8) .eq. 'DIATOMIC') then
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  if (n12(ia).eq.1 .and. n12(ib).eq.1) then
                     if (use(ia) .or. use(ib)) then
                        nrat = nrat + 1
                        irat(1,nrat) = ia
                        irat(2,nrat) = ib
                        krat(nrat) = bl(i)
                     end if
                  end if
               end do
c
c     fix bonds and angle in triatomics to give a rigid molecule
c
            else if (rattyp(1:9) .eq. 'TRIATOMIC') then
               do i = 1, nangle
                  ia = iang(1,i)
                  ib = iang(2,i)
                  ic = iang(3,i)
                  if (n12(ia)+n12(ib)+n12(ic) .eq. 4) then
                     rab = bl(bndlist(1,ia))
                     rbc = bl(bndlist(1,ic))
                     cosine = cos(anat(i)/radian)
                     rac = sqrt(rab**2+rbc**2-2.0d0*rab*rbc*cosine)
                     if (use(ia) .or. use(ib)) then
                        nrat = nrat + 1
                        irat(1,nrat) = ia
                        irat(2,nrat) = ib
                        krat(nrat) = rab
                     end if
                     if (use(ib) .or. use(ic)) then
                        nrat = nrat + 1
                        irat(1,nrat) = ib
                        irat(2,nrat) = ic
                        krat(nrat) = rbc
                     end if
                     if (use(ia) .or. use(ib) .or. use(ic)) then
                        nrat = nrat + 1
                        irat(1,nrat) = ia
                        irat(2,nrat) = ic
                        krat(nrat) = rac
                     end if
                  end if
               end do
c
c     fix bonds and angles of each water to give a rigid molecule
c
            else if (rattyp(1:5) .eq. 'WATER') then
               do i = 1, n
                  nh = 0
                  if (atomic(i) .eq. 8) then
                     do j = 1, n12(i)
                        if (atomic(i12(j,i)) .eq. 1)  nh = nh + 1
                     end do
                  end if
                  if (nh .ge. 2) then
                     do j = 1, n12(i)
                        ilist = bndlist(j,i)
                        ia = ibnd(1,ilist)
                        ib = ibnd(2,ilist)
                        if (use(ia) .or. use(ib)) then
                           nrat = nrat + 1
                           irat(1,nrat) = ia
                           irat(2,nrat) = ib
                           krat(nrat) = bl(ilist)
                        end if
                     end do
                     do j = 1, 2*n12(i)-3
                        ilist = anglist(j,i)
                        ia = iang(1,ilist)
                        ib = iang(2,ilist)
                        ic = iang(3,ilist)
                        if (use(ia) .or. use(ib) .or. use(ic)) then
                           do m = 1, n12(ib)
                              if (i12(m,ib) .eq. ia) then
                                 rab = bl(bndlist(m,ib))
                              else if (i12(m,ib) .eq. ic) then
                                 rbc = bl(bndlist(m,ib))
                              end if
                           end do
                           cosine = cos(anat(ilist)/radian)
                           rac = sqrt(rab*rab+rbc*rbc
     &                                   -2.0d0*rab*rbc*cosine)
                           nrat = nrat + 1
                           irat(1,nrat) = ia
                           irat(2,nrat) = ic
                           krat(nrat) = rac
                           call chkangle (ia,ib,ic)
                        end if
                     end do
                  end if
               end do
c
c     fix all bonds to hydrogen atoms at their ideal length
c
            else
               do i = 1, nbond
                  ia = ibnd(1,i)
                  ib = ibnd(2,i)
                  if (atomic(ia).eq.1 .or. atomic(ib).eq.1) then
                     if (use(ia) .or. use(ib)) then
                        nrat = nrat + 1
                        irat(1,nrat) = ia
                        irat(2,nrat) = ib
                        krat(nrat) = bl(i)
                     end if
                  end if
               end do
            end if
         end if
      end do
c
c     process keywords containing specific distance constraints
c
      do k = 1, nkey
         next = 1
         record = keyline(k)
         call upcase (record)
         call gettext (record,keyword,next)
         if (keyword(1:16) .eq. 'RATTLE-DISTANCE ') then
            call getnumb (record,ia,next)
            call getnumb (record,ib,next)
            rab = 0.0d0
            string = record(next:120)
            read (string,*,err=10,end=10)  rab
   10       continue
            if (rab .eq. 0.0d0) then
               do i = 1, n12(ia)
                  if (i12(i,ia) .eq. ib) then
                     rab = bl(bndlist(i,ia))
                  end if
               end do
            end if
            if (rab .eq. 0.0d0) then
               rab = sqrt((x(ia)-x(ib))**2 + (y(ia)-y(ib))**2
     &                           + (z(ia)-z(ib))**2)
            end if
            done = .false.
            do j = 1, nrat
               ja = irat(1,j)
               jb = irat(2,j)
               if ((ia.eq.ja .and. ib.eq.jb) .or.
     &             (ia.eq.jb .and. ib.eq.ja)) then
                  done = .true.
                  krat(j) = rab
               end if
            end do
            if (.not. done) then
               nrat = nrat + 1
               irat(1,nrat) = ia
               irat(2,nrat) = ib
               krat(nrat) = rab
            end if
c
c     process keywords containing atom group spatial constraints
c
         else if (keyword(1:13) .eq. 'RATTLE-PLANE ') then
            call getnumb (record,ia,next)
            nratx = nratx + 1
            iratx(nratx) = ia
            kratx(nratx) = 1
         else if (keyword(1:12) .eq. 'RATTLE-LINE ') then
            call getnumb (record,ia,next)
            nratx = nratx + 1
            iratx(nratx) = ia
            kratx(nratx) = 2
         else if (keyword(1:14) .eq. 'RATTLE-ORIGIN ') then
            call getnumb (record,ia,next)
            nratx = nratx + 1
            iratx(nratx) = ia
            kratx(nratx) = 3
         end if
      end do
c
c     find and remove any duplicate distance constraints
c
      do i = 1, nrat-1
         ia = irat(1,i)
         ib = irat(2,i)
         do j = i+1, nrat
            ja = irat(1,j)
            jb = irat(2,j)
            if ((ia.eq.ja .and. ib.eq.jb) .or.
     &          (ia.eq.jb .and. ib.eq.ja))  krat(j) = -1.0d0
         end do
      end do
      k = nrat
      do i = k, 1, -1
         if (krat(i) .lt. 0.0d0) then
            do j = i, k-1
               irat(1,j) = irat(1,j+1)
               irat(2,j) = irat(2,j+1)
               krat(j) = krat(j+1)
            end do
            nrat = nrat - 1
         end if
      end do
c
c     set flag to apply minimum image to intermolecular constraints
c
      do i = 1, nrat
         ia = irat(1,i)
         ib = irat(2,i)
         if (use_bounds .and. (molcule(ia).ne.molcule(ib))) then
            ratimage(i) = .true.
         else if (use_polymer) then
            ratimage(i) = .true.
         else
            ratimage(i) = .false.
         end if
      end do
c
c     if no rattle constraints are used, turn off its use
c
      if (nrat.eq.0 .and. nratx.eq.0)  use_rattle = .false.
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine chkangle  --  eliminate redundant constraints  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "chkangle" tests angles to be constrained for their presence
c     in small rings and removes constraints that are redundant
c
c     note this version correctly handles isolated small rings,
c     but should remove one additional redundant constraint for
c     each ring fusion
c
c
      subroutine chkangle (ia,ib,ic)
      implicit none
      include 'sizes.i'
      include 'couple.i'
      include 'ring.i'
      include 'shake.i'
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,imin
      logical remove
c
c
c     all internal angles of 3-membered rings are redundant
c
      remove = .false.
      if (nring3 .ne. 0) then
         do i = 1, n12(ia)
            j = i12(i,ia)
            if (j .eq. ic)  remove = .true.
         end do
      end if
c
c     for 4-membered rings, two internal angles are redundant
c
      if (nring4 .ne. 0) then
         do i = 1, n12(ia)
            id = i12(i,ia)
            if (id .ne. ib) then
               do j = 1, n12(id)
                  k = i12(j,id)
                  if (k .eq. ic) then
                     imin = min(ia,ib,ic,id)
                     if (ib .eq. imin)  remove = .true.
                     if (id .eq. imin)  remove = .true.
                  end if
               end do
            end if
         end do
      end if
c
c     for 5-membered rings, one internal angle is redundant
c
      if (nring5 .ne. 0) then
         do i = 1, n12(ia)
            id = i12(i,ia)
            if (id.ne.ib .and. id.ne.ic) then
               do j = 1, n12(ic)
                  ie = i12(j,ic)
                  if (ie.ne.ib .and. ie.ne.ia) then
                     do k = 1, n12(id)
                        if (i12(k,id) .eq. ie) then
                           imin = min(ia,ib,ic,id,ie)
                           if (ib .eq. imin)  remove = .true.
                        end if
                     end do
                  end if
               end do
            end if
         end do
      end if
c
c     if the constraint if redundant, remove it from the list
c
      if (remove)  nrat = nrat - 1
      return
      end
