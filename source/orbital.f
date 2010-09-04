c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine orbital  --  setup for pisystem calculation  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "orbital" finds and organizes lists of atoms in a pisystem,
c     bonds connecting pisystem atoms and torsions whose two
c     central atoms are both pisystem atoms
c
c
      subroutine orbital
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bond.i'
      include 'couple.i'
      include 'iounit.i'
      include 'keys.i'
      include 'piorbs.i'
      include 'potent.i'
      include 'tors.i'
      integer i,j,k,m,ib,ic
      integer iorb,jorb,next
      integer piatoms(maxpi)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set the default values for the pisystem variables
c
      do i = 1, maxpi
         piatoms(i) = 0
      end do
      reorbit = 1
c
c     check the keywords for any lists of pisystem atoms
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'PISYSTEM ') then
            string = record(next:120)
            read (string,*,err=10,end=10)  (piatoms(k),k=1,maxpi)
         end if
   10    continue
      end do
c
c     quit if no pisystem was found for consideration
c
      if (piatoms(1) .eq. 0) then
         use_orbit = .false.
         return
      else
         use_orbit = .true.
      end if
c
c     organize and make lists of the pisystem atoms
c
      do i = 1, n
         listpi(i) = .false.
      end do
      i = 1
      do while (piatoms(i) .ne. 0)
         if (piatoms(i) .gt. 0) then
            listpi(piatoms(i)) = .true.
            i = i + 1
         else
            do j = -piatoms(i), piatoms(i+1)
               listpi(j) = .true.
            end do
            i = i + 2
         end if
      end do
      norbit = 0
      do i = 1, n
         if (listpi(i)) then
            norbit = norbit + 1
            iorbit(norbit) = i
         end if
      end do
c
c     quit if the molecule contains too many piorbitals
c
      if (norbit .gt. maxpi) then
         write (iout,20)
   20    format (' ORBITAL  --  Too many Pi-Orbitals;',
     &           ' Increase MAXPI')
         call fatal
      end if
c
c     find three atoms which define a plane
c     perpendicular to each orbital
c
      call piplane
c
c     find and store the pisystem bonds
c
      nbpi = 0
      do i = 1, norbit-1
         iorb = iorbit(i)
         do j = i, norbit
            jorb = iorbit(j)
            do k = 1, n12(iorb)
               if (i12(k,iorb) .eq. jorb) then
                  nbpi = nbpi + 1
                  do m = 1, nbond
                     if (iorb.eq.ibnd(1,m) .and.
     &                   jorb.eq.ibnd(2,m)) then
                        ibpi(1,nbpi) = m
                        ibpi(2,nbpi) = i
                        ibpi(3,nbpi) = j
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end do
      end do
c
c     find and store the pisystem torsions
c
      ntpi = 0
      do i = 1, ntors
         ib = itors(2,i)
         ic = itors(3,i)
         if (listpi(ib) .and. listpi(ic)) then
            do j = 1, nbpi
               k = ibpi(1,j)
               if (ib.eq.ibnd(1,k).and.ic.eq.ibnd(2,k) .or.
     &             ib.eq.ibnd(2,k).and.ic.eq.ibnd(1,k)) then
                  ntpi = ntpi + 1
                  itpi(1,ntpi) = i
                  itpi(2,ntpi) = j
                  goto 40
               end if
            end do
   40       continue
         end if
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine piplane  --  plane perpendicular to orbital  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "piplane" selects the three atoms which specify the plane
c     perpendicular to each p-orbital; the current version will
c     fail in certain situations, including ketenes, allenes,
c     and isolated or adjacent triple bonds
c
c
      subroutine piplane
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'iounit.i'
      include 'piorbs.i'
      integer i,j,iorb
      integer atmnum,trial
      integer alpha,beta,gamma
      integer attach
      logical done
c
c
c     for each pisystem atom, find a set of atoms which define
c     the p-orbital's plane based on piatom's atomic number and
c     the number and type of attached atoms
c
      do i = 1, norbit
         iorb = iorbit(i)
         attach = n12(iorb)
         atmnum = atomic(iorb)
         done = .false.
c
c     most common case of an atom bonded to three atoms
c
         if (attach .eq. 3) then
            piperp(1,i) = i12(1,iorb)
            piperp(2,i) = i12(2,iorb)
            piperp(3,i) = i12(3,iorb)
            done = .true.
c
c     any non-alkyne atom bonded to exactly two atoms
c
         else if (attach.eq.2 .and. atmnum.ne.6) then
            piperp(1,i) = iorb
            piperp(2,i) = i12(1,iorb)
            piperp(3,i) = i12(2,iorb)
            done = .true.
c
c     atom bonded to four different atoms (usually two lone
c     pairs and two "real" atoms); use the "real" atoms
c
         else if (attach .eq. 4) then
            piperp(1,i) = iorb
            do j = 1, n12(iorb)
               trial = i12(j,iorb)
               if (atomic(trial) .ne. 0) then
                  if (piperp(2,i) .eq. 0) then
                     piperp(2,i) = trial
                  else
                     piperp(3,i) = trial
                     done = .true.
                  end if
               end if
            end do
c
c     "carbonyl"-type oxygen atom, third atom is any atom
c     attached to the "carbonyl carbon"; fails for ketenes
c
         else if (attach.eq.1 .and. atmnum.eq.8) then
            alpha = i12(1,iorb)
            beta = i12(1,alpha)
            if (beta .eq. iorb)  beta = i12(2,alpha)
            piperp(1,i) = iorb
            piperp(2,i) = alpha
            piperp(3,i) = beta
            done = .true.
c
c     an sp nitrogen atom, third atom must be a gamma atom
c
         else if (attach.eq.1 .and. atmnum.eq.7) then
            alpha = i12(1,iorb)
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.iorb .and. listpi(trial) .and.
     &             n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma .eq. alpha)  gamma = i12(2,beta)
            piperp(1,i) = iorb
            piperp(2,i) = alpha
            piperp(3,i) = gamma
c
c     an sp carbon atom; third atom must be an atom attached
c     to the non-sp piatom bonded to the original carbon
c
         else if (attach.eq.2 .and. atmnum.eq.6) then
            alpha = i12(1,iorb)
            if ((n12(alpha).eq.2 .and. atomic(alpha).eq.6) .or.
     &          (n12(alpha).eq.1 .and. atomic(alpha).eq.7))
     &         alpha = i12(2,iorb)
            do j = 1, n12(iorb)
               trial = i12(j,iorb)
               if (trial.ne.iorb .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            do j = 1, n12(alpha)
               trial = i12(j,alpha)
               if (trial.ne.iorb .and. trial.ne.alpha .and.
     &             listpi(trial) .and. n12(trial).eq.3) then
                  beta = trial
                  done = .true.
               end if
            end do
            gamma = i12(1,beta)
            if (gamma.eq.iorb .or. gamma.eq.alpha)  gamma = i12(2,beta)
            piperp(1,i) = iorb
            piperp(2,i) = alpha
            piperp(3,i) = gamma
         end if
c
c     quit if the p-orbital plane remains undefined
c
         if (.not. done) then
            write (iout,10)  iorb
   10       format(/,' PIPLANE  --  Failure to Define a',
     &                ' p-Orbital Plane for Atom',i6)
            call fatal
         end if
      end do
      return
      end
