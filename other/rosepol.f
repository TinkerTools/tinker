c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2021  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program rosepol  --  manipulate atomic multipole values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "rosepol" generates initial sets of charge penetration and
c     polarizability parameters from a minimal coordinates file
c
c
      program rosepol
      use iounit
      use potent
      implicit none
      integer mode,idma
      integer freeunit
      logical exist,query
      character*240 string
c
c
c     get the desired type of coordinate file modification
c
      call initial
c
c     perform the desired multipole manipulation operation
c
      call getxyz
      call attach
      call field
      call katom
      call kmpole
      call kpolar
      call setpolar
      call setpgrp
      call avgpole
      call prtpole
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine setpolar  --  define the polarization model  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "setpolar" assigns atomic polarizabilities, Thole damping or
c     charge penetration parameters, and allows user modification
c
c     note this routine contains directly coded scale factors, and
c     Thole and charge penetration values for atom types that should
c     be updated if the default force field values are modified
c
c
      subroutine setpolar
      use atomid
      use atoms
      use chgpen
      use couple
      use fields
      use iounit
      use mplpot
      use mpole
      use polar
      use polpot
      implicit none
      integer i,j,k,m
      integer jj,ia,ib
      integer atn,next
      real*8 pol,thl
      real*8 pel,pal
      real*8 sixth
      logical exist,query
      logical change
      logical aromatic
      logical chkarom
      character*1 answer
      character*240 record
      character*240 string
      external chkarom
c
c
c     allow the user to select the polarization model
c
      forcefield = 'AMOEBA'
      use_thole = .true.
      use_dirdamp = .false.
      use_chgpen = .false.
      dpequal = .false.
      query = .true.
      answer = ' '
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  answer
         call upcase (answer)
         if (answer.eq.'A' .or. answer.eq.'H')  query = .false.
      end if
   10 continue
      if (query) then
         answer = 'A'
         write (iout,20)
   20    format (/,' Choose Either the AMOEBA or HIPPO Polarization',
     &              ' Model [A] :  ',$)
         read (input,30)  record
   30    format (a240)
         next = 1
         call gettext (record,answer,next)
         call upcase (answer)
      end if
      if (answer .eq. 'H') then
         forcefield = 'HIPPO'
         use_thole = .false.
         use_chgpen = .true.
         dpequal = .true.
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(polarity))  allocate (polarity(n))
      if (.not. allocated(thole))  allocate (thole(n))
      if (.not. allocated(dirdamp))  allocate (dirdamp(n))
      if (.not. allocated(pdamp))  allocate (pdamp(n))
      if (.not. allocated(pcore))  allocate (pcore(n))
      if (.not. allocated(pval))  allocate (pval(n))
      if (.not. allocated(palpha))  allocate (palpha(n))
c
c     zero out the polarization and charge penetration values
c
      do i = 1, n
         polarity(i) = 0.0d0
         thole(i) = 0.0d0
         dirdamp(i) = 0.0d0
         pdamp(i) = 0.0d0
         pcore(i) = 0.0d0
         pval(i) = 0.0d0
         palpha(i) = 0.0d0
      end do
c
c     set multipole and polarization scale factors for AMOEBA
c
      if (forcefield .eq. 'AMOEBA') then
         m2scale = 0.0d0
         m3scale = 0.0d0
         m4scale = 0.4d0
         m5scale = 0.8d0
         p2scale = 0.0d0
         p3scale = 0.0d0
         p4scale = 1.0d0
         p5scale = 1.0d0
         p2iscale = 0.0d0
         p3iscale = 0.0d0
         p4iscale = 0.5d0
         p5iscale = 1.0d0
         d1scale = 0.0d0
         d2scale = 1.0d0
         d3scale = 1.0d0
         d4scale = 1.0d0
         u1scale = 1.0d0
         u2scale = 1.0d0
         u3scale = 1.0d0
         u4scale = 1.0d0
      end if
c
c     set multipole and polarization scale factors for HIPPO
c
      if (forcefield .eq. 'HIPPO') then
         m2scale = 0.0d0
         m3scale = 0.0d0
         m4scale = 0.4d0
         m5scale = 0.8d0
         p2scale = 0.0d0
         p3scale = 0.5d0
         p4scale = 1.0d0
         p5scale = 1.0d0
         p2iscale = 0.0d0
         p3iscale = 0.0d0
         p4iscale = 0.0d0
         p5iscale = 0.5d0
         d1scale = 0.0d0
         d2scale = 1.0d0
         d3scale = 1.0d0
         d4scale = 1.0d0
         u1scale = 1.0d0
         u2scale = 1.0d0
         u3scale = 1.0d0
         u4scale = 1.0d0
         w2scale = 0.2d0
         w3scale = 1.0d0
         w4scale = 1.0d0
         w5scale = 1.0d0
      end if
c
c     assign default atomic polarizabilities for AMOEBA model
c
      if (forcefield .eq. 'AMOEBA') then
         do i = 1, n
            thole(i) = 0.39d0
            atn = atomic(i)
            if (atn .eq. 1) then
               polarity(i) = 0.496d0
            else if (atn .eq. 5) then
               polarity(i) = 1.600d0
            else if (atn .eq. 6) then
               polarity(i) = 1.334d0
            else if (atn .eq. 7) then
               polarity(i) = 1.073d0
            else if (atn .eq. 8) then
               polarity(i) = 0.837d0
            else if (atn .eq. 9) then
               polarity(i) = 0.507d0
            else if (atn .eq. 14) then
               polarity(i) = 3.640d0
            else if (atn .eq. 15) then
               polarity(i) = 1.828d0
            else if (atn .eq. 16) then
               polarity(i) = 3.300d0
            else if (atn .eq. 17) then
               polarity(i) = 2.500d0
            else if (atn .eq. 35) then
               polarity(i) = 3.595d0
            else if (atn .eq. 53) then
               polarity(i) = 5.705d0
            end if
         end do
c
c     alter polarizabilities for alkene/aromatic carbon and hydrogen
c
         do i = 1, n
            atn = atomic(i)
            if (atn .eq. 1) then
               j = i12(1,i)
               if (atomic(j).eq.6 .and. n12(j).eq.3) then
                  polarity(i) = 0.696d0
                  do k = 1, n12(j)
                     m = i12(k,j)
                     if (atomic(m).eq.8 .and. n12(m).eq.1) then
                        polarity(i) = 0.494d0
                     end if
                  end do
               end if
            else if (atn .eq. 6) then
               if (n12(i) .eq. 3) then
                  polarity(i) = 1.75d0
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.8 .and. n12(k).eq.1) then
                        polarity(i) = 1.334d0
                     end if
                  end do
               end if
            end if
         end do
      end if
c
c     assign default atom-based parameters for HIPPO model
c
      if (forcefield .eq. 'HIPPO') then
         do i = 1, n
            atn = atomic(i)
            if (atn .eq. 1) then
               pcore(i) = 1.0d0
               polarity(i) = 0.373d0
               palpha(i) = 4.3225d0
               k = atomic(i12(1,i))
               if (k .eq. 6) then
                  do j = 1, n13(i)
                     m = atomic(i13(j,i))
                     if ((atomic(m).ne.6.or.n12(m).ne.4)
     &                     .and. atomic(m).ne.1)  goto 40
                  end do
                  do j = 1, n14(i)
                     m = i14(j,i)
                     if ((atomic(m).ne.6.or.n12(m).ne.4)
     &                     .and. atomic(m).ne.1)  goto 40
                  end do
                  polarity(i) = 0.504d0
                  palpha(i) = 4.9530d0
   40             continue
                  aromatic = chkarom (k)
                  if (aromatic) then
                     polarity(i) = 0.1106d0
                     palpha(i) = 4.9530d0
                  end if
               else if (k .eq. 7) then
                  polarity(i) = 0.005d0
                  palpha(i) = 5.5155d0
               else if (k .eq. 8) then
                  polarity(i) = 0.3698d0
                  palpha(i) = 4.7441d0
               else if (k .eq. 16) then
                  polarity(i) = 0.2093d0
                  palpha(i) = 4.3952d0
               end if
            else if (atn .eq. 5) then
               pcore(i) = 3.0d0
               polarity(i) = 1.6d0
               palpha(i) = 0.0d0      !! missing parameter
            else if (atn .eq. 6) then
               pcore(i) = 4.0d0
               polarity(i) = 0.9354d0
               palpha(i) = 4.5439d0
               do j = 1, n12(i)
                  k = i12(j,i)
                  if ((atomic(k).ne.6.or.n12(k).ne.4)
     &                  .and. atomic(k).ne.1)  goto 50
               end do
               do j = 1, n13(i)
                  k = atomic(i13(j,i))
                  if ((atomic(k).ne.6.or.n12(k).ne.4)
     &                  .and. atomic(k).ne.1)  goto 50
               end do
               polarity(i) = 0.755d0
               palpha(i) = 4.2998d0
   50          continue
               if (n12(i) .eq. 3) then
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.6 .and. n12(k).eq.3) then
                        polarity(i) = 1.9384d0
                        palpha(i) = 3.5491d0
                     end if
                  end do
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.8 .and. n12(k).eq.1) then
                        polarity(i) = 0.6577d0
                        palpha(i) = 5.9682d0
                     end if
                  end do
               end if
               if (chkarom(i)) then
                  polarity(i) = 1.5624d0
                  palpha(i) = 3.8056d0
                  do j = 1, n12(i)
                     k = atomic(i12(j,i))
                     if (k.ne.6 .and. k.ne.1) then
                        polarity(i) = 1.2811d0
                        palpha(i) = 3.8066d0
                     end if
                  end do
               end if
               if (n12(i) .eq. 2) then
                  polarity(i) = 1.0d0    !! missing parameter
                  palpha(i) = 0.0d0      !! missing parameter
               end if
            else if (atn .eq. 7) then
               pcore(i) = 5.0d0
               polarity(i) = 1.4289d0
               palpha(i) = 3.9882d0
               if (n12(i) .eq. 3) then
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.6 .and. n12(k).eq.3) then
                        polarity(i) = 1.4545d0
                        palpha(i) = 3.9413d0
                     end if
                  end do
               end if
               if (chkarom(i)) then
                  polarity(i) = 1.3037d0
                  palpha(i) = 3.9434d0
               end if
            else if (atn .eq. 8) then
               pcore(i) = 6.0d0
               polarity(i) = 0.6645d0
               palpha(i) = 4.7004d0
               if (n12(i) .eq. 1) then
                  k = i12(1,i)
                  if (atomic(k).eq.6 .and. n12(k).eq.3) then
                     polarity(i) = 1.4266d0
                     palpha(i) = 4.2263d0
                     do j = 1, n13(i)
                        m = i13(j,i)
                        if (atomic(m).eq.8 .and. n12(m).eq.1) then
                           polarity(i) = 1.8809d0
                           palpha(i) = 4.0355d0
                        end if
                     end do
                  end if
                  if (atomic(k) .eq. 15) then
                     jj = 0
                     do j = 1, n12(k)
                        m = i12(j,k)
                        if (atomic(m).eq.8 .and. n12(m).eq.1) then
                           jj = jj + 1
                        end if
                     end do
                     if (jj .eq. 1) then
                        polarity(i) = 1.0d0
                        palpha(i) = 4.3312d0
                     else
                        polarity(i) = 1.0d0
                        palpha(i) = 4.4574d0
                     end if
                  end if
               end if
            else if (atn .eq. 9) then
               pcore(i) = 7.0d0
               polarity(i) = 0.5d0
               palpha(i) = 5.5080d0
            else if (atn .eq. 15) then
               pcore(i) = 5.0d0
               polarity(i) = 1.8d0
               palpha(i) = 2.8130d0
            else if (atn .eq. 16) then
               pcore(i) = 6.0d0
               polarity(i) = 3.1967d0
               palpha(i) = 3.3620d0
               if (n12(i) .gt. 2) then
                  polarity(i) = 2.458d0
                  palpha(i) = 2.7272d0
               end if
            else if (atn .eq. 17) then
               pcore(i) = 7.0d0
               polarity(i) = 2.366d0
               palpha(i) = 3.6316d0
            else if (atn .eq. 35) then
               pcore(i) = 7.0d0
               polarity(i) = 3.4458d0
               palpha(i) = 3.2008d0
            else if (atn .eq. 53) then
               pcore(i) = 7.0d0
               polarity(i) = 5.5d0
               palpha(i) = 0.0d0      !! missing parameter
            end if
         end do
      end if
c
c     set valence electrons from number of core electrons
c
      do i = 1, n
         pval(i) = pole(1,i) - pcore(i)
      end do
c
c     compute the Thole polarizability damping values
c
      sixth = 1.0d0 / 6.0d0
      do i = 1, n
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     allow the user to manually alter polarizability values
c
      change = .false.
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=60,end=60)  i
         if (i .eq. 0)  query = .false.
      end if
   60 continue
      do while (query)
         i = 0
         k = 0
         if (use_thole) then
            pol = 0.0d0
            thl = 0.39d0
            write (iout,70)
   70       format (/,' Enter Atom Number, Polarizability & Thole',
     &                 ' Values :  ',$)
            read (input,80)  record
   80       format (a240)
            read (record,*,err=90,end=90)  k,pol,thl
   90       continue
            if (k .ne. 0) then
               i = pollist(k)
               if (pol .eq. 0.0d0)  pol = polarity(i)
               if (thl .eq. 0.0d0)  thl = thole(i)
            end if
         else if (use_chgpen) then
            pol = 0.0d0
            pel = 0.0d0
            pal = 0.0d0
            write (iout,100)
  100       format (/,' Enter Atom Number, Polarize, Core & Damp',
     &                 ' Values :  ',$)
            read (input,110)  record
  110       format (a240)
            read (record,*,err=120,end=120)  k,pol,pel,pal
  120       continue
            if (k .ne. 0) then
               i = pollist(k)
               if (pol .eq. 0.0d0)  pol = polarity(i)
               if (pel .eq. 0.0d0)  pel = pcore(i)
               if (pal .eq. 0.0d0)  pal = palpha(i)
            end if
         end if
         if (k .eq. 0) then
            query = .false.
         else
            i = pollist(k)
         end if
         if (i .ne. 0) then
            change = .true.
            polarity(i) = pol
            if (use_thole) then
               thole(i) = thl
               pdamp(i) = polarity(i)**sixth
            else if (use_chgpen) then
               pcore(i) = pel
               palpha(i) = pal
               pval(i) = pole(1,i) - pcore(i)
            end if
         end if
      end do
c
c     repeat polarizability values if parameters were altered
c
      if (change) then
         write (iout,130)
  130    format (/,' Atomic Polarizabilities for Multipole Sites :')
         if (use_thole) then
            write (iout,140)
  140       format (/,5x,'Atom',5x,'Name',7x,'Polarize',10x,'Thole',/)
         else if (use_chgpen) then
            write (iout,150)
  150       format (/,5x,'Atom',5x,'Name',7x,'Polarize',4x,'Core Chg',
     &                 8x,'Damp',/)
         end if
         do k = 1, n
            i = pollist(k)
            if (use_thole) then
               if (i .eq. 0) then
                  write (iout,160)  k,name(k)
  160             format (i8,6x,a3,12x,'--',13x,'--')
               else
                  write (iout,170)  k,name(k),polarity(i),thole(i)
  170             format (i8,6x,a3,4x,f12.4,3x,f12.4)
               end if
            else if (use_chgpen) then
               if (i .eq. 0) then
                  write (iout,180)  k,name(k)
  180             format (i8,6x,a3,12x,'--',13x,'--',10x,'--')
               else
                  write (iout,190)  k,name(k),polarity(i),pcore(i),
     &                              palpha(i)
  190             format (i8,6x,a3,4x,f12.4,3x,2f12.4)
               end if
            end if
         end do
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine setpgrp  --  define the polarization groups  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "setpgrp" chooses the polarization groups as defined by bonds
c     separating groups, and allows user modification of the groups
c
c
      subroutine setpgrp
      use atomid
      use atoms
      use bndstr
      use couple
      use iounit
      use kpolr
      use ring
      implicit none
      integer i,j,k,m
      integer mode
      integer ia,ib,ic
      integer ita,itb,itc
      integer ata,atb
      integer n12a,n12b
      logical exist,query
      logical chkarom,split
      logical aroma,aromb
      character*240 record
      character*240 string
c
c
c     get the desired type of coordinate file modification
c
      mode = -1
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         if (mode.ge.1 .and. mode.le.2)  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Choose Method for Division into Polarization',
     &              ' Groups :',
     &           //,4x,'(1) Put All Atoms in One Polarization Group',
     &           /,4x,'(2) Separate into Groups at Rotatable Bonds',
     &           /,4x,'(3) Manual Entry of Bonds Separating Groups')
         do while (mode.lt.1 .or. mode.gt.3)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
            if (mode .le. 0)  mode = 1
   50       continue
         end do
      end if
c
c     initialize by placing all atoms in one polarization group
c
      do i = 1, n
         do j = 1, n12(i)
            pgrp(j,i) = i12(j,i)
         end do
      end do
c
c     separate into polarization groups at rotatable bonds
c
      if (mode .eq. 2) then
         call bonds
         do k = 1, nbond
            ia = ibnd(1,k)
            ib = ibnd(2,k)
            n12a = n12(ia)
            n12b = n12(ib)
            ata = atomic(ia)
            atb = atomic(ib)
            ita = 10*ata + n12a
            itb = 10*atb + n12b
            aroma = chkarom(ia)
            aromb = chkarom(ib)
            split = .true.
c
c     remove bonds involving univalent atoms
c
            if (min(n12a,n12b) .le. 1)  split = .false.
c
c     remove bonds internal to aromatic ring
c
            if (aroma .and. aromb) then
               do i = 1, nring5
                  m = 0
                  do j = 1, 5
                     if (iring5(j,i) .eq. ia)  m = m + 1
                     if (iring5(j,i) .eq. ib)  m = m + 1
                  end do
                  if (m .eq. 2)  split = .false.
               end do
               do i = 1, nring6
                  m = 0
                  do j = 1, 6
                     if (iring6(j,i) .eq. ia)  m = m + 1
                     if (iring6(j,i) .eq. ib)  m = m + 1
                  end do
                  if (m .eq. 2)  split = .false.
               end do
            end if
c
c     remove bonds with sp-hybridized carbon atom
c
            if (ita.eq.62 .or. itb.eq.62)  split = .false.
c
c     remove the C=C bond of terminal alkene
c
            if (ita.eq.63 .and. .not.aroma .and.
     &             itb.eq.63 .and. .not.aromb) then
               split = .false.
               do i = 1, n12a
                  ic = i12(i,ia)
                  if (ic .ne. ib) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc.eq.63 .or. itc.eq.73 .or.
     &                   itc.eq.72 .or. itc.eq.81) then
                        split = .true.
                     end if
                  end if
               end do
               if (split) then
                  split = .false.
                  do i = 1, n12b
                     ic = i12(i,ib)
                     if (ic .ne. ia) then
                        itc = 10*atomic(ic) + n12(ic)
                        if (itc.eq.63 .or. itc.eq.72 .or.
     &                      itc.eq.73 .or. itc.eq.81) then
                           split = .true.
                        end if
                     end if
                  end do
               end if
            end if
c
c     remove the C-O bonds of alcohol and ether
c
            if (ita.eq.82 .and. itb.eq.64) then
               do i = 1, n12a
                  ic = i12(i,ia)
                  if (ic .ne. ib) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc.eq.11 .or. itc.eq.64)  split = .false.
                  end if
               end do
            else if (itb.eq.82 .and. ita.eq.64) then
               do i = 1, n12b
                  ic = i12(i,ib)
                  if (ic .ne. ia) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc.eq.11 .or. itc.eq.64)  split = .false.
                  end if
               end do
            end if
c
c     remove the C-O bonds of carboxylic acid and ester
c
            if (ita.eq.82 .and. itb.eq.63) then
               do i = 1, n12b
                  ic = i12(i,ib)
                  itc = 10*atomic(ic) + n12(ic)
                  if (itc .eq. 81)  split = .false.
               end do
            else if (itb.eq.82 .and. ita.eq.63) then
               do i = 1, n12a
                  ic = i12(i,ia)
                  itc = 10*atomic(ic) + n12(ic)
                  if (itc .eq. 81)  split = .false.
               end do
            end if
c
c     remove the C-N bonds of alkyl amine
c
            if (ita.eq.73 .and. itb.eq.64) then
               m = 0
               do i = 1, n12a
                  ic = i12(i,ia)
                  if (ic .ne. ib) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc.eq.11 .or. itc.eq.64)  m = m + 1
                  end if
               end do
               if (m .eq. 2)  split = .false.
            else if (itb.eq.73 .and. ita.eq.64) then
               m = 0
               do i = 1, n12b
                  ic = i12(i,ib)
                  if (ic .ne. ia) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc.eq.11 .or. itc.eq.64)  m = m + 1
                  end if
               end do
               if (m .eq. 2)  split = .false.
            end if
c
c     remove the C-N bonds of amide, urea, amidine and guanidinium
c
            if (ita.eq.73 .and. itb.eq.63) then
               do i = 1, n12b
                  ic = i12(i,ib)
                  if (ic .ne. ia) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc .eq. 81)  split = .false.
                     if (itc .eq. 73)  split = .false.
                  end if
               end do
            else if (itb.eq.73 .and. ita.eq.63) then
               do i = 1, n12a
                  ic = i12(i,ia)
                  if (ic .ne. ib) then
                     itc = 10*atomic(ic) + n12(ic)
                     if (itc .eq. 81)  split = .false.
                     if (itc .eq. 73)  split = .false.
                  end if
               end do
            end if
c
c     remove any P-X and S-X bonds with X = N or O
c
            if (ata.eq.15 .or. ata.eq.16) then
               if (atb.eq.7 .or. atb.eq.8)  split = .false.
            else if (atb.eq.15 .or. atb.eq.16) then
               if (ata.eq.7 .or. ata.eq.8)  split = .false.
            end if
c
c     modify membership to split groups at allowed bonds
c
            if (split) then
               do i = 1, n12a
                  if (pgrp(i,ia) .eq. ib) then
                     do j = i+1, n12a
                        pgrp(j-1,ia) = pgrp(j,ia)
                     end do
                     pgrp(n12a,ia) = 0
                  end if
               end do
               do i = 1, n12b
                  if (pgrp(i,ib) .eq. ia) then
                     do j = i+1, n12b
                        pgrp(j-1,ib) = pgrp(j,ib)
                     end do
                     pgrp(n12b,ib) = 0
                  end if
               end do
            end if
         end do
c
c     allow modification of polarization group one bond at a time
c
      else if (mode .eq. 3) then
         write (iout,60)
   60    format (/,' All atoms are placed initially into one',
     &              ' polarization group;',
     &           /,' This can be modified by entering a series',
     &              ' of bonded atom pairs',
     &           /,' that separate the molecule into distinct',
     &              ' polarization groups')
c
c     get the bonds that separate the polarization groups
c
         query = .true.
         i = -1
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=70,end=70)  i
            if (i .eq. 0)  query = .false.
         end if
   70    continue
         do while (query)
            ia = 0
            ib = 0
            write (iout,80)
   80       format (/,' Enter a Bond between Polarization Groups',
     &                 ' [<Enter>=Exit] :  ',$)
            read (input,90)  record
   90       format (a240)
            read (record,*,err=100,end=100)  ia,ib
  100       continue
            if (ia.eq.0 .or. ib.eq.0) then
               query = .false.
            else
               do i = 1, n12(ia)
                  if (pgrp(i,ia) .eq. ib) then
                     do j = i+1, n12(ia)
                        pgrp(j-1,ia) = pgrp(j,ia)
                     end do
                     pgrp(n12(ia),ia) = 0
                  end if
               end do
               do i = 1, n12(ib)
                  if (pgrp(i,ib) .eq. ia) then
                     do j = i+1, n12(ib)
                        pgrp(j-1,ib) = pgrp(j,ib)
                     end do
                     pgrp(n12(ib),ib) = 0
                  end if
               end do
            end if
         end do
      end if
c
c     find the polarization groups and their connectivities
c
      call polargrp
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function chkarom  --  check for atom in aromatic ring  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "chkatom" tests for the presence of a specified atom as a
c     member of an aromatic ring
c
c
      function chkarom (iatom)
      use atomid
      use couple
      use ring
      implicit none
      integer i,j,k
      integer iatom
      logical chkarom
      logical member
      logical trigonal
c
c
c     determine membership in 5-membered aromatic ring
c
      chkarom = .false.
      do i = 1, nring5
         trigonal = .true.
         member = .false.
         do j = 1, 5
            k = iring5(j,i)
            if (k .eq. iatom)  member = .true.
            if (atomic(k).eq.6 .and. n12(k).ne.3)  trigonal = .false.
            if (atomic(k).eq.7 .and. n12(k).eq.4)  trigonal = .false.
         end do
         if (member .and. trigonal)  chkarom = .true.
      end do
c
c     determine membership in 6-membered aromatic ring
c
      do i = 1, nring6
         trigonal = .true.
         member = .false.
         do j = 1, 6
            k = iring6(j,i)
            if (k .eq. iatom)  member = .true.
            if (atomic(k).eq.6 .and. n12(k).ne.3)  trigonal = .false.
            if (atomic(k).eq.7 .and. n12(k).eq.4)  trigonal = .false.
         end do
         if (member .and. trigonal)  chkarom = .true.
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine avgpole  --  condense multipole atom types  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "avgpole" condenses the number of multipole atom types based
c     upon atoms with equivalent attachments and additional user
c     specified sets of equivalent atoms
c
c
      subroutine avgpole
      use atomid
      use atoms
      use couple
      use iounit
      use kpolr
      use mpole
      use sizes
      implicit none
      integer i,j,k,m
      integer it,mt
      integer ix,iy,iz
      integer in,jn,ni,nj
      integer size,numtyp
      integer priority
      integer nlist,nave
      integer tmin,tmax
      integer xaxe,yaxe,zaxe
      integer indx(4)
      integer, allocatable :: ci(:)
      integer, allocatable :: cj(:)
      integer, allocatable :: list(:)
      integer, allocatable :: tsort(:)
      integer, allocatable :: pkey(:)
      integer, allocatable :: pgrt(:,:)
      real*8 pave(13)
      logical done,repeat
      logical header,exist
      logical query,condense
      logical match,diff
      logical useframe
      logical symm
      logical yzero,xyzero
      character*1 answer
      character*4 pa,pb,pc,pd
      character*16 ptlast
      character*16, allocatable :: pt(:)
      character*240 record
      character*240 string
c
c
c     check for user requested reduction of equivalent types
c
      condense = .true.
      answer = 'Y'
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  answer
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Condense Symmetric Atoms to Equivalent Types',
     &              ' [Y] :  ',$)
         read (input,30)  answer
   30    format (a1)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  condense = .false.
c
c     perform dynamic allocation of some local arrays
c
      if (condense) then
         allocate (ci(n))
         allocate (cj(n))
         size = 40
         allocate (list(max(n,size)))
c
c     condense groups of equivalent atoms to the same atom type
c
         header = .true.
         do i = 1, n
            list(i) = 0
         end do
         do i = 1, n-1
            ni = n12(i) + n13(i) + n14(i) + n15(i)
            m = 0
            do k = 1, n12(i)
               m = m + 1
               in = i12(k,i)
               ci(m) = 2000 + 10*atomic(in) + n12(in)
            end do
            do k = 1, n13(i)
               m = m + 1
               in = i13(k,i)
               ci(m) = 3000 + 10*atomic(in) + n12(in)
            end do
            do k = 1, n14(i)
               m = m + 1
               in = i14(k,i)
               ci(m) = 4000 + 10*atomic(in) + n12(in)
            end do
            do k = 1, n15(i)
               m = m + 1
               in = i15(k,i)
               ci(m) = 5000 + 10*atomic(in) + n12(in)
            end do
            call sort (ni,ci)
            do j = i+1, n
               if (atomic(i) .eq. atomic(j)) then
                  nj = n12(j) + n13(j) + n14(j) + n15(j)
                  if (nj .eq. ni) then
                     m = 0
                     do k = 1, n12(j)
                        m = m + 1
                        jn = i12(k,j)
                        cj(m) = 2000 + 10*atomic(jn) + n12(jn)
                     end do
                     do k = 1, n13(j)
                        m = m + 1
                        jn = i13(k,j)
                        cj(m) = 3000 + 10*atomic(jn) + n12(jn)
                     end do
                     do k = 1, n14(j)
                        m = m + 1
                        jn = i14(k,j)
                        cj(m) = 4000 + 10*atomic(jn) + n12(jn)
                     end do
                     do k = 1, n15(j)
                        m = m + 1
                        jn = i15(k,j)
                        cj(m) = 5000 + 10*atomic(jn) + n12(jn)
                     end do
                     call sort (nj,cj)
                     match = .true.
                     do k = 1, ni
                        if (ci(k) .ne. cj(k)) then
                           match = .false.
                           goto 40
                        end if
   40                   continue
                     end do
                     if (match) then
                        tmin = min(type(i),type(j))
                        tmax = max(type(i),type(j))
                        do k = 1, n
                           if (type(k) .eq. tmax)  type(k) = tmin
                        end do
                        if (list(i).eq.0 .or. list(j).eq.0) then
                           list(i) = 1
                           list(j) = 1
                        end if
                     end if
                  end if
               end if
            end do
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (ci)
         deallocate (cj)
c
c     perform dynamic allocation of some local arrays
c
         allocate (tsort(maxtyp))
         allocate (pkey(maxtyp))
         allocate (pt(maxtyp))
c
c     identify atoms with the same atom type number, or find
c     atoms with equivalent local frame defining atom types
c
         useframe = .false.
         do i = 1, npole
            it = type(i)
            zaxe = 0
            xaxe = 0
            yaxe = 0
            if (useframe) then
               if (zaxis(i) .ne. 0)  zaxe = type(zaxis(i))
               if (xaxis(i) .ne. 0)  xaxe = type(xaxis(i))
               if (yaxis(i) .ne. 0)  yaxe = type(yaxis(i))
            end if
            size = 4
            call numeral (it,pa,size)
            call numeral (zaxe,pb,size)
            call numeral (xaxe,pc,size)
            call numeral (yaxe,pd,size)
            pt(i) = pa//pb//pc//pd
         end do
         call sort7 (npole,pt,pkey)
c
c     average the multipole values at equivalent atom sites
c
         nave = 1
         ptlast = '                '
         do i = 1, npole
            it = pkey(i)
            if (pt(i) .eq. ptlast) then
               nave = nave + 1
               do j = 1, 13
                  pave(j) = pave(j) + pole(j,it)
               end do
               if (i .eq. npole) then
                  do j = 1, 13
                     pave(j) = pave(j) / dble(nave)
                  end do
                  do m = 1, nave
                     mt = pkey(i-m+1)
                     do j = 1, 13
                        pole(j,mt) = pave(j)
                     end do
                  end do
               end if
            else
               if (nave .ne. 1) then
                  do j = 1, 13
                     pave(j) = pave(j) / dble(nave)
                  end do
                  do m = 1, nave
                     mt = pkey(i-m)
                     do j = 1, 13
                        pole(j,mt) = pave(j)
                     end do
                  end do
               end if
               nave = 1
               do j = 1, 13
                  pave(j) = pole(j,it)
               end do
               ptlast = pt(i)
            end if
         end do
c
c     perform deallocation of some local arrays
c
         deallocate (list)
         deallocate (tsort)
         deallocate (pkey)
         deallocate (pt)
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pgrt(maxval,maxtyp))
c
c     set polarization groups over the condensed atom types
c
      do i = 1, n
         do j = 1, maxval
            pgrt(j,i) = 0
         end do
      end do
      do i = 1, npole
         it = type(i)
         k = 0
         do j = 1, maxval
            if (pgrp(j,i) .ne. 0) then
               k = j
               pgrp(j,it) = type(pgrp(j,i))
            else
               pgrp(j,it) = 0
            end if
         end do
         call sort8 (k,pgrp(1,it))
         do j = 1, k
            do m = 1, maxval
               if (pgrt(m,it) .eq. 0) then
                  pgrt(m,it) = pgrp(j,it)
                  goto 50
               end if
            end do
   50       continue
         end do
      end do
      do i = 1, npole
         it = type(i)
         do j = 1, maxval
            pgrp(j,it) = pgrt(j,it)
            if (pgrp(j,it) .ne. 0)  k = j
         end do
         call sort8 (k,pgrp(1,it))
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pgrt)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine prtpole  --  create file with final multipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtpole" creates a coordinates file, and a key file with
c     atomic multipoles corrected for intergroup polarization
c
c
      subroutine prtpole
      use atoms
      use atomid
      use chgpen
      use files
      use keys
      use kpolr
      use mplpot
      use mpole
      use polar
      use polpot
      use sizes
      use units
      implicit none
      integer i,j,k,it
      integer ixyz,ikey
      integer size
      integer xaxe
      integer yaxe
      integer zaxe
      integer freeunit
      integer trimtext
      character*4 pa,pb,pc,pd
      character*16 ptlast
      character*240 keyfile
      character*240 record
      logical, allocatable :: prttyp(:)
c
c
c     create a file with coordinates and connectivities
c
      ixyz = freeunit ()
      call prtxyz (ixyz)
c
c     output some definitions and parameters to a keyfile
c
      ikey = freeunit ()
      keyfile = filename(1:leng)//'.key'
      call version (keyfile,'new')
      open (unit=ikey,file=keyfile,status='new')
c
c     copy the contents of any previously existing keyfile
c
      do i = 1, nkey
         record = keyline(i)
         size = trimtext (record)
         write (ikey,10)  record(1:size)
   10    format (a)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (prttyp(maxtyp))
      do i = 1, maxtyp
         prttyp(i) = .true.
      end do
c
c     output any charge penetration parameters to the keyfile
c
      do i = 1, maxtyp
         prttyp(i) = .true.
      end do
      if (use_chgpen) then
         if (n .ne. 0) then
            write (ikey,20)
   20       format ()
         end if
         do i = 1, npole
            it = type(i)
            if (prttyp(it)) then
               write (ikey,30)  it,pcore(i),palpha(i)
   30          format ('chgpen',9x,i5,5x,2f11.4)
            end if
            prttyp(it) = .false.
         end do
      end if
c
c     output the polarizability parameters to the keyfile
c
      do i = 1, maxtyp
         prttyp(i) = .true.
      end do
      if (n .ne. 0) then
         write (ikey,40)
   40    format ()
      end if
      do i = 1, npole
         it = type(i)
         if (prttyp(it)) then
            k = 0
            do j = 1, maxval
               if (pgrp(j,it) .ne. 0)  k = j
            end do
            call sort8 (k,pgrp(1,it))
            if (use_thole) then
               write (ikey,50)  it,polarity(i),thole(i),
     &                           (pgrp(j,it),j=1,k)
   50          format ('polarize',7x,i5,5x,2f11.4,2x,20i5)
            else if (use_chgpen) then
               write (ikey,60)  it,polarity(i),(pgrp(j,it),j=1,k)
   60          format ('polarize',7x,i5,5x,f11.4,6x,20i7)
            end if
         end if
         prttyp(it) = .false.
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (prttyp)
      close (unit=ikey)
      return
      end
