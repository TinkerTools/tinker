c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program prtmol2  --  output Tripos MOL2 structure file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtmol2" writes out a set of coordinates in Tripos MOL2
c     format to an external disk file
c
c
      subroutine prtmol2 (imol2)
      use atoms
      use bndstr
      use files
      use iounit
      use titles
      implicit none
      integer i,j,imol2
      integer subnum
      real*8, allocatable :: atmchg(:)
      logical opened
      character*2, allocatable :: bndtyp(:)
      character*3 subnam
      character*5, allocatable :: atmtyp(:)
      character*8, allocatable :: atmnam(:)
      character*240 mol2file
c
c
c     open output unit if not already done
c
      inquire (unit=imol2,opened=opened)
      if (.not. opened) then
         mol2file = filename(1:leng)//'.mol2'
         call version (mol2file,'new')
         open (unit=imol2,file=mol2file,status='new')
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (atmnam(n))
      allocate (atmtyp(n))
      allocate (atmchg(n))
      allocate (bndtyp(nbond))
c
c     write the MOLECULE record type indicator
c
      write (imol2,10)
   10 format ('@<TRIPOS>MOLECULE')
      if (ltitle .eq. 0) then
         write (imol2,20)
   20    format ('****')
      else
         write (imol2,30)  title(1:ltitle)
   30    format (a)
      end if
      write (imol2,40)  n,nbond,1
   40 format (3i7)
      write (imol2,50)
   50 format ('SMALL')
      write (imol2,60)
   60 format ('USER_CHARGES')
c
c     determine MOL2 atom names/types and bond types
c
      call setmol2 (atmnam,atmtyp,atmchg,bndtyp)
c
c     write the ATOM record type indicator
c
      write (imol2,70)
   70 format (/,'@<TRIPOS>ATOM')
      do i = 1, n
         subnum = 1
         subnam = '<1>'
         write (imol2,80)  i,atmnam(i),x(i),y(i),z(i),atmtyp(i),
     &                     subnum,subnam,atmchg(i)
   80    format (i7,2x,a8,3f12.6,2x,a5,i4,2x,a3,f7.2)
      end do
c
c     write the BOND record type indicator
c
      write (imol2,90)
   90 format (/,'@<TRIPOS>BOND')
      do i = 1, nbond
         write (imol2,100)  i,(ibnd(j,i),j=1,2),bndtyp(i)
  100    format (3i7,2x,a2)
      end do
c
c     write the SUBSTRUCTURE record type indicator
c
      write (imol2,110)
  110 format (/,'@<TRIPOS>SUBSTRUCTURE')
      write (imol2,120)  1,'****',1,'TEMP',0,'****','****',0,'ROOT'
  120 format (i7,2x,a4,i7,2x,a4,i7,2x,a4,2x,a4,i7,2x,a4)
c
c     perform deallocation of some local arrays
c
      deallocate (atmnam)
      deallocate (atmtyp)
      deallocate (bndtyp)
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=imol2)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  program setmol2  --  set MOL2 atom, charge & bond values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setmol2" assigns MOL2 atom names/types/charges and bond types
c     based upon atomic numbers and connectivity
c
c
      subroutine setmol2 (atmnam,atmtyp,atmchg,bndtyp)
      use angbnd
      use atmlst
      use atomid
      use atoms
      use bndstr
      use couple
      use ring
      use tors
      implicit none
      integer i,j,k,m
      integer ia,ib,ic,id
      integer ka,kb,kc
      integer tenthou,thousand
      integer hundred,tens,ones
      integer atmnum,it,nlist
      integer list(12)
      real*8 atmchg(*)
      logical aromat
      logical done,proceed
      logical terma,termd
      character*1 ta,tb,tc,td
      character*1 digit(0:9)
      character*2 bndtyp(*)
      character*5 number
      character*5 atmtyp(*)
      character*8 atmnam(*)
      data digit  / '0','1','2','3','4','5','6','7','8','9' /
c
c
c     initialize atom_names, atom_types, charges and bond_types
c
      do i = 1, n
         atmnam(i) = '        '
         atmtyp(i) = '     '
         atmchg(i) = 0.0d0
         bndtyp(i) = '  '
      end do
c
c     construct the generic MOL2 atom_name for each atom
c
      do i = 1, n
         tenthou = 1 / 10000
         thousand = (i - 10000*tenthou) / 1000
         hundred = (i - 10000*tenthou - 1000*thousand) / 100
         tens = (i - 10000*tenthou - 1000*thousand - 100*hundred) / 10
         ones = i - 10000*tenthou - 1000*thousand
     &             - 100*hundred - 10*tens
         number(1:1) = digit(tenthou)
         number(2:2) = digit(thousand)
         number(3:3) = digit(hundred)
         number(4:4) = digit(tens)
         number(5:5) = digit(ones)
         if (number(1:1) .eq. '0')  number(1:1) = ' '
         if (number(2:2).eq.'0' .and. number(1:1).eq.' ') then
            number(2:2) = ' '
         end if
         if (number(3:3).eq.'0' .and. number(2:2).eq.' ') then
            number(3:3) = ' '
         end if
         if (number(4:4).eq.'0' .and. number(3:3).eq.' ') then
            number(4:4) = ' '
         end if
         atmnam(i) = name(i)//number
         do j = 1, 7
            do while (atmnam(i)(j:j) .eq. ' ')
               do k = j, 7
                  atmnam(i)(k:k) = atmnam(i)(k+1:k+1)
               end do
               atmnam(i)(8:8) = '*'
            end do
         end do
         do j = 1, 8
            if (atmnam(i)(j:j) .eq. '*')  atmnam(i)(j:j) = ' '
         end do
      end do
c
c     determine the element types based upon atom names
c
      do i = 1, n
         it = n12(i)
         atomic(i) = 0
         call upcase (name)
         if (name(i)(1:1) .eq. 'H')  atomic(i) = 1
         if (name(i)(1:2).eq.'LI' .and. it.eq.0)  atomic(i) = 3
         if (name(i).eq.'F  ' .or. name(i).eq.'F- ')  atomic(i) = 9
         if (name(i)(1:2).eq.'NA' .and. it.eq.0)  atomic(i) = 11
         if (name(i)(1:2).eq.'MG' .and. it.eq.0)  atomic(i) = 12
         if (name(i)(1:2).eq.'AL' .and. it.eq.0)  atomic(i) = 13
         if (name(i)(1:2) .eq. 'SI')  atomic(i) = 14
         if (name(i).eq.'CL ' .or. name(i).eq.'CL-')  atomic(i) = 17
         if (name(i)(1:1).eq.'K' .and. it.eq.0)  atomic(i) = 19
         if (name(i)(1:2).eq.'CA' .and. it.eq.0)  atomic(i) = 20
         if (name(i)(1:2).eq.'CR' .and. it.eq.0)  atomic(i) = 24
         if (name(i)(1:2) .eq. 'MN')  atomic(i) = 25
         if (name(i)(1:2) .eq. 'FE')  atomic(i) = 26
         if (name(i)(1:2).eq.'CO' .and. ic.eq.0)  atomic(i) = 27
         if (name(i)(1:2) .eq. 'CU')  atomic(i) = 29
         if (name(i)(1:2) .eq. 'ZN')  atomic(i) = 30
         if (name(i)(1:2) .eq. 'SE')  atomic(i) = 34
         if (name(i).eq.'BR ' .or. name(i).eq.'BR-')  atomic(i) = 35
         if (name(i)(1:2) .eq. 'MO')  atomic(i) = 42
         if (name(i)(1:2) .eq. 'SN')  atomic(i) = 50
         if (name(i).eq.'I  ' .or. name(i).eq.'I- ')  atomic(i) = 53
         if (atomic(i) .eq. 0) then
            if (name(i)(1:1) .eq. 'C')  atomic(i) = 6
            if (name(i)(1:1) .eq. 'N')  atomic(i) = 7
            if (name(i)(1:1) .eq. 'O')  atomic(i) = 8
            if (name(i)(1:1) .eq. 'S')  atomic(i) = 16
         end if
      end do
c
c     assign the generic MOL2 atom_type for each atom
c
      do i = 1, n
         it = 0
         do k = 1, n12(i)
            if (atomic(i12(k,i)) .ne. 0)  it = it + 1
         end do
         atmnum = atomic(i)
         if (atmnum .eq. 0) then
            if (name(i) .eq. 'LP ')  atmtyp(i) = 'LP   '
            if (name(i) .eq. 'DU ')  atmtyp(i) = 'Du   '
         else if (atmnum .eq. 1) then
            atmtyp(i) = 'H    '
         else if (atmnum .eq. 3) then
            atmtyp(i) = 'Li   '
         else if (atmnum .eq. 6) then
            if (it .eq. 4)  atmtyp(i) = 'C.3  '
            if (it .eq. 3)  atmtyp(i) = 'C.2  '
            if (it .eq. 2)  atmtyp(i) = 'C.1  '
         else if (atmnum .eq. 7) then
            if (it .eq. 4)  atmtyp(i) = 'N.4  '
            if (it .eq. 3)  atmtyp(i) = 'N.3  '
            if (it .eq. 2)  atmtyp(i) = 'N.2  '
            if (it .eq. 1)  atmtyp(i) = 'N.1  '
         else if (atmnum .eq. 8) then
            if (it .ge. 2)  atmtyp(i) = 'O.3  '
            if (it .eq. 1)  atmtyp(i) = 'O.2  '
         else if (atmnum .eq. 9) then
            atmtyp(i) = 'F    '
         else if (atmnum .eq. 11) then
            atmtyp(i) = 'Na   '
         else if (atmnum .eq. 12) then
            atmtyp(i) = 'Mg   '
         else if (atmnum .eq. 13) then
            atmtyp(i) = 'Al   '
         else if (atmnum .eq. 14) then
            atmtyp(i) = 'Si   '
         else if (atmnum .eq. 15) then
            atmtyp(i) = 'P.3  '
         else if (atmnum .eq. 16) then
            if (it .ge. 2)  atmtyp(i) = 'S.3  '
            if (it .le. 1)  atmtyp(i) = 'S.2  '
         else if (atmnum .eq. 17) then
            atmtyp(i) = 'Cl   '
         else if (atmnum .eq. 19) then
            atmtyp(i) = 'K    '
         else if (atmnum .eq. 20) then
            atmtyp(i) = 'Ca   '
         else if (atmnum .eq. 24) then
            atmtyp(i) = 'Cr.oh'
         else if (atmnum .eq. 25) then
            atmtyp(i) = 'Mn   '
         else if (atmnum .eq. 26) then
            atmtyp(i) = 'Fe   '
         else if (atmnum .eq. 27) then
            atmtyp(i) = 'Co.oh'
         else if (atmnum .eq. 29) then
            atmtyp(i) = 'Cu   '
         else if (atmnum .eq. 30) then
            atmtyp(i) = 'Zn   '
         else if (atmnum .eq. 35) then
            atmtyp(i) = 'Br   '
         else if (atmnum .eq. 42) then
            atmtyp(i) = 'Mo   '
         else if (atmnum .eq. 50) then
            atmtyp(i) = 'Sn   '
         else if (atmnum .eq. 53) then
            atmtyp(i) = 'I    '
         end if
      end do
c
c     handle 5-membered rings for MOL2 atom_type assignment
c
      do i = 1, nring5
         aromat = .true.
         do j = 1, 5
            k = iring5(j,i)
            if (atomic(k).ne.6 .and. atomic(k).ne.7)  aromat = .false.
            if (atomic(k).eq.6 .and. n12(k).ne.3)  aromat = .false.
            if (atomic(k).eq.7 .and. n12(k).eq.4)  aromat = .false.
            if (atomic(k).eq.6 .and. n12(i).eq.3) then
               do m = 1, n12(k)
                  kc = i12(m,k)
                  if (atomic(kc).eq.8 .and. n12(kc).eq.1) then
                     aromat = .false.
                  end if
               end do
            end if
         end do
         if (aromat) then
            do j = 1, 5
               k = iring5(j,i)
c              atmtyp(k) = atmtyp(k)(1:1)//'.ar '
               if (atomic(k).eq.7 .and. n12(k).eq.3)
     &            atmtyp(k) = 'N.pl3'
            end do
         end if
      end do
c
c     handle 6-membered rings for MOL2 atom_type assignment
c
      do i = 1, nring6
         aromat = .true.
         do j = 1, 6
            k = iring6(j,i)
            if (atomic(k).ne.6 .and. atomic(k).ne.7)  aromat = .false.
            if (atomic(k).eq.6 .and. n12(k).ne.3)  aromat = .false.
            if (atomic(k).eq.7 .and. n12(k).eq.4)  aromat = .false.
            if (atomic(k).eq.6 .and. n12(i).eq.3) then
               do m = 1, n12(k)
                  kc = i12(m,k)
                  if (atomic(kc).eq.8 .and. n12(kc).eq.1) then
                     aromat = .false.
                  end if
               end do
            end if
         end do
         if (aromat) then
            do j = 1, 6
               k = iring6(j,i)
               atmtyp(k) = atmtyp(k)(1:1)//'.ar '
            end do
         end if
      end do
c
c     handle 7-membered rings for MOL2 atom_type assignment
c
      do i = 1, nring7
         aromat = .true.
         do j = 1, 7
            m = iring7(j,i)
            if (atomic(m).ne.6 .or. n12(m).ne.3)  aromat = .false.
         end do
         if (aromat) then
            do j = 1, 7
               list(j) = iring7(j,i)
            end do
            do k = 1, nring5
               do j = 1, 5
                  list(j+7) = iring5(j,k)
               end do
               nlist = 12
               call sort8 (nlist,list)
               if (nlist .eq. 10) then
                  aromat = .true.
                  do j = 1, 5
                     m = iring5(j,i)
                     if (atomic(m).ne.6 .or. n12(m).ne.3)
     &                  aromat = .false.
                  end do
                  if (aromat) then
                     do j = 1, 7
                        atmtyp(iring7(j,i)) = 'C.ar '
                     end do
                     do j = 1, 5
                        atmtyp(iring5(j,k)) = 'C.ar '
                     end do
                  end if
               end if
            end do
         end if
      end do
c
c     handle amide nitrogens for MOL2 atom_type assignment
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if (atomic(ib).eq.6 .and. n12(ib).eq.3) then
            if (atomic(ia).eq.8 .and. n12(ia).eq.1 .and.
     &          atomic(ic).eq.7 .and. n12(ic).eq.3) then
               atmtyp(ic) = 'N.am '
            end if
            if (atomic(ic).eq.8 .and. n12(ic).eq.1 .and.
     &          atomic(ia).eq.7 .and. n12(ia).eq.3) then
               atmtyp(ia) = 'N.am '
            end if
         end if
      end do
c
c     handle guanidinium carbons for MOL2 atom_type assignment
c
      do i = 1, n
         if (atomic(i).eq.6 .and. n12(i).eq.3) then
            k = 0
            do m = 1, n12(i)
               kc = i12(m,i)
               if (atomic(kc).eq.7 .and. n12(kc).eq.3)  k = k + 1
            end do
            if (k .eq. 3)  atmtyp(i) = 'C.cat'
         end if
      end do
c
c     handle carboxylate oxygens for MOL2 atom_type assignment
c
      do i = 1, nangle
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         if ((atomic(ia).eq.8.and.n12(ia).eq.1) .and.
     &       (atomic(ic).eq.8.and.n12(ic).eq.1)) then
            if (atomic(ib).eq.6.and.n12(ib).eq.3) then
               atmtyp(ia) = 'O.co2'
               atmtyp(ic) = 'O.co2'
            end if
         end if
      end do
c
c     handle phosphate oxygens for MOL2 atom_type assignment
c
      do i = 1, n
         if (atomic(i).eq.8 .and. n12(i).eq.1) then
            if (atomic(i12(1,i)) .eq. 15)  atmtyp(i) = 'O.co2'
         end if
      end do
c
c     handle sulfoxide sulfurs for MOL2 atom_type assignment
c
      do i = 1, n
         if (atomic(i).eq.16 .and. n12(i).eq.3) then
            k = 0
            do m = 1, n12(i)
               kc = i12(m,i)
               if (atomic(kc).eq.8 .and. n12(kc).eq.1)  k = k + 1
            end do
            if (k .eq. 1)  atmtyp(i) = 'S.o  '
         end if
      end do
c
c     handle sulfone sulfurs for MOL2 atom_type assignment
c
      do i = 1, n
         if (atomic(i).eq.16 .and. n12(i).eq.4) then
            k = 0
            do m = 1, n12(i)
               kc = i12(m,i)
               if (atomic(kc).eq.8 .and. n12(kc).eq.1)  k = k + 1
            end do
            if (k .ge. 2)  atmtyp(i) = 'S.o2 '
         end if
      end do
c
c     assign the generic MOL2 bond_type for each bond
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         bndtyp(i) = '1 '
         if (atmtyp(ia)(3:3).eq.'2' .and.
     &       atmtyp(ib)(3:3).eq.'2') then
            bndtyp(i) = '2 '
         end if
         if (atmtyp(ia)(3:3).eq.'2' .and.
     &       atmtyp(ib)(3:3).eq.'1') then
            bndtyp(i) = '2 '
         end if
         if (atmtyp(ia)(3:3).eq.'1' .and.
     &       atmtyp(ib)(3:3).eq.'2') then
            bndtyp(i) = '2 '
         end if
         if (atmtyp(ia)(3:3).eq.'1' .and.
     &       atmtyp(ib)(3:3).eq.'1') then
            bndtyp(i) = '3 '
         end if
      end do
c
c     handle aromaticity for MOL2 bond_type assignment
c
      do i = 1, nbond
         ia = ibnd(1,i)
         ib = ibnd(2,i)
         if (atmtyp(ia)(3:4).eq.'ar' .and.
     &       atmtyp(ib)(3:4).eq.'ar') then
            do k = 1, nring5
               m = 0
               do j = 1, 5
                  kc = iring5(j,k)
                  if (kc.eq.ia .or. kc.eq.ib)  m = m + 1
               end do
               if (m .eq. 2)  bndtyp(i) = 'ar'
            end do
            do k = 1, nring6
               m = 0
               do j = 1, 6
                  kc = iring6(j,k)
                  if (kc.eq.ia .or. kc.eq.ib)  m = m + 1
               end do
               if (m .eq. 2)  bndtyp(i) = 'ar'
            end do
            do k = 1, nring7
               m = 0
               do j = 1, 7
                  kc = iring7(j,k)
                  if (kc.eq.ia .or. kc.eq.ib)  m = m + 1
               end do
               if (m .eq. 2)  bndtyp(i) = 'ar'
            end do
         end if
      end do
c
c     handle conjugation for MOL2 bond_type assignment
c
      done = .false.
      dowhile (.not. done)
         done = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            proceed = .false.
            do k = 1, n12(ib)
               j = bndlist(k,ib)
               ka = ibnd(1,j)
               kb = ibnd(2,j)
               if (ka.eq.ic .or. kb.eq.ic) then
                  if (bndtyp(j).eq.'2 ' .or.
     &                bndtyp(j).eq.'3 ') then
                     m = j
                     proceed = .true.
                  end if
               end if
            end do
            if (proceed) then
               ta = atmtyp(ia)(3:3)
               tb = atmtyp(ib)(3:3)
               tc = atmtyp(ic)(3:3)
               td = atmtyp(id)(3:3)
               if ((ta.eq.'2' .or. ta.eq.'1') .and.
     &             (tb.eq.'2' .or. tb.eq.'1') .and.
     &             (tc.eq.'2' .or. tc.eq.'1') .and.
     &             (td.eq.'2' .or. td.eq.'1')) then
                  terma = .true.
                  do k = 1, n12(ia)
                     j = bndlist(k,ia)
                     ka = ibnd(1,j)
                     kb = ibnd(2,j)
                     if (ka.ne.ib .and. kb.ne.ib) then
                        if (bndtyp(j).eq.'2 ' .or.
     &                      bndtyp(j).eq.'3 ') then
                           terma = .false.
                        end if
                     end if
                  end do
                  termd = .true.
                  do k = 1, n12(id)
                     j = bndlist(k,id)
                     ka = ibnd(1,j)
                     kb = ibnd(2,j)
                     if (ka.ne.ic .and. kb.ne.ic) then
                        if (bndtyp(j).eq.'2 ' .or.
     &                      bndtyp(j).eq.'3 ') then
                           termd = .false.
                        end if
                     end if
                  end do
                  if (terma .or. termd) then
                     bndtyp(m) = '1 '
                     done = .false.
                  end if
               end if
            end if
         end do
      end do
c
c     assign the generic MOL2 charge for each atom
c
      do i = 1, n
         it = 0
         do k = 1, n12(i)
            if (atomic(i12(k,i)) .ne. 0)  it = it + 1
         end do
         atmnum = atomic(i)
         if (atmnum.eq.3 .and. it.eq.0)  atmchg(i) = 1.0d0
         if (atmnum.eq.9 .and. it.eq.0)  atmchg(i) = -1.0d0
         if (atmnum.eq.11 .and. it.eq.0)  atmchg(i) = 1.0d0
         if (atmnum.eq.12 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.13 .and. it.eq.0)  atmchg(i) = 3.0d0
         if (atmnum.eq.17 .and. it.eq.0)  atmchg(i) = -1.0d0
         if (atmnum.eq.19 .and. it.eq.0)  atmchg(i) = 1.0d0
         if (atmnum.eq.20 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.24 .and. it.eq.0)  atmchg(i) = 3.0d0
         if (atmnum.eq.25 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.26 .and. it.eq.0)  atmchg(i) = 3.0d0
         if (atmnum.eq.27 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.29 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.30 .and. it.eq.0)  atmchg(i) = 2.0d0
         if (atmnum.eq.35 .and. it.eq.0)  atmchg(i) = -1.0d0
         if (atmnum.eq.42 .and. it.eq.0)  atmchg(i) = 4.0d0
         if (atmnum.eq.53 .and. it.eq.0)  atmchg(i) = -1.0d0
c
c     handle ammonium nitrogens for MOL2 charge assignment
c
         if (atmtyp(i) .eq. 'N.4  ') then
            atmchg(i) = 1.0d0
         end if
c
c     handle pyridinium nitrogens for MOL2 charge assignment
c
         if (atmtyp(i).eq.'N.ar ' .and. it.eq.3) then
            atmchg(i) = 1.0d0
         end if
c
c     handle carboxylate oxygens for MOL2 charge assignment
c
         if (atmtyp(i) .eq. 'O.co2') then
            kc = i12(1,i)
            if (atomic(kc) .eq. 6) then
               atmchg(i) = -1.0d0
               atmchg(kc) = 1.0d0
            end if
         end if
c
c     handle phosphate groups for MOL2 charge assignment
c
         if (atmtyp(i) .eq. 'P.3  ') then
            atmchg(i) = 1.0d0
            do m = 1, n12(i)
               kc = i12(m,i)
               if (atmtyp(kc) .eq. 'O.co2') then
                  atmchg(kc) = -1.0d0
               end if
            end do
         end if
c
c     handle sulfonate groups for MOL2 charge assignment
c
         if (atmtyp(i) .eq. 'S.o2 ') then
            k = 0
            do m = 1, n12(i)
               kc = i12(m,i)
               if (atmtyp(kc) .eq. 'O.2  ')  k = k + 1
            end do
            if (k .ge. 3) then
               if (k .eq. 3)  atmchg(i) = 0.5d0
               do m = 1, n12(i)
                  kc = i12(m,i)
                  if (atmtyp(kc) .eq. 'O.2  ') then
                     atmchg(kc) = -0.5d0
                  end if
               end do
            end if
         end if
      end do
      return
      end
