c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kpolar  --  assign polarizability parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kpolar" assigns atomic dipole polarizabilities to the atoms
c     within the structure and processes any new or changed values
c
c     literature references:
c
c     A. C. Simmonett, F. C. Pickard IV, J. W. Ponder and B. R. Brooks,
c     "An Empirical Extrapolation Scheme for Efficient Treatment of
c     Induced Dipoles", Journal of Chemical Physics, 145, 164101 (2016)
c     [OPT method]
c
c     F. Aviat, L. Lagardere and J.-P. Piquemal, "The Truncated
c     Conjugate Gradient (TCG), a Non-Iterative/Fixed-Cost Strategy for
c     Computing Polarization in Molecular Dynamics: Fast Evaluation of
c     Analytical Forces", Journal of Chemical Physics, 147, 161724
c     (2018)  [TCG method]
c
c
      subroutine kpolar
      use atoms
      use chgpen
      use expol
      use inform
      use iounit
      use keys
      use kpolpr
      use kpolr
      use mplpot
      use mpole
      use polar
      use polopt
      use polpot
      use polpcg
      use poltcg
      use potent
      implicit none
      integer i,j,k
      integer ii,kk
      integer ia,ib,it
      integer next,size
      integer nlist,npg
      integer number
      integer pg(maxval)
      integer, allocatable :: list(:)
      integer, allocatable :: rlist(:)
      real*8 pol,thl,thd
      real*8 sixth
      logical header
      character*4 pa,pb
      character*8 blank,pt
      character*20 keyword
      character*20 text
      character*240 record
      character*240 string
c
c
c     set the default values for polarization variables
c
      polprt = .false.
c
c     set defaults for PCG induced dipole parameters
c
      pcgprec = .true.
      pcgguess = .true.
      pcgpeek = 1.0d0
c
c     set defaults for TCG induced dipole parameters
c
      tcgorder = 0
      tcgguess = .true.
      tcgpeek = 1.0d0
      if (poltyp .eq. 'TCG   ')  poltyp = 'TCG2  '
      if (poltyp .eq. 'TCG0  ') then
         poltyp = 'DIRECT'
      else if (poltyp .eq. 'TCG1  ') then
         poltyp = 'TCG   '
         tcgorder = 1
      else if (poltyp(1:3) .eq. 'TCG') then
         poltyp = 'TCG   '
         tcgorder = 2
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(copt))  deallocate (copt)
      if (allocated(copm))  deallocate (copm)
      allocate (copt(0:maxopt))
      allocate (copm(0:maxopt))
c
c     set defaults for OPT induced dipole coefficients
c
      optorder = 0
      do i = 0, maxopt
         copt(i) = 0.0d0
         copm(i) = 0.0d0
      end do
      if (poltyp .eq. 'OPT   ')  poltyp = 'OPT4  '
      if (poltyp .eq. 'OPT1  ') then
         copt(0) = 0.530d0
         copt(1) = 0.604d0
      else if (poltyp .eq. 'OPT2  ') then
         copt(0) = 0.042d0
         copt(1) = 0.635d0
         copt(2) = 0.414d0
      else if (poltyp .eq. 'OPT3  ') then
         copt(0) = -0.132d0
         copt(1) = 0.218d0
         copt(2) = 0.637d0
         copt(3) = 0.293d0
      else if (poltyp .eq. 'OPT4  ') then
         copt(0) = -0.071d0
         copt(1) = -0.096d0
         copt(2) = 0.358d0
         copt(3) = 0.587d0
         copt(4) = 0.216d0
      else if (poltyp .eq. 'OPT5  ') then
         copt(0) = -0.005d0
         copt(1) = -0.129d0
         copt(2) = -0.026d0
         copt(3) = 0.465d0
         copt(4) = 0.528d0
         copt(5) = 0.161d0
      else if (poltyp .eq. 'OPT6  ') then
         copt(0) = 0.014d0
         copt(1) = -0.041d0
         copt(2) = -0.172d0
         copt(3) = 0.073d0
         copt(4) = 0.535d0
         copt(5) = 0.467d0
         copt(6) = 0.122d0
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
c
c     set defaults for numbers and lists of polarizable atoms
c
      nlist = 0
      do i = 1, n
         list(i) = 0
      end do
c
c     get keywords containing polarization-related options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'POLARIZABLE ') then
            read (string,*,err=10,end=10)  (list(j),j=nlist+1,n)
   10       continue
            do while (list(nlist+1) .ne. 0)
               nlist = nlist + 1
            end do
         else if (keyword(1:12) .eq. 'POLAR-PRINT ') then
            polprt = .true.
         else if (keyword(1:12) .eq. 'PCG-PRECOND ') then
            pcgprec = .true.
         else if (keyword(1:14) .eq. 'PCG-NOPRECOND ') then
            pcgprec = .false.
         else if (keyword(1:10) .eq. 'PCG-GUESS ') then
            pcgguess = .true.
         else if (keyword(1:12) .eq. 'PCG-NOGUESS ') then
            pcgguess = .false.
         else if (keyword(1:9) .eq. 'PCG-PEEK ') then
            read (string,*,err=20,end=20)  pcgpeek
         else if (keyword(1:10) .eq. 'TCG-GUESS ') then
            tcgguess = .true.
         else if (keyword(1:12) .eq. 'TCG-NOGUESS ') then
            tcgguess = .false.
         else if (keyword(1:9) .eq. 'TCG-PEEK ') then
            read (string,*,err=20,end=20)  tcgpeek
         else if (keyword(1:10) .eq. 'OPT-COEFF ') then
            do j = 0, maxopt
               copt(j) = 0.0d0
            end do
            read (string,*,err=20,end=20)  (copt(j),j=0,maxopt)
         end if
   20    continue
      end do
c
c     get maximum coefficient order for OPT induced dipoles
c
      if (poltyp(1:3) .eq. 'OPT') then
         poltyp = 'OPT   '
         do i = 1, maxopt
            if (copt(i) .ne. 0.0d0)  optorder = max(i,optorder)
         end do
         do i = 0, optorder
            do j = optorder, i, -1
               copm(i) = copm(i) + copt(j)
            end do
         end do
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(ipolar))  deallocate (ipolar)
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(thole))  deallocate (thole)
      if (allocated(tholed))  deallocate (tholed)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(udir))  deallocate (udir)
      if (allocated(udirp))  deallocate (udirp)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(douind))  deallocate (douind)
      allocate (ipolar(n))
      allocate (polarity(n))
      allocate (thole(n))
      allocate (tholed(n))
      allocate (pdamp(n))
      allocate (udir(3,n))
      allocate (udirp(3,n))
      allocate (uind(3,n))
      allocate (uinp(3,n))
      allocate (douind(n))
      if (allocated(uopt))  deallocate (uopt)
      if (allocated(uoptp))  deallocate (uoptp)
      if (allocated(fopt))  deallocate (fopt)
      if (allocated(foptp))  deallocate (foptp)
      if (poltyp .eq. 'OPT') then
         allocate (uopt(0:optorder,3,n))
         allocate (uoptp(0:optorder,3,n))
         allocate (fopt(0:optorder,10,n))
         allocate (foptp(0:optorder,10,n))
      end if
c
c     set the atoms allowed to have nonzero induced dipoles
c
      do i = 1, n
         douind(i) = .true.
      end do
      i = 1
      do while (list(i) .ne. 0)
         if (i .eq. 1) then
            do j = 1, n
               douind(j) = .false.
            end do
         end if
         if (list(i).gt.0 .and. list(i).le.n) then
            j = list(i)
            if (.not. douind(j)) then
               douind(j) = .true.
            end if
         else if (list(i).lt.0 .and. list(i).ge.-n) then
            do j = abs(list(i)), abs(list(i+1))
               if (.not. douind(j)) then
                  douind(j) = .true.
               end if
            end do
            i = i + 1
         end if
         i = i + 1
      end do
c
c     perform dynamic allocation of some local arrays
c
      deallocate (list)
c
c     process keywords containing polarizability parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:9) .eq. 'POLARIZE ') then
            k = 0
            pol = 0.0d0
            thl = -1.0d0
            thd = -1.0d0
            do j = 1, maxval
               pg(j) = 0
            end do
            call getnumb (record,k,next)
            call gettext (record,text,next)
            read (text,*,err=30,end=30)  pol
            call gettext (record,text,next)
            j = 1
            call getnumb (text,pg(1),j)
            if (pg(1) .eq. 0) then
               read (text,*,err=30,end=30)  thl
               call gettext (record,text,next)
               j = 1
               call getnumb (text,pg(1),j)
               string = record(next:240)
               if (pg(1) .eq. 0) then
                  read (text,*,err=30,end=30)  thd
                  read (string,*,err=30,end=30)  (pg(j),j=1,maxval)
               else
                  read (string,*,err=30,end=30)  (pg(j),j=2,maxval)
               end if
            else
               string = record(next:240)
               read (string,*,err=30,end=30)  (pg(j),j=2,maxval)
            end if
   30       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,40)
   40             format (/,' Additional Atomic Dipole',
     &                       ' Polarizability Parameters :')
                  if (thd .ge. 0.0d0) then
                     write (iout,50)
   50                format (/,5x,'Atom Type',11x,'Alpha',7x,
     &                          'Thole',6x,'TholeD',5x,
     &                          'Group Atom Types',/)
                  else if (thl .ge. 0.0d0) then
                     write (iout,60)
   60                format (/,5x,'Atom Type',11x,'Alpha',7x,
     &                          'Thole',5x,'Group Atom Types',/)
                  else
                     write (iout,70)
   70                format (/,5x,'Atom Type',11x,'Alpha',5x,
     &                          'Group Atom Types',/)
                  end if
               end if
               if (k .le. maxtyp) then
                  polr(k) = pol
                  athl(k) = max(0.0d0,thl)
                  dthl(k) = max(0.0d0,thd)
                  do j = 1, maxval
                     pgrp(j,k) = pg(j)
                     if (pg(j) .eq. 0) then
                        npg = j - 1
                        goto 80
                     end if
                  end do
   80             continue
                  if (.not. silent) then
                     if (thd .ge. 0.0d0) then
                        write (iout,90)  k,pol,thl,thd,(pg(j),j=1,npg)
   90                   format (4x,i8,8x,f10.3,2x,f10.3,2x,f10.3,
     &                             7x,20i5)
                     else if (thl .ge. 0.0d0) then
                        write (iout,100)  k,pol,thl,(pg(j),j=1,npg)
  100                   format (4x,i8,8x,f10.3,2x,f10.3,7x,20i5)
                     else
                        write (iout,110)  k,pol,(pg(j),j=1,npg)
  110                   format (4x,i8,8x,f10.3,7x,20i5)
                     end if
                  end if
               else
                  write (iout,120)
  120             format (/,' KPOLAR  --  Too many Dipole',
     &                       ' Polarizability Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     process keywords with specific pair polarization values
c
      blank = '        '
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'POLPAIR ') then
            ia = 0
            ib = 0
            thl = -1.0d0
            thd = -1.0d0
            string = record(next:240)
            read (string,*,err=130,end=130)  ia,ib,thl,thd
  130       continue
            if (header .and. .not.silent) then
               header = .false.
               write (iout,140)
  140          format (/,' Additional Polarization Parameters',
     &                    ' for Specific Pairs :')
               if (thd .ge. 0.0d0) then
                  write (iout,150)
  150             format (/,5x,'Atom Types',14x,'Thole',
     &                       9x,'TholeD',/)
               else if (thl .ge. 0.0d0) then
                  write (iout,160)
  160             format (/,5x,'Atom Types',14x,'Thole',/)
               end if
            end if
            if (thd.ge.0.0d0 .and. .not.silent) then
               write (iout,170)  ia,ib,thl,thd
  170          format (6x,2i4,5x,2f15.4)
            else if (thl.ge.0.0d0 .and. .not.silent) then
               write (iout,180)  ia,ib,thl
  180          format (6x,2i4,5x,f15.4)
            end if
            size = 4
            call numeral (ia,pa,size)
            call numeral (ib,pb,size)
            if (ia .le. ib) then
               pt = pa//pb
            else
               pt = pb//pa
            end if
            do k = 1, maxnpp
               if (kppr(k).eq.blank .or. kppr(k).eq.pt) then
                  kppr(k) = pt
                  thlpr(k) = max(thl,0.0d0)
                  thdpr(k) = max(thd,0.0d0)
                  goto 200
               end if
            end do
            write (iout,190)
  190       format (/,' KPOLAR  --  Too many Special Pair',
     &                 ' Thole Parameters')
            abort = .true.
  200       continue
         end if
      end do
c
c     find and store the atomic dipole polarizability parameters
c
      sixth = 1.0d0 / 6.0d0
      npolar = n
      do i = 1, n
         polarity(i) = 0.0d0
         thole(i) = 0.0d0
         tholed(i) = 0.0d0
         pdamp(i) = 0.0d0
         it = type(i)
         if (it .ne. 0) then
            polarity(i) = polr(it)
            thole(i) = athl(it)
            tholed(i) = dthl(it)
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(jpolar))  deallocate (jpolar)
      allocate (jpolar(n))
c
c     perform dynamic allocation of some local arrays
c
      allocate (list(n))
      allocate (rlist(maxtyp))
c
c     set atom type index into condensed pair Thole matrices
c
      nlist = n
      do i = 1, n
         list(i) = type(i)
         jpolar(i) = list(i)
      end do
      call sort8 (nlist,list)
      do i = 1, maxtyp
         rlist(i) = 0
      end do
      do i = 1, n
         j = jpolar(i)
         if (rlist(j) .eq. 0) then
            do k = 1, nlist
               if (list(k) .eq. j)  rlist(j) = k
            end do
         end if
      end do
      do i = 1, n
         jpolar(i) = rlist(type(i))
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(thlval))  deallocate (thlval)
      if (allocated(thdval))  deallocate (thdval)
      allocate (thlval(nlist,nlist))
      allocate (thdval(nlist,nlist))
c
c     use combination rules for pairwise Thole damping values
c
      do ii = 1, nlist
         i = list(ii)
         do kk = ii, nlist
            k = list(kk)
            thl = min(athl(i),athl(k))
            if (thl .eq. 0.0d0)  thl = max(athl(i),athl(k))
            thd = min(dthl(i),dthl(k))
            if (thd .eq. 0.0d0)  thd = max(dthl(i),dthl(k))
            thlval(ii,kk) = thl
            thlval(kk,ii) = thl
            thdval(ii,kk) = thd
            thdval(kk,ii) = thd
         end do
      end do
c
c     apply Thole damping values for special atom type pairs
c
      do i = 1, maxnpp
         if (kppr(i) .eq. blank)  goto 210
         ia = rlist(number(kppr(i)(1:4)))
         ib = rlist(number(kppr(i)(5:8)))
         if (ia.ne.0 .and. ib.ne.0) then
            thlval(ia,ib) = thlpr(i)
            thlval(ib,ia) = thlpr(i)
            thdval(ia,ib) = thdpr(i)
            thdval(ib,ia) = thdpr(i)
         end if
      end do
  210 continue
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (rlist)
c
c     setup exchange polarization via variable polarizability
c
      call kexpol
c
c     remove zero or undefined electrostatic sites from the list
c
      if ((use_polar .or. use_repel .or. use_xrepel .or. use_solv)
     &       .and. .not.use_chgtrn) then
         npole = 0
         ncp = 0
         npolar = 0
         nexpol = 0
         do i = 1, n
            if (polarity(i) .eq. 0.0d0)  douind(i) = .false.
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               mono0(i) = pole(1,i)
               if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
               if (polarity(i) .ne. 0.0d0) then
                  npolar = npolar + 1
                  ipolar(npolar) = i
                  douind(i) = .true.
               end if
               if (tholed(i) .ne. 0.0d0)  use_tholed = .true.
               if (kpep(i) .ne. 0.0d0)  nexpol = nexpol + 1
            end if
         end do
      end if
c
c     test multipoles at chiral sites and invert if necessary
c
      if (use_polar .and. .not.use_chgtrn)  call chkpole
c
c     assign polarization group connectivity of each atom
c
      call polargrp
c
c     turn off polarizable multipole potentials if not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (ncp .ne. 0)  use_chgpen = .true.
      if (npolar .eq. 0)  use_polar = .false.
      if (ncp .ne. 0)  use_thole = .false.
      if (use_tholed)  use_thole = .true.
      if (nexpol .ne. 0)  use_expol = .true.
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine polargrp  --  polarization group connectivity  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "polargrp" generates members of the polarization group of
c     each atom and separate lists of the 1-2, 1-3 and 1-4 group
c     connectivities
c
c
      subroutine polargrp
      use atoms
      use couple
      use iounit
      use kpolr
      use polgrp
      implicit none
      integer i,j,k,m
      integer it,jt
      integer jj,kk
      integer start,stop
      integer nkeep,nlist
      integer maxkeep,maxlist
      integer, allocatable :: keep(:)
      integer, allocatable :: list(:)
      integer, allocatable :: mask(:)
      logical done,abort
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(np11))  deallocate (np11)
      if (allocated(np12))  deallocate (np12)
      if (allocated(np13))  deallocate (np13)
      if (allocated(np14))  deallocate (np14)
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
      allocate (np11(n))
      allocate (np12(n))
      allocate (np13(n))
      allocate (np14(n))
      allocate (ip11(maxp11,n))
      allocate (ip12(maxp12,n))
      allocate (ip13(maxp13,n))
      allocate (ip14(maxp14,n))
c
c     initialize size and connectivity of polarization groups
c
      do i = 1, n
         np11(i) = 1
         ip11(1,i) = i
         np12(i) = 0
         np13(i) = 0
         np14(i) = 0
      end do
c
c     set termination flag and temporary group storage
c
      abort = .false.
      maxkeep = 100
      maxlist = 10000
c
c     find the directly connected group members for each atom
c
      do i = 1, n
         it = type(i)
         if (it .ne. 0) then
            do j = 1, n12(i)
               jj = i12(j,i)
               jt = type(jj)
               do k = 1, maxval
                  kk = pgrp(k,it)
                  if (kk .eq. 0)  goto 20
                  if (pgrp(k,it) .eq. jt) then
                     if (np11(i) .lt. maxp11) then
                        np11(i) = np11(i) + 1
                        ip11(np11(i),i) = jj
                     else
                        write (iout,10)
   10                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 30
                     end if
                  end if
               end do
   20          continue
            end do
         end if
      end do
   30 continue
c
c     make sure all connected group members are bidirectional
c
      do i = 1, n
         do j = 1, np11(i)
            k = ip11(j,i)
            do m = 1, np11(k)
               if (ip11(m,k) .eq. i)  goto 50
            end do
            write (iout,40)  min(i,k),max(i,k)
   40       format (/,' POLARGRP  --  Check Polarization Groups for',
     &                 ' Atoms',i9,' and',i9)
            abort = .true.
   50       continue
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (keep(maxkeep))
      allocate (list(maxlist))
      allocate (mask(n))
c
c     find any other group members for each atom in turn
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         done = .false.
         start = 1
         stop = np11(i)
         do j = start, stop
            jj = ip11(j,i)
            if (jj .lt. i) then
               done = .true.
               np11(i) = np11(jj)
               do k = 1, np11(i)
                  ip11(k,i) = ip11(k,jj)
               end do
            else
               mask(jj) = i
            end if
         end do
         do while (.not. done)
            done = .true.
            do j = start, stop
               jj = ip11(j,i)
               do k = 1, np11(jj)
                  kk = ip11(k,jj)
                  if (mask(kk) .ne. i) then
                     if (np11(i) .lt. maxp11) then
                        np11(i) = np11(i) + 1
                        ip11(np11(i),i) = kk
                     else
                        write (iout,60)
   60                   format (/,' POLARGRP  --  Too many Atoms',
     &                             ' in Polarization Group')
                        abort = .true.
                        goto 70
                     end if
                     mask(kk) = i
                  end if
               end do
            end do
            if (np11(i) .ne. stop) then
               done = .false.
               start = stop + 1
               stop = np11(i)
            end if
         end do
         call sort (np11(i),ip11(1,i))
      end do
   70 continue
      if (abort)  call fatal
c
c     loop over atoms finding all the 1-2 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         nkeep = 0
         do j = 1, np11(i)
            jj = ip11(j,i)
            do k = 1, n12(jj)
               kk = i12(k,jj)
               if (mask(kk) .ne. i) then
                  nkeep = nkeep + 1
                  keep(nkeep) = kk
               end if
            end do
         end do
         nlist = 0
         do j = 1, nkeep
            jj = keep(j)
            do k = 1, np11(jj)
               kk = ip11(k,jj)
               nlist = nlist + 1
               list(nlist) = kk
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp12) then
            np12(i) = nlist
            do j = 1, nlist
               ip12(j,i) = list(j)
            end do
         else
            write (iout,80)
   80       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-2 Polarization Group')
            abort = .true.
            goto 90
         end if
      end do
   90 continue
c
c     loop over atoms finding all the 1-3 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np12(i)
            jj = ip12(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp13) then
            np13(i) = nlist
            do j = 1, nlist
               ip13(j,i) = list(j)
            end do
         else
            write (iout,100)
  100       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-3 Polarization Group')
            abort = .true.
            goto 110
         end if
      end do
  110 continue
c
c     loop over atoms finding all the 1-4 group relationships
c
      do i = 1, n
         mask(i) = 0
      end do
      do i = 1, n
         do j = 1, np11(i)
            jj = ip11(j,i)
            mask(jj) = i
         end do
         do j = 1, np12(i)
            jj = ip12(j,i)
            mask(jj) = i
         end do
         do j = 1, np13(i)
            jj = ip13(j,i)
            mask(jj) = i
         end do
         nlist = 0
         do j = 1, np13(i)
            jj = ip13(j,i)
            do k = 1, np12(jj)
               kk = ip12(k,jj)
               if (mask(kk) .ne. i) then
                  nlist = nlist + 1
                  list(nlist) = kk
               end if
            end do
         end do
         call sort8 (nlist,list)
         if (nlist .le. maxp14) then
            np14(i) = nlist
            do j = 1, nlist
               ip14(j,i) = list(j)
            end do
         else
            write (iout,120)
  120       format (/,' POLARGRP  --  Too many Atoms',
     &                 ' in 1-4 Polarization Group')
            abort = .true.
            goto 130
         end if
      end do
  130 continue
      if (abort)  call fatal
c
c     perform deallocation of some local arrays
c
      deallocate (keep)
      deallocate (list)
      deallocate (mask)
      return
      end
