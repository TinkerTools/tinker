c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program crystal  --  fractional coordinate manipulations  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "crystal" is a utility which converts between fractional and
c     Cartesian coordinates, and can generate full unit cells from
c     asymmetric units
c
c
      program crystal
      use atoms
      use bound
      use boxes
      use iounit
      use files
      use math
      implicit none
      integer maxspace
      parameter (maxspace=84)
      integer i,ixyz,mode
      integer na,nb,nc
      integer next,freeunit
      logical exist,query
      character*1 answer
      character*8 sgroup(maxspace)
      character*240 xyzfile
      character*240 record
      character*240 string
      data sgroup / 'P1      ', 'P2      ', 'P1(-)   ', 'P21     ',
     &              'C2      ', 'Pm      ', 'Pc      ', 'Cm      ',
     &              'Cc      ', 'P2/m    ', 'P21/m   ', 'C2/m    ',
     &              'P2/c    ', 'P21/c   ', 'P21/n   ', 'P21/a   ',
     &              'C2/c    ', 'P21212  ', 'P212121 ', 'C2221   ',
     &              'Pca21   ', 'Pmn21   ', 'Pna21   ', 'Pn21a   ',
     &              'Cmc21   ', 'Aba2    ', 'Fdd2    ', 'Iba2    ',
     &              'Pnna    ', 'Pmna    ', 'Pcca    ', 'Pbam    ',
     &              'Pccn    ', 'Pbcm    ', 'Pnnm    ', 'Pmmn    ',
     &              'Pbcn    ', 'Pbca    ', 'Pnma    ', 'Cmcm    ',
     &              'Cmca    ', 'Fddd    ', 'Ibam    ', 'P41     ',
     &              'P43     ', 'I4(-)   ', 'P4/n    ', 'P42/n   ',
     &              'I4/m    ', 'I41/a   ', 'P41212  ', 'P43212  ',
     &              'P4(-)21m', 'P4(-)21c', 'P4(-)m2 ', 'I41/amd ',
     &              'P31     ', 'P32     ', 'R3      ', 'P3(-)   ',
     &              'R3(-)   ', 'P3121   ', 'P3221   ', 'R3m     ',
     &              'R3c     ', 'R3(-)c  ', 'P63     ', 'P63/m   ',
     &              'P63/mmc ', 'Pa3(-)  ', 'P43m    ', 'I4(-)3m ',
     &              'P4(-)3n ', 'I4(-)3d ', 'Pm3(-)m ', 'Pn3(-)n ',
     &              'Pm3(-)n ', 'Pn3(-)m ', 'Fm3(-)m ', 'Fm3(-)c ',
     &              'Fd3(-)m ', 'Fm3(-)c ', 'Im3(-)m ', 'Im3(-)d '/
c
c
c     get and read the Cartesian coordinates file
c
      call initial
      call getxyz
c
c     find out which unit cell manipulation to perform
c
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         if (mode.ge.1 .and. mode.le.5)  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Tinker Crystal Structure Utility Can :',
     &           //,4x,'(1) Convert Fractional to Cartesian Coords',
     &           /,4x,'(2) Convert Cartesian to Fractional Coords',
     &           /,4x,'(3) Move Any Stray Molecules into Unit Cell',
     &           /,4x,'(4) Make a Unit Cell from an Asymmetric Unit',
     &           /,4x,'(5) Make a Big Block from a Single Unit Cell')
         do while (mode.lt.1 .or. mode.gt.5)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     get any cell dimensions found in the keyword list
c
      call unitcell
c
c     determine the space group if it will be needed later
c
      if (mode .eq. 4) then
         do i = 1, maxspace
            if (spacegrp .eq. sgroup(i))  goto 120
         end do
   60    continue
         write (iout,70)  (sgroup(i),i=1,32)
   70    format (/,' Available Crystallographic Space Groups :',/,
     &           /,3x,'(1) ',a8,5x,'(2) ',a8,5x,'(3) ',a8,
     &              5x,'(4) ',a8,
     &           /,3x,'(5) ',a8,5x,'(6) ',a8,5x,'(7) ',a8,
     &              5x,'(8) ',a8,
     &           /,3x,'(9) ',a8,4x,'(10) ',a8,4x,'(11) ',a8,
     &              4x,'(12) ',a8,
     &           /,2x,'(13) ',a8,4x,'(14) ',a8,4x,'(15) ',a8,
     &              4x,'(16) ',a8,
     &           /,2x,'(17) ',a8,4x,'(18) ',a8,4x,'(19) ',a8,
     &              4x,'(20) ',a8,
     &           /,2x,'(21) ',a8,4x,'(22) ',a8,4x,'(23) ',a8,
     &              4x,'(24) ',a8,
     &           /,2x,'(25) ',a8,4x,'(26) ',a8,4x,'(27) ',a8,
     &              4x,'(28) ',a8,
     &           /,2x,'(29) ',a8,4x,'(30) ',a8,4x,'(31) ',a8,
     &              4x,'(32) ',a8)
         write (iout,80)  (sgroup(i),i=33,64)
   80    format (2x,'(33) ',a8,4x,'(34) ',a8,4x,'(35) ',a8,
     &              4x,'(36) ',a8,
     &           /,2x,'(37) ',a8,4x,'(38) ',a8,4x,'(39) ',a8,
     &              4x,'(40) ',a8,
     &           /,2x,'(41) ',a8,4x,'(42) ',a8,4x,'(43) ',a8,
     &              4x,'(44) ',a8,
     &           /,2x,'(45) ',a8,4x,'(46) ',a8,4x,'(47) ',a8,
     &              4x,'(48) ',a8,
     &           /,2x,'(49) ',a8,4x,'(50) ',a8,4x,'(51) ',a8,
     &              4x,'(52) ',a8,
     &           /,2x,'(53) ',a8,4x,'(54) ',a8,4x,'(55) ',a8,
     &              4x,'(56) ',a8,
     &           /,2x,'(57) ',a8,4x,'(58) ',a8,4x,'(59) ',a8,
     &              4x,'(60) ',a8,
     &           /,2x,'(61) ',a8,4x,'(62) ',a8,4x,'(63) ',a8,
     &              4x,'(64) ',a8)
         write (iout,90)  (sgroup(i),i=65,maxspace)
   90    format (2x,'(76) ',a8,4x,'(66) ',a8,4x,'(67) ',a8,
     &              4x,'(68) ',a8,
     &           /,2x,'(69) ',a8,4x,'(70) ',a8,4x,'(71) ',a8,
     &              4x,'(72) ',a8,
     &           /,2x,'(73) ',a8,4x,'(74) ',a8,4x,'(75) ',a8,
     &              4x,'(76) ',a8,
     &           /,2x,'(77) ',a8,4x,'(78) ',a8,4x,'(79) ',a8,
     &              4x,'(80) ',a8,
     &           /,2x,'(81) ',a8,4x,'(82) ',a8,4x,'(83) ',a8,
     &              4x,'(84) ',a8)
         write (iout,100)
  100    format (/,' Enter the Number of the Desired Choice :  ',$)
         read (input,110)  i
  110    format (i10)
         if (i.lt.1 .or. i.gt.maxspace)  goto 60
         spacegrp = sgroup(i)
  120    continue
      end if
c
c     if not in keyfile, get the unit cell axis lengths
c
      do while (xbox .eq. 0.0d0)
         write (iout,130)
  130    format (/,' Enter Unit Cell Axis Lengths :  ',$)
         read (input,140)  record
  140    format (a240)
         read (record,*,err=150,end=150)  xbox,ybox,zbox
  150    continue
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
         use_bounds = .true.
      end do
c
c     if not in keyfile, get the unit cell angle values
c
      do while (alpha .eq. 0.0d0)
         write (iout,160)
  160    format (/,' Enter Unit Cell Axis Angles :   ',$)
         read (input,170)  record
  170    format (a240)
         read (record,*,err=180,end=180)  alpha,beta,gamma
  180    continue
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = alpha
         if (gamma .eq. 0.0d0)  gamma = alpha
         if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &          .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end do
c
c     find constants for coordinate interconversion
c
      call lattice
c
c     print out the initial cell dimensions to be used
c
      write (iout,190)  xbox,ybox,zbox,alpha,beta,gamma
  190 format (/,' Unit Cell Dimensions :       a     =',f12.4,
     &        /,'                              b     =',f12.4,
     &        /,'                              c     =',f12.4,
     &        /,'                             Alpha  =',f12.4,
     &        /,'                             Beta   =',f12.4,
     &        /,'                             Gamma  =',f12.4)
c
c     convert Cartesian to fractional coordinates
c
      if (mode.ne.1 .and. mode.ne.3) then
         do i = 1, n
            z(i) = (z(i)/gamma_term) / zbox
            y(i) = ((y(i)-z(i)*zbox*beta_term)/gamma_sin) / ybox
            x(i) = (x(i)-y(i)*ybox*gamma_cos-z(i)*zbox*beta_cos) / xbox
         end do
      end if
c
c     apply the appropriate space group symmetry operators
c
      if (mode .eq. 4) then
         write (iout,200)  spacegrp
  200    format (/,' Space Group Symbol :',8x,a8)
         call symmetry (spacegrp)
      end if
c
c     replicate the unit cell to make a block of unit cells
c
      if (mode .eq. 5) then
         na = 0
         nb = 0
         nc = 0
         write (iout,210)
  210    format (/,' Enter Number of Replicates along a-, b- and',
     &              ' c-Axes [1 1 1] :   ',$)
         read (input,220)  record
  220    format (a240)
         read (record,*,err=230,end=230)  na,nb,nc
  230    continue
         if (na .eq. 0)  na = 1
         if (nb .eq. 0)  nb = na
         if (nc .eq. 0)  nc = na
         if (na*nb*nc*n .gt. maxatm) then
            write (iout,240)  maxatm
  240       format (/,' CRYSTAL  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
         call bigblock (na,nb,nc)
         write (iout,250)  na,nb,nc,xbox,ybox,zbox
  250    format (/,' Dimensions of the',i3,' x',i3,' x',i3,
     &              ' Cell Block :',
     &           //,' New Cell Dimensions :       a    =',f10.4,
     &            /,'                             b    =',f10.4,
     &            /,'                             c    =',f10.4)
      end if
c
c     convert fractional to Cartesian coordinates
c
      if (mode.ne.2 .and. mode.ne.3) then
         do i = 1, n
            x(i) = x(i)*xbox + y(i)*ybox*gamma_cos + z(i)*zbox*beta_cos
            y(i) = y(i)*ybox*gamma_sin + z(i)*zbox*beta_term
            z(i) = z(i)*zbox*gamma_term
         end do
      end if
c
c     merge fragments to form complete connected molecules
c
      if (mode .eq. 4) then
         answer = ' '
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,260)
  260       format (/,' Attempt to Merge Fragments to Form Full',
     &                 ' Molecules [N] :   ',$)
            read (input,270)  record
  270       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  call molmerge
      end if
c
c     translate any stray molecules back into the unit cell
c
      if (mode .eq. 3) then
         call field
         call katom
         call molecule
         call bounds
      else if (mode .eq. 4) then
         answer = ' '
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,280)
  280       format (/,' Move Any Stray Molecules into Unit Cell',
     &                 ' [N] :   ',$)
            read (input,290)  record
  290       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y') then
            call field
            call katom
            call molecule
            call bounds
         end if
      end if
c
c     optionally move unit cell center to coordinate origin
c
      if (mode .eq. 4) then
         answer = ' '
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,300)
  300       format (/,' Locate Center of Unit Cell at Coordinate',
     &                 ' Origin [N] :   ',$)
            read (input,310)  record
  310       format (a240)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y') then
            do i = 1, n
               z(i) = (z(i)/gamma_term) / zbox
               y(i) = ((y(i)-z(i)*zbox*beta_term)/gamma_sin) / ybox
               x(i) = (x(i)-y(i)*ybox*gamma_cos-z(i)*zbox*beta_cos)
     &                                  / xbox
            end do
            do i = 1, n
               x(i) = x(i) - 0.5d0
               y(i) = y(i) - 0.5d0
               z(i) = z(i) - 0.5d0
            end do
            do i = 1, n
               x(i) = x(i)*xbox + y(i)*ybox*gamma_cos
     &                   + z(i)*zbox*beta_cos
               y(i) = y(i)*ybox*gamma_sin + z(i)*zbox*beta_term
               z(i) = z(i)*zbox*gamma_term
            end do
         end if
      end if
c
c     write out the new coordinates to a file
c
      ixyz = freeunit ()
      if (mode .eq. 2) then
         xyzfile = filename(1:leng)//'.frac'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      else
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
      call prtxyz (ixyz)
      close (unit=ixyz)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine molmerge  --  connect fragments into molecules  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "molmerge" connects fragments and removes duplicate atoms
c     during generation of a unit cell from an asymmetric unit
c
c
      subroutine molmerge
      use atomid
      use atoms
      use couple
      use molcul
      implicit none
      integer i,j,k,m,h
      integer im,km,in,kn
      real*8 r,eps
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      logical ih,kh
      logical merge,join
      logical, allocatable :: omit(:)
      logical, allocatable :: hydro(:)
c
c
c     parse the system to find molecules and fragments
c
      call molecule
c
c     perform dynamic allocation of some local arrays
c
      allocate (omit(n))
c
c     zero out the list of atoms to be deleted
c
      do i = 1, n
         omit(i) = .false.
      end do
c
c     first pass tests all pairs for duplicate atoms
c
      eps = 0.01d0
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            km = molcule(k)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            merge = .false.
            if (r .lt. eps)  merge = .true.
c
c    translate molecular fragment to the closest image
c
            if (merge) then
               xr = xr - x(k) + xi
               yr = yr - y(k) + yi
               zr = zr - z(k) + zi
               do j = imol(1,km), imol(2,km)
                  m = kmol(j)
                  x(m) = x(m) + xr
                  y(m) = y(m) + yr
                  z(m) = z(m) + zr
               end do
c
c    connections between partially duplicated fragments
c
               omit(k) = .true.
               do j = 1, n12(k)
                  m = i12(j,k)
                  join = .true.
                  do h = 1, n12(m)
                     if (i12(h,m) .eq. i)  join = .false.
                  end do
                  if (join) then
                     n12(m) = n12(m) + 1
                     i12(n12(m),m) = i
                  end if
                  join = .true.
                  do h = 1, n12(i)
                     if (i12(h,i) .eq. m)  join = .false.
                  end do
                  if (join) then
                     n12(i) = n12(i) + 1
                     i12(n12(i),i) = m
                  end if
               end do
            end if
         end do
      end do
c
c     delete any duplicated atoms identical by symmetry
c
      j = n
      do i = j, 1, -1
         if (omit(i))  call delete (i)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (omit)
c
c     parse the system to find molecules and fragments
c
      call molecule
c
c     perform dynamic allocation of some local arrays
c
      allocate (hydro(n))
c
c     find hydrogen atoms for use in connectivity assignment
c
      do i = 1, n
         hydro(i) = .false.
         if (atomic(i) .eq. 1)  hydro(i) = .true.
         if (name(i)(1:1) .eq. 'H')  hydro(i) = .true.
      end do
c
c     second pass tests all pairs for atoms to be bonded
c
      do i = 1, n-1
         im = molcule(i)
         in = n12(i)
         ih = hydro(i)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         do k = i+1, n
            km = molcule(k)
            kn = n12(k)
            kh = hydro(k)
            xr = x(k) - xi
            yr = y(k) - yi
            zr = z(k) - zi
            call image (xr,yr,zr)
            r = sqrt(xr*xr + yr*yr + zr*zr)
            merge = .false.
            if (im .ne. km) then
               eps = 2.0d0
               if (ih .or. kh)  eps = 1.6d0
               if (ih .and. kh)  eps = 0.0d0
               if (r .lt. eps)  merge = .true.
            end if
c
c    translate molecular fragment to the closest image
c
            if (merge) then
               xr = xr - x(k) + xi
               yr = yr - y(k) + yi
               zr = zr - z(k) + zi
               do j = imol(1,km), imol(2,km)
                  m = kmol(j)
                  x(m) = x(m) + xr
                  y(m) = y(m) + yr
                  z(m) = z(m) + zr
               end do
c
c     connection between bonded atoms in different fragments
c
               n12(i) = n12(i) + 1
               i12(n12(i),i) = k
               n12(k) = n12(k) + 1
               i12(n12(k),k) = i
            end if
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (hydro)
c
c     sort the connected atom lists into ascending order
c
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine cellatom  --  add new atom to the unit cell  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "cellatom" completes the addition of a symmetry related atom
c     to a unit cell by updating the atom type and attachment arrays
c
c
      subroutine cellatom (jj,j)
      use atomid
      use atoms
      use couple
      implicit none
      integer i,j,jj,delta
c
c
c     attachments of replicated atom are analogous to base atom
c
      delta = jj - j
      n12(jj) = n12(j)
      do i = 1, n12(j)
         i12(i,jj) = i12(i,j) + delta
      end do
      type(jj) = type(j)
      name(jj) = name(j)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine bigblock  --  create a block of unit cells  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "bigblock" replicates the coordinates of a single unit cell
c     to give a larger unit cell as a block of repeated units
c
c
      subroutine bigblock (na,nb,nc)
      use atoms
      use boxes
      implicit none
      integer i,j,k
      integer ii,jj,nsym
      integer na,nb,nc
      real*8, allocatable :: trans(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      nsym = na * nb * nc
      allocate (trans(3,nsym))
c
c     construct translation offsets for the replicated cells
c
      nsym = 0
      do i = (1-na)/2, na/2
         do j = (1-nb)/2, nb/2
            do k = (1-nc)/2, nc/2
               nsym = nsym + 1
               trans(1,nsym) = i
               trans(2,nsym) = j
               trans(3,nsym) = k
            end do
         end do
      end do
c
c     put the original cell at the top of the replica list
c
      do i = 1, nsym
         if (trans(1,i).eq.0 .and. trans(2,i).eq.0
     &           .and. trans(3,i).eq.0)  k = i
      end do
      do i = k, 2, -1
         trans(1,i) = trans(1,i-1)
         trans(2,i) = trans(2,i-1)
         trans(3,i) = trans(3,i-1)
      end do
      trans(1,1) = 0
      trans(2,1) = 0
      trans(3,1) = 0
c
c     translate the original unit cell to make a block of cells
c
      do i = 2, nsym
         ii = (i-1) * n
         do j = 1, n
            jj = j + ii
            x(jj) = x(j) + trans(1,i)
            y(jj) = y(j) + trans(2,i)
            z(jj) = z(j) + trans(3,i)
            call cellatom (jj,j)
         end do
      end do
      n = nsym * n
c
c     update the cell dimensions and fractional coordinates
c
      xbox = xbox * dble(na)
      ybox = ybox * dble(nb)
      zbox = zbox * dble(nc)
      do i = 1, n
         x(i) = x(i) / dble(na)
         y(i) = y(i) / dble(nb)
         z(i) = z(i) / dble(nc)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (trans)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine symmetry  --  apply space group symmetry  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "symmetry" applies symmetry operators to the fractional
c     coordinates of the asymmetric unit in order to generate
c     the symmetry related atoms of the full unit cell
c
c
      subroutine symmetry (spacegrp)
      use atoms
      implicit none
      integer i,j,k
      integer ii,jj,kk
      integer nsym,noff
      real*8 one3,two3
      real*8 xoff,yoff,zoff
      character*10 spacegrp
c
c
c     P1 space group  (International Tables 1)
c
      if (spacegrp .eq. 'P1       ') then
         nsym = 1
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P1(-) space group  (International Tables 2)
c
      else if (spacegrp .eq. 'P1(-)   ') then
         nsym = 2
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P2 space group  (International Tables 3)
c
      else if (spacegrp .eq. 'P2      ') then
         nsym = 2
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21 space group  (International Tables 4)
c
      else if (spacegrp .eq. 'P21     ') then
         nsym = 2
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     C2 space group  (International Tables 5)
c
      else if (spacegrp .eq. 'C2      ') then
         nsym = 2
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Pm space group  (International Tables 6)
c
      else if (spacegrp .eq. 'Pm      ') then
         nsym = 2
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pc space group  (International Tables 7)
c
      else if (spacegrp .eq. 'Pc      ') then
         nsym = 2
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Cm space group  (International Tables 8)
c
      else if (spacegrp .eq. 'Cm      ') then
         nsym = 2
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Cc space group  (International Tables 9)
c
      else if (spacegrp .eq. 'Cc      ') then
         nsym = 2
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + 0.5d0 + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P2/m space group  (International Tables 10)
c
      else if (spacegrp .eq. 'P2/m    ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/m space group  (International Tables 11)
c
      else if (spacegrp .eq. 'P21/m   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     C2/m space group  (International Tables 12)
c
      else if (spacegrp .eq. 'C2/m    ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P2/c space group  (International Tables 13)
c
      else if (spacegrp .eq. 'P2/c    ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/c space group  (International Tables 14)
c
      else if (spacegrp .eq. 'P21/c   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/n space group  (International Tables 14)
c
      else if (spacegrp .eq. 'P21/n   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P21/a space group  (International Tables 14)
c
      else if (spacegrp .eq. 'P21/a   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     C2/c space group  (International Tables 15)
c
      else if (spacegrp .eq. 'C2/c    ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P21212 space group  (International Tables 18)
c
      else if (spacegrp .eq. 'P21212  ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P212121 space group  (International Tables 19)
c
      else if (spacegrp .eq. 'P212121 ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     C2221 space group  (International Tables 20)
c
      else if (spacegrp .eq. 'C2221   ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Pca21 space group  (International Tables 29)
c
      else if (spacegrp .eq. 'Pca21   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pmn21 space group  (International Tables 31)
c
      else if (spacegrp .eq. 'Pmn21   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pna21 space group  (International Tables 33)
c
      else if (spacegrp .eq. 'Pna21   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pn21a space group  (International Tables 33)
c
      else if (spacegrp .eq. 'Pn21a   ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Cmc21 space group  (International Tables 36)
c
      else if (spacegrp .eq. 'Cmc21   ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Aba2 space group  (International Tables 41)
c
      else if (spacegrp .eq. 'Aba2    ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Fdd2 space group  (International Tables 43)
c
      else if (spacegrp .eq. 'Fdd2    ') then
         nsym = 4
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Iba2 space group  (International Tables 45)
c
      else if (spacegrp .eq. 'Iba2    ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Pnna space group  (International Tables 52)
c
      else if (spacegrp .eq. 'Pnna    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pmna space group  (International Tables 53)
c
      else if (spacegrp .eq. 'Pmna    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pcca space group  (International Tables 54)
c
      else if (spacegrp .eq. 'Pcca    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbam space group  (International Tables 55)
c
      else if (spacegrp .eq. 'Pbam    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pccn space group  (International Tables 56)
c
      else if (spacegrp .eq. 'Pccn    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbcm space group  (International Tables 57)
c
      else if (spacegrp .eq. 'Pbcm    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pnnm space group  (International Tables 58)
c
      else if (spacegrp .eq. 'Pnnm    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pmmn space group  (International Tables 59)
c
      else if (spacegrp .eq. 'Pmmn    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbcn space group  (International Tables 60)
c
      else if (spacegrp .eq. 'Pbcn    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pbca space group  (International Tables 61)
c
      else if (spacegrp .eq. 'Pbca    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pnma space group  (International Tables 62)
c
      else if (spacegrp .eq. 'Pnma    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Cmcm space group  (International Tables 63)
c
      else if (spacegrp .eq. 'Cmcm    ') then
         nsym = 8
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Cmca space group  (International Tables 64)
c
      else if (spacegrp .eq. 'Cmca    ') then
         nsym = 8
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Fddd space group  (International Tables 70)
c
      else if (spacegrp .eq. 'Fddd    ') then
         nsym = 8
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Ibam space group  (International Tables 72)
c
      else if (spacegrp .eq. 'Ibam    ') then
         nsym = 8
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P41 space group  (International Tables 76)
c
      else if (spacegrp .eq. 'P41     ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.75d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P43 space group  (International Tables 78)
c
      else if (spacegrp .eq. 'P43     ') then
         nsym = 4
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = 0.75d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.25d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I4(-) space group  (International Tables 82)
c
      else if (spacegrp .eq. 'I4(-)   ') then
         nsym = 4
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P4/n space group  (International Tables 85)
c
      else if (spacegrp .eq. 'P4/n    ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P42/n space group  (International Tables 86)
c
      else if (spacegrp .eq. 'P42/n   ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I4/m space group  (International Tables 87)
c
      else if (spacegrp .eq. 'I4/m    ') then
         nsym = 8
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     I41/a space group  (International Tables 88)
c
      else if (spacegrp .eq. 'I41/a   ') then
         nsym = 8
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P41212 space group  (International Tables 92)
c
      else if (spacegrp .eq. 'P41212  ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.75d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.25d0 - z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.75d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P43212 space group  (International Tables 96)
c
      else if (spacegrp .eq. 'P43212  ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.75d0 + z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.25d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.75d0 - z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.25d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)21m space group  (International Tables 113)
c
      else if (spacegrp .eq. 'P4(-)21m') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)21c space group  (International Tables 114)
c
      else if (spacegrp .eq. 'P4(-)21c') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)m2 space group  (International Tables 115)
c
      else if (spacegrp .eq. 'P4(-)m2 ') then
         nsym = 8
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 5) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I41/amd space group  (Intl. Tables 141, origin at center)
c
      else if (spacegrp .eq. 'I41/amd ') then
         nsym = 16
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P31 space group  (International Tables 144)
c
      else if (spacegrp .eq. 'P31     ') then
         nsym = 3
         noff = 1
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = one3 + z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = two3 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P32 space group  (International Tables 145)
c
      else if (spacegrp .eq. 'P32     ') then
         nsym = 3
         noff = 1
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = two3 + z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = one3 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     R3 space group  (International Tables 146)
c
      else if (spacegrp .eq. 'R3      ') then
         nsym = 3
         noff = 3
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = two3
                  yoff = one3
                  zoff = one3
               else if (k .eq. 3) then
                  xoff = one3
                  yoff = two3
                  zoff = two3
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P3(-) space group  (International Tables 147)
c
      else if (spacegrp .eq. 'P3(-)   ') then
         nsym = 6
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = -z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     R3(-) space group  (International Tables 148)
c
      else if (spacegrp .eq. 'R3(-)   ') then
         nsym = 6
         noff = 3
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = two3
                  yoff = one3
                  zoff = one3
               else if (k .eq. 3) then
                  xoff = one3
                  yoff = two3
                  zoff = two3
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = y(j) + xoff
                     y(jj) = y(j) - x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) - y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P3121 space group  (International Tables 152)
c
      else if (spacegrp .eq. 'P3121   ') then
         nsym = 6
         noff = 1
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = one3 + z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = two3 + z(j)
               else if (i .eq. 4) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = x(j) - y(j)
                  y(jj) = -y(j)
                  z(jj) = two3 - z(j)
               else if (i .eq. 6) then
                  x(jj) = -x(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = one3 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P3221 space group  (International Tables 154)
c
      else if (spacegrp .eq. 'P3221   ') then
         nsym = 6
         noff = 1
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = two3 + z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = one3 + z(j)
               else if (i .eq. 4) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = x(j) - y(j)
                  y(jj) = -y(j)
                  z(jj) = one3 - z(j)
               else if (i .eq. 6) then
                  x(jj) = -x(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = two3 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     R3m space group  (International Tables 160)
c
      else if (spacegrp .eq. 'R3m     ') then
         nsym = 6
         noff = 3
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = two3
                  yoff = one3
                  zoff = one3
               else if (k .eq. 3) then
                  xoff = one3
                  yoff = two3
                  zoff = two3
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     R3c space group  (International Tables 161)
c
      else if (spacegrp .eq. 'R3c     ') then
         nsym = 6
         noff = 3
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = two3
                  yoff = one3
                  zoff = one3
               else if (k .eq. 3) then
                  xoff = one3
                  yoff = two3
                  zoff = two3
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = x(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     R3(-)c space group  (International Tables 167)
c
      else if (spacegrp .eq. 'R3(-)c  ') then
         nsym = 12
         noff = 3
         one3 = 1.0d0 / 3.0d0
         two3 = 2.0d0 / 3.0d0
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = two3
                  yoff = one3
                  zoff = one3
               else if (k .eq. 3) then
                  xoff = one3
                  yoff = two3
                  zoff = two3
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = x(j) - y(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) - x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = y(j) + xoff
                     y(jj) = y(j) - x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = x(j) - y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = y(j) - x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = x(j) + xoff
                     y(jj) = x(j) - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P63 space group  (International Tables 173)
c
      else if (spacegrp .eq. 'P63     ') then
         nsym = 6
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P63/m space group  (International Tables 176)
c
      else if (spacegrp .eq. 'P63/m   ') then
         nsym = 12
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 7) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = -z(j)
               else if (i .eq. 9) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 10) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 11) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 12) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 - z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P6(3)/mmc space group  (Intl. Tables 194, Hexagonal Close Packed)
c
      else if (spacegrp .eq. 'P63/mmc ') then
         nsym = 24
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 4) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 5) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 6) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 7) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 8) then
                  x(jj) = x(j) - y(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 9) then
                  x(jj) = -x(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = -z(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 11) then
                  x(jj) = y(j) - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 12) then
                  x(jj) = x(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 13) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 14) then
                  x(jj) = y(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = -z(j)
               else if (i .eq. 15) then
                  x(jj) = x(j) - y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 16) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 17) then
                  x(jj) = -y(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 18) then
                  x(jj) = y(j) - x(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 19) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 20) then
                  x(jj) = y(j) - x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 21) then
                  x(jj) = x(j)
                  y(jj) = x(j) - y(j)
                  z(jj) = z(j)
               else if (i .eq. 22) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 23) then
                  x(jj) = x(j) - y(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 24) then
                  x(jj) = -x(j)
                  y(jj) = y(j) - x(j)
                  z(jj) = 0.5d0 + z(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pa3(-) space group  (International Tables 205)
c
      else if (spacegrp .eq. 'Pa3(-)  ') then
         nsym = 24
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = -y(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 8) then
                  x(jj) = -z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 11) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = -x(j)
               else if (i .eq. 12) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = -z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 13) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 14) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 15) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 16) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = z(j)
               else if (i .eq. 17) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 18) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = y(j)
               else if (i .eq. 19) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 20) then
                  x(jj) = z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 21) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 22) then
                  x(jj) = y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 23) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = x(j)
               else if (i .eq. 24) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = z(j)
                  z(jj) = 0.5d0 - x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     P4(-)3m space group  (International Tables 215)
c
      else if (spacegrp .eq. 'P4(-)3m ') then
         nsym = 24
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 7) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 12) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = z(j)
               else if (i .eq. 14) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 15) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 16) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 17) then
                  x(jj) = x(j)
                  y(jj) = z(j)
                  z(jj) = y(j)
               else if (i .eq. 18) then
                  x(jj) = -x(j)
                  y(jj) = z(j)
                  z(jj) = -y(j)
               else if (i .eq. 19) then
                  x(jj) = -x(j)
                  y(jj) = -z(j)
                  z(jj) = y(j)
               else if (i .eq. 20) then
                  x(jj) = x(j)
                  y(jj) = -z(j)
                  z(jj) = -y(j)
               else if (i .eq. 21) then
                  x(jj) = z(j)
                  y(jj) = y(j)
                  z(jj) = x(j)
               else if (i .eq. 22) then
                  x(jj) = z(j)
                  y(jj) = -y(j)
                  z(jj) = -x(j)
               else if (i .eq. 23) then
                  x(jj) = -z(j)
                  y(jj) = y(j)
                  z(jj) = -x(j)
               else if (i .eq. 24) then
                  x(jj) = -z(j)
                  y(jj) = -y(j)
                  z(jj) = x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I4(-)3m space group  (Intl. Tables 217, Body Centered Cubic)
c
      else if (spacegrp .eq. 'I4(-)3m ') then
         nsym = 24
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = -x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = -z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     P4(-)3n space group  (International Tables 218)
c
      else if (spacegrp .eq. 'P4(-)3n ') then
         nsym = 24
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 7) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 12) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 14) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 15) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 16) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 17) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 18) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 19) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 20) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 21) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 22) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 23) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 24) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     I4(-)3d space group  (International Tables 220)
c
      else if (spacegrp .eq. 'I4(-)3d ') then
         nsym = 24
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Pm3(-)m space group  (International Tables 221)
c
      else if (spacegrp .eq. 'Pm3(-)m ') then
         nsym = 48
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 7) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 12) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 14) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 15) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 16) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = z(j)
               else if (i .eq. 17) then
                  x(jj) = x(j)
                  y(jj) = z(j)
                  z(jj) = -y(j)
               else if (i .eq. 18) then
                  x(jj) = -x(j)
                  y(jj) = z(j)
                  z(jj) = y(j)
               else if (i .eq. 19) then
                  x(jj) = -x(j)
                  y(jj) = -z(j)
                  z(jj) = -y(j)
               else if (i .eq. 20) then
                  x(jj) = x(j)
                  y(jj) = -z(j)
                  z(jj) = y(j)
               else if (i .eq. 21) then
                  x(jj) = z(j)
                  y(jj) = y(j)
                  z(jj) = -x(j)
               else if (i .eq. 22) then
                  x(jj) = z(j)
                  y(jj) = -y(j)
                  z(jj) = x(j)
               else if (i .eq. 23) then
                  x(jj) = -z(j)
                  y(jj) = y(j)
                  z(jj) = x(j)
               else if (i .eq. 24) then
                  x(jj) = -z(j)
                  y(jj) = -y(j)
                  z(jj) = -x(j)
               else if (i .eq. 25) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 26) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 27) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 28) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 29) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 30) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 31) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 32) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 33) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 34) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 35) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 36) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 37) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = z(j)
               else if (i .eq. 38) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = z(j)
               else if (i .eq. 39) then
                  x(jj) = -y(j)
                  y(jj) = x(j)
                  z(jj) = -z(j)
               else if (i .eq. 40) then
                  x(jj) = y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 41) then
                  x(jj) = -x(j)
                  y(jj) = -z(j)
                  z(jj) = y(j)
               else if (i .eq. 42) then
                  x(jj) = x(j)
                  y(jj) = -z(j)
                  z(jj) = -y(j)
               else if (i .eq. 43) then
                  x(jj) = x(j)
                  y(jj) = z(j)
                  z(jj) = y(j)
               else if (i .eq. 44) then
                  x(jj) = -x(j)
                  y(jj) = z(j)
                  z(jj) = -y(j)
               else if (i .eq. 45) then
                  x(jj) = -z(j)
                  y(jj) = -y(j)
                  z(jj) = x(j)
               else if (i .eq. 46) then
                  x(jj) = -z(j)
                  y(jj) = y(j)
                  z(jj) = -x(j)
               else if (i .eq. 47) then
                  x(jj) = z(j)
                  y(jj) = -y(j)
                  z(jj) = -x(j)
               else if (i .eq. 48) then
                  x(jj) = z(j)
                  y(jj) = y(j)
                  z(jj) = x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pn3(-)n space group  (Intl. Tables 222, origin at center)
c
      else if (spacegrp .eq. 'Pn3(-)n ') then
         nsym = 48
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = z(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 12) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 14) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 15) then
                  x(jj) = y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = z(j)
               else if (i .eq. 16) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = x(j)
                  z(jj) = z(j)
               else if (i .eq. 17) then
                  x(jj) = x(j)
                  y(jj) = z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 18) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = z(j)
                  z(jj) = y(j)
               else if (i .eq. 19) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 20) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = y(j)
               else if (i .eq. 21) then
                  x(jj) = z(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 22) then
                  x(jj) = z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = x(j)
               else if (i .eq. 23) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = y(j)
                  z(jj) = x(j)
               else if (i .eq. 24) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 25) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 26) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 27) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 28) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 29) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 30) then
                  x(jj) = -z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 31) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = -y(j)
               else if (i .eq. 32) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 33) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 34) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = -z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 35) then
                  x(jj) = -y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 36) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = -x(j)
               else if (i .eq. 37) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 38) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 39) then
                  x(jj) = -y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = -z(j)
               else if (i .eq. 40) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 41) then
                  x(jj) = -x(j)
                  y(jj) = -z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 42) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -z(j)
                  z(jj) = -y(j)
               else if (i .eq. 43) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 44) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = -y(j)
               else if (i .eq. 45) then
                  x(jj) = -z(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 46) then
                  x(jj) = -z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -x(j)
               else if (i .eq. 47) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = -y(j)
                  z(jj) = -x(j)
               else if (i .eq. 48) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pm3(-)n space group  (International Tables 223)
c
      else if (spacegrp .eq. 'Pm3(-)n ') then
         nsym = 48
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 7) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 12) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 14) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 15) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 16) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 17) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 18) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 19) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 20) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 21) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 22) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 23) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 24) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 25) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 26) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = -z(j)
               else if (i .eq. 27) then
                  x(jj) = x(j)
                  y(jj) = -y(j)
                  z(jj) = z(j)
               else if (i .eq. 28) then
                  x(jj) = -x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 29) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 30) then
                  x(jj) = -z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 31) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = -y(j)
               else if (i .eq. 32) then
                  x(jj) = z(j)
                  y(jj) = -x(j)
                  z(jj) = y(j)
               else if (i .eq. 33) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 34) then
                  x(jj) = y(j)
                  y(jj) = -z(j)
                  z(jj) = x(j)
               else if (i .eq. 35) then
                  x(jj) = -y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 36) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = -x(j)
               else if (i .eq. 37) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 38) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 39) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 40) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 41) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 42) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 43) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 44) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 45) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 46) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 47) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 48) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Pn3(-)m space group  (Intl. Tables 224, origin at center)
c
      else if (spacegrp .eq. 'Pn3(-)m ') then
         nsym = 48
         noff = 1
         do i = 1, nsym
            ii = (i-1) * n
            do j = 1, n
               jj = j + ii
               if (i .eq. 1) then
                  x(jj) = x(j)
                  y(jj) = y(j)
                  z(jj) = z(j)
               else if (i .eq. 2) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = z(j)
               else if (i .eq. 3) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 4) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 5) then
                  x(jj) = z(j)
                  y(jj) = x(j)
                  z(jj) = y(j)
               else if (i .eq. 6) then
                  x(jj) = z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 7) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = y(j)
               else if (i .eq. 8) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 9) then
                  x(jj) = y(j)
                  y(jj) = z(j)
                  z(jj) = x(j)
               else if (i .eq. 10) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = z(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 11) then
                  x(jj) = y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 12) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = x(j)
               else if (i .eq. 13) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = -z(j)
               else if (i .eq. 14) then
                  x(jj) = -y(j)
                  y(jj) = -x(j)
                  z(jj) = -z(j)
               else if (i .eq. 15) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 16) then
                  x(jj) = -y(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 17) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = -y(j)
               else if (i .eq. 18) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 19) then
                  x(jj) = -x(j)
                  y(jj) = -z(j)
                  z(jj) = -y(j)
               else if (i .eq. 20) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -z(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 21) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -x(j)
               else if (i .eq. 22) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 23) then
                  x(jj) = -z(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 24) then
                  x(jj) = -z(j)
                  y(jj) = -y(j)
                  z(jj) = -x(j)
               else if (i .eq. 25) then
                  x(jj) = -x(j)
                  y(jj) = -y(j)
                  z(jj) = -z(j)
               else if (i .eq. 26) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = -z(j)
               else if (i .eq. 27) then
                  x(jj) = 0.5d0 + x(j)
                  y(jj) = -y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 28) then
                  x(jj) = -x(j)
                  y(jj) = 0.5d0 + y(j)
                  z(jj) = 0.5d0 + z(j)
               else if (i .eq. 29) then
                  x(jj) = -z(j)
                  y(jj) = -x(j)
                  z(jj) = -y(j)
               else if (i .eq. 30) then
                  x(jj) = -z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 31) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = 0.5d0 + x(j)
                  z(jj) = -y(j)
               else if (i .eq. 32) then
                  x(jj) = 0.5d0 + z(j)
                  y(jj) = -x(j)
                  z(jj) = 0.5d0 + y(j)
               else if (i .eq. 33) then
                  x(jj) = -y(j)
                  y(jj) = -z(j)
                  z(jj) = -x(j)
               else if (i .eq. 34) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = -z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 35) then
                  x(jj) = -y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = 0.5d0 + x(j)
               else if (i .eq. 36) then
                  x(jj) = 0.5d0 + y(j)
                  y(jj) = 0.5d0 + z(j)
                  z(jj) = -x(j)
               else if (i .eq. 37) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = z(j)
               else if (i .eq. 38) then
                  x(jj) = y(j)
                  y(jj) = x(j)
                  z(jj) = z(j)
               else if (i .eq. 39) then
                  x(jj) = 0.5d0 - y(j)
                  y(jj) = x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 40) then
                  x(jj) = y(j)
                  y(jj) = 0.5d0 - x(j)
                  z(jj) = 0.5d0 - z(j)
               else if (i .eq. 41) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = y(j)
               else if (i .eq. 42) then
                  x(jj) = x(j)
                  y(jj) = 0.5d0 - z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 43) then
                  x(jj) = x(j)
                  y(jj) = z(j)
                  z(jj) = y(j)
               else if (i .eq. 44) then
                  x(jj) = 0.5d0 - x(j)
                  y(jj) = z(j)
                  z(jj) = 0.5d0 - y(j)
               else if (i .eq. 45) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = x(j)
               else if (i .eq. 46) then
                  x(jj) = 0.5d0 - z(j)
                  y(jj) = y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 47) then
                  x(jj) = z(j)
                  y(jj) = 0.5d0 - y(j)
                  z(jj) = 0.5d0 - x(j)
               else if (i .eq. 48) then
                  x(jj) = z(j)
                  y(jj) = y(j)
                  z(jj) = x(j)
               end if
               call cellatom (jj,j)
            end do
         end do
c
c     Fm3(-)m space group  (Intl. Tables 225, Face Centered Cubic)
c
      else if (spacegrp .eq. 'Fm3(-)m ') then
         nsym = 48
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = -x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = -z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = -x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = -z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Fm3(-)c space group  (International Tables 226)
c
      else if (spacegrp .eq. 'Fm3(-)c ') then
         nsym = 48
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Fd3(-)m space group  (Intl. Tables 227, origin at center)
c
      else if (spacegrp .eq. 'Fd3(-)m ') then
         nsym = 48
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Fd3(-)c space group  (Intl. Tables 228, origin at center)
c
      else if (spacegrp .eq. 'Fd3(-)c ') then
         nsym = 48
         noff = 4
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.0d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               else if (k .eq. 3) then
                  xoff = 0.5d0
                  yoff = 0.0d0
                  zoff = 0.5d0
               else if (k .eq. 4) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.0d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = -y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = -z(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = x(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = z(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Im3(-)m space group  (Intl. Tables 229, Body Centered Cubic)
c
      else if (spacegrp .eq. 'Im3(-)m ') then
         nsym = 48
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = -x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = -z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = -x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = -z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = -y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = -y(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = y(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = x(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = -x(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = -z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = z(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = z(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
c
c     Ia3(-)d space group  (International Tables 230)
c
      else if (spacegrp .eq. 'Ia3(-)d ') then
         nsym = 48
         noff = 2
         do i = 1, nsym
            ii = (i-1) * noff * n
            do k = 1, noff
               kk = ii + (k-1)*n
               if (k .eq. 1) then
                  xoff = 0.0d0
                  yoff = 0.0d0
                  zoff = 0.0d0
               else if (k .eq. 2) then
                  xoff = 0.5d0
                  yoff = 0.5d0
                  zoff = 0.5d0
               end if
               do j = 1, n
                  jj = j + kk
                  if (i .eq. 1) then
                     x(jj) = x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 2) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 3) then
                     x(jj) = -x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 4) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 5) then
                     x(jj) = z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 6) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 7) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 8) then
                     x(jj) = -z(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 9) then
                     x(jj) = y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 10) then
                     x(jj) = -y(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 11) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 12) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 13) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 14) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 15) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 16) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 17) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 18) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 19) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 20) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 21) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 22) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 23) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  else if (i .eq. 24) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 25) then
                     x(jj) = -x(j) + xoff
                     y(jj) = -y(j) + yoff
                     z(jj) = -z(j) + zoff
                  else if (i .eq. 26) then
                     x(jj) = 0.5d0 + x(j) + xoff
                     y(jj) = y(j) + yoff
                     z(jj) = 0.5d0 - z(j) + zoff
                  else if (i .eq. 27) then
                     x(jj) = x(j) + xoff
                     y(jj) = 0.5d0 - y(j) + yoff
                     z(jj) = 0.5d0 + z(j) + zoff
                  else if (i .eq. 28) then
                     x(jj) = 0.5d0 - x(j) + xoff
                     y(jj) = 0.5d0 + y(j) + yoff
                     z(jj) = z(j) + zoff
                  else if (i .eq. 29) then
                     x(jj) = -z(j) + xoff
                     y(jj) = -x(j) + yoff
                     z(jj) = -y(j) + zoff
                  else if (i .eq. 30) then
                     x(jj) = 0.5d0 - z(j) + xoff
                     y(jj) = 0.5d0 + x(j) + yoff
                     z(jj) = y(j) + zoff
                  else if (i .eq. 31) then
                     x(jj) = 0.5d0 + z(j) + xoff
                     y(jj) = x(j) + yoff
                     z(jj) = 0.5d0 - y(j) + zoff
                  else if (i .eq. 32) then
                     x(jj) = z(j) + xoff
                     y(jj) = 0.5d0 - x(j) + yoff
                     z(jj) = 0.5d0 + y(j) + zoff
                  else if (i .eq. 33) then
                     x(jj) = -y(j) + xoff
                     y(jj) = -z(j) + yoff
                     z(jj) = -x(j) + zoff
                  else if (i .eq. 34) then
                     x(jj) = y(j) + xoff
                     y(jj) = 0.5d0 - z(j) + yoff
                     z(jj) = 0.5d0 + x(j) + zoff
                  else if (i .eq. 35) then
                     x(jj) = 0.5d0 - y(j) + xoff
                     y(jj) = 0.5d0 + z(j) + yoff
                     z(jj) = x(j) + zoff
                  else if (i .eq. 36) then
                     x(jj) = 0.5d0 + y(j) + xoff
                     y(jj) = z(j) + yoff
                     z(jj) = 0.5d0 - x(j) + zoff
                  else if (i .eq. 37) then
                     x(jj) = 0.25d0 - y(j) + xoff
                     y(jj) = 0.75d0 - x(j) + yoff
                     z(jj) = 0.75d0 + z(j) + zoff
                  else if (i .eq. 38) then
                     x(jj) = 0.25d0 + y(j) + xoff
                     y(jj) = 0.25d0 + x(j) + yoff
                     z(jj) = 0.25d0 + z(j) + zoff
                  else if (i .eq. 39) then
                     x(jj) = 0.75d0 - y(j) + xoff
                     y(jj) = 0.75d0 + x(j) + yoff
                     z(jj) = 0.25d0 - z(j) + zoff
                  else if (i .eq. 40) then
                     x(jj) = 0.75d0 + y(j) + xoff
                     y(jj) = 0.25d0 - x(j) + yoff
                     z(jj) = 0.75d0 - z(j) + zoff
                  else if (i .eq. 41) then
                     x(jj) = 0.25d0 - x(j) + xoff
                     y(jj) = 0.75d0 - z(j) + yoff
                     z(jj) = 0.75d0 + y(j) + zoff
                  else if (i .eq. 42) then
                     x(jj) = 0.75d0 + x(j) + xoff
                     y(jj) = 0.25d0 - z(j) + yoff
                     z(jj) = 0.75d0 - y(j) + zoff
                  else if (i .eq. 43) then
                     x(jj) = 0.25d0 + x(j) + xoff
                     y(jj) = 0.25d0 + z(j) + yoff
                     z(jj) = 0.25d0 + y(j) + zoff
                  else if (i .eq. 44) then
                     x(jj) = 0.75d0 - x(j) + xoff
                     y(jj) = 0.75d0 + z(j) + yoff
                     z(jj) = 0.25d0 - y(j) + zoff
                  else if (i .eq. 45) then
                     x(jj) = 0.25d0 - z(j) + xoff
                     y(jj) = 0.75d0 - y(j) + yoff
                     z(jj) = 0.75d0 + x(j) + zoff
                  else if (i .eq. 46) then
                     x(jj) = 0.75d0 - z(j) + xoff
                     y(jj) = 0.75d0 + y(j) + yoff
                     z(jj) = 0.25d0 - x(j) + zoff
                  else if (i .eq. 47) then
                     x(jj) = 0.75d0 + z(j) + xoff
                     y(jj) = 0.25d0 - y(j) + yoff
                     z(jj) = 0.75d0 - x(j) + zoff
                  else if (i .eq. 48) then
                     x(jj) = 0.25d0 + z(j) + xoff
                     y(jj) = 0.25d0 + y(j) + yoff
                     z(jj) = 0.25d0 + x(j) + zoff
                  end if
                  call cellatom (jj,j)
               end do
            end do
         end do
      end if
c
c     set the total number of atoms in the full unit cell
c
      n = nsym * noff * n
      return
      end
