c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2000 by P. Bagossi, P. Ren & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  program poledit  --  manipulate atomic multipole values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "poledit" provides for modification and manipulation of
c     the atomic multipole electrostatic models used in Tinker
c
c
      program poledit
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
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Tinker Multipole Editing Utility Can :',
     &           //,4x,'(1) Multipole Parameters from GDMA Output',
     &           /,4x,'(2) Alter Local Coordinate Frame Definitions',
     &           /,4x,'(3) Removal of Intramolecular Polarization')
         do while (mode.lt.1 .or. mode.gt.3)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     perform the desired multipole manipulation operation
c
      if (mode .eq. 1) then
         use_mpole = .true.
         use_polar = .true.
         idma = freeunit ()
         call readgdma (idma)
         call initprm
         call molsetup
         call setframe
         call rotframe
         call setpolar
         call alterpol
         call avgpole
         call prtpole
      else if (mode .eq. 2) then
         call getxyz
         call attach
         call field
         call katom
         call kmpole
         call kpolar
         call fixframe
         call prtpole
      else if (mode .eq. 3) then
         call getxyz
         call attach
         call field
         call katom
         call kmpole
         call kpolar
         call alterpol
         call avgpole
         call prtpole
      end if
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine molsetup  --  set molecule for polarization  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "molsetup" generates trial parameters needed to perform
c     polarizable multipole calculations on a structure read
c     from distributed multipole analysis output
c
c
      subroutine molsetup
      use atomid
      use atoms
      use couple
      use files
      use mpole
      use polar
      use ptable
      implicit none
      integer i,j
      integer atn,size
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 ri,rij,dij
      real*8, allocatable :: rad(:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (rad(n))
c
c     set base atomic radii from covalent radius values
c
      do i = 1, n
         rad(i) = 0.76d0
         atn = atomic(i)
         if (atn .ne. 0)  rad(i) = covrad(atn)
         rad(i) = 1.1d0 * rad(i)
      end do
c
c     assign atom connectivities based on interatomic distances
c
      do i = 1, n
         n12(i) = 0
      end do
      do i = 1, n-1
         xi = x(i)
         yi = y(i)
         zi = z(i)
         ri = rad(i)
         do j = i+1, n
            xr = x(j) - xi
            yr = y(j) - yi
            zr = z(j) - zi
            rij = ri + rad(j)
            dij = sqrt(xr*xr + yr*yr + zr*zr)
            if (dij .lt. rij) then
               n12(i) = n12(i) + 1
               i12(n12(i),i) = j
               n12(j) = n12(j) + 1
               i12(n12(j),j) = i
            end if
         end do
      end do
      do i = 1, n
         call sort (n12(i),i12(1,i))
      end do
c
c     find the bonds, angles, torsions and small rings
c
      call attach
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     perform deallocation of some local arrays
c
      deallocate (rad)
c
c     assign unique atom types and set the valence values
c
      size = min(20,leng)
      do i = 1, n
         type(i) = i
         valence(i) = n12(i)
         story(i) = filename(1:size)
      end do
c
c     assign the standard atomic weight by atomic number
c
      do i = 1, n
         mass(i) = 1.0d0
         atn = atomic(i)
         if (atn .ne. 0)  mass(i) = atmass(atn)
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(ipole))  allocate (ipole(n))
      if (.not. allocated(polsiz))  allocate (polsiz(n))
      if (.not. allocated(pollist))  allocate (pollist(n))
c
c     set atomic multipole sites and polarizability indices
c
      npole = n
      npolar = n
      do i = 1, n
         ipole(i) = i
         polsiz(i) = 13
         pollist(i) = i
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setframe  --  define local coordinate frames  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setframe" assigns a local coordinate frame at each atomic
c     multipole site using high priority connected atoms along axes
c
c
      subroutine setframe
      use atomid
      use atoms
      use couple
      use iounit
      use mpole
      implicit none
      integer i,j,k,m,kb
      integer ia,ib,ic,id
      integer kab,kac,kad
      integer kbc,kbd,kcd
      integer priority
      logical exist,query
      logical change
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(zaxis))  allocate (zaxis(npole))
      if (.not. allocated(xaxis))  allocate (xaxis(npole))
      if (.not. allocated(yaxis))  allocate (yaxis(npole))
      if (.not. allocated(polaxe))  allocate (polaxe(npole))
c
c     initialize the local frame type and defining atoms
c
      do i = 1, npole
         polaxe(i) = 'None'
         zaxis(i) = 0
         xaxis(i) = 0
         yaxis(i) = 0
      end do
c
c     assign the local frame definition for an isolated atom
c
      do i = 1, npole
         j = n12(i)
         if (j .eq. 0) then
            polaxe(i) = 'None'
            zaxis(i) = 0
            xaxis(i) = 0
            yaxis(i) = 0
c
c     assign the local frame definition for a monovalent atom
c
         else if (j .eq. 1) then
            ia = i12(1,i)
            if (n12(ia) .eq. 1) then
               polaxe(i) = 'Z-Only'
               zaxis(i) = ia
               xaxis(i) = 0
               yaxis(i) = 0
            else
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               m = 0
               do k = 1, n12(ia)
                  ib = i12(k,ia)
                  kb = atomic(ib)
                  if (kb.gt.m .and. ib.ne.i) then
                     xaxis(i) = ib
                     m = kb
                  end if
               end do
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a divalent atom
c
         else if (j .eq. 2) then
            ia = i12(1,i)
            ib = i12(2,i)
            kab = priority (i,ia,ib)
            if (kab .eq. ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab .eq. ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               yaxis(i) = 0
            else
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a trivalent atom
c
         else if (j .eq. 3) then
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i12(3,i)
            kab = priority (i,ia,ib)
            kac = priority (i,ia,ic)
            kbc = priority (i,ib,ic)
            if (kab.eq.0 .and. kac.eq.0) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.ia .and. kac.eq.ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               if (kbc .eq. ic)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kab.eq.ib .and. kbc.eq.ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               if (kac .eq. ic)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kac.eq.ic .and. kbc.eq.ic) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ic
               xaxis(i) = ia
               if (kab .eq. ib)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kab .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kac .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kbc .eq. 0) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = 0
            end if
c
c     assign the local frame definition for a tetravalent atom
c
         else if (j .eq. 4) then
            ia = i12(1,i)
            ib = i12(2,i)
            ic = i12(3,i)
            id = i12(4,i)
            kab = priority (i,ia,ib)
            kac = priority (i,ia,ic)
            kad = priority (i,ia,id)
            kbc = priority (i,ib,ic)
            kbd = priority (i,ib,id)
            kcd = priority (i,ic,id)
            if (kab.eq.0 .and. kac.eq.0 .and. kad.eq.0) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.ia .and. kac.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ia
               xaxis(i) = ib
               if (kbc.eq.ic .and. kcd.eq.ic)  xaxis(i) = ic
               if (kbd.eq.id .and. kcd.eq.id)  xaxis(i) = id
               if (kbc.eq.ic .and. kcd.eq.0)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kab.eq.ib .and. kbc.eq.ib .and. kbd.eq.ib) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ib
               xaxis(i) = ia
               if (kac.eq.ic .and. kcd.eq.ic)  xaxis(i) = ic
               if (kad.eq.id .and. kcd.eq.id)  xaxis(i) = id
               if (kac.eq.ic .and. kcd.eq.0)  xaxis(i) = ic
               yaxis(i) = 0
            else if (kac.eq.ic .and. kbc.eq.ic .and. kcd.eq.ic) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = ic
               xaxis(i) = ia
               if (kab.eq.ib .and. kbd.eq.ib)  xaxis(i) = ib
               if (kad.eq.id .and. kbd.eq.id)  xaxis(i) = id
               if (kab.eq.ib .and. kbd.eq.0)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kad.eq.id .and. kbd.eq.id .and. kcd.eq.id) then
               polaxe(i) = 'Z-then-X'
               zaxis(i) = id
               xaxis(i) = ia
               if (kab.eq.ib .and. kbc.eq.ib)  xaxis(i) = ib
               if (kac.eq.ic .and. kbc.eq.ic)  xaxis(i) = ic
               if (kab.eq.ib .and. kbc.eq.0)  xaxis(i) = ib
               yaxis(i) = 0
            else if (kab.eq.0 .and. kac.eq.0 .and. kbc.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = ic
            else if (kab.eq.0 .and. kad.eq.0 .and. kbd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = id
            else if (kac.eq.0 .and. kad.eq.0 .and. kcd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = id
            else if (kbc.eq.0 .and. kbd.eq.0 .and. kcd.eq.0) then
               polaxe(i) = 'Z-Bisect'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = id
            else if (kab.eq.0 .and. kac.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ib
               yaxis(i) = 0
            else if (kac.eq.0 .and. kab.eq.ia .and. kad.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kad.eq.0 .and. kab.eq.ia .and. kac.eq.ia) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ia
               xaxis(i) = id
               yaxis(i) = 0
            else if (kbc.eq.0 .and. kab.eq.ib .and. kbd.eq.ib) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = ic
               yaxis(i) = 0
            else if (kbd.eq.0 .and. kab.eq.ib .and. kbc.eq.ib) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ib
               xaxis(i) = id
               yaxis(i) = 0
            else if (kcd.eq.0 .and. kac.eq.ic .and. kbc.eq.ic) then
               polaxe(i) = 'Bisector'
               zaxis(i) = ic
               xaxis(i) = id
               yaxis(i) = 0
            end if
         end if
      end do
c
c     list the local frame definition for each multipole site
c
      write (iout,10)
   10 format (/,' Local Frame Definition for Multipole Sites :')
      write (iout,20)
   20 format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',
     &           2x,'X Axis',2x,'Y Axis',/)
      do i = 1, npole
         write (iout,30)  i,name(i),polaxe(i),zaxis(i),
     &                    xaxis(i),yaxis(i)
   30    format (i8,6x,a3,7x,a8,2x,3i8)
      end do
c
c     allow the user to manually alter local coordinate frames
c
      change = .false.
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=40,end=40)  i
         if (i .eq. 0)  query = .false.
      end if
   40 continue
      do while (query)
         i = 0
         ia = 0
         ib = 0
         ic = 0
         write (iout,50)
   50    format (/,' Enter Altered Local Frame Definition',
     &              ' [<Enter>=Exit] :  ',$)
         read (input,60)  record
   60    format (a240)
         read (record,*,err=70,end=70)  i,ia,ib,ic
   70    continue
         if (i .eq. 0) then
            query = .false.
         else
            change = .true.
            if (ia .eq. 0)  polaxe(i)= 'None'
            if (ia.ne.0 .and. ib.eq.0)  polaxe(i) = 'Z-Only'
            if (ia.gt.0 .and. ib.gt.0)  polaxe(i) = 'Z-then-X'
            if (ia.lt.0 .or. ib.lt.0)  polaxe(i) = 'Bisector'
            if (ib.lt.0 .and. ic.lt.0)  polaxe(i) = 'Z-Bisect'
            if (max(ia,ib,ic) .lt. 0)  polaxe(i) = '3-Fold'
            zaxis(i) = abs(ia)
            xaxis(i) = abs(ib)
            yaxis(i) = abs(ic)
         end if
      end do
c
c     repeat local frame list if definitions were altered
c
      if (change) then
         write (iout,80)
   80    format (/,' Local Frame Definition for Multipole Sites :')
         write (iout,90)
   90    format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',
     &              2x,'X Axis',2x,'Y Axis',/)
         do i = 1, npole
            write (iout,100)  i,name(i),polaxe(i),zaxis(i),
     &                        xaxis(i),yaxis(i)
  100       format (i8,6x,a3,7x,a8,2x,3i8)
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function priority  --  atom priority for axis assignment  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "priority" decides which of two connected atoms should be
c     preferred during construction of a local coordinate frame
c     and returns that atom number; if the two atoms have equal
c     priority then a zero is returned
c
c
      function priority (i,ia,ib)
      use atomid
      use couple
      implicit none
      integer i,k,m
      integer ia,ib
      integer ka,kb
      integer ka1,kb1
      integer ka2,kb2
      integer priority
c
c
c     get priority based on atomic number and connected atoms
c
      ka = atomic(ia)
      kb = atomic(ib)
      if (ka .gt. kb) then
         priority = ia
      else if (kb .gt. ka) then
         priority = ib
      else
         ka1 = 0
         ka2 = 0
         do k = 1, n12(ia)
            m = i12(k,ia)
            if (i .ne. m) then
               m = atomic(m)
               if (m .ge. ka1) then
                  ka2 = ka1
                  ka1 = m
               else if (m .gt. ka2) then
                  ka2 = m
               end if
            end if
         end do
         kb1 = 0
         kb2 = 0
         do k = 1, n12(ib)
            m = i12(k,ib)
            if (i .ne. m) then
               m = atomic(m)
               if (m .gt. kb1) then
                  kb2 = kb1
                  kb1 = m
               else if (m .gt. kb2) then
                  kb2 = m
               end if
            end if
         end do
         if (n12(ia) .lt. n12(ib)) then
            priority = ia
         else if (n12(ib) .lt. n12(ia)) then
            priority = ib
         else if (ka1 .gt. kb1) then
            priority = ia
         else if (kb1 .gt. ka1) then
            priority = ib
         else if (ka2 .gt. kb2) then
            priority = ia
         else if (kb2 .gt. ka2) then
            priority = ib
         else
            priority = 0
         end if
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine rotframe  --  convert multipoles to local frame  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "rotframe" takes the global multipole moments and rotates them
c     into the local coordinate frame defined at each atomic site
c
c
      subroutine rotframe
      use atomid
      use atoms
      use iounit
      use mpole
      use units
      implicit none
      integer i,j,k
      integer xaxe
      integer yaxe
      integer zaxe
      real*8 a(3,3)
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(pole))  allocate (pole(maxpole,n))
      if (.not. allocated(rpole))  allocate (rpole(maxpole,n))
      if (.not. allocated(spole))  allocate (spole(maxpole,n))
      if (.not. allocated(srpole))  allocate (srpole(maxpole,n))
c
c     store the global multipoles in the local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         rpole(1,i) = pole(1,i)
         do j = 2, 4
            rpole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            rpole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     print the local frame Cartesian atomic multipoles
c
      write (iout,10)
   10 format (/,' Local Frame Cartesian Multipole Moments :')
      do i = 1, n
         k = pollist(i)
         if (k .eq. 0) then
            write (iout,20)  i,name(i),atomic(i)
   20       format (/,' Atom:',i8,9x,'Name:',3x,a3,7x,
     &                 'Atomic Number:',i8)
            write (iout,30)
   30       format (/,' No Atomic Multipole Moments for this Site')
         else
            zaxe = zaxis(k)
            xaxe = xaxis(k)
            yaxe = yaxis(k)
            if (yaxe .lt. 0)  yaxe = -yaxe
            write (iout,40)  i,name(i),atomic(i)
   40       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,50)  polaxe(k),zaxe,xaxe,yaxe
   50       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,60)  rpole(1,k)
   60       format (/,' Charge:',10x,f15.5)
            write (iout,70)  rpole(2,k),rpole(3,k),rpole(4,k)
   70       format (' Dipole:',10x,3f15.5)
            write (iout,80)  rpole(5,k)
   80       format (' Quadrupole:',6x,f15.5)
            write (iout,90)  rpole(8,k),rpole(9,k)
   90       format (18x,2f15.5)
            write (iout,100)  rpole(11,k),rpole(12,k),rpole(13,k)
  100       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine fixframe  --  alter the local frame definition  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "fixframe" is a service routine that alters the local frame
c     definition for specified atoms
c
c
      subroutine fixframe
      use atomid
      use atoms
      use couple
      use files
      use keys
      use kpolr
      use iounit
      use mpole
      use polar
      use units
      implicit none
      integer i,j,k
      integer ia,ib,ic
      integer xaxe
      integer yaxe
      integer zaxe
      real*8 eps,ci,cj
      real*8 big,sum
      real*8 a(3,3)
      logical query,change
      character*240 record
c
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     list the local frame definition for each multipole site
c
      write (iout,10)
   10 format (/,' Local Frame Definition for Multipole Sites :')
      write (iout,20)
   20 format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',2x,
     &          'X Axis',2x,'Y Axis',/)
      do i = 1, n
         k = pollist(i)
         if (k .eq. 0) then
            write (iout,30)  i,name(i)
   30       format (i8,6x,a3,10x,'--',11x,'--',6x,'--',6x,'--')
         else
            zaxe = zaxis(k)
            xaxe = xaxis(k)
            yaxe = yaxis(k)
            if (yaxe .lt. 0)  yaxe = -yaxe
            write (iout,40)  i,name(i),polaxe(k),zaxe,xaxe,yaxe
   40       format (i8,6x,a3,7x,a8,2x,3i8)
         end if
      end do
c
c     allow the user to manually alter local coordinate frames
c
      query = .true.
      change = .false.
      do while (query)
         i = 0
         k = 0
         ia = 0
         ib = 0
         ic = 0
         write (iout,50)
   50    format (/,' Enter Altered Local Frame Definition',
     &              ' [<Enter>=Exit] :  ',$)
         read (input,60)  record
   60    format (a240)
         read (record,*,err=70,end=70)  k,ia,ib,ic
   70    continue
         if (k .eq. 0) then
            query = .false.
         else
            i = pollist(k)
         end if
         if (i .ne. 0) then
            change = .true.
            if (ia .eq. 0)  polaxe(i) = 'None'
            if (ia.ne.0 .and. ib.eq.0)  polaxe(i) = 'Z-Only'
            if (ia.gt.0 .and. ib.gt.0)  polaxe(i) = 'Z-then-X'
            if (ia.lt.0  .or. ib.lt.0)  polaxe(i) = 'Bisector'
            if (ib.lt.0 .and. ic.lt.0)  polaxe(i) = 'Z-Bisect'
            if (max(ia,ib,ic)  .lt. 0)  polaxe(i) = '3-Fold'
            zaxis(i) = abs(ia)
            xaxis(i) = abs(ib)
            yaxis(i) = abs(ic)
         end if
      end do
c
c     repeat local frame list if definitions were altered
c
      if (change) then
         write (iout,80)
   80    format (/,' Local Frame Definition for Multipole Sites :')
         write (iout,90)
   90    format (/,5x,'Atom',5x,'Name',6x,'Axis Type',5x,'Z Axis',2x,
     &             'X Axis',2x,'Y Axis',/)
         do k = 1, npole
            i = pollist(k)
            if (i .eq. 0) then
               write (iout,100)  k,name(k)
  100          format (i8,6x,a3,10x,'--',11x,'--',6x,'--',6x,'--')
            else
               zaxe = zaxis(i)
               xaxe = xaxis(i)
               yaxe = yaxis(i)
               if (yaxe .lt. 0)  yaxe = -yaxe
               write (iout,110)  k,name(k),polaxe(i),zaxe,xaxe,yaxe
  110          format (i8,6x,a3,7x,a8,2x,3i8)
            end if
         end do
      end if
c
c     store the global multipoles in the local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         pole(1,i) = pole(1,i)
         do j = 2, 4
            pole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     regularize the multipole moments to desired precision
c
      eps = 0.00001d0
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = dble(nint(pole(j,i)/eps)) * eps
         end do
      end do
c
c     maintain integer net charge over atomic multipoles
c
      k = 0
      big = 0.0d0
      sum = 0.0d0
      do i = 1, npole
         sum = sum + pole(1,i)
         ci = abs(pole(1,i))
         if (ci .gt. big) then
            do j = 1, npole
               cj = abs(pole(1,j))
               if (i.ne.j .and. ci.eq.cj)  goto 120
            end do
            k = i
            big = ci
  120       continue
         end if
      end do
      sum = sum - dble(nint(sum))
      if (k .ne. 0)  pole(1,k) = pole(1,k) - sum
c
c     maintain traceless quadrupole at each multipole site
c
      do i = 1, npole
         sum = pole(5,i) + pole(9,i) + pole(13,i)
         big = max(abs(pole(5,i)),abs(pole(9,i)),abs(pole(13,i)))
         k = 0
         if (big .eq. abs(pole(5,i)))  k = 5
         if (big .eq. abs(pole(9,i)))  k = 9
         if (big .eq. abs(pole(13,i)))  k = 13
         if (pole(9,i) .eq. pole(13,i))  k = 5
         if (pole(5,i) .eq. pole(13,i))  k = 9
         if (pole(5,i) .eq. pole(9,i))  k = 13
         if (k .ne. 0)  pole(k,i) = pole(k,i) - sum
      end do
c
c     print the altered local frame atomic multipole values
c
      write (iout,130)
  130 format (/,' Multipoles With Altered Local Frame Definition :')
      do i = 1, n
         k = pollist(i)
         if (i .eq. 0) then
            write (iout,140)  i,name(i),atomic(i)
  140       format (/,' Atom:',i8,9x,'Name:',3x,a3,7x,
     &                 'Atomic Number:',i8)
            write (iout,150)
  150       format (/,' No Atomic Multipole Moments for this Site')
         else
            zaxe = zaxis(k)
            xaxe = xaxis(k)
            yaxe = yaxis(k)
            if (yaxe .lt. 0)  yaxe = -yaxe
            write (iout,160)  i,name(i),atomic(i)
  160       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,170)  polaxe(k),zaxe,xaxe,yaxe
  170       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,180)  pole(1,k)
  180       format (/,' Charge:',10x,f15.5)
            write (iout,190)  pole(2,k),pole(3,k),pole(4,k)
  190       format (' Dipole:',10x,3f15.5)
            write (iout,200)  pole(5,k)
  200       format (' Quadrupole:',6x,f15.5)
            write (iout,210)  pole(8,k),pole(9,k)
  210       format (18x,2f15.5)
            write (iout,220)  pole(11,k),pole(12,k),pole(13,k)
  220       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setpolar  --  define polarization and groups  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setpolar" assigns atomic polarizabilities, damping parameters
c     and polarization groups with user modification of these values
c
c
      subroutine setpolar
      use atomid
      use atoms
      use couple
      use iounit
      use kpolr
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k,m
      integer jj,ia,ib
      integer atn,next
      real*8 pol,thl
      real*8 sixth
      logical exist,query
      logical change
      logical aromatic
      logical chkarom
      character*1 answer
      character*6 pmodel
      character*240 record
      character*240 string
      external chkarom
c
c
c     allow the user to select the polarization model
c
      pmodel = 'AMOEBA'
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
      if (answer .eq. 'H')  pmodel = 'HIPPO '
c
c     set polarization damping scheme based on model selected
c
      if (pmodel .eq. 'AMOEBA') then
         use_thole = .true.
         use_chgpen = .false.
      else if (pmodel .eq. 'HIPPO ') then
         use_thole = .false.
         use_chgpen = .true.
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(polarity))  allocate (polarity(n))
      if (.not. allocated(thole))  allocate (thole(n))
      if (.not. allocated(pdamp))  allocate (pdamp(n))
c
c     assign default atomic polarizabilities for Thole model
c
      if (pmodel .eq. 'AMOEBA') then
         do i = 1, n
            polarity(i) = 0.0d0
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
c     alter Thole values for alkene/aromatic carbon and hydrogen
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
c
c     assign default atomic polarizabilities for HIPPO model
c
      else if (pmodel .eq. 'HIPPO ') then
         do i = 1, n
            polarity(i) = 0.0d0
            thole(i) = 0.0d0
            atn = atomic(i)
            if (atn .eq. 1) then
               polarity(i) = 0.373d0
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
   40             continue
                  aromatic = chkarom (k)
                  if (aromatic)  polarity(i) = 0.1106d0
               end if
               if (k .eq. 7)  polarity(i) = 0.005d0
               if (k .eq. 8)  polarity(i) = 0.3698d0
               if (k .eq. 16)  polarity(i) = 0.2093d0
            else if (atn .eq. 5) then
               polarity(i) = 1.6d0
            else if (atn .eq. 6) then
               polarity(i) = 0.755d0
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
               polarity(i) = 0.9354d0
   50          continue
               if (n12(i) .eq. 3) then
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.6 .and. n12(k).eq.3) then
                        polarity(i) = 1.9384d0
                     end if
                  end do
                  do j = 1, n12(i)
                     k = i12(j,i)
                     if (atomic(k).eq.8 .and. n12(k).eq.1) then
                        polarity(i) = 0.6577d0
                     end if
                  end do
               end if
               if (chkarom(i)) then
                  polarity(i) = 1.5624d0
                  do j = 1, n12(i)
                     k = atomic(i12(j,i))
                     if (k.ne.6 .and. k.ne.1)  polarity(i) = 1.2811d0
                  end do
               end if
               if (n12(i) .eq. 2)  polarity(i) = 1.0d0
            else if (atn .eq. 7) then
               polarity(i) = 1.4289d0
               if (n12(i) .eq. 3)  polarity(i) = 1.4545d0
               if (chkarom(i))  polarity(i) = 1.3037d0
            else if (atn .eq. 8) then
               polarity(i) = 0.6645d0
               if (n12(i) .eq. 1) then
                  k = i12(1,i)
                  if (atomic(k).eq.6 .and. n12(k).eq.3) then
                     polarity(i) = 1.4266d0
                     do j = 1, n13(i)
                        m = i13(j,i)
                        if (atomic(m).eq.8 .and. n12(m).eq.1) then
                           polarity(i) = 1.8809d0
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
                     if (jj .eq. 1)  polarity(i) = 1.0d0
                     if (jj .eq. 2)  polarity(i) = 1.0d0
                     if (jj .eq. 3)  polarity(i) = 1.0d0
                     if (jj .eq. 4)  polarity(i) = 1.0d0
                  end if
               end if
            else if (atn .eq. 9) then
               polarity(i) = 0.5d0
            else if (atn .eq. 15) then
               polarity(i) = 1.8d0
            else if (atn .eq. 16) then
               polarity(i) = 3.1967d0
               if (n12(i) .eq. 4)  polarity(i) = 2.458d0
            else if (atn .eq. 17) then
               polarity(i) = 2.366d0
            else if (atn .eq. 35) then
               polarity(i) = 3.4458d0
            else if (atn .eq. 35) then
               polarity(i) = 5.5d0
            end if
         end do
      end if
c
c     compute the Thole polarizability damping values
c
      sixth = 1.0d0 / 6.0d0
      do i = 1, npole
         if (thole(i) .eq. 0.0d0) then
            pdamp(i) = 0.0d0
         else
            pdamp(i) = polarity(i)**sixth
         end if
      end do
c
c     list the polariability values for each multipole site
c
      write (iout,60)
   60 format (/,' Atomic Polarizabilities for Multipole Sites :')
      if (pmodel .eq. 'AMOEBA') then
         write (iout,70)
   70    format (/,5x,'Atom',5x,'Name',7x,'Polarize',10x,'Thole',/)
      else
         write (iout,80)
   80    format (/,5x,'Atom',5x,'Name',7x,'Polarize',/)
      end if
      do k = 1, n
         i = pollist(k)
         if (pmodel .eq. 'AMOEBA') then
            if (i .eq. 0) then
               write (iout,90)  k,name(k)
   90          format (i8,6x,a3,12x,'--',13x,'--')
            else
               write (iout,100)  k,name(k),polarity(i),thole(i)
  100          format (i8,6x,a3,4x,f12.4,3x,f12.4)
            end if
         else
            if (i .eq. 0) then
               write (iout,110)  k,name(k)
  110          format (i8,6x,a3,12x,'--')
            else
               write (iout,120)  k,name(k),polarity(i)
  120          format (i8,6x,a3,4x,f12.4)
            end if
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
         read (string,*,err=130,end=130)  i
         if (i .eq. 0)  query = .false.
      end if
  130 continue
      do while (query)
         i = 0
         k = 0
         pol = 0.0d0
         thl = 0.39d0
         write (iout,140)
  140    format (/,' Enter Atom Number & Polarizability Values',
     &              ' [<Enter>=Exit] :  ',$)
         read (input,150)  record
  150    format (a240)
         read (record,*,err=160,end=160)  k,pol,thl
  160    continue
         if (k .eq. 0) then
            query = .false.
         else
            i = pollist(k)
         end if
         if (i .ne. 0) then
            change = .true.
            polarity(i) = pol
            thole(i) = thl
         end if
      end do
c
c     repeat polarizability values if parameters were altered
c
      if (change) then
         write (iout,170)
  170    format (/,' Atomic Polarizabilities for Multipole Sites :')
         if (pmodel .eq. 'AMOEBA') then
            write (iout,180)
  180       format (/,5x,'Atom',5x,'Name',7x,'Polarize',10x,'Thole',/)
         else
            write (iout,190)
  190       format (/,5x,'Atom',5x,'Name',7x,'Polarize',/)
         end if
         do k = 1, n
            i = pollist(k)
            if (pmodel .eq. 'AMOEBA') then
               if (i .eq. 0) then
                  write (iout,200)  k,name(k)
  200             format (i8,6x,a3,12x,'--',13x,'--')
               else
                  write (iout,210)  k,name(k),polarity(i),thole(i)
  210             format (i8,6x,a3,4x,f12.4,3x,f12.4)
               end if
            else
               if (i .eq. 0) then
                  write (iout,220)  k,name(k)
  220             format (i8,6x,a3,12x,'--')
               else
                  write (iout,230)  k,name(k),polarity(i)
  230             format (i8,6x,a3,4x,f12.4)
               end if
            end if
         end do
      end if
c
c     use bonded atoms as initial guess at polarization groups
c
      write (iout,240)
  240 format (/,' The default is to place all atoms into one',
     &           ' polarization group;',
     &        /,' This can be altered by entering a series of',
     &           ' bonded atom pairs',
     &        /,' that separate the molecule into distinct',
     &           ' polarization groups')
      do i = 1, n
         do j = 1, n12(i)
            pgrp(j,i) = i12(j,i)
         end do
      end do
c
c     get the bonds that separate the polarization groups
c
      query = .true.
      i = -1
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=250,end=250)  i
         if (i .eq. 0)  query = .false.
      end if
  250 continue
      do while (query)
         ia = 0
         ib = 0
         write (iout,260)
  260    format (/,' Enter a Bond between Polarization Groups',
     &              ' [<Enter>=Exit] :  ',$)
         read (input,270)  record
  270    format (a240)
         read (record,*,err=280,end=280)  ia,ib
  280    continue
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
      call polargrp
c
c     list the polarization group for each multipole site
c
      write (iout,290)
  290 format (/,' Polarization Groups for Multipole Sites :')
      write (iout,300)
  300 format (/,5x,'Atom',5x,'Name',7x,'Polarization Group',
     &           ' Definition',/)
      do i = 1, n
         k = 0
         do j = 1, maxval
            if (pgrp(j,i) .ne. 0)  k = j
         end do
         write (iout,310)  i,name(i),(pgrp(j,i),j=1,k)
  310    format (i8,6x,a3,8x,20i6)
      end do
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
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine alterpol  --  alter multipoles for polarization  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "alterpol" finds an output set of atomic multipole parameters
c     which when used with an intergroup polarization model will
c     give the same electrostatic potential around the molecule as
c     the input set of multipole parameters with all atoms in one
c     polarization group
c
c     for example, the input parameters could be from a distributed
c     multipole analysis of a molecular wavefunction and the output
c     will be the parameter values that achieve the same potential
c     in the presence of intergroup (intramolecular) polarization
c
c
      subroutine alterpol
      use atomid
      use atoms
      use iounit
      use mpole
      use polar
      use units
      implicit none
      integer i,j,k
      integer xaxe
      integer yaxe
      integer zaxe
      real*8 a(3,3)
c
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute induced dipoles to be removed from QM multipoles
c
      call interpol
c
c     remove induced dipoles from global frame multipoles
c
      do i = 1, npole
         rpole(2,i) = rpole(2,i) - uind(1,i)
         rpole(3,i) = rpole(3,i) - uind(2,i)
         rpole(4,i) = rpole(4,i) - uind(3,i)
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     rotate the multipoles from global frame to local frame
c
      do i = 1, npole
         call rotmat (i,a)
         call invert (3,a)
         call rotsite (i,a)
      end do
c
c     copy the rotated multipoles back to local frame array
c
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = rpole(j,i)
         end do
      end do
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         rpole(1,i) = pole(1,i)
         do j = 2, 4
            rpole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            rpole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     print multipoles with intergroup polarization removed
c
      write (iout,10)
   10 format (/,' Multipoles after Removal of Intergroup',
     &           ' Polarization :')
      do i = 1, n
         k = pollist(i)
         if (i .eq. 0) then
            write (iout,20)  i,name(i),atomic(i)
   20       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,30)
   30       format (/,' No Atomic Multipole Moments for this Site')
         else
            zaxe = zaxis(k)
            xaxe = xaxis(k)
            yaxe = yaxis(k)
            if (yaxe .lt. 0)  yaxe = -yaxe
            write (iout,40)  i,name(i),atomic(i)
   40       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,50)  polaxe(k),zaxe,xaxe,yaxe
   50       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,60)  rpole(1,k)
   60       format (/,' Charge:',10x,f15.5)
            write (iout,70)  rpole(2,k),rpole(3,k),rpole(4,k)
   70       format (' Dipole:',10x,3f15.5)
            write (iout,80)  rpole(5,k)
   80       format (' Quadrupole:',6x,f15.5)
            write (iout,90)  rpole(8,k),rpole(9,k)
   90       format (18x,2f15.5)
            write (iout,100)  rpole(11,k),rpole(12,k),rpole(13,k)
  100       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine interpol  --  get intergroup induced dipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "interpol" is computes induced dipole moments for use during
c     removal of intergroup polarization
c
c
      subroutine interpol
      use atoms
      use iounit
      use mpole
      use polar
      use polpot
      use units
      implicit none
      integer i,j,k,iter
      integer maxiter
      real*8 eps,epsold
      real*8 polmin,norm
      real*8 a,b,sum
      real*8, allocatable :: poli(:)
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: rsd(:,:)
      real*8, allocatable :: zrsd(:,:)
      real*8, allocatable :: conj(:,:)
      real*8, allocatable :: vec(:,:)
      logical done
c
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(uind))  allocate (uind(3,npole))
c
c     perform dynamic allocation of some local arrays
c
      allocate (poli(npole))
      allocate (field(3,npole))
      allocate (rsd(3,npole))
      allocate (zrsd(3,npole))
      allocate (conj(3,npole))
      allocate (vec(3,npole))
c
c     compute induced dipoles as polarizability times field
c
      call dfieldi (field)
      do i = 1, npole
         do j = 1, 3
            uind(j,i) = polarity(i) * field(j,i)
         end do
      end do
c
c     for direct-only models set mutual scale factors to zero
c
      if (poltyp .eq. 'DIRECT') then
         u1scale = 0.0d0
         u2scale = 0.0d0
         u3scale = 0.0d0
         u4scale = 0.0d0
      end if
c
c     compute intergroup induced dipole moments via CG algorithm
c
      done = .false.
      maxiter = 500
      iter = 0
      eps = 100.0d0
      polmin = 0.00000001d0
      call ufieldi (field)
      do i = 1, npole
         poli(i) = max(polmin,polarity(i))
         do j = 1, 3
            rsd(j,i) = field(j,i)
            zrsd(j,i) = rsd(j,i) * poli(i)
            conj(j,i) = zrsd(j,i)
         end do
      end do
c
c     iterate the intergroup induced dipoles and check convergence
c
      do while (.not. done)
         iter = iter + 1
         do i = 1, npole
            do j = 1, 3
               vec(j,i) = uind(j,i)
               uind(j,i) = conj(j,i)
            end do
         end do
         call ufieldi (field)
         do i = 1, npole
            do j = 1, 3
               uind(j,i) = vec(j,i)
               vec(j,i) = conj(j,i)/poli(i) - field(j,i)
            end do
         end do
         a = 0.0d0
         sum = 0.0d0
         do i = 1, npole
            do j = 1, 3
               a = a + conj(j,i)*vec(j,i)
               sum = sum + rsd(j,i)*zrsd(j,i)
            end do
         end do
         if (a .ne. 0.0d0)  a = sum / a
         do i = 1, npole
            do j = 1, 3
               uind(j,i) = uind(j,i) + a*conj(j,i)
               rsd(j,i) = rsd(j,i) - a*vec(j,i)
            end do
         end do
         b = 0.0d0
         do i = 1, npole
            do j = 1, 3
               zrsd(j,i) = rsd(j,i) * poli(i)
               b = b + rsd(j,i)*zrsd(j,i)
            end do
         end do
         if (sum .ne. 0.0d0)  b = b / sum
         eps = 0.0d0
         do i = 1, npole
            do j = 1, 3
               conj(j,i) = zrsd(j,i) + b*conj(j,i)
               eps = eps + rsd(j,i)*rsd(j,i)
            end do
         end do
         eps = debye * sqrt(eps/dble(npolar))
         epsold = eps
         if (iter .eq. 1) then
            write (iout,10)
   10       format (/,' Determination of Intergroup Induced',
     &                 ' Dipoles :',
     &              //,4x,'Iter',8x,'RMS Change (Debye)',/)
         end if
         write (iout,20)  iter,eps
   20    format (i8,7x,f16.10)
         if (eps .lt. poleps)  done = .true.
         if (eps .gt. epsold)  done = .true.
         if (iter .ge. maxiter)  done = .true.
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (poli)
      deallocate (field)
      deallocate (rsd)
      deallocate (zrsd)
      deallocate (conj)
      deallocate (vec)
c
c     terminate the calculation if dipoles failed to converge
c
      if (eps .gt. poleps) then
         write (iout,30)
   30    format (/,' INTERPOL  --  Warning, Induced Dipoles',
     &              ' are not Converged')
         call prterr
         call fatal
      end if
c
c     print out a list of the final induced dipole moments
c
      write (iout,40)
   40 format (/,' Intergroup Induced Dipoles to be Removed',
     &           ' (Debye) :')
      write (iout,50)
   50 format (/,4x,'Atom',14x,'X',11x,'Y',11x,'Z',
     &           9x,'Total'/)
      do i = 1, npole
         k = ipole(i)
         norm = sqrt(uind(1,i)**2+uind(2,i)**2+uind(3,i)**2)
         write (iout,60)  k,(debye*uind(j,i),j=1,3),debye*norm
   60    format (i8,5x,4f12.4)
      end do
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine dfieldi  --  find permanent multipole field  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "dfieldi" computes the electrostatic field due to permanent
c     multipole moments
c
c
      subroutine dfieldi (field)
      use atoms
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5,rr7
      real*8 rr3i,rr5i,rr7i
      real*8 rr3k,rr5k,rr7k
      real*8 ci,dix,diy,diz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 dmpi(7),dmpk(7)
      real*8, allocatable :: dscale(:)
      real*8 field(3,*)
c
c
c     zero out the induced dipole and the field at each site
c
      do ii = 1, npole
         do j = 1, 3
            uind(j,ii) = 0.0d0
            field(j,ii) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
c
c     set array needed to scale atom and group interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
      end do
c
c     find the electrostatic field due to permanent multipoles
c
      do ii = 1, npole-1
         i = ipole(ii)
         ci = rpole(1,ii)
         dix = rpole(2,ii)
         diy = rpole(3,ii)
         diz = rpole(4,ii)
         qixx = rpole(5,ii)
         qixy = rpole(6,ii)
         qixz = rpole(7,ii)
         qiyy = rpole(9,ii)
         qiyz = rpole(10,ii)
         qizz = rpole(13,ii)
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               dscale(i12(j,i)) = p2scale
               do k = 1, np11(i)
                  if (i12(j,i) .eq. ip11(k,i))
     &               dscale(i12(j,i)) = p2iscale
               end do
            end do
            do j = 1, n13(i)
               dscale(i13(j,i)) = p3scale
               do k = 1, np11(i)
                  if (i13(j,i) .eq. ip11(k,i))
     &               dscale(i13(j,i)) = p3iscale
               end do
            end do
            do j = 1, n14(i)
               dscale(i14(j,i)) = p4scale
               do k = 1, np11(i)
                  if (i14(j,i) .eq. ip11(k,i))
     &               dscale(i14(j,i)) = p4iscale
               end do
            end do
            do j = 1, n15(i)
               dscale(i15(j,i)) = p5scale
               do k = 1, np11(i)
                  if (i15(j,i) .eq. ip11(k,i))
     &               dscale(i15(j,i)) = p5iscale
               end do
            end do
         else
            do j = 1, np11(i)
               dscale(ip11(j,i)) = d1scale
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = d2scale
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = d3scale
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = d4scale
            end do
         end if
c
c     evaluate higher-numbered sites with the original site
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            r2 = xr*xr + yr* yr + zr*zr
            r = sqrt(r2)
            ck = rpole(1,kk)
            dkx = rpole(2,kk)
            dky = rpole(3,kk)
            dkz = rpole(4,kk)
            qkxx = rpole(5,kk)
            qkxy = rpole(6,kk)
            qkxz = rpole(7,kk)
            qkyy = rpole(9,kk)
            qkyz = rpole(10,kk)
            qkzz = rpole(13,kk)
c
c     intermediates involving moments and separation distance
c
            dir = dix*xr + diy*yr + diz*zr
            qix = qixx*xr + qixy*yr + qixz*zr
            qiy = qixy*xr + qiyy*yr + qiyz*zr
            qiz = qixz*xr + qiyz*yr + qizz*zr
            qir = qix*xr + qiy*yr + qiz*zr
            dkr = dkx*xr + dky*yr + dkz*zr
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
            qkr = qkx*xr + qky*yr + qkz*zr
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               damp = pdi * pdamp(kk)
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(kk))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7 = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)
                  end if
               end if
               rr3 = scale3 / (r*r2)
               rr5 = 3.0d0 * scale5 / (r*r2*r2)
               rr7 = 15.0d0 * scale7 / (r*r2*r2*r2)
               fid(1) = -xr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dkx + 2.0d0*rr5*qkx
               fid(2) = -yr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dky + 2.0d0*rr5*qky
               fid(3) = -zr*(rr3*ck-rr5*dkr+rr7*qkr)
     &                     - rr3*dkz + 2.0d0*rr5*qkz
               fkd(1) = xr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*dix - 2.0d0*rr5*qix
               fkd(2) = yr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*diy - 2.0d0*rr5*qiy
               fkd(3) = zr*(rr3*ci+rr5*dir+rr7*qir)
     &                     - rr3*diz - 2.0d0*rr5*qiz
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kk)
               valk = pval(kk)
               alphak = palpha(kk)
               call dampdir (r,alphai,alphak,dmpi,dmpk)
               rr3 = 1.0d0 / (r*r2)
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr3i = dmpi(3) * rr3
               rr5i = dmpi(5) * rr5
               rr7i = dmpi(7) * rr7
               rr3k = dmpk(3) * rr3
               rr5k = dmpk(5) * rr5
               rr7k = dmpk(7) * rr7
               fid(1) = -xr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkx + 2.0d0*rr5k*qkx
               fid(2) = -yr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dky + 2.0d0*rr5k*qky
               fid(3) = -zr*(rr3*corek + rr3k*valk
     &                     - rr5k*dkr + rr7k*qkr)
     &                     - rr3k*dkz + 2.0d0*rr5k*qkz
               fkd(1) = xr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*dix - 2.0d0*rr5i*qix
               fkd(2) = yr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diy - 2.0d0*rr5i*qiy
               fkd(3) = zr*(rr3*corei + rr3i*vali
     &                     + rr5i*dir + rr7i*qir)
     &                     - rr3i*diz - 2.0d0*rr5i*qiz
            end if
            do j = 1, 3
               field(j,ii) = field(j,ii) + fid(j)*dscale(k)
               field(j,kk) = field(j,kk) + fkd(j)*dscale(k)
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            do j = 1, n12(i)
               dscale(i12(j,i)) = 1.0d0
            end do
            do j = 1, n13(i)
               dscale(i13(j,i)) = 1.0d0
            end do
            do j = 1, n14(i)
               dscale(i14(j,i)) = 1.0d0
            end do
            do j = 1, n15(i)
               dscale(i15(j,i)) = 1.0d0
            end do
         else
            do j = 1, np11(i)
               dscale(ip11(j,i)) = 1.0d0
            end do
            do j = 1, np12(i)
               dscale(ip12(j,i)) = 1.0d0
            end do
            do j = 1, np13(i)
               dscale(ip13(j,i)) = 1.0d0
            end do
            do j = 1, np14(i)
               dscale(ip14(j,i)) = 1.0d0
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine ufieldi  --  find induced intergroup field  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "ufieldi" computes the electrostatic field due to intergroup
c     induced dipole moments
c
c     literature reference:
c
c     P. Ren and J. W. Ponder, "Consistent Treatment of Inter- and
c     Intramolecular Polarization in Molecular Mechanics Calculations",
c     Journal of Computational Chemistry, 23, 1497-1506 (2002)
c
c
      subroutine ufieldi (field)
      use atoms
      use chgpen
      use couple
      use mplpot
      use mpole
      use polar
      use polgrp
      use polpot
      implicit none
      integer i,j,k
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 uix,uiy,uiz
      real*8 ukx,uky,ukz
      real*8 uir,ukr
      real*8 corei,corek
      real*8 vali,valk
      real*8 alphai,alphak
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 pdi,pti,pgamma
      real*8 fiu(3),fku(3)
      real*8 dmpik(5)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: uscale(:)
      real*8, allocatable :: wscale(:)
      real*8, allocatable :: gscale(:)
      real*8 field(3,*)
c
c
c     zero out the value of the field at each site
c
      do i = 1, npole
         do j = 1, 3
            field(j,i) = 0.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale(n))
      allocate (pscale(n))
      allocate (uscale(n))
      allocate (wscale(n))
      allocate (gscale(n))
c
c     set arrays needed to scale atom and group interactions
c
      do i = 1, n
         dscale(i) = 1.0d0
         pscale(i) = 1.0d0
         uscale(i) = 1.0d0
         wscale(i) = 1.0d0
         gscale(i) = 0.0d0
      end do
c
c     find the electrostatic field due to induced dipoles
c
      do ii = 1, npole-1
         i = ipole(ii)
         uix = uind(1,ii)
         uiy = uind(2,ii)
         uiz = uind(3,ii)
         pdi = pdamp(ii)
         pti = thole(ii)
         if (use_chgpen) then
            corei = pcore(ii)
            vali = pval(ii)
            alphai = palpha(ii)
         end if
c
c     set exclusion coefficients for connected atoms
c
         if (dpequal) then
            if (use_chgpen) then
               do j = 1, n12(i)
                  gscale(i12(j,i)) = w2scale - p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  gscale(i12(j,i)) = w2scale - p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  gscale(i13(j,i)) = w3scale - p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  gscale(i13(j,i)) = w3scale - p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  gscale(i14(j,i)) = w4scale - p4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  gscale(i14(j,i)) = w4scale - p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  gscale(i15(j,i)) = w5scale - p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  gscale(i15(j,i)) = w5scale - p5iscale
                  end do
               end do
            else
               do j = 1, n12(i)
                  pscale(i12(j,i)) = p2scale
                  do k = 1, np11(i)
                     if (i12(j,i) .eq. ip11(k,i))
     &                  pscale(i12(j,i)) = p2iscale
                  end do
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = p3scale
                  do k = 1, np11(i)
                     if (i13(j,i) .eq. ip11(k,i))
     &                  pscale(i13(j,i)) = p3iscale
                  end do
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = w4scale
                  do k = 1, np11(i)
                     if (i14(j,i) .eq. ip11(k,i))
     &                  pscale(i14(j,i)) = p4iscale
                  end do
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = p5scale
                  do k = 1, np11(i)
                     if (i15(j,i) .eq. ip11(k,i))
     &                  pscale(i15(j,i)) = p5iscale
                  end do
               end do
               do j = 1, np11(i)
                  uscale(ip11(j,i)) = u1scale
               end do
               do j = 1, np12(i)
                  uscale(ip12(j,i)) = u2scale
               end do
               do j = 1, np13(i)
                  uscale(ip13(j,i)) = u3scale
               end do
               do j = 1, np14(i)
                  uscale(ip14(j,i)) = u4scale
               end do
               do j = ii+1, npole
                  k = ipole(j)
                  gscale(k) = uscale(k) - pscale(k)
               end do
            end if
         else
            if (use_chgpen) then
               do j = 1, n12(i)
                  wscale(i12(j,i)) = w2scale
               end do
               do j = 1, n13(i)
                  wscale(i13(j,i)) = w3scale
               end do
               do j = 1, n14(i)
                  wscale(i14(j,i)) = w4scale
               end do
               do j = 1, n15(i)
                  wscale(i15(j,i)) = w5scale
               end do
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = d1scale
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = d2scale
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = d3scale
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = d4scale
               end do
               do j = ii+1, npole
                  k = ipole(j)
                  gscale(k) = wscale(k) - dscale(k)
               end do
            else
               do j = 1, np11(i)
                  gscale(ip11(j,i)) = u1scale - d1scale
               end do
               do j = 1, np12(i)
                  gscale(ip12(j,i)) = u2scale - d2scale
               end do
               do j = 1, np13(i)
                  gscale(ip13(j,i)) = u3scale - d3scale
               end do
               do j = 1, np14(i)
                  gscale(ip14(j,i)) = u4scale - d4scale
               end do
            end if
         end if
c
c     evaluate higher-numbered sites with the original site
c
         do kk = ii+1, npole
            k = ipole(kk)
            xr = x(k) - x(i)
            yr = y(k) - y(i)
            zr = z(k) - z(i)
            r2 = xr*xr + yr* yr + zr*zr
            r = sqrt(r2)
            ukx = uind(1,kk)
            uky = uind(2,kk)
            ukz = uind(3,kk)
c
c     intermediates involving moments and separation distance
c
            uir = xr*uix + yr*uiy + zr*uiz
            ukr = xr*ukx + yr*uky + zr*ukz
c
c     find the field components for Thole polarization damping
c
            if (use_thole) then
               scale3 = 1.0d0
               scale5 = 1.0d0
               damp = pdi * pdamp(kk)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(kk))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                  end if
               end if
c
c     find the field components for charge penetration damping
c
            else if (use_chgpen) then
               corek = pcore(kk)
               valk = pval(kk)
               alphak = palpha(kk)
               call dampmut (r,alphai,alphak,dmpik)
               scale3 = dmpik(3)
               scale5 = dmpik(5)
            end if
            rr3 = -scale3 / (r*r2)
            rr5 = 3.0d0 * scale5 / (r*r2*r2)
            fiu(1) = rr3*ukx + rr5*ukr*xr
            fiu(2) = rr3*uky + rr5*ukr*yr
            fiu(3) = rr3*ukz + rr5*ukr*zr
            fku(1) = rr3*uix + rr5*uir*xr
            fku(2) = rr3*uiy + rr5*uir*yr
            fku(3) = rr3*uiz + rr5*uir*zr
            do j = 1, 3
               field(j,ii) = field(j,ii) + fiu(j)*gscale(k)
               field(j,kk) = field(j,kk) + fku(j)*gscale(k)
            end do
         end do
c
c     reset exclusion coefficients for connected atoms
c
         if (dpequal) then
            if (use_chgpen) then
               do j = 1, n12(i)
                  gscale(i12(j,i)) = 0.0d0
               end do
               do j = 1, n13(i)
                  gscale(i13(j,i)) = 0.0d0
               end do
               do j = 1, n14(i)
                  gscale(i14(j,i)) = 0.0d0
               end do
               do j = 1, n15(i)
                  gscale(i15(j,i)) = 0.0d0
               end do
            else
               do j = 1, np11(i)
                  uscale(ip11(j,i)) = 1.0d0
                  gscale(ip11(j,i)) = 0.0d0
               end do
               do j = 1, np12(i)
                  uscale(ip12(j,i)) = 1.0d0
                  gscale(ip12(j,i)) = 0.0d0
               end do
               do j = 1, np13(i)
                  uscale(ip13(j,i)) = 1.0d0
                  gscale(ip13(j,i)) = 0.0d0
               end do
               do j = 1, np14(i)
                  uscale(ip14(j,i)) = 1.0d0
                  gscale(ip14(j,i)) = 0.0d0
               end do
               do j = 1, n12(i)
                  pscale(i12(j,i)) = 1.0d0
                  gscale(i12(j,i)) = 0.0d0
               end do
               do j = 1, n13(i)
                  pscale(i13(j,i)) = 1.0d0
                  gscale(i13(j,i)) = 0.0d0
               end do
               do j = 1, n14(i)
                  pscale(i14(j,i)) = 1.0d0
                  gscale(i14(j,i)) = 0.0d0
               end do
               do j = 1, n15(i)
                  pscale(i15(j,i)) = 1.0d0
                  gscale(i15(j,i)) = 0.0d0
               end do
            end if
         else
            if (use_chgpen) then
               do j = 1, np11(i)
                  dscale(ip11(j,i)) = 1.0d0
                  gscale(ip11(j,i)) = 0.0d0
               end do
               do j = 1, np12(i)
                  dscale(ip12(j,i)) = 1.0d0
                  gscale(ip12(j,i)) = 0.0d0
               end do
               do j = 1, np13(i)
                  dscale(ip13(j,i)) = 1.0d0
                  gscale(ip13(j,i)) = 0.0d0
               end do
               do j = 1, np14(i)
                  dscale(ip14(j,i)) = 1.0d0
                  gscale(ip14(j,i)) = 0.0d0
               end do
               do j = 1, n12(i)
                  wscale(i12(j,i)) = 1.0d0
                  gscale(i12(j,i)) = 0.0d0
               end do
               do j = 1, n13(i)
                  wscale(i13(j,i)) = 1.0d0
                  gscale(i13(j,i)) = 0.0d0
               end do
               do j = 1, n14(i)
                  wscale(i14(j,i)) = 1.0d0
                  gscale(i14(j,i)) = 0.0d0
               end do
               do j = 1, n15(i)
                  wscale(i15(j,i)) = 1.0d0
                  gscale(i15(j,i)) = 0.0d0
               end do
            else
               do j = 1, np11(i)
                  gscale(ip11(j,i)) = 0.0d0
               end do
               do j = 1, np12(i)
                  gscale(ip12(j,i)) = 0.0d0
               end do
               do j = 1, np13(i)
                  gscale(ip13(j,i)) = 0.0d0
               end do
               do j = 1, np14(i)
                  gscale(ip14(j,i)) = 0.0d0
               end do
            end if
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (dscale)
      deallocate (pscale)
      deallocate (uscale)
      deallocate (wscale)
      deallocate (gscale)
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
c     on equivalent attached atoms and additional user specified
c     sets of equivalent atoms
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
      use units
      implicit none
      integer i,j,k,m
      integer it,jt,kt,mt
      integer ia,ja,ka,ma
      integer in,jn
      integer mintyp
      integer size,nsing
      integer nsame,nave
      integer xaxe,yaxe
      integer zaxe
      integer list(20)
      integer, allocatable :: ising(:)
      integer, allocatable :: jsing(:)
      integer, allocatable :: isame(:)
      integer, allocatable :: tsort(:)
      integer, allocatable :: tkey(:)
      integer, allocatable :: pkey(:)
      integer, allocatable :: pgrt(:,:)
      real*8 pave(13)
      logical done,repeat
      logical header
      character*4 pa,pb,pc,pd
      character*16 ptlast
      character*16, allocatable :: pt(:)
      character*240 record
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (ising(n))
      allocate (jsing(n))
      allocate (isame(n))
      allocate (tsort(n))
      allocate (tkey(n))
      allocate (pkey(n))
      allocate (pt(n))
c
c     store the atomic numbers of atoms attached to each atom
c
      do i = 1, n
         do j = 1, 4
            list(j) = 0
         end do
         do j = 1, n12(i)
            list(j) = atomic(i12(j,i))
         end do
         call sort (n12(i),list)
         size = 4
         call numeral (list(1),pa,size)
         call numeral (list(2),pb,size)
         call numeral (list(3),pc,size)
         call numeral (list(4),pd,size)
         pt(i) = pa//pb//pc//pd
      end do
c
c     partially automated determination of equivalent atoms
c
      repeat = .true.
      dowhile (repeat)
         repeat = .false.
c
c     find the monovalent and singly attached atom groups
c
         nsing = 0
         do i = 1, n
            ising(i) = 0
            jsing(i) = 0
            k = 0
            m = 0
            do j = 1, n12(i)
               if (n12(i12(j,i)) .gt. 1) then
                  k = k + 1
                  m = i12(j,i)
               end if
            end do
            if (k .eq. 1) then
               nsing = nsing + 1
               ising(nsing) = i
               jsing(nsing) = m
            end if
         end do
c
c     condense equivalent attached atom groups to same type
c
         header = .true.
         do i = 1, nsing-1
            ia = ising(i)
            it = atomic(ia)
            in = type(jsing(i))
            do j = i+1, nsing
               ja = ising(j)
               jt = atomic(ja)
               jn = type(jsing(j))
               if (it.eq.jt .and. in.eq.jn
     &                .and. pt(ia).eq.pt(ja)) then
                  if (type(ia) .ne. type(ja)) then
                     mintyp = min(type(ia),type(ja))
                     type(ia) = mintyp
                     type(ja) = mintyp
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Equivalent Atom Pairs Condensed',
     &                             ' to the Same Atom Type :',/)
                     end if
                     write (iout,20)  ia,ja
   20                format (' Atoms',i6,2x,'and',i6,2x,'Set to',
     &                          ' Equivalent Types')
                  end if
                  do k = 1, n12(ia)
                     ka = i12(k,ia)
                     kt = atomic(ka)
                     do m = 1, n12(ja)
                        ma = i12(m,ja)
                        mt = atomic(ma)
                        if (kt .eq. mt) then
                           if (type(ka) .ne. type(ma)) then
                              mintyp = min(type(ka),type(ma))
                              type(ka) = mintyp
                              type(ma) = mintyp
                              if (header) then
                                 header = .false.
                                 write (iout,30)
   30                            format (/,' Equivalent Atoms to Set',
     &                                      ' to Same Atom Type :',/)
                              end if
                              write (iout,40)  ka,ma
   40                         format (' Atoms',i6,2x,'and',i6,2x,
     &                                   'are Equivalent')
                           end if
                        end if
                     end do
                  end do
               end if
            end do
         end do
c
c     query for sets of atoms to condense to a single type
c
         done = .false.
         dowhile (.not. done)
            do i = 1, 20
               list(i) = 0
            end do
            write (iout,50)
   50       format (/,' Enter Further Sets of Equivalent Atoms',
     &                 ' [<Enter>=Exit] :  ',$)
            read (input,60)  record
   60       format (a240)
            read (record,*,err=70,end=70)  (list(i),i=1,20)
   70       continue
c
c     process the input groups to a list of equivalent atoms
c
            nsame = 0
            do i = 1, n
               isame(i) = 0
            end do
            i = 1
            do while (list(i) .ne. 0)
               list(i) = max(-n,min(n,list(i)))
               if (list(i) .gt. 0) then
                  k = list(i)
                  nsame = nsame + 1
                  isame(nsame) = k
                  i = i + 1
               else
                  list(i+1) = max(-n,min(n,list(i+1)))
                  do k = abs(list(i)), abs(list(i+1))
                     nsame = nsame + 1
                     isame(nsame) = k
                  end do
                  i = i + 2
               end if
            end do
            if (nsame .eq. 0)  done = .true.
c
c     assign equivalent atoms to a common atom type number
c
            if (nsame .ne. 0) then
               repeat = .true.
               call sort8 (nsame,isame)
               k = type(isame(1))
               do i = 1, nsame
                  type(isame(i)) = k
               end do
            end if
         end do
c
c     renumber the atom types to remove deleted type numbers
c
         do i = 1, n
            tsort(i) = type(i)
         end do
         call sort3 (n,tsort,tkey)
         k = 0
         m = 0
         do i = 1, n
            if (tsort(i) .ne. k) then
               m = m + 1
               k = tsort(i)
            end if
            type(tkey(i)) = m
         end do
c
c     print the atoms, atom types and local frame definitions
c
         write (iout,80)
   80    format (/,' Atom Type and Local Frame Definition',
     &              ' for Each Atom :',
     &           //,5x,'Atom',4x,'Type',6x,'Local Frame',10x,
     &              'Frame Defining Atoms',/)
         do i = 1, n
            k = pollist(i)
            write (iout,90)  i,type(i),polaxe(k),zaxis(k),
     &                       xaxis(k),yaxis(k)
   90       format (2i8,9x,a8,6x,3i8)
         end do
      end do
c
c     locate the equivalently defined multipole sites
c
      do i = 1, npole
         k = ipole(i)
         it = type(k)
         zaxe = 0
         xaxe = 0
         yaxe = 0
         if (zaxis(i) .ne. 0)  zaxe = type(zaxis(i))
         if (xaxis(i) .ne. 0)  xaxe = type(xaxis(i))
         if (yaxis(i) .ne. 0)  yaxe = type(yaxis(i))
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
      deallocate (ising)
      deallocate (jsing)
      deallocate (isame)
      deallocate (tsort)
      deallocate (tkey)
      deallocate (pkey)
      deallocate (pt)
c
c     perform dynamic allocation of some local arrays
c
      allocate (pgrt(maxval,n))
c
c     set polarization groups over the condensed atom types
c
      do i = 1, n
         do j = 1, maxval
            pgrt(j,i) = 0
         end do
      end do
      do i = 1, npole
         it = type(ipole(i))
         do j = 1, maxval
            if (pgrp(j,i) .ne. 0)  k = j
            pgrp(j,it) = type(pgrp(j,i))
         end do
         call sort8 (k,pgrp(1,it))
         do j = 1, k
            do m = 1, maxval
               if (pgrt(m,it) .eq. 0) then
                  pgrt(m,it) = pgrp(j,it)
                  goto 100
               end if
            end do
  100       continue
         end do
      end do
      do i = 1, npole
         it = type(ipole(i))
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
c
c     convert dipole and quadrupole moments back to atomic units
c
      do i = 1, npole
         pole(1,i) = pole(1,i)
         do j = 2, 4
            pole(j,i) = pole(j,i) / bohr
         end do
         do j = 5, 13
            pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
         end do
      end do
c
c     regularize the multipole moments to standardized values
c
      call fixpole
c
c     print the final multipole values for force field use
c
      write (iout,110)
  110 format (/,' Multipole Moments for Polarizable Force Field :')
      do i = 1, n
         k = pollist(i)
         if (k .eq. 0) then
            write (iout,120)  i,name(i),atomic(i)
  120       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,130)
  130       format (/,' No Atomic Multipole Moments for this Site')
         else
            zaxe = zaxis(k)
            xaxe = xaxis(k)
            yaxe = yaxis(k)
            if (yaxe .lt. 0)  yaxe = -yaxe
            write (iout,140)  i,name(i),atomic(i)
  140       format (/,' Atom:',i8,9x,'Name:',3x,a3,
     &                 7x,'Atomic Number:',i8)
            write (iout,150)  polaxe(k),zaxe,xaxe,yaxe
  150       format (/,' Local Frame:',12x,a8,6x,3i8)
            write (iout,160)  pole(1,k)
  160       format (/,' Charge:',10x,f15.5)
            write (iout,170)  pole(2,k),pole(3,k),pole(4,k)
  170       format (' Dipole:',10x,3f15.5)
            write (iout,180)  pole(5,k)
  180       format (' Quadrupole:',6x,f15.5)
            write (iout,190)  pole(8,k),pole(9,k)
  190       format (18x,2f15.5)
            write (iout,200)  pole(11,k),pole(12,k),pole(13,k)
  200       format (18x,3f15.5)
         end if
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fixpole  --  regularize the multipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fixpole" removes multipole values that are zero by symmetry,
c     rounds moments to desired precision, and enforces integer net
c     charge and traceless quadrupoles
c
c
      subroutine fixpole
      use atoms
      use mpole
      implicit none
      integer i,j,k
      real*8 eps,big,sum
      real*8 ci,cj
      logical yzero
c
c
c     remove multipole components that are zero by symmetry
c
      do i = 1, npole
         yzero = .false.
         if (yaxis(i) .eq. 0)  yzero = .true.
         if (polaxe(i) .eq. 'Bisector')  yzero = .true.
         if (polaxe(i) .eq. 'Z-Bisect')  yzero = .true.
         if (zaxis(i).eq.0 .or. zaxis(i).gt.n) then
            pole(13,i) = 0.0d0
         end if
         if (xaxis(i).eq.0 .or. xaxis(i).gt.n) then
            pole(2,i) = 0.0d0
            pole(5,i) = -0.5d0 * pole(13,i)
            pole(7,i) = 0.0d0
            pole(9,i) = pole(5,i)
            pole(11,i) = 0.0d0
         end if
         if (yzero) then
            pole(3,i) = 0.0d0
            pole(6,i) = 0.0d0
            pole(8,i) = 0.0d0
            pole(10,i) = 0.0d0
            pole(12,i) = 0.0d0
         end if
      end do
c
c     regularize multipole moments to desired precision
c
      eps = 0.00001d0
      do i = 1, npole
         do j = 1, 13
            pole(j,i) = dble(nint(pole(j,i)/eps)) * eps
         end do
      end do
c
c     enforce integer net charge over atomic multipoles
c
      k = 0
      big = 0.0d0
      sum = 0.0d0
      do i = 1, npole
         sum = sum + pole(1,i)
         ci = abs(pole(1,i))
         if (ci .gt. big) then
            do j = 1, n
               cj = abs(pole(1,j))
               if (i.ne.j .and. ci.eq.cj)  goto 10
            end do
            k = i
            big = ci
   10       continue
         end if
      end do
      sum = sum - dble(nint(sum))
      if (k .ne. 0)  pole(1,k) = pole(1,k) - sum
c
c     enforce traceless quadrupole at each multipole site
c
      do i = 1, npole
         sum = pole(5,i) + pole(9,i) + pole(13,i)
         big = max(abs(pole(5,i)),abs(pole(9,i)),abs(pole(13,i)))
         k = 0
         if (big .eq. abs(pole(5,i)))  k = 5
         if (big .eq. abs(pole(9,i)))  k = 9
         if (big .eq. abs(pole(13,i)))  k = 13
         if (k .ne. 0)  pole(k,i) = pole(k,i) - sum
      end do
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
      use files
      use keys
      use kpolr
      use mpole
      use polar
      use sizes
      use units
      implicit none
      integer i,j,k,m
      integer it,jt,kt,mt
      integer ixyz,ikey
      integer size,tmax
      integer xaxe
      integer yaxe
      integer zaxe
      integer freeunit
      integer trimtext
      integer, allocatable :: pkey(:)
      character*4 pa,pb,pc,pd
      character*16 ptlast
      character*16, allocatable :: pt(:)
      character*240 keyfile
      character*240 record
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
      if (nkey .ne. 0) then
         write (ikey,20)
   20    format ()
      end if
c
c     output the atom definitions to the keyfile as appropriate
c
      tmax = 0
      do i = 1, n
         it = type(ipole(i))
         if (it .gt. tmax) then
            write (ikey,30)  type(i),type(i),name(i),story(i),
     &                       atomic(i),mass(i),valence(i)
   30       format ('atom',6x,2i5,4x,a3,3x,'"',a20,'"',i10,f10.3,i5)
         end if
         tmax = max(it,tmax)
      end do
      if (n .ne. 0) then
         write (ikey,40)
   40    format ()
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (pt(npole))
      allocate (pkey(npole))
c
c     locate the equivalently defined multipole sites
c
      do i = 1, npole
         k = ipole(i)
         it = type(k)
         zaxe = 0
         xaxe = 0
         yaxe = 0
         if (zaxis(i) .ne. 0)  zaxe = type(zaxis(i))
         if (xaxis(i) .ne. 0)  xaxe = type(xaxis(i))
         if (yaxis(i) .ne. 0)  yaxe = type(yaxis(i))
         size = 4
         call numeral (it,pa,size)
         call numeral (zaxe,pb,size)
         call numeral (xaxe,pc,size)
         call numeral (yaxe,pd,size)
         pt(i) = pa//pb//pc//pd
      end do
      call sort7 (npole,pt,pkey)
c
c     output the local frame multipole values to the keyfile
c
      ptlast = '                '
      do k = 1, npole
         i = pkey(k)
         it = type(ipole(i))
         if (pt(k) .ne. ptlast) then
            ptlast = pt(k)
            zaxe = type(zaxis(i))
            xaxe = type(xaxis(i))
            yaxe = type(yaxis(i))
            if (yaxe .lt. 0)  yaxe = -yaxe
            if (polaxe(i) .eq. 'None') then
               write (ikey,50)  it,pole(1,i)
   50          format ('multipole',1x,i5,21x,f11.5)
            else if (polaxe(i) .eq. 'Z-Only') then
               write (ikey,60)  it,zaxe,pole(1,i)
   60          format ('multipole',1x,2i5,16x,f11.5)
            else if (polaxe(i) .eq. 'Z-then-X') then
               if (yaxis(i) .eq. 0) then
                  write (ikey,70)  it,zaxe,xaxe,pole(1,i)
   70             format ('multipole',1x,3i5,11x,f11.5)
               else
                  write (ikey,80)  it,zaxe,xaxe,yaxe,pole(1,i)
   80             format ('multipole',1x,4i5,6x,f11.5)
               end if
            else if (polaxe(i) .eq. 'Bisector') then
               if (yaxis(i) .eq. 0) then
                  write (ikey,90)  it,-zaxe,-xaxe,pole(1,i)
   90             format ('multipole',1x,3i5,11x,f11.5)
               else
                  write (ikey,100)  it,-zaxe,-xaxe,yaxe,pole(1,i)
  100             format ('multipole',1x,4i5,6x,f11.5)
               end if
            else if (polaxe(i) .eq. 'Z-Bisect') then
               write (ikey,110)  it,zaxe,-xaxe,-yaxe,pole(1,i)
  110          format ('multipole',1x,4i5,6x,f11.5)
            else if (polaxe(i) .eq. '3-Fold') then
               write (ikey,120)  it,-zaxe,-xaxe,-yaxe,pole(1,i)
  120          format ('multipole',1x,4i5,6x,f11.5)
            end if
            write (ikey,130)  pole(2,i),pole(3,i),pole(4,i)
  130       format (36x,3f11.5)
            write (ikey,140)  pole(5,i)
  140       format (36x,f11.5)
            write (ikey,150)  pole(8,i),pole(9,i)
  150       format (36x,2f11.5)
            write (ikey,160)  pole(11,i),pole(12,i),pole(13,i)
  160       format (36x,3f11.5)
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (pt)
      deallocate (pkey)
c
c     output the polarizability parameters to the keyfile
c
      if (n .ne. 0) then
         write (ikey,170)
  170    format ()
      end if
      tmax = 0
      do i = 1, npole
         it = type(ipole(i))
         if (it .gt. tmax) then
            k = 0
            do j = 1, maxval
               if (pgrp(j,it) .ne. 0)  k = j
            end do
            call sort8 (k,pgrp(1,it))
            write (ikey,180)  it,polarity(i),thole(i),
     &                        (pgrp(j,it),j=1,k)
  180       format ('polarize',2x,i5,5x,2f11.4,2x,20i5)
         end if
         tmax = max(it,tmax)
      end do
      close (unit=ikey)
      return
      end
