c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program xyzedit  --  editing of Cartesian coordinates  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "xyzedit" provides for modification and manipulation
c     of the contents of Cartesian coordinates files
c
c
      program xyzedit
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use iounit
      use limits
      use math
      use molcul
      use ptable
      use titles
      use units
      use usage
      implicit none
      integer i,j,k,it
      integer ixyz,imod
      integer init,stop
      integer nmode,mode
      integer natom,atmnum
      integer nlist
      integer offset,origin
      integer oldtype,newtype
      integer freeunit
      integer trimtext
      integer, allocatable :: list(:)
      integer, allocatable :: keep(:)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 xcm,ycm,zcm
      real*8 xnew,ynew,znew
      real*8 xorig,yorig,zorig
      real*8 ri,rij,dij
      real*8 phi,theta,psi
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 dist2,cut2
      real*8 random,norm,weigh
      real*8, allocatable :: rad(:)
      real*8 a(3,3)
      logical exist,query
      logical opened,multi
      logical append
      character*240 xyzfile
      character*240 modfile
      character*240 record
      character*240 string
      external random,merge
c
c
c     initialize various constants and the output flags
c
      call initial
      opened = .false.
      multi = .false.
      nmode = 23
      offset = 0
c
c     try to get a filename from the command line arguments
c
      call nextarg (xyzfile,exist)
      if (exist) then
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Cartesian Coordinate File Name :  ',$)
         read (input,20)  xyzfile
   20    format (a240)
         call basefile (xyzfile)
         call suffix (xyzfile,'xyz','old')
         inquire (file=xyzfile,exist=exist)
      end do
c
c     open and then read the Cartesian coordinates file
c
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     get the force field definition and assign atom types
c
      call attach
      call field
      call katom
c
c     present a list of possible coordinate modifications
c
      write (iout,30)
   30 format (/,' The Tinker XYZ File Editing Utility Can :',
     &        //,4x,'(1) Offset the Numbers of the Current Atoms',
     &        /,4x,'(2) Remove User Specified Individual Atoms',
     &        /,4x,'(3) Remove User Specified Types of Atoms',
     &        /,4x,'(4) Delete Inactive Atoms Beyond Cutoff Range',
     &        /,4x,'(5) Insertion of Individual Specified Atoms',
     &        /,4x,'(6) Replace Old Atom Type with a New Type',
     &        /,4x,'(7) Assign Connectivities for Linear Chain',
     &        /,4x,'(8) Assign Connectivities Based on Distance',
     &        /,4x,'(9) Convert Units from Bohrs to Angstroms',
     &        /,3x,'(10) Invert thru Origin to Give Mirror Image',
     &        /,3x,'(11) Translate All Atoms by an X,Y,Z-Vector',
     &        /,3x,'(12) Translate Center of Mass to the Origin',
     &        /,3x,'(13) Translate a Specified Atom to the Origin',
     &        /,3x,'(14) Translate and Rotate to Inertial Frame',
     &        /,3x,'(15) Move to Specified Rigid Body Coordinates',
     &        /,3x,'(16) Move Stray Molecules into Periodic Box',
     &        /,3x,'(17) Trim a Periodic Box to a Smaller Size',
     &        /,3x,'(18) Make Truncated Octahedron from Cubic Box',
     &        /,3x,'(19) Make Rhombic Dodecahedron from Cubic Box',
     &        /,3x,'(20) Append a Second XYZ File to Current One',
     &        /,3x,'(21) Create and Fill a Periodic Boundary Box',
     &        /,3x,'(22) Soak Current Molecule in Box of Solvent',
     &        /,3x,'(23) Place Monoatomic Ions around a Solute')
c
c     get the desired type of coordinate file modification
c
   40 continue
      abort = .false.
      mode = -1
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=50,end=50)  mode
         if (mode.ge.0 .and. mode.le.nmode)  query = .false.
      end if
   50 continue
      if (query) then
         do while (mode.lt.0 .or. mode.gt.nmode)
            mode = 0
            write (iout,60)
   60       format (/,' Number of the Desired Choice [<Enter>=Exit]',
     &                 ' :  ',$)
            read (input,70,err=40,end=80)  mode
   70       format (i10)
   80       continue
         end do
      end if
c
c     open the file to be used for the output coordinates
c
      if (mode.gt.0 .and. .not.opened) then
         opened = .true.
         imod = freeunit ()
         modfile = filename(1:leng)//'.xyz'
         call version (modfile,'new')
         open (unit=imod,file=modfile,status='new')
      end if
c
c     get the offset value to be used in atom renumbering
c
      if (mode .eq. 1) then
   90    continue
         offset = 0
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=100,end=100)  offset
            query = .false.
         end if
  100    continue
         if (query) then
            write (iout,110)
  110       format (/,' Offset used to Renumber the Atoms [0] :  ',$)
            read (input,120,err=90)  offset
  120       format (i10)
         end if
         do while (.not. abort)
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     remove a specified list of individual atoms
c
      if (mode .eq. 2) then
         allocate (list(n))
         nlist = 0
         do i = 1, n
            list(i) = 0
         end do
         write (iout,130)
  130    format (/,' Numbers of the Atoms to be Removed :  ',$)
         read (input,140)  record
  140    format (a240)
         read (record,*,err=150,end=150)  (list(i),i=1,n)
  150    continue
         do while (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         do i = 1, nlist
            if (list(i) .gt. n)  list(i) = n
            if (list(i) .lt. -n)  list(i) = -n
         end do
         call sort4 (nlist,list)
         do while (.not. abort)
            do i = nlist, 1, -1
               if (i .gt. 1) then
                  if (list(i-1) .lt. 0) then
                     do j = abs(list(i)), abs(list(i-1)), -1
                        call delete (j)
                     end do
                  else if (list(i) .gt. 0) then
                     call delete (list(i))
                  end if
               else if (list(i) .gt. 0) then
                  call delete (list(i))
               end if
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     remove atoms with any of a specified list of atom types
c
      if (mode .eq. 3) then
         allocate (list(n))
         nlist = 0
         do i = 1, n
            list(i) = 0
         end do
         write (iout,160)
  160    format (/,' Atom Types to be Removed :  ',$)
         read (input,170)  record
  170    format (a240)
         read (record,*,err=180,end=180)  (list(i),i=1,n)
  180    continue
         do while (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         natom = n
         do while (.not. abort)
            do i = natom, 1, -1
               it = type(i)
               do j = 1, nlist
                  if (list(j) .eq. it) then
                     call delete (i)
                     goto 190
                  end if
               end do
  190          continue
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     remove atoms that are inactive and lie outside all cutoffs
c
      if (mode .eq. 4) then
         call active
         call cutoffs
         cut2 = 0.0d0
         if (vdwcut .le. 1000.0d0)  cut2 = max(vdwcut**2,cut2)
         if (chgcut .le. 1000.0d0)  cut2 = max(chgcut**2,cut2)
         if (dplcut .le. 1000.0d0)  cut2 = max(dplcut**2,cut2)
         if (mpolecut .le. 1000.0d0)  cut2 = max(mpolecut**2,cut2)
         if (cut2 .eq. 0.0d0)  cut2 = 1.0d16
         allocate (list(n))
         allocate (keep(n))
         do while (.not. abort)
            nlist = 0
            do i = 1, n
               keep(i) = 0
            end do
            do i = 1, n
               if (.not. use(i)) then
                  do j = 1, n12(i)
                     keep(i12(j,i)) = i
                  end do
                  do j = 1, n13(i)
                     keep(i13(j,i)) = i
                  end do
                  do j = 1, n14(i)
                     keep(i14(j,i)) = i
                  end do
                  xi = x(i)
                  yi = y(i)
                  zi = z(i)
                  do j = 1, n
                     if (use(j)) then
                        if (keep(j) .eq. i)  goto 200
                        dist2 = (x(j)-xi)**2+(y(j)-yi)**2+(z(j)-zi)**2
                        if (dist2 .le. cut2)  goto 200
                     end if
                  end do
                  nlist = nlist + 1
                  list(nlist) = i
  200             continue
               end if
            end do
            do i = nlist, 1, -1
               call delete (list(i))
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         deallocate (keep)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     insert a specified list of individual atoms
c
      if (mode .eq. 5) then
         allocate (list(n))
         nlist = 0
         do i = 1, n
            list(i) = 0
         end do
         write (iout,210)
  210    format (/,' Numbers of the Atoms to be Inserted :  ',$)
         read (input,220)  record
  220    format (a240)
         read (record,*,err=230,end=230)  (list(i),i=1,n)
  230    continue
         do while (list(nlist+1) .ne. 0)
            nlist = nlist + 1
         end do
         call sort4 (nlist,list)
         do while (.not. abort)
            do i = nlist, 1, -1
               if (i .gt. 1) then
                  if (list(i-1) .lt. 0) then
                     do j = abs(list(i-1)), abs(list(i))
                        call insert (j)
                     end do
                  else if (list(i) .gt. 0) then
                     call insert (list(i))
                  end if
               else if (list(i) .gt. 0) then
                  call insert (list(i))
               end if
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     get an old atom type and new atom type for replacement
c
      if (mode .eq. 6) then
  240    continue
         oldtype = 0
         newtype = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  oldtype
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  newtype
  250    continue
         if (oldtype.eq.0 .or. newtype.eq.0) then
            write (iout,260)
  260       format (/,' Numbers of the Old and New Atom Types :  ',$)
            read (input,270)  record
  270       format (a240)
         end if
         read (record,*,err=240,end=240)  oldtype,newtype
         do while (.not. abort)
            do i = 1, n
               if (type(i) .eq. oldtype) then
                  type(i) = newtype
               end if
            end do
            call katom
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     assign atom connectivities to produce a linear chain
c
      if (mode .eq. 7) then
         do while (.not. abort)
            do i = 1, n
               n12(i) = 0
               if (i .ne. 1) then
                  n12(i) = n12(i) + 1
                  i12(n12(i),i) = i - 1
               end if
               if (i .ne. n) then
                  n12(i) = n12(i) + 1
                  i12(n12(i),i) = i + 1
               end if
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     assign atom connectivities based on interatomic distances
c
      if (mode .eq. 8) then
         allocate (rad(n))
         do while (.not. abort)
            call unitcell
            call lattice
            do i = 1, n
               rad(i) = 0.76d0
               atmnum = atomic(i)
               if (atmnum .ne. 0)  rad(i) = covrad(atmnum)
               rad(i) = 1.15d0 * rad(i)
            end do
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
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (rad)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     convert the coordinate units from Bohrs to Angstroms
c
      if (mode .eq. 9) then
         do while (.not. abort)
            do i = 1, n
               x(i) = x(i) * bohr
               y(i) = y(i) * bohr
               z(i) = z(i) * bohr
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     get mirror image by inverting coordinates through origin
c
      if (mode .eq. 10) then
         do while (.not. abort)
            do i = 1, n
               x(i) = -x(i)
               y(i) = -y(i)
               z(i) = -z(i)
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     translate the entire system by a specified x,y,z-vector
c
      if (mode .eq. 11) then
         xr = 0.0d0
         yr = 0.0d0
         zr = 0.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=280,end=280)  xr
         call nextarg (string,exist)
         if (exist)  read (string,*,err=280,end=280)  yr
         call nextarg (string,exist)
         if (exist)  read (string,*,err=280,end=280)  zr
  280    continue
         if (xr.eq.0.0d0 .and. yr.eq.0.0d0 .and. zr.eq.0.0d0) then
            write (iout,290)
  290       format (/,' Enter Translation Vector Components :  ',$)
            read (input,300)  record
  300       format (a240)
            read (record,*,err=310,end=310)  xr,yr,zr
  310       continue
         end if
         do while (.not. abort)
            do i = 1, n
               x(i) = x(i) + xr
               y(i) = y(i) + yr
               z(i) = z(i) + zr
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     translate the center of mass to the coordinate origin
c
      if (mode .eq. 12) then
         do while (.not. abort)
            xcm = 0.0d0
            ycm = 0.0d0
            zcm = 0.0d0
            norm = 0.0d0
            do i = 1, n
               weigh = mass(i)
               xcm = xcm + x(i)*weigh
               ycm = ycm + y(i)*weigh
               zcm = zcm + z(i)*weigh
               norm = norm + weigh
            end do
            xcm = xcm / norm
            ycm = ycm / norm
            zcm = zcm / norm
            do i = 1, n
               x(i) = x(i) - xcm
               y(i) = y(i) - ycm
               z(i) = z(i) - zcm
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     translate to place a specified atom at the origin
c
      if (mode .eq. 13) then
         origin = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=320,end=320)  origin
  320    continue
         if (origin .eq. 0) then
            write (iout,330)
  330       format (/,' Number of the Atom to Move to the Origin',
     &                 ' :  ',$)
            read (input,340)  origin
  340       format (i10)
         end if
         do while (.not. abort)
            xorig = x(origin)
            yorig = y(origin)
            zorig = z(origin)
            do i = 1, n
               x(i) = x(i) - xorig
               y(i) = y(i) - yorig
               z(i) = z(i) - zorig
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     translate and rotate into standard orientation
c
      if (mode .eq. 14) then
         do while (.not. abort)
            call inertia (2)
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     translate and rotate to specified rigid body coordinates
c
      if (mode .eq. 15) then
         xcm = 0.0d0
         ycm = 0.0d0
         zcm = 0.0d0
         phi = 0.0d0
         theta = 0.0d0
         psi = 0.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  xcm
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  ycm
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  zcm
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  phi
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  theta
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  psi
  350    continue
         if (min(xcm,ycm,zcm,phi,theta,psi).eq.0.0d0 .and.
     &       max(xcm,ycm,zcm,phi,theta,psi).eq.0.0d0) then
            write (iout,360)
  360       format (/,' Enter Rigid Body Coordinates :  ',$)
            read (input,370)  record
  370       format (a240)
            read (record,*,err=380,end=380)  xcm,ycm,zcm,phi,theta,psi
  380       continue
         end if
         call inertia (2)
         phi = phi / radian
         theta = theta / radian
         psi = psi / radian
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
         do while (.not. abort)
            do i = 1, n
               xorig = x(i)
               yorig = y(i)
               zorig = z(i)
               x(i) = a(1,1)*xorig + a(2,1)*yorig + a(3,1)*zorig + xcm
               y(i) = a(1,2)*xorig + a(2,2)*yorig + a(3,2)*zorig + ycm
               z(i) = a(1,3)*xorig + a(2,3)*yorig + a(3,3)*zorig + zcm
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     move stray molecules back into original periodic box
c
      if (mode .eq. 16) then
         do while (.not. abort)
            call unitcell
            if (use_bounds) then
               call lattice
               call molecule
               call bounds
            end if
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     remove molecules to trim periodic box to smaller size
c
      if (mode .eq. 17) then
         xnew = 0.0d0
         ynew = 0.0d0
         znew = 0.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=390,end=390)  xnew
         call nextarg (string,exist)
         if (exist)  read (string,*,err=390,end=390)  ynew
         call nextarg (string,exist)
         if (exist)  read (string,*,err=390,end=390)  znew
  390    continue
         do while (xnew .eq. 0.0d0)
            write (iout,400)
  400       format (/,' Enter Periodic Box Dimensions (X,Y,Z) :  ',$)
            read (input,410)  record
  410       format (a240)
            read (record,*,err=420,end=420)  xnew,ynew,znew
  420       continue
         end do
         if (ynew .eq. 0.0d0)  ynew = xnew
         if (znew .eq. 0.0d0)  znew = xnew
         xbox = xnew
         ybox = ynew
         zbox = znew
         call lattice
         call molecule
         allocate (list(n))
         allocate (keep(n))
         do while (.not. abort)
            do i = 1, n
               list(i) = 1
            end do
            do i = 1, nmol
               init = imol(1,i)
               stop = imol(2,i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = init, stop
                  k = kmol(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               weigh = molmass(i)
               xcm = xcm / weigh
               ycm = ycm / weigh
               zcm = zcm / weigh
               if (abs(xcm).gt.xbox2 .or. abs(ycm).gt.ybox2
     &                   .or. abs(zcm).gt.zbox2) then
                  do j = init, stop
                     k = kmol(j)
                     list(k) = 0
                  end do
               end if
            end do
            k = 0
            do i = 1, n
               if (list(i) .ne. 0) then
                  k = k + 1
                  keep(k) = i
                  list(i) = k
               end if
            end do
            n = k
            do k = 1, n
               i = keep(k)
               name(k) = name(i)
               x(k) = x(i)
               y(k) = y(i)
               z(k) = z(i)
               type(k) = type(i)
               n12(k) = n12(i)
               do j = 1, n12(k)
                  i12(j,k) = list(i12(j,i))
               end do
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort) then
               multi = .true.
               xbox = xnew
               ybox = ynew
               zbox = znew
               call lattice
               call molecule
            end if
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         deallocate (keep)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     trim cube to truncated octahedron or rhombic dodecahedron
c
      if (mode.eq.18 .or. mode.eq.19) then
         call unitcell
         dowhile (xbox .eq. 0.0d0)
            write (iout,430)
  430       format (/,' Enter Edge Length of Cubic Periodic Box :  ',$)
            read (input,440)  record
  440       format (a240)
            read (record,*,err=450,end=450)  xbox
  450       continue
         end do
         if (mode .eq. 18)  octahedron = .true.
         if (mode .eq. 19)  dodecadron = .true.
         call molecule
         allocate (list(n))
         allocate (keep(n))
         do while (.not. abort)
            do i = 1, n
               list(i) = 1
            end do
            do i = 1, nmol
               init = imol(1,i)
               stop = imol(2,i)
               xcm = 0.0d0
               ycm = 0.0d0
               zcm = 0.0d0
               do j = init, stop
                  k = kmol(j)
                  weigh = mass(k)
                  xcm = xcm + x(k)*weigh
                  ycm = ycm + y(k)*weigh
                  zcm = zcm + z(k)*weigh
               end do
               weigh = molmass(i)
               xcm = xcm / weigh
               ycm = ycm / weigh
               zcm = zcm / weigh
               if (octahedron) then
                  xcm = xcm - xbox*nint(xcm/xbox)
                  ycm = ycm - ybox*nint(ycm/ybox)
                  zcm = zcm - zbox*nint(zcm/zbox)
                  if (abs(xcm)+abs(ycm)+abs(zcm) .gt. 0.75*xbox) then
                     do j = init, stop
                        k = kmol(j)
                        list(k) = 0
                     end do
                  end if
               else if (dodecadron) then
                  xcm = xcm - xbox*nint(xcm/xbox)
                  ycm = ycm - ybox*nint(ycm/ybox)
                  zcm = zcm - root2*zbox*nint(zcm/(zbox*root2))
                  if (abs(xcm)+abs(ycm)+abs(root2*zcm) .gt. xbox) then
                     do j = init, stop
                        k = kmol(j)
                        list(k) = 0
                     end do
                  end if
               end if
            end do
            k = 0
            do i = 1, n
               if (list(i) .ne. 0) then
                  k = k + 1
                  keep(k) = i
                  list(i) = k
               end if
            end do
            n = k
            do k = 1, n
               i = keep(k)
               name(k) = name(i)
               x(k) = x(i)
               y(k) = y(i)
               z(k) = z(i)
               type(k) = type(i)
               n12(k) = n12(i)
               do j = 1, n12(k)
                  i12(j,k) = list(i12(j,i))
               end do
            end do
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         deallocate (list)
         deallocate (keep)
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     append a second file to the current coordinates file
c
      if (mode .eq. 20) then
         append = .false.
         do while (.not. abort)
            call makeref (1)
            if (append) then
               call getref (3)
            else
               call getxyz
               call makeref (3)
               append = .true.
            end if
            call merge (1)
            call makeref (1)
            call readxyz (ixyz)
            if (.not. abort)  multi = .true.
            if (multi) then
               call makeref (2)
               call getref (1)
               call prtmod (imod,offset)
               call getref (2)
            end if
         end do
         if (.not. multi) then
            call getref (1)
            goto 40
         end if
      end if
c
c     create random box full of the current coordinates file
c
      if (mode .eq. 21) then
         call makebox
      end if
c
c     solvate the current system by insertion into a solvent box
c
      if (mode .eq. 22) then
         call soak
      end if
c
c     replace random solvent molecules outside solute with ions
c
      if (mode .eq. 23) then
         call molecule
         call addions
      end if
c
c     output final coordinates for single frame and print info
c
      if (opened .and. .not.multi) then
         call prtmod (imod,offset)
      end if
      if (opened) then
         close (unit=imod)
         write (iout,460)  modfile(1:trimtext(modfile))
  460    format (/,' New Coordinates Written To :  ',a)
      end if
      close (unit=ixyz)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine prtmod  --  Cartesian coordinates with offset  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "prtmod" writes out a set of modified Cartesian coordinates
c     with an optional atom number offset to an external disk file
c
c
      subroutine prtmod (imod,offset)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use inform
      use titles
      implicit none
      integer i,k,imod
      integer offset
      integer size,crdsiz
      real*8 crdmin,crdmax
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
c
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n+offset .ge. 100000)  atmc = 'i7'
      if (n+offset .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0d0)  crdsiz = 7
      if (crdmax .ge. 10000.0d0)  crdsiz = 7
      if (crdmin .le. -10000.0d0)  crdsiz = 8
      if (crdmax .ge. 100000.0d0)  crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral (crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (imod,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (imod,fstr(1:9))  n,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (imod,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the coordinate line for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         write (imod,fstr)  i+offset,name(i),x(i),y(i),z(i),type(i),
     &                      (i12(k,i)+offset,k=1,n12(i))
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine makebox  --  build periodic box from monomers  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "makebox" builds a periodic box of a desired size by randomly
c     copying a specified number of monomers into a target box size,
c     followed by optional excluded volume refinement
c
c
      subroutine makebox
      use atomid
      use atoms
      use boxes
      use couple
      use iounit
      implicit none
      integer i,j,k,m
      integer ncopy
      integer offset
      real*8 xcm,ycm,zcm
      real*8 phi,theta,psi
      real*8 cphi,ctheta,cpsi
      real*8 sphi,stheta,spsi
      real*8 random,reduce
      real*8 norm,weigh
      real*8, allocatable :: x0(:)
      real*8, allocatable :: y0(:)
      real*8, allocatable :: z0(:)
      real*8 a(3,3)
      logical exist,query
      logical refine
      character*1 answer
      character*240 record
      character*240 string
c
c
c     get the number of copies of the monomer to be used
c
      ncopy = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  ncopy
   10 continue
      if (ncopy .eq. 0)  then
         write (iout,20)
   20    format (/,' Enter Number of Copies to Put in Box :  ',$)
         read (input,30)  ncopy
   30    format (i10)
      end if
c
c     find the size of the periodic box to be constructed
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  xbox
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  ybox
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  zbox
   40 continue
      do while (xbox .eq. 0.0d0)
         write (iout,50)
   50    format (/,' Enter Periodic Box Dimensions (X,Y,Z) :  ',$)
         read (input,60)  record
   60    format (a240)
         read (record,*,err=70,end=70)  xbox,ybox,zbox
   70    continue
      end do
      if (ybox .eq. 0.0d0)  ybox = xbox
      if (zbox .eq. 0.0d0)  zbox = xbox
      orthogonal = .true.
c
c     decide whether to use excluded volume refinement
c
      refine = .true.
      answer = 'Y'
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=80,end=80)  answer
         query = .false.
      end if
   80 continue
      if (query) then
         write (iout,90)
   90    format (/,' Refine the Periodic Box Configuration',
     &              ' [Y] :  ',$)
         read (input,100)  answer
  100    format (a1)
      end if
      call upcase (answer)
      if (answer .eq. 'N')  refine = .false.
c
c     center the monomer and reduce its size to avoid overlap
c
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      norm = 0.0d0
      do i = 1, n
         weigh = mass(i)
         xcm = xcm + x(i)*weigh
         ycm = ycm + y(i)*weigh
         zcm = zcm + z(i)*weigh
         norm = norm + weigh
      end do
      xcm = xcm / norm
      ycm = ycm / norm
      zcm = zcm / norm
      allocate (x0(n))
      allocate (y0(n))
      allocate (z0(n))
      reduce = 0.001d0
      do i = 1, n
         x(i) = x(i) - xcm
         y(i) = y(i) - ycm
         z(i) = z(i) - zcm
         if (refine) then
            x(i) = reduce * x(i)
            y(i) = reduce * y(i)
            z(i) = reduce * z(i)
         end if
         x0(i) = x(i)
         y0(i) = y(i)
         z0(i) = z(i)
      end do
c
c     randomly place monomer copies in the periodic box
c
      do k = 1, ncopy
         offset = (k-1) * n
         xcm = xbox * (random()-0.5d0)
         ycm = ybox * (random()-0.5d0)
         zcm = zbox * (random()-0.5d0)
         phi = 360.0d0 * random ()
         theta = 360.0d0 * random ()
         psi = 360.0d0 * random ()
         cphi = cos(phi)
         sphi = sin(phi)
         ctheta = cos(theta)
         stheta = sin(theta)
         cpsi = cos(psi)
         spsi = sin(psi)
         a(1,1) = ctheta * cphi
         a(2,1) = spsi*stheta*cphi - cpsi*sphi
         a(3,1) = cpsi*stheta*cphi + spsi*sphi
         a(1,2) = ctheta * sphi
         a(2,2) = spsi*stheta*sphi + cpsi*cphi
         a(3,2) = cpsi*stheta*sphi - spsi*cphi
         a(1,3) = -stheta
         a(2,3) = ctheta * spsi
         a(3,3) = ctheta * cpsi
         do i = 1, n
            j = i + offset
            name(j) = name(i)
            type(j) = type(i)
            mass(j) = mass(i)
            x(j) = a(1,1)*x0(i) + a(2,1)*y0(i) + a(3,1)*z0(i) + xcm
            y(j) = a(1,2)*x0(i) + a(2,2)*y0(i) + a(3,2)*z0(i) + ycm
            z(j) = a(1,3)*x0(i) + a(2,3)*y0(i) + a(3,3)*z0(i) + zcm
            n12(j) = n12(i)
            do m = 1, n12(i)
               i12(m,j) = i12(m,i) + offset
            end do
         end do
      end do
      deallocate (x0)
      deallocate (y0)
      deallocate (z0)
      offset = 0
      n = ncopy * n
c
c     optionally perform excluded volume coordinate refinement
c
      if (refine) then
         call boxmin
         call bounds
      else
         call lattice
         call molecule
         call bounds
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine boxmin  --  expand molecules into periodic box  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "boxmin" uses minimization of valence and vdw potential energy
c     to expand and refine a collection of solvent molecules in a
c     periodic box
c
c
      subroutine boxmin
      use atomid
      use atoms
      use boxes
      use inform
      use limits
      use linmin
      use minima
      use neigh
      use output
      use potent
      use scales
      use vdw
      use vdwpot
      implicit none
      integer i,j,k,nvar
      integer ii,kk
      real*8 minimum
      real*8 boxmin1
      real*8 grdmin
      real*8 boxsiz
      real*8, allocatable :: xx(:)
      external boxmin1
      external optsave
c
c
c     setup for minimization with only valence and vdw terms
c
      call mechanic
      call potoff
      use_bond = .true.
      use_angle = .true.
      use_opbend = .true.
      use_opdist = .true.
      use_improp = .true.
      use_imptor = .true.
      use_tors = .true.
      use_vdw = .true.
c
c     set artificial Lennard-Jones vdw values for the system
c
      vdwtyp = 'LENNARD-JONES'
      nvdw = n
      do i = 1, n
         ivdw(i) = i
         jvdw(i) = class(i)
         ired(i) = i
      end do
      do i = 1, n-1
         ii = jvdw(i)
         do k = i+1, n
            kk = jvdw(k)
            if (atomic(i).eq.1 .and. atomic(k).eq.1) then
               radmin(ii,kk) = 2.90d0
               epsilon(ii,kk) = 0.016d0
               radmin4(ii,kk) = 2.90d0
               epsilon4(ii,kk) = 0.016d0
            else if (atomic(i).eq.1 .or. atomic(k).eq.1) then
               radmin(ii,kk) = 3.35d0
               epsilon(ii,kk) = 0.040d0
               radmin4(ii,kk) = 3.35d0
               epsilon4(ii,kk) = 0.040d0
            else
               radmin(ii,kk) = 3.80d0
               epsilon(ii,kk) = 0.100d0
               radmin4(ii,kk) = 3.80d0
               epsilon4(ii,kk) = 0.100d0
            end if
         end do
      end do
c
c     cutoff values and neighbor lists for vdw interactions
c
      use_list = .false.
      use_vlist = .false.
      vdwcut = 5.0d0
      vdwtaper = 4.5d0
      lbuffer = 1.0d0
      boxsiz = min(xbox,ybox,zbox)
      if (boxsiz .gt. 2.0d0*(vdwcut+lbuffer)) then
         use_list = .true.
         use_vlist = .true.
         dovlst = .true.
         lbuf2 = (0.5d0*lbuffer)**2
         vbuf2 = (vdwcut+lbuffer)**2
         vbufx = (vdwcut+2.0d0*lbuffer)**2
         maxvlst = int(sqrt(vbuf2)**3) + 100
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (use_vlist) then
         if (allocated(nvlst))  deallocate (nvlst)
         if (allocated(vlst))  deallocate (vlst)
         if (allocated(xvold))  deallocate (xvold)
         if (allocated(yvold))  deallocate (yvold)
         if (allocated(zvold))  deallocate (zvold)
         allocate (nvlst(n))
         allocate (vlst(maxvlst,n))
         allocate (xvold(n))
         allocate (yvold(n))
         allocate (zvold(n))
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (.not. allocated(scale))  allocate (scale(3*n))
c
c     mark for use of all atoms, and set scale factors
c
      nvar = 0
      do i = 1, n
         do j = 1, 3
            nvar = nvar + 1
            scale(nvar) = 12.0d0
         end do
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (xx(nvar))
c
c     scale the coordinates of each active atom
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         xx(nvar) = x(i) * scale(nvar)
         nvar = nvar + 1
         xx(nvar) = y(i) * scale(nvar)
         nvar = nvar + 1
         xx(nvar) = z(i) * scale(nvar)
      end do
c
c     make the call to the optimization routine
c
      iprint = 100
      maxiter = 10000
      stpmax = 10.0d0
      grdmin = 1.0d0
      coordtype = 'NONE'
      call lbfgs (nvar,xx,minimum,grdmin,boxmin1,optsave)
c
c     unscale the final coordinates for active atoms
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar) / scale(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar) / scale(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar) / scale(nvar)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xx)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  function boxmin1  --  energy and gradient for boxmin  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "boxmin1" is a service routine that computes the energy and
c     gradient during refinement of a periodic box
c
c
      function boxmin1 (xx,g)
      use atoms

      use energi
      use potent
      use repel

      use scales
      implicit none
      integer i,nvar
      real*8 e,boxmin1
      real*8 xx(*)
      real*8 g(*)
      real*8, allocatable :: derivs(:,:)
c
c
c     convert optimization parameters to atomic coordinates
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         x(i) = xx(nvar) / scale(nvar)
         nvar = nvar + 1
         y(i) = xx(nvar) / scale(nvar)
         nvar = nvar + 1
         z(i) = xx(nvar) / scale(nvar)
      end do
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     compute and store the energy and gradient
c
      call gradient (e,derivs)
      boxmin1 = e
c
c     convert gradient components to optimization parameters
c
      nvar = 0
      do i = 1, n
         nvar = nvar + 1
         g(nvar) = derivs(1,i) / scale(nvar)
         nvar = nvar + 1
         g(nvar) = derivs(2,i) / scale(nvar)
         nvar = nvar + 1
         g(nvar) = derivs(3,i) / scale(nvar)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine soak  --  place a solute into a solvent box  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "soak" takes a currently defined solute system and places
c     it into a solvent box, with removal of any solvent molecules
c     that overlap the solute
c
c
      subroutine soak
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use iounit
      use molcul
      use refer
      implicit none
      integer i,j,k
      integer ii,jj
      integer n12i,n12k
      integer isolv,icount
      integer ntot,freeunit
      integer, allocatable :: map(:)
      real*8 xi,yi,zi
      real*8 xr,yr,zr
      real*8 rik2,close2
      real*8 dxx,dxx2
      real*8 dxh,dxh2
      real*8 dhh,dhh2
      logical exist,header
      logical, allocatable :: remove(:)
      character*240 solvfile
      external merge
c
c
c     make a copy of the solute coordinates and connectivities
c
      call makeref (1)
c
c     get the filename for the solvent box coordinates
c
      call nextarg (solvfile,exist)
      if (exist) then
         call basefile (solvfile)
         call suffix (solvfile,'xyz','old')
         inquire (file=solvfile,exist=exist)
      end if
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Name of Solvent Box Coordinates :  ',$)
         read (input,20)  solvfile
   20    format (a240)
         call basefile (solvfile)
         call suffix (solvfile,'xyz','old')
         inquire (file=solvfile,exist=exist)
      end do
c
c     read the coordinate file containing the solvent atoms
c
      isolv = freeunit ()
      open (unit=isolv,file=solvfile,status='old')
      rewind (unit=isolv)
      call readxyz (isolv)
      close (unit=isolv)
c
c     combine solute and solvent into a single coordinate set
c
      call merge (1)
      call basefile (solvfile)
      call getkey
c
c     reset the default values for unitcell dimensions
c
      xbox = 0.0d0
      ybox = 0.0d0
      zbox = 0.0d0
      alpha = 0.0d0
      beta = 0.0d0
      gamma = 0.0d0
c
c     count number of molecules and set lattice parameters
c
      call molecule
      call unitcell
      call lattice
c
c     set distance cutoffs for solute-solvent close contacts
c
      dxx = 2.40d0
      dxh = 2.19d0
      dhh = 1.82d0
      dxx2 = dxx * dxx
      dxh2 = dxh * dxh
      dhh2 = dhh * dhh
c
c     perform dynamic allocation of some local arrays
c
      allocate (map(n))
      allocate (remove(nmol))
c
c     initialize the list of solvent molecules to be deleted
c
      do i = 1, nmol
         remove(i) = .false.
      end do
c
c     print header information when processing large systems
c
      icount = 0
      header = .true.
      if (n-nref(1) .ge. 10000) then
         write (iout,30)
   30    format (/,' Scan for Solvent Molecules to be Removed :')
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(nref,n,x,y,z,n12,molcule,dxx2,dxh2,dhh2,remove,
!$OMP& header,icount,iout)
!$OMP DO schedule(guided)
c
c     search for close contacts between solute and solvent
c
      do i = nref(1)+1, n
         if (.not. remove(molcule(i))) then
            xi = x(i)
            yi = y(i)
            zi = z(i)
            n12i = n12(i)
            do k = 1, nref(1)
               n12k = n12(k)
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               call imagen (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
               if (n12i.gt.1 .and. n12k.gt.1) then
                  close2 = dxx2
               else if (n12i.gt.1 .or. n12k.gt.1) then
                  close2 = dxh2
               else
                  close2 = dhh2
               end if
               if (rik2 .lt. close2) then
                  remove(molcule(i)) = .true.
                  goto 40
               end if
            end do
   40       continue
         end if
         icount = icount + 1
         if (mod(icount,10000) .eq. 0) then
            if (header) then
               header = .false.
               write (iout,50)
   50          format ()
            end if
            write (iout,60)  10000*(icount/10000)
   60       format (' Solvent Atoms Processed',i15)
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     print final status when processing large systems
c
      icount = n - nref(1)
      if (mod(icount,10000).ne.0 .and. icount.gt.10000) then
         write (iout,70)  icount
   70    format (' Solvent Atoms Processed',i15)
      end if
c
c     delete solvent molecules that are too close to the solute
c
      k = nref(1)
      ntot = k
      do i = nref(1)+1, n
         map(i) = 0
         if (.not. remove(molcule(i))) then
            k = k + 1
            map(i) = k
            ntot = k
         end if
      end do
      do i = nref(1)+1, n
         ii = map(i)
         if (ii .ne. 0) then
            x(ii) = x(i)
            y(ii) = y(i)
            z(ii) = z(i)
            name(ii) = name(i)
            type(ii) = type(i)
            k = 0
            do j = 1, n12(i)
               jj = map(i12(j,i))
               if (jj .ne. 0) then
                  k = k + 1
                  i12(k,ii) = jj
               end if
            end do
            n12(ii) = k
         end if
      end do
      n = ntot
c
c     perform deallocation of some local arrays
c
      deallocate (map)
      deallocate (remove)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine addions  --  placement of ions around solute  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "addions" takes a currently defined solvated system and
c     places ions, with removal of solvent molecules
c
c
      subroutine addions
      use atomid
      use atoms
      use couple
      use iounit
      use katoms
      use molcul
      implicit none
      integer i,j,k
      integer nsolute,size
      integer start,stop
      integer icount,iontyp
      integer ncopy,ranatm
      integer, allocatable :: list(:)
      integer, allocatable :: isolute(:)
      real*8 xi,yi,zi
      real*8 xr,yr,zr,rik2
      real*8 close,close2
      real*8 rand,random
      real*8 weigh,xmid,ymid,zmid
      real*8, allocatable :: xion(:)
      real*8, allocatable :: yion(:)
      real*8, allocatable :: zion(:)
      logical exist,header,done
      logical, allocatable :: remove(:)
      character*240 record
      character*240 string
      external random
c
c
c     perform dynamic allocation of some local arrays
c
      size = 40
      allocate (list(size))
      allocate (isolute(n))
c
c     get the range atoms numbers constituting the solute
c
      do i = 1, size
         list(i) = 0
      end do
      i = 0
      do while (exist)
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=10,end=10)  list(i+1)
            i = i + 1
         end if
      end do
   10 continue
      if (i .eq. 0) then
         write (iout,20)
   20    format (/,' Enter Atom Numbers in Solute Molecules :  ',$)
         read (input,30)  record
   30    format (a240)
         read (record,*,err=40,end=40)  (list(i),i=1,size)
   40    continue
      end if
      i = 1
      nsolute = 0
      do while (list(i) .ne. 0)
         list(i) = max(-n,min(n,list(i)))
         if (list(i) .gt. 0) then
            k = list(i)
            nsolute = nsolute + 1
            isolute(nsolute) = k
            i = i + 1
         else
            list(i+1) = max(-n,min(n,list(i+1)))
            do k = abs(list(i)), abs(list(i+1))
               nsolute = nsolute + 1
               isolute(nsolute) = k
            end do
            i = i + 2
         end if
      end do
c
c     get the atom type of ion to be added and number of copies
c
   50 continue
      iontyp = 0
      ncopy = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  iontyp
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  ncopy
   60 continue
      if (iontyp.eq.0 .or. ncopy.eq.0) then
         write (iout,70)
   70    format (/,' Enter Ion Atom Type Number & Copies to Add :  ',$)
         read (input,80)  record
   80    format (a240)
      end if
      read (record,*,err=50,end=50)  iontyp,ncopy
c
c     set minimum distance cutoff for solute-ion contacts
c
      close = 6.0d0
      close2 = close * close
c
c     perform dynamic allocation of some local arrays
c
      allocate (remove(nmol))
      allocate (xion(ncopy))
      allocate (yion(ncopy))
      allocate (zion(ncopy))
c
c     initialize the list of solvent molecules to be deleted
c
      do i = 1, nmol
         remove(i) = .false.
      end do
c
c     print header information when processing large systems
c
      icount = 0
      header = .true.
      if (n .ge. 10000) then
         write (iout,90)
   90    format (/,' Scan for Available Locations to Place Ions :')
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(n,x,y,z,molcule,close2,remove,header,nsolute,
!$OMP& isolute,icount,iout)
!$OMP DO schedule(guided)
c
c     search for short distance between solute and solvent
c
      do i = 1, n
         if (.not. remove(molcule(i))) then
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = 1, nsolute
               j = isolute(k)
               xr = x(j) - xi
               yr = y(j) - yi
               zr = z(j) - zi
               call imagen (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
               if (rik2 .lt. close2) then
                  remove(molcule(i)) = .true.
                  goto 100
               end if
            end do
  100       continue
         end if
         icount = icount + 1
         if (mod(icount,10000) .eq. 0) then
            if (header) then
               header = .false.
               write (iout,110)
  110          format ()
            end if
            write (iout,120)  10000*(icount/10000)
  120       format (' Solvent Atoms Processed',i15)
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (isolute)
c
c     print final status when processing large systems
c
      if (mod(n,10000).ne.0 .and. n.gt.10000) then
         write (iout,130)  n
  130    format (' Solvent Atoms Processed',i15)
      end if
c
c     randomly replace the solvent molecules with ions
c
      do i = 1, ncopy
         done = .false.
         do while (.not. done)
            rand = random ()
            ranatm = int(rand*dble(n)) + 1
c
c     check solute distance, then delete polyatomic molecule
c
            if (.not. remove(molcule(ranatm))) then
               start = imol(1,molcule(ranatm))
               stop = imol(2,molcule(ranatm))
               if (start .eq. stop) then
                  done = .false.
               else
                  done = .true.
                  xmid = 0.0d0
                  ymid = 0.0d0
                  zmid = 0.0d0
                  do k = stop, start, -1
                     weigh = mass(k)
                     xmid = xmid + x(k)*weigh
                     ymid = ymid + y(k)*weigh
                     zmid = zmid + z(k)*weigh
                     call delete (k)
                  end do
                  weigh = molmass(molcule(ranatm))
                  xion(i) = xmid / weigh
                  yion(i) = ymid / weigh
                  zion(i) = zmid / weigh
               end if
            end if
         end do
      end do
c
c     insert new monoatomic ions at saved centers of mass
c
      do i = 1, ncopy
         n = n + 1
         name(n) = symbol(iontyp)
         x(n) = xion(i)
         y(n) = yion(i)
         z(n) = zion(i)
         type(n) = iontyp
         n12(n) = 0
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (remove)
      deallocate (xion)
      deallocate (yion)
      deallocate (zion)
      return
      end
