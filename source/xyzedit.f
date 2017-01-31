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
      use sizes
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
      integer i,j,k,m,it
      integer ixyz,imod
      integer init,stop
      integer nmode,mode
      integer natom,atmnum
      integer nlist,ncopy
      integer offset,origin
      integer oldtype,newtype
      integer freeunit
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
      real*8, allocatable :: x0(:)
      real*8, allocatable :: y0(:)
      real*8, allocatable :: z0(:)
      real*8 a(3,3)
      logical exist,query
      logical opened,multi
      logical append
      character*240 xyzfile
      character*240 modfile
      character*240 record
      character*240 string
      external merge
c
c
c     initialize various constants and the output flags
c
      call initial
      opened = .false.
      multi = .false.
      nmode = 20
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
      call field
      call katom
c
c     present a list of possible coordinate modifications
c
      write (iout,30)
   30 format (/,' The TINKER XYZ File Editing Utility Can :',
     &        //,4x,'(1) Offset the Numbers of the Current Atoms',
     &        /,4x,'(2) Deletion of Individual Specified Atoms',
     &        /,4x,'(3) Deletion of Specified Types of Atoms',
     &        /,4x,'(4) Deletion of Atoms outside Cutoff Range',
     &        /,4x,'(5) Insertion of Individual Specified Atoms',
     &        /,4x,'(6) Replace Old Atom Type with a New Type',
     &        /,4x,'(7) Assign Connectivities for Linear Chain',
     &        /,4x,'(8) Assign Connectivities based on Distance',
     &        /,4x,'(9) Convert Units from Bohrs to Angstroms',
     &        /,3x,'(10) Invert thru Origin to give Mirror Image',
     &        /,3x,'(11) Translate All Atoms by an X,Y,Z-Vector',
     &        /,3x,'(12) Translate Center of Mass to the Origin',
     &        /,3x,'(13) Translate a Specified Atom to the Origin',
     &        /,3x,'(14) Translate and Rotate to Inertial Frame',
     &        /,3x,'(15) Move to Specified Rigid Body Coordinates',
     &        /,3x,'(16) Move Stray Molecules into Periodic Box',
     &        /,3x,'(17) Delete Molecules outside of Periodic Box',
     &        /,3x,'(18) Append a Second XYZ File to Current One',
     &        /,3x,'(19) Create and Fill a Periodic Boundary Box',
     &        /,3x,'(20) Soak Current Molecule in Box of Solvent')
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
   60       format (/,' Number of the Desired Choice [<CR>=Exit] :  ',$)
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
         write (iout,100)
  100    format (/,' Offset used to Renumber the Atoms [0] :  ',$)
         read (input,110,err=90)  offset
  110    format (i10)
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
         write (iout,120)
  120    format (/,' Numbers of the Atoms to be Removed :  ',$)
         read (input,130)  record
  130    format (a240)
         read (record,*,err=140,end=140)  (list(i),i=1,n)
  140    continue
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
         write (iout,150)
  150    format (/,' Atom Types to be Removed :  ',$)
         read (input,160)  record
  160    format (a240)
         read (record,*,err=170,end=170)  (list(i),i=1,n)
  170    continue
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
                     goto 180
                  end if
               end do
  180          continue
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
                        if (keep(j) .eq. i)  goto 190
                        dist2 = (x(j)-xi)**2+(y(j)-yi)**2+(z(j)-zi)**2
                        if (dist2 .le. cut2)  goto 190
                     end if
                  end do
                  nlist = nlist + 1
                  list(nlist) = i
  190             continue
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
         write (iout,200)
  200    format (/,' Numbers of the Atoms to be Inserted :  ',$)
         read (input,210)  record
  210    format (a240)
         read (record,*,err=220,end=220)  (list(i),i=1,n)
  220    continue
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
  230    continue
         oldtype = 0
         newtype = 0
         write (iout,240)
  240    format (/,' Numbers of the Old and New Atom Types :  ',$)
         read (input,250)  record
  250    format (a240)
         read (record,*,err=230,end=230)  oldtype,newtype
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
         write (iout,260)
  260    format (/,' Enter Translation Vector Components :  ',$)
         read (input,270)  record
  270    format (a240)
         read (record,*,err=280,end=280)  xr,yr,zr
  280    continue
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
         write (iout,290)
  290    format (/,' Number of the Atom to Move to the Origin :  ',$)
         read (input,300)  origin
  300    format (i10)
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
         write (iout,310)
  310    format (/,' Enter Rigid Body Coordinates :  ',$)
         read (input,320)  record
  320    format (a240)
         read (record,*,err=330,end=330)  xcm,ycm,zcm,phi,theta,psi
  330    continue
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
         do while (xnew .eq. 0.0d0)
            write (iout,340)
  340       format (/,' Enter Periodic Box Dimensions (X,Y,Z) :  ',$)
            read (input,350)  record
  350       format (a240)
            read (record,*,err=360,end=360)  xnew,ynew,znew
  360       continue
            if (ynew .eq. 0.0d0)  ynew = xnew
            if (znew .eq. 0.0d0)  znew = xnew
         end do
         allocate (list(n))
         allocate (keep(n))
         do while (.not. abort)
            xbox = xnew
            ybox = ynew
            zbox = znew
            call lattice
            call molecule
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
      if (mode .eq. 18) then
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
      if (mode .eq. 19) then
         write (iout,370)
  370    format (/,' Enter Number of Copies to Put in Box :  ',$)
         read (input,380)  ncopy
  380    format (i10)
         xbox = 0.0d0
         ybox = 0.0d0
         zbox = 0.0d0
         do while (xbox .eq. 0.0d0)
            write (iout,390)
  390       format (/,' Enter Periodic Box Dimensions (X,Y,Z) :  ',$)
            read (input,400)  record
  400       format (a240)
            read (record,*,err=410,end=410)  xbox,ybox,zbox
  410       continue
            if (ybox .eq. 0.0d0)  ybox = xbox
            if (zbox .eq. 0.0d0)  zbox = xbox
         end do
         orthogonal = .true.
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
         do i = 1, n
            x(i) = x(i) - xcm
            y(i) = y(i) - ycm
            z(i) = z(i) - zcm
            x0(i) = x(i)
            y0(i) = y(i)
            z0(i) = z(i)
         end do
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
         call lattice
         call molecule
         call bounds
      end if
c
c     solvate the current system by insertion into a solvent box
c
      if (mode .eq. 20) then
         call soak
      end if
c
c     output final coordinates for single frame input file
c
      if (opened .and. .not.multi) then
         call prtmod (imod,offset)
      end if
c
c     perform any final tasks before program exit
c
      if (opened) then
         close (unit=imod)
         write (iout,420)  modfile
  420    format (/,' New Coordinates File Written To :  ',a)
      end if
      close (unit=ixyz)
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
      use sizes
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
      use sizes
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
      integer isolv,icount
      integer ntot,freeunit
      integer, allocatable :: map(:)
      real*8 xi,yi,zi
      real*8 xr,yr,zr,rik2
      real*8 close,close2
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
c     set distance cutoff for solute-solvent close contacts
c
      close = 1.5d0
      close2 = close * close
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
   30    format (/,' Scan for Solvent Molecules to be Removed')
      end if
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(nref,n,x,y,z,molcule,close2,remove,header,icount)
!$OMP DO schedule(guided)
c
c     search for close contacts between solute and solvent
c
      do i = nref(1)+1, n
         if (.not. remove(molcule(i))) then
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do k = 1, nref(1)
               xr = x(k) - xi
               yr = y(k) - yi
               zr = z(k) - zi
               call imagen (xr,yr,zr)
               rik2 = xr*xr + yr*yr + zr*zr
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
