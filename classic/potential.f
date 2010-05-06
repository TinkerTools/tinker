c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2008 by Chuanjie Wu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program potential  --  compute electrostatic potential  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potential" calculates the electrostatic potential for a
c     molecule at a set of grid points; optionally compares to a
c     target potential or optimizes electrostatic parameters
c
c
      program potential
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'minima.i'
      include 'output.i'
      include 'potent.i'
      include 'potfit.i'
      include 'charge.i'
      include 'mpole.i'
      include 'titles.i'
      include 'units.i'
      integer i,j,k
      integer ixyz,ipot,igrd
      integer next,mode
      integer nmodel,nvar
      integer nglist,nflist
      integer freeunit
      integer trimtext
      integer glist(maxatm)
      integer flist(maxatm)
      real*8 xi,yi,zi,pot
      real*8 x0,y0,z0
      real*8 xx0,xy0,xz0
      real*8 yy0,yz0,zz0
      real*8 minimum,grdmin
      real*8 potfit1
      real*8 xx(maxvar)
      logical exist,query
      logical dogrid,docube
      logical domodel,dopair
      logical dotarget,dofit
      logical dofull
      character*1 answer
      character*20 keyword
      character*120 record
      character*120 string
      character*120 xyzfile
      character*120 potfile
      character*120 gridfile
      external potfit1
      external optsave
c
c
c     setup the computation and assign some default values
c
      call initial
      nmodel = 1
      nconf = 0
      dogrid = .false.
      docube = .false.
      domodel = .false.
      dopair = .false.
      dotarget = .false.
      dofit = .false.
      fit_mpl = .true.
      fit_dpl = .true.
      fit_qdp = .true.
c
c     initialize target molecular dipole and quadrupole values
c
      use_dpl = .false.
      use_qdp = .false.
      do i = 1, maxref
         xdpl0(i) = 0.0d0
         ydpl0(i) = 0.0d0
         zdpl0(i) = 0.0d0
         xxqdp0(i) = 0.0d0
         xyqdp0(i) = 0.0d0
         xzqdp0(i) = 0.0d0
         yyqdp0(i) = 0.0d0
         yzqdp0(i) = 0.0d0
         zzqdp0(i) = 0.0d0
      end do
c
c     find electrostatic potential manipulation to be performed
c
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
   20    format (/,' The TINKER Electrostatic Potential Facility Can :',
     &           //,4x,'(1) Create an Input File for Gaussian CUBEGEN',
     &           /,4x,'(2) Get QM Potential from a Gaussian CUBE File',
     &           /,4x,'(3) Calculate the Model Potential for a System',
     &           /,4x,'(4) Compare the Model Potentials of Two Systems',
     &           /,4x,'(5) Compare a Model Potential to a Target Grid',
     &           /,4x,'(6) Fit Electrostatic Parameters to Target Grid')
         dowhile (mode.lt.1 .or. mode.gt.6)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
      if (mode .eq. 1) then
         dogrid = .true.
      else if (mode .eq. 2) then
         docube = .true.
      else if (mode .eq. 3) then
         domodel = .true.
      else if (mode .eq. 4) then
         nmodel = 2
         dopair = .true.
      else if (mode .eq. 5) then
         dotarget = .true.
      else if (mode .eq. 6) then
         dotarget = .true.
         dofit = .true.
      end if
c
c     read the electrostatic potential from a Gaussian CUBE file
c
      if (docube) then
         call nextarg (potfile,exist)
         if (exist) then
            call basefile (potfile)
            call suffix (potfile,'pot')
            call version (potfile,'old')
            inquire (file=potfile,exist=exist)
         end if
         dowhile (.not. exist)
            write (iout,60)
   60       format (/,' Enter the Gaussian CUBE File Name :  ',$)
            read (input,70)  potfile
   70       format (a120)
            call basefile (potfile)
            call suffix (potfile,'cube')
            call version (potfile,'old')
            inquire (file=potfile,exist=exist)
         end do
         ipot = freeunit ()
         open (unit=ipot,file=potfile,status ='old')
         rewind (unit=ipot)
         read (ipot,80)  title
   80    format (1x,a120)
         ltitle = trimtext (title)
         read (ipot,90)
   90    format ()
         read (ipot,100)  n
  100    format (i5)
         read (ipot,110)  npgrid(1)
  110    format (i5)
         do i = 1, n+2
            read (ipot,120)
  120       format ()
         end do
         do i = 1, npgrid(1)
            read (ipot,130)  record
  130       format (a120)
            read (record,*)  xi,yi,zi,pot
            pgrid(1,i,1) = xi
            pgrid(2,i,1) = yi
            pgrid(3,i,1) = zi
            epot(1,i,1) = hartree * pot
         end do
         close (unit=ipot)
         potfile = filename
         call suffix (potfile,'pot')
         call version (potfile,'new')
         open (unit=ipot,file=potfile,status ='new')
         rewind (unit=ipot)
         write (ipot,140)  npgrid(1),title(1:ltitle)
  140    format (i8,2x,a)
         do i = 1, npgrid(1)
            xi = pgrid(1,i,1)
            yi = pgrid(2,i,1)
            zi = pgrid(3,i,1)
            pot = epot(1,i,1)
            write (ipot,150)  i,xi,yi,zi,pot
  150       format (i8,3x,3f12.6,2x,f12.4)
         end do
         close (unit=ipot)
         write (iout,160)  potfile(1:trimtext(potfile))
  160    format (/,' Electrostatic Potential Written to File :  ',a)
         goto 380
      end if
c
c     read the first structure and get connectivity info
c
      call getxyz
      call attach
      call active
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     reopen the structure file and get all conformations
c
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz')
      call version (xyzfile,'old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
      dowhile (.not. abort)
         nconf = nconf + 1
         call makeref (nconf)
         call readxyz (ixyz)
      end do
      close (unit=ixyz)
      call getref (1)
      if (nconf .gt. 1) then
         write (iout,170)  nconf
  170    format (/,' Conformations for Potential Analysis :',i6)
      end if
c
c     get the force field parameters and setup electrostatics
c
      call cutoffs
      call field
      call katom
      if (use_charge .or. use_chgdpl)  call kcharge
      if (use_dipole .or. use_chgdpl)  call kdipole
      if (use_mpole .or. use_polar)  call kmpole
      if (use_polar)  call kpolar
c
c     set defaults for the active grid atoms and fit atoms
c
      nglist = 0
      nflist = 0
      ngatm = n
      nfatm = n
      do i = 1, n
         glist(i) = 0
         flist(i) = 0
         gatm(i) = .true.
         fatm(i) = .true.
      end do
c
c     get control parameters and target values from keyfile
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:16) .eq. 'POTENTIAL-ATOMS ') then
            read (string,*,err=180,end=180)  (glist(k),k=nglist+1,n)
  180       continue
            dowhile (glist(nglist+1) .ne. 0)
               nglist = nglist + 1
               glist(nglist) = max(-n,min(n,glist(nglist)))
            end do
         else if (keyword(1:14) .eq. 'POTENTIAL-FIT ') then
            read (string,*,err=190,end=190)  (flist(k),k=nflist+1,n)
  190       continue
            dowhile (flist(nflist+1) .ne. 0)
               nflist = nflist + 1
               flist(nflist) = max(-n,min(n,flist(nflist)))
            end do
         else if (keyword(1:13) .eq. 'FIX-MONOPOLE ') then
            fit_mpl = .false.
         else if (keyword(1:11) .eq. 'FIX-DIPOLE ') then
            fit_dpl = .false.
         else if (keyword(1:15) .eq. 'FIX-QUADRUPOLE ') then
            fit_qdp = .false.
         else if (keyword(1:14) .eq. 'TARGET-DIPOLE ') then
            use_dpl = .true.
            k = 1
            read (string,*,err=200,end=200)  x0,y0,z0,k
  200       continue
            xdpl0(k) = x0
            ydpl0(k) = y0
            zdpl0(k) = z0
         else if (keyword(1:18) .eq. 'TARGET-QUADRUPOLE ') then
            use_qdp = .true.
            k = 1
            read (string,*,err=210,end=210)  xx0,xy0,xz0,yy0,yz0,zz0,k
  210       continue
            xxqdp0(k) = xx0
            xyqdp0(k) = xy0
            xzqdp0(k) = xz0
            yyqdp0(k) = yy0
            yzqdp0(k) = yz0
            zzqdp0(k) = zz0
         end if
      end do
c
c     set active grid atoms to only those marked for use
c
      i = 1
      dowhile (glist(i) .ne. 0)
         if (i .eq. 1) then
            ngatm = 0
            do k = 1, n
               gatm(k) = .false.
            end do
         end if
         if (glist(i) .gt. 0) then
            k = glist(i)
            if (.not. gatm(k)) then
               gatm(k) = .true.
               ngatm = ngatm + 1
            end if
            i = i + 1
         else
            do k = abs(glist(i)), abs(glist(i+1))
               if (.not. gatm(k)) then
                  gatm(k) = .true.
                  ngatm = ngatm + 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     set active fitting atoms to only those marked for use
c
      i = 1
      dowhile (flist(i) .ne. 0)
         if (i .eq. 1) then
            nfatm = 0
            do k = 1, n
               fatm(k) = .false.
            end do
         end if
         if (flist(i) .gt. 0) then
            k = flist(i)
            if (.not. fatm(k)) then
               fatm(k) = .true.
               nfatm = nfatm + 1
            end if
            i = i + 1
         else
            do k = abs(flist(i)), abs(flist(i+1))
               if (.not. fatm(k)) then
                  fatm(k) = .true.
                  nfatm = nfatm + 1
               end if
            end do
            i = i + 2
         end if
      end do
c
c     generate potential grid based on the molecular surface
c
      if (.not. dotarget) then
         do i = 1, nconf
            call getref (i)
            call potgrid (i)
         end do
      end if
c
c     get name of optional second structure for comparison
c
      if (dopair) then
         call nextarg (xyzfile,exist)
         if (exist) then
            call basefile (xyzfile)
            call suffix (xyzfile,'xyz')
            call version (xyzfile,'old')
            inquire (file=xyzfile,exist=exist)
         end if
         dowhile (.not. exist)
            write (iout,220)
  220       format (/,' Enter Name of Second Coordinate File :  ',$)
            read (input,230)  xyzfile
  230       format (a120)
            call basefile (xyzfile)
            call suffix (xyzfile,'xyz')
            call version (xyzfile,'old')
            inquire (file=xyzfile,exist=exist)
         end do
      end if
c
c     get optional file with grid points and target potential
c
      if (dotarget) then
         call nextarg (potfile,exist)
         if (exist) then
            call basefile (potfile)
            call suffix (potfile,'pot')
            call version (potfile,'old')
            inquire (file=potfile,exist=exist)
         end if
         dowhile (.not. exist)
            write (iout,240)
  240       format (/,' Enter Target Grid/Potential File Name :  ',$)
            read (input,250)  potfile
  250       format (a120)
            call basefile (potfile)
            call suffix (potfile,'pot')
            call version (potfile,'old')
            inquire (file=potfile,exist=exist)
         end do
      end if
c
c     decide whether to output potential at each grid point
c
      dofull = .false.
      if (domodel .or. dopair .or. dotarget) then
         call nextarg (answer,exist)
         if (.not. exist) then
            write (iout,260)
  260       format (/,' Output Potential Value at Each Grid Point',
     &                 ' [N] :  ',$)
            read (input,270)  record
  270       format (a120)
            next = 1
            call gettext (record,answer,next)
         end if
         call upcase (answer)
         if (answer .eq. 'Y')  dofull = .true.
      end if
c
c     read the grid points where potential will be computed
c
      if (dotarget) then
         ipot = freeunit ()
         open (unit=ipot,file=potfile,status='old')
         rewind (unit=ipot)
         do i = 1, nconf
            call readpot (ipot,i)
         end do
         close (unit=ipot)
      end if
c
c     output the number of potential grid points to be used
c
      do i = 1, nconf
         if (i .eq. 1) then
            write (iout,280)
  280       format ()
         end if
         if (npgrid(i) .gt. maxpgrd) then
            write (iout,290)
  290       format (' POTENTIAL  --  Too many Grid Points;',
     &                 ' Increase MAXGRID')
            call fatal
         else if (nconf .eq. 1) then
            write (iout,300)  npgrid(1)
  300       format (' Electrostatic Potential Grid Points :',6x,i10)
         else
            write (iout,310)  i,npgrid(i)
  310       format (' Potential Grid Points for Conformer',i4,' :',
     &                 2x,i10)
         end if
      end do
c
c     output grid points at which to compute QM potential
c
      if (dogrid) then
         igrd = freeunit ()
         gridfile = filename
         call suffix (gridfile,'grid')
         call version (gridfile,'new')
         open (unit=igrd,file=gridfile,status='new')
         do j = 1, nconf
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               write (igrd,320)  xi,yi,zi
  320          format (3f15.8)
            end do
         end do
         close (unit=igrd)
         write (iout,330)  gridfile(1:trimtext(gridfile))
  330    format (/,' Gaussian CUBEGEN Input Written to File :   ',a)
         write (iout,340)
  340    format (/,' Next, run the Gaussian CUBEGEN program; for',
     &              ' example:',
     &           /,5x,'cubegen 0 potential=MP2 xxx.fchk',
     &              ' xxx.cube -5 h < xxx.grid',
     &           //,' See the Gaussian documentation for additional',
     &              ' details;',
     &           /,' After CUBEGEN, rerun TINKER POTENTIAL using',
     &              ' Option 2')
      end if
c
c     get termination criterion for fitting as RMS gradient
c
      if (dofit) then
         grdmin = -1.0d0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=350,end=350)  grdmin
  350    continue
         if (grdmin .le. 0.0d0) then
            write (iout,360)
  360       format (/,' Enter RMS Gradient Termination Criterion',
     &                 ' [0.001] :  ',$)
            read (input,370)  grdmin
  370       format (f20.0)
         end if
         if (grdmin .le. 0.0d0)  grdmin = 0.001d0
      end if
c
c     setup the potential computation for each input structure
c
      if (.not. dogrid) then
         do k = 1, nmodel
            if (k .eq. 2) then
               ixyz = freeunit ()
               open (unit=ixyz,file=xyzfile,status='old')
               rewind (unit=ixyz)
               do j = 1, nconf
                  call readxyz (ixyz)
                  call makeref (j)
               end do
               close (unit=ixyz)
               call cutoffs
               call field
               call katom
               if (use_charge .or. use_chgdpl)  call kcharge
               if (use_dipole .or. use_chgdpl)  call kdipole
               if (use_mpole .or. use_polar)  call kmpole
               if (use_polar)  call kpolar
            end if
c
c     compute the value of the potential at each grid point
c
            do j = 1, nconf
               call getref (j)
               if (use_mpole) then
                  call chkpole
                  call rotpole
               end if
               if (use_polar) then
                  call induce
               end if
               do i = 1, npgrid(j)
                  xi = pgrid(1,i,j)
                  yi = pgrid(2,i,j)
                  zi = pgrid(3,i,j)
                  call potpoint (xi,yi,zi,pot)
                  epot(k,i,j) = pot
               end do
            end do
         end do
c
c     output overall statistics for electrostatic potential
c
         call potstat (dofull,domodel,dopair,dotarget)
      end if
c
c     set initial electrostatic parameters and do optimization
c
      if (dofit) then
         call prmvar (nvar,xx)
         hguess = 1.0d-4
         coordtype = 'NONE'
         call ocvm (nvar,xx,minimum,grdmin,potfit1,optsave)
c
c     move optimized values back to parameters and output
c
         call varprm (nvar,xx,0,0.0d0)
         call prmvar (nvar,xx)
c
c     evaluate the potential value at each grid point
c
         do j = 1, nconf
            call getref (j)
            if (use_mpole) then
               call chkpole
               call rotpole
            end if
            if (use_polar) then
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
            end do
         end do
c
c     output overall statistics and optimized parameters
c
         call potstat (dofull,domodel,dopair,dotarget)
         call prtfit
      end if
  380 continue
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine readpot  --  get and assign potential grid  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "readpot" gets a set of grid points and target electrostatic
c     potential values from an external disk file
c
c
      subroutine readpot (ipot,iconf)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'potfit.i'
      integer i,j,k
      integer ipot,iconf
      integer npoint,atmnum
      real*8 xi,yi,zi
      real*8 big,small
      real*8 r2,dist
      real*8 rad(maxatm)
      character*120 record
c
c
c     read the grid points and target potential from a file
c
      npoint = 0
      read (ipot,10,err=20,end=20)  record
   10 format (a120)
      read (record,*,err=20,end=20)  npoint
   20 continue
      do i = 1, npoint
         pgrid(1,i,iconf) = 0.0d0
         pgrid(2,i,iconf) = 0.0d0
         pgrid(3,i,iconf) = 0.0d0
         epot(2,i,iconf) = 0.0d0
         read (ipot,30,err=40,end=40)  record
   30    format (a120)
         read (record,*,err=40,end=40)  k,(pgrid(j,i,iconf),j=1,3),
     &                                  epot(2,i,iconf)
   40    continue
      end do
c
c     set base atomic radii from traditional Bondi values
c
      do i = 1, n
         rad(i) = 1.70d0
         atmnum = atomic(i)
         if (atmnum .eq. 0)  rad(i) = 0.00d0
         if (atmnum .eq. 1)  rad(i) = 1.20d0
         if (atmnum .eq. 2)  rad(i) = 1.40d0
         if (atmnum .eq. 6)  rad(i) = 1.70d0
         if (atmnum .eq. 7)  rad(i) = 1.55d0
         if (atmnum .eq. 8)  rad(i) = 1.52d0
         if (atmnum .eq. 9)  rad(i) = 1.47d0
         if (atmnum .eq. 10)  rad(i) = 1.54d0
         if (atmnum .eq. 14)  rad(i) = 2.10d0
         if (atmnum .eq. 15)  rad(i) = 1.80d0
         if (atmnum .eq. 16)  rad(i) = 1.80d0
         if (atmnum .eq. 17)  rad(i) = 1.75d0
         if (atmnum .eq. 18)  rad(i) = 1.88d0
         if (atmnum .eq. 35)  rad(i) = 1.85d0
         if (atmnum .eq. 36)  rad(i) = 2.02d0
         if (atmnum .eq. 53)  rad(i) = 1.98d0
         if (atmnum .eq. 54)  rad(i) = 2.16d0
      end do
c
c     assign each grid point to atom on molecular surface
c
      big = 1000.0d0
      do i = 1, npoint
         small = big
         xi = pgrid(1,i,iconf)
         yi = pgrid(2,i,iconf)
         zi = pgrid(3,i,iconf)
         do k = 1, n
            r2 = (xi-x(k))**2 + (yi-y(k))**2 + (zi-z(k))**2
            dist = sqrt(r2) - rad(k)
            if (dist .lt. small) then
               small = dist
               ipgrid(i,iconf) = k
            end if
         end do
      end do
c
c     use potential grid points only for active grid atoms
c
      k = npoint
      npoint = 0
      do i = 1, k
         if (gatm(ipgrid(i,iconf))) then
            npoint = npoint + 1
            ipgrid(npoint,iconf) = ipgrid(i,iconf)
            pgrid(1,npoint,iconf) = pgrid(1,i,iconf)
            pgrid(2,npoint,iconf) = pgrid(2,i,iconf)
            pgrid(3,npoint,iconf) = pgrid(3,i,iconf)
            epot(2,npoint,iconf) = epot(2,i,iconf)
         end if
      end do
      npgrid(iconf) = npoint
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine potgrid  --  generate shells of grid points  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potgrid" generates electrostatic potential grid points in
c     radially distributed shells outside the molecular surface
c
c
      subroutine potgrid (iconf)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'keys.i'
      include 'math.i'
      include 'potfit.i'
      integer maxdot
      parameter (maxdot=50000)
      integer i,j,k,m
      integer iconf,next
      integer npoint,nshell
      integer ndot,atmnum
      real*8 r2,roffset
      real*8 spacing
      real*8 density
      real*8 round
      real*8 xi,yi,zi
      real*8 xj,yj,zj
      real*8 xr,yr,zr
      real*8 rad(maxatm)
      real*8 rad2(maxatm)
      real*8 dot(3,maxdot)
      character*20 keyword
      character*120 record
      character*120 string
c
c
c     set default values for grid point generation parameters
c
      npoint = 0
      nshell = 4
      spacing = 0.35d0
      density = 4.0d0 * pi / spacing**2
      roffset = 1.0d0
      round = 0.000001d0
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:17) .eq. 'POTENTIAL-SHELLS ') then
            read (string,*,err=10,end=10)  nshell
         else if (keyword(1:18) .eq. 'POTENTIAL-SPACING ') then
            read (string,*,err=10,end=10)  spacing
            density = 4.0d0 * pi / spacing**2
         else if (keyword(1:17) .eq. 'POTENTIAL-OFFSET ') then
            read (string,*,err=10,end=10)  roffset
         end if
   10    continue
      end do
c
c     set base atomic radii from traditional Bondi values
c
      do i = 1, n
         rad(i) = 1.70d0
         atmnum = atomic(i)
         if (atmnum .eq. 0)  rad(i) = 0.00d0
         if (atmnum .eq. 1)  rad(i) = 1.20d0
         if (atmnum .eq. 2)  rad(i) = 1.40d0
         if (atmnum .eq. 6)  rad(i) = 1.70d0
         if (atmnum .eq. 7)  rad(i) = 1.55d0
         if (atmnum .eq. 8)  rad(i) = 1.52d0
         if (atmnum .eq. 9)  rad(i) = 1.47d0
         if (atmnum .eq. 10)  rad(i) = 1.54d0
         if (atmnum .eq. 14)  rad(i) = 2.10d0
         if (atmnum .eq. 15)  rad(i) = 1.80d0
         if (atmnum .eq. 16)  rad(i) = 1.80d0
         if (atmnum .eq. 17)  rad(i) = 1.75d0
         if (atmnum .eq. 18)  rad(i) = 1.88d0
         if (atmnum .eq. 35)  rad(i) = 1.85d0
         if (atmnum .eq. 36)  rad(i) = 2.02d0
         if (atmnum .eq. 53)  rad(i) = 1.98d0
         if (atmnum .eq. 54)  rad(i) = 2.16d0
         rad(i) = rad(i) + roffset
         rad2(i) = rad(i)**2
      end do
c
c     find points on each of the molecular surface shells
c
      do m = 1, nshell
         if (m .ne. 1) then
            do i = 1, n
               rad(i) = rad(i) + spacing
               rad2(i) = rad(i)**2
            end do
         end if
         do i = 1, n
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ndot = int(density*rad2(i))
            if (ndot .gt. maxdot) then
               write (iout,20)
   20          format (/,' POTGRID  --  Too many Surface Grid Points;',
     &                    ' Increase MAXDOT')
               call fatal
            end if
            call sphere (ndot,dot)
            do j = 1, ndot
               xj = xi + rad(i)*dot(1,j)
               yj = yi + rad(i)*dot(2,j)
               zj = zi + rad(i)*dot(3,j)
               xj = dble(nint(xj/round)) * round
               yj = dble(nint(yj/round)) * round
               zj = dble(nint(zj/round)) * round
               do k = 1, i-1
                  xr = xj - x(k)
                  yr = yj - y(k)
                  zr = zj - z(k)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .lt. rad2(k))  goto 30
               end do
               do k = i+1, n
                  xr = xj - x(k)
                  yr = yj - y(k)
                  zr = zj - z(k)
                  r2 = xr*xr + yr*yr + zr*zr
                  if (r2 .lt. rad2(k))  goto 30
               end do
               npoint = npoint + 1
               ipgrid(npoint,iconf) = i
               pgrid(1,npoint,iconf) = xj
               pgrid(2,npoint,iconf) = yj
               pgrid(3,npoint,iconf) = zj
   30          continue
            end do
         end do
      end do
c
c     use potential grid points only for active grid atoms
c
      k = npoint
      npoint = 0
      do i = 1, k
         if (gatm(ipgrid(i,iconf))) then
            npoint = npoint + 1
            ipgrid(npoint,iconf) = ipgrid(i,iconf)
            pgrid(1,npoint,iconf) = pgrid(1,i,iconf)
            pgrid(2,npoint,iconf) = pgrid(2,i,iconf)
            pgrid(3,npoint,iconf) = pgrid(3,i,iconf)
         end if
      end do
      npgrid(iconf) = npoint
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine potpoint  --  electrostatic potential at point  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "potpoint" calculates the electrostatic potential at a grid
c     point "i" as the total electrostatic interaction energy of
c     a molecule with a positive charge located at the grid point
c
c
      subroutine potpoint (xi,yi,zi,pot)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'chgpot.i'
      include 'dipole.i'
      include 'mpole.i'
      include 'polar.i'
      include 'potent.i'
      include 'units.i'
      integer k,kk,k1,k2
      real*8 e,ei,pot
      real*8 ec,ed,em,ep
      real*8 xi,yi,zi
      real*8 xk,yk,zk
      real*8 xr,yr,zr
      real*8 r,r2
      real*8 rk2,rkr3,dotk
      real*8 rr1,rr3,rr5
      real*8 f,fi,ci,ck
      real*8 dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qkx,qky,qkz
      real*8 scd,scq,scu
c
c
c     zero out charge, dipole and polarizable multipole energies
c
      ec = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
c
c     set charge of probe site and electrostatic constants
c
      f = electric / dielec
      ci = 1.0d0
      fi = f * ci
c
c     calculate the charge contribution to the potential
c
      if (use_charge) then
         do k = 1, nion
            kk = iion(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            r2 = xr*xr + yr* yr + zr*zr
            r = sqrt(r2)
            e = fi * pchg(k) / r
            ec = ec + e
         end do
      end if
c
c     calculate the bond dipole contribution to the potential
c
      if (use_dipole) then
         do k = 1, ndipole
            k1 = idpl(1,k)
            k2 = idpl(2,k)
            xk = x(k2) - x(k1)
            yk = y(k2) - y(k1)
            zk = z(k2) - z(k1)
            xr = x(k1) + xk*sdpl(k) - xi
            yr = y(k1) + yk*sdpl(k) - yi
            zr = z(k1) + zk*sdpl(k) - zi
            r2 = xr*xr + yr* yr + zr*zr
            rk2 = xk*xk + yk*yk + zk*zk
            rkr3 = sqrt(rk2*r2) * r2
            dotk = xk*xr + yk*yr + zk*zr
            e = (fi/debye) * bdpl(k) * dotk / rkr3
            ed = ed + e
         end do
      end if
c
c     calculate the multipole contribution to the potential
c
      if (use_mpole .or. use_polar) then
         do k = 1, npole
            kk = ipole(k)
            xr = x(kk) - xi
            yr = y(kk) - yi
            zr = z(kk) - zi
            r2 = xr*xr + yr* yr + zr*zr
            r = sqrt(r2)
            ck = rpole(1,k)
            dkx = rpole(2,k)
            dky = rpole(3,k)
            dkz = rpole(4,k)
            qkxx = rpole(5,k)
            qkxy = rpole(6,k)
            qkxz = rpole(7,k)
            qkyy = rpole(9,k)
            qkyz = rpole(10,k)
            qkzz = rpole(13,k)
            ukx = uind(1,k)
            uky = uind(2,k)
            ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
            qkx = qkxx*xr + qkxy*yr + qkxz*zr
            qky = qkxy*xr + qkyy*yr + qkyz*zr
            qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate scalar products for permanent and induced
c
            scd = dkx*xr + dky*yr + dkz*zr
            scq = qkx*xr + qky*yr + qkz*zr
            scu = ukx*xr + uky*yr + ukz*zr
c
c     compute the energy contributions for this interaction
c
            rr1 = 1.0d0 / r
            rr3 = rr1 / r2
            rr5 = 3.0d0 * rr3 / r2
            e = ck*rr1 - scd*rr3 + scq*rr5
            ei = -scu * rr3
c
c     increment the overall multipole and polarization energies
c
            e = fi * e
            ei = fi * ei
            em = em + e
            ep = ep + ei
         end do
      end if
c
c     potential is sum of all interactions with probe site
c
      pot = ec + ed + em + ep
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function potfit1  --  potential fit error and gradient  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "potfit1" is a service routine that computes the RMS error
c     and gradient for electrostatic parameters fit to a potential
c
c
      function potfit1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'moment.i'
      include 'potent.i'
      include 'potfit.i'
      integer i,j,k,m
      integer nvar,npoint
      real*8 potfit1
      real*8 pot,eps
      real*8 e,e0
      real*8 er,ec,et
      real*8 xi,yi,zi
      real*8 cscale,tscale
      real*8 xx(maxvar)
      real*8 g(maxvar)
c
c
c     initialize scaling factors for error and gradient
c
      npoint = 0
      do j = 1, nconf
         npoint = npoint + npgrid(j)
      end do
      cscale = 100000000.0d0 / dble(nconf)
      tscale = 10000.0d0 / dble(nconf)
      eps = 0.000001d0
c
c     copy optimization values to electrostatic parameters
c
      call varprm (nvar,xx,0,0.0d0)
c
c     find total error by cycling over all conformations
c
      er = 0.0d0
      ec = 0.0d0
      et = 0.0d0
      do j = 1, nconf
         call getref (j)
         if (use_mpole) then
            call chkpole
            call rotpole
         end if
         if (use_polar) then
            call induce
         end if
c
c     get the RMS potential error summed over grid points
c
         do i = 1, npgrid(j)
            xi = pgrid(1,i,j)
            yi = pgrid(2,i,j)
            zi = pgrid(3,i,j)
            call potpoint (xi,yi,zi,pot)
            epot(1,i,j) = pot
            er = er + (epot(1,i,j)-epot(2,i,j))**2
         end do
c
c     get deviation from integral net molecular charge
c
         call momfull
         ec = ec + cscale*(netchg-dble(nint(netchg)))**2
c
c     get deviation from dipole and quadrupole targets
c
         if (use_dpl) then
            et = et + tscale*(xdpl-xdpl0(j))**2
            et = et + tscale*(ydpl-ydpl0(j))**2
            et = et + tscale*(zdpl-zdpl0(j))**2
         end if
         if (use_qdp) then
            et = et + tscale*(xxqdp-xxqdp0(j))**2
            et = et + tscale*(xyqdp-xyqdp0(j))**2
            et = et + tscale*(xzqdp-xzqdp0(j))**2
            et = et + tscale*(yyqdp-yyqdp0(j))**2
            et = et + tscale*(yzqdp-yzqdp0(j))**2
            et = et + tscale*(zzqdp-zzqdp0(j))**2
         end if
      end do
      er = sqrt(er/dble(npoint))
      potfit1 = er + ec + et
c
c     compute numerical gradient for electrostatic parameters
c
      m = nvar
      do k = 1, m
         call varprm (nvar,xx,k,-0.5d0*eps)
         er = 0.0d0
         ec = 0.0d0
         et = 0.0d0
         do j = 1, nconf
            call getref (j)
            if (use_mpole) then
               call chkpole
               call rotpole
            end if
            if (use_polar) then
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
               er = er + (epot(1,i,j)-epot(2,i,j))**2
            end do
            call momfull
            ec = ec + cscale*(netchg-dble(nint(netchg)))**2
            if (use_dpl) then
               et = et + tscale*(xdpl-xdpl0(j))**2
               et = et + tscale*(ydpl-ydpl0(j))**2
               et = et + tscale*(zdpl-zdpl0(j))**2
            end if
            if (use_qdp) then
               et = et + tscale*(xxqdp-xxqdp0(j))**2
               et = et + tscale*(xyqdp-xyqdp0(j))**2
               et = et + tscale*(xzqdp-xzqdp0(j))**2
               et = et + tscale*(yyqdp-yyqdp0(j))**2
               et = et + tscale*(yzqdp-yzqdp0(j))**2
               et = et + tscale*(zzqdp-zzqdp0(j))**2
            end if
         end do
         er = sqrt(er/dble(npoint))
         e0 = er + ec + et
         call varprm (nvar,xx,k,0.5d0*eps)
         er = 0.0d0
         ec = 0.0d0
         et = 0.0d0
         do j = 1, nconf
            call getref (j)
            if (use_mpole) then
               call chkpole
               call rotpole
            end if
            if (use_polar) then
               call induce
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               call potpoint (xi,yi,zi,pot)
               epot(1,i,j) = pot
               er = er + (epot(1,i,j)-epot(2,i,j))**2
            end do
            call momfull
            ec = ec + cscale*(netchg-dble(nint(netchg)))**2
            if (use_dpl) then
               et = et + tscale*(xdpl-xdpl0(j))**2
               et = et + tscale*(ydpl-ydpl0(j))**2
               et = et + tscale*(zdpl-zdpl0(j))**2
            end if
            if (use_qdp) then
               et = et + tscale*(xxqdp-xxqdp0(j))**2
               et = et + tscale*(xyqdp-xyqdp0(j))**2
               et = et + tscale*(xzqdp-xzqdp0(j))**2
               et = et + tscale*(yyqdp-yyqdp0(j))**2
               et = et + tscale*(yzqdp-yzqdp0(j))**2
               et = et + tscale*(zzqdp-zzqdp0(j))**2
            end if
         end do
         er = sqrt(er/dble(npoint))
         e = er + ec + et
         g(k) = (e-e0) / eps
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prmvar  --  electrostatics to optimization  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prmvar" determines the optimization values from the
c     corresponding electrostatic potential energy parameters
c
c
      subroutine prmvar (nvar,xx)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'charge.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'potfit.i'
      include 'units.i'
      integer i,ii,it
      integer k,kk,kt
      integer nvar
      real*8 dterm,qterm
      real*8 xx(maxvar)
      logical done
      character*17 prmtyp
c
c
c     zero out the total number of optimization parameters
c
      nvar = 0
      dterm = 1.0d0 / bohr
      qterm = 3.0d0 / bohr**2
c
c     list active atoms when not all are used in optimization
c
      if (nfatm .ne. n) then
         write (iout,10)
   10    format (/,' Atomic Parameters Included in Potential Fitting :',
     &           //,3x,'Atom',10x,'Atom Name',9x,'Atom Type',
     &              9x,'Parameters',/)
         do i = 1, nion
            ii = iion(i)
            if (fatm(ii)) then
               it = type(ii)
               prmtyp = 'Partial Charge'
               write (iout,20)  ii,name(ii),it,prmtyp
   20          format (i6,15x,a3,10x,i6,13x,a)
            end if
         end do
         do i = 1, npole
            ii = ipole(i)
            if (fatm(ii)) then
               it = type(ii)
               prmtyp = 'Atomic Multipoles'
               write (iout,30)  ii,name(ii),it,prmtyp
   30          format (i6,15x,a3,10x,i6,13x,a)
            end if
         end do
      end if
c
c     print header information for electrostatic parameters
c
      write (iout,40)
   40 format (/,' Potential Fitting of Electrostatic Parameters :',
     &        //,1x,'Parameter',6x,'Atom Type',9x,'Category',
     &           12x,'Value',8x,'Fixed',/)
c
c     get optimization parameters from partial charge values
c
      do i = 1, nion
         done = .true.
         ii = iion(i)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            it = type(ii)
            do k = 1, i-1
               kk = iion(k)
               kt = type(kk)
               if (kt.eq.it .and. fatm(kk))  done = .true.
            end do
         end if
         if (.not. done) then
            if (pchg(i) .ne. 0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pchg(i)
               write (iout,50)  nvar,it,'Charge  ',pchg(i)
   50          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,60)  it,'Charge  ',pchg(i)
   60          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
         end if
      end do
c
c     get optimization parameters from atomic multipole values
c
      do i = 1, npole
         done = .true.
         ii = ipole(i)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            it = type(ii)
            do k = 1, i-1
               kk = ipole(k)
               kt = type(kk)
               if (kt.eq.it .and. fatm(kk))  done = .true.
            end do
         end if
         if (.not. done) then
            if (fit_mpl .and. pole(1,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(1,i)
               write (iout,70)  nvar,it,'Monopole',pole(1,i)
   70          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,80)  it,'Monopole',pole(1,i)
   80          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(2,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(2,i)
               write (iout,90)  nvar,it,'X-Dipole',dterm*pole(2,i)
   90          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,100)  it,'X-Dipole',dterm*pole(2,i)
  100          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(3,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(3,i)
               write (iout,110)  nvar,it,'Y-Dipole',dterm*pole(3,i)
  110          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,120)  it,'Y-Dipole',dterm*pole(3,i)
  120          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_dpl .and. pole(4,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(4,i)
               write (iout,130)  nvar,it,'Z-Dipole',dterm*pole(4,i)
  130          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,140)  it,'Z-Dipole',dterm*pole(4,i)
  140          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(5,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(5,i)
               write (iout,150)  nvar,it,'XX-Quad ',qterm*pole(5,i)
  150          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,160)  it,'XX-Quad ',qterm*pole(5,i)
  160          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(6,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(6,i)
               write (iout,170)  nvar,it,'XY-Quad ',qterm*pole(6,i)
  170          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,180)  it,'XY-Quad ',qterm*pole(6,i)
  180          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(7,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(7,i)
               write (iout,190)  nvar,it,'XZ-Quad ',qterm*pole(7,i)
  190          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,200)  it,'XZ-Quad ',qterm*pole(7,i)
  200          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(9,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(9,i)
               write (iout,210)  nvar,it,'YY-Quad ',qterm*pole(9,i)
  210          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,220)  it,'YY-Quad ',qterm*pole(9,i)
  220          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(10,i).ne.0.0d0) then
               nvar = nvar + 1
               xx(nvar) = pole(10,i)
               write (iout,230)  nvar,it,'YZ-Quad ',qterm*pole(10,i)
  230          format (i6,7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,240)  it,'YZ-Quad ',qterm*pole(10,i)
  240          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
            if (fit_qdp .and. pole(13,i).ne.0.0d0) then
               write (iout,250)  it,'ZZ-Quad ',qterm*pole(13,i)
  250          format (4x,'--',7x,i8,13x,a8,5x,f12.5)
            else
               write (iout,260)  it,'ZZ-Quad ',qterm*pole(13,i)
  260          format (4x,'--',7x,i8,13x,a8,5x,f12.5,10x,'X')
            end if
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine varprm  --  optimization to electrostatics  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "varprm" copies the current optimization values into the
c     corresponding electrostatic potential energy parameters
c
c
      subroutine varprm (nvar,xx,ivar,eps)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'charge.i'
      include 'mpole.i'
      include 'potent.i'
      include 'potfit.i'
      integer i,j,k
      integer ii,it
      integer kk,kt
      integer nvar,ivar
      real*8 eps
      real*8 xx(maxvar)
      logical done
c
c
c     zero out the total number of optimization parameters
c
      nvar = 0
c
c     translate optimization values back to partial charges
c
      do i = 1, nion
         done = .true.
         ii = iion(i)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            it = type(ii)
            do k = 1, i-1
               kk = iion(k)
               kt = type(kk)
               if (kt.eq.it .and. fatm(kk))  done = .true.
            end do
         end if
         if (.not. done) then
            if (pchg(i) .ne. 0.0d0) then
               nvar = nvar + 1
               pchg(i) = xx(nvar)
               if (ivar .eq. nvar)  pchg(i) = pchg(i) + eps
               do k = i+1, nion
                  kk = iion(k)
                  kt = type(kk)
                  if (kt.eq.it .and. fatm(kk)) then
                     pchg(k) = pchg(i)
                  end if
               end do
            end if
         end if
      end do
c
c     translate optimization values back to atomic multipoles
c
      do i = 1, npole
         done = .true.
         ii = ipole(i)
         if (fatm(ii))  done = .false.
         if (.not. done) then
            it = type(ii)
            do k = 1, i-1
               kk = ipole(k)
               kt = type(kk)
               if (kt.eq.it .and. fatm(kk))  done = .true.
            end do
         end if
         if (.not. done) then
            if (fit_mpl .and. pole(1,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(1,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(1,i) = pole(1,i) + eps
            end if
            if (fit_dpl .and. pole(2,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(2,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(2,i) = pole(2,i) + eps
            end if
            if (fit_dpl .and. pole(3,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(3,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(3,i) = pole(3,i) + eps
            end if
            if (fit_dpl .and. pole(4,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(4,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(4,i) = pole(4,i) + eps
            end if
            if (fit_qdp .and. pole(5,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(5,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(5,i) = pole(5,i) + eps
            end if
            if (fit_qdp .and. pole(6,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(6,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(6,i) = pole(6,i) + eps
               pole(8,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(8,i) = pole(8,i) + eps
            end if
            if (fit_qdp .and. pole(7,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(7,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(7,i) = pole(7,i) + eps
               pole(11,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(11,i) = pole(11,i) + eps
            end if
            if (fit_qdp .and. pole(9,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(9,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(9,i) = pole(9,i) + eps
            end if
            if (fit_qdp .and. pole(10,i).ne.0.0d0) then
               nvar = nvar + 1
               pole(10,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(10,i) = pole(10,i) + eps
               pole(12,i) = xx(nvar)
               if (ivar .eq. nvar)  pole(12,i) = pole(12,i) + eps
            end if
            if (fit_qdp .and. pole(13,i).ne.0.0d0) then
               pole(13,i) = -pole(5,i) - pole(9,i)
            end if
            do k = i+1, npole
               kk = ipole(k)
               kt = type(kk)
               if (kt.eq.it .and. fatm(kk)) then
                  do j = 1, 13
                     pole(j,k) = pole(j,i)
                  end do
               end if
            end do
         end if
      end do
c
c     check chiral multipoles and rotate into global frame
c
      if (use_mpole) then
         call chkpole
         call rotpole
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine potstat  --  electrostatic potential statistics  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "potstat" computes and prints statistics for the electrostatic
c     potential over a set of grid points
c
c
      subroutine potstat (dofull,domodel,dopair,dotarget)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'potfit.i'
      include 'titles.i'
      integer i,j,k,ipot
      integer npoint
      integer freeunit
      integer trimtext
      integer natm(maxatm)
      real*8 xi,yi,zi
      real*8 pave1,pave2
      real*8 tave,uave,rmsd
      real*8 patm1(maxatm)
      real*8 patm2(maxatm)
      real*8 rmsa(maxatm)
      logical dofull,domodel
      logical dopair,dotarget
      character*120 potfile
c
c
c     output potential values for each model at each point
c
      if (dofull) then
         if (domodel) then
            ipot = freeunit ()
            potfile = filename(1:leng)//'.pot'
            call version (potfile,'new')
            open (unit=ipot,file=potfile,status='new')
         end if
         do j = 1, nconf
            if (nconf .eq. 1) then
               write (iout,10)
   10          format (/,' Electrostatic Potential at Each Grid',
     &                    ' Point :',
     &                 /,8x,'(Kcal/mole per unit charge)')
            else
               write (iout,20)  j
   20          format (/,' Electrostatic Potential at Grid Points',
     &                    ' for Conformer',i4,' :',
     &                 /,12x,'(Kcal/mole per unit charge)')
            end if
            if (dotarget) then
               write (iout,30)
   30          format (/,3x,'Point',15x,'XYZ-Coordinates',15x,
     &                    'Potential',5x,'Target',/)
            else if (dopair) then
               write (iout,40)
   40          format (/,3x,'Point',15x,'XYZ-Coordinates',13x,
     &                    'Potential 1',3x,'Potential 2',/)
            else if (domodel) then
               write (iout,50)
   50          format (/,3x,'Point',15x,'XYZ-Coordinates',14x,
     &                    'Potential',/)
               write (ipot,60)  npgrid(j),title(1:ltitle)
   60          format (i8,2x,a)
            end if
            do i = 1, npgrid(j)
               xi = pgrid(1,i,j)
               yi = pgrid(2,i,j)
               zi = pgrid(3,i,j)
               if (dotarget .or. dopair) then
                  write (iout,70)  i,xi,yi,zi,epot(1,i,j),epot(2,i,j)
   70             format (i8,3x,3f12.6,2x,2f12.4)
               else if (domodel) then
                  write (iout,80)  i,xi,yi,zi,epot(1,i,j)
   80             format (i8,3x,3f12.6,2x,f12.4)
                  write (ipot,90)  i,xi,yi,zi,epot(1,i,j)
   90             format (i8,3x,3f12.6,2x,f12.4)
               end if
            end do
         end do
         if (domodel) then
            close (unit=ipot)
            write (iout,100)  potfile(1:trimtext(potfile))
  100       format (/,' Electrostatic Potential Written to File :  ',a)
         end if
      end if
c
c     find average electrostatic potential around each atom
c
      write (iout,110)
  110 format (/,' Average Electrostatic Potential over Atoms :',
     &        /,6x,'(Kcal/mole per unit charge)')
      if (dotarget) then
         write (iout,120)
  120    format (/,4x,'Atom',9x,'Points',8x,'Potential',
     &              8x,'Target',8x,'RMS Diff',/)
      else if (dopair) then
         write (iout,130)
  130    format (/,4x,'Atom',9x,'Points',7x,'Potential 1',
     &              4x,'Potential 2',6x,'RMS Diff',/)
      else if (domodel) then
         write (iout,140)
  140    format (/,4x,'Atom',9x,'Points',8x,'Potential',/)
      end if
      do i = 1, n
         natm(i) = 0
         patm1(i) = 0.0d0
         patm2(i) = 0.0d0
         rmsa(i) = 0.0d0
      end do
      do j = 1, nconf
         do i = 1, npgrid(j)
            k = ipgrid(i,j)
            natm(k) = natm(k) + 1
            patm1(k) = patm1(k) + epot(1,i,j)
            patm2(k) = patm2(k) + epot(2,i,j)
            rmsa(k) = rmsa(k) + (epot(1,i,j)-epot(2,i,j))**2
         end do
      end do
      do i = 1, n
         if (natm(i) .ne. 0) then
            patm1(i) = patm1(i) / dble(natm(i))
            patm2(i) = patm2(i) / dble(natm(i))
            rmsa(i) = sqrt(rmsa(i)/dble(natm(i)))
         end if
         if (gatm(i)) then
            if (dotarget .or. dopair) then
               write (iout,150)  i,natm(i),patm1(i),patm2(i),rmsa(i)
  150          format (i8,6x,i8,5x,f12.4,3x,f12.4,3x,f12.4)
            else if (domodel) then
               write (iout,160)  i,natm(i),patm1(i)
  160          format (i8,6x,i8,5x,f12.4)
            end if
         end if
      end do
c
c     overall averages for the sets of electrostatic potentials
c
      npoint = 0
      pave1 = 0.0d0
      pave2 = 0.0d0
      tave = 0.0d0
      uave = 0.0d0
      rmsd = 0.0d0
      do j = 1, nconf
         npoint = npoint + npgrid(j)
         do i = 1, npgrid(j)
            pave1 = pave1 + abs(epot(1,i,j))
            pave2 = pave2 + abs(epot(2,i,j))
            tave = tave + epot(1,i,j) - epot(2,i,j)
            uave = uave + abs(epot(1,i,j)-epot(2,i,j))
            rmsd = rmsd + (epot(1,i,j)-epot(2,i,j))**2
         end do
      end do
      pave1 = pave1 / dble(npoint)
      pave2 = pave2 / dble(npoint)
      tave = tave / dble(npoint)
      uave = uave / dble(npoint)
      rmsd = sqrt(rmsd/dble(npoint))
      if (dopair) then
         write (iout,170)  pave1
  170    format (/,' Electrostatic Potential over all Grid Points :',
     &           //,' Average Magnitude for Potential 1 :',6x,f12.4)
      else
         write (iout,180)  pave1
  180    format (/,' Electrostatic Potential over all Grid Points :',
     &           //,' Average Magnitude for Potential :',8x,f12.4)
      end if
      if (dotarget) then
         write (iout,190)  pave2,tave,uave,rmsd
  190    format (' Average Magnitude for Target :',11x,f12.4,
     &           /,' Average Signed Potential Difference :',4x,f12.4,
     &           /,' Average Unsigned Potential Difference :',2x,f12.4,
     &           /,' Root Mean Square Potential Difference :',2x,f12.4)
      else if (dopair) then
         write (iout,200)  pave2,tave,uave,rmsd
  200    format (' Average Magnitude for Potential 2 :',6x,f12.4,
     &           /,' Average Signed Potential Difference :',4x,f12.4,
     &           /,' Average Unsigned Potential Difference :',2x,f12.4,
     &           /,' Root Mean Square Potential Difference :',2x,f12.4)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine prtfit  --  create file with optimal parameters  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "prtfit" makes a key file containing results from fitting a
c     charge or multipole model to an electrostatic potential grid
c
c
      subroutine prtfit
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'charge.i'
      include 'files.i'
      include 'keys.i'
      include 'mpole.i'
      include 'potfit.i'
      include 'units.i'
      integer i,j,k
      integer ii,kk
      integer it,kt
      integer ix,iy,iz
      integer ikey,size
      integer freeunit
      integer trimtext
      real*8 dterm,qterm
      real*8 eps,sum,big
      logical doprint,header
      character*120 keyfile
      character*120 record
c
c
c     convert dipole and quadrupole moments to atomic units
c
      dterm = 1.0d0 / bohr
      qterm = 3.0d0 / bohr**2
      do i = 1, npole
         do j = 2, 4
            pole(j,i) = dterm * pole(j,i)
         end do
         do j = 5, 13
            pole(j,i) = qterm * pole(j,i)
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
c     maintain traceless quadrupole at each multipole site
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
c
c     open a new keyfile to contain the optimized parameters
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
c     print a header for the fitted multipole parameters
c
      write (ikey,20)
   20 format (/,'#',/,'# Results of Electrostatic Potential Fitting',
     &           /,'#')
c
c     output the optimized charge values to the keyfile
c
      header = .true.
      do i = 1, nion
         ii = iion(i)
         it = type(ii)
         doprint = .false.
         if (fatm(ii)) then
            doprint = .true.
            do k = 1, i-1
               kk = iion(k)
               kt = type(kk)
               if (fatm(kk) .and. it.eq.kt)  doprint = .false.
            end do
         end if
         if (doprint) then
            if (header) then
               header = .false.
               write (ikey,30)
   30          format ()
            end if
            write (ikey,40)  it,pchg(i)
   40       format ('charge',4x,i5,10x,f11.4)
         end if
      end do
c
c     output the optimized multipole values to the keyfile
c
      header = .true.
      do i = 1, npole
         ii = ipole(i)
         it = type(ii)
         doprint = .false.
         if (fatm(ii)) then
            doprint = .true.
            do k = 1, i-1
               kk = ipole(k)
               kt = type(kk)
               if (fatm(kk) .and. it.eq.kt)  doprint = .false.
            end do
         end if
         if (doprint) then
            if (header) then
               header = .false.
               write (ikey,50)
   50          format ()
            end if
            iz = zaxis(i)
            ix = xaxis(i)
            iy = yaxis(i)
            if (iz .ne. 0)  iz = type(iz)
            if (ix .ne. 0)  ix = type(ix)
            if (iy .ne. 0)  iy = type(iy)
            if (polaxe(i) .eq. 'Z-then-X') then
               if (yaxis(i) .eq. 0) then
                  write (ikey,60)  it,iz,ix,pole(1,i)
   60             format ('multipole',1x,3i5,11x,f11.5)
               else
                  write (ikey,70)  it,iz,ix,iy,pole(1,i)
   70             format ('multipole',1x,4i5,6x,f11.5)
               end if
            else if (polaxe(i) .eq. 'Bisector') then
               if (yaxis(i) .eq. 0) then
                  write (ikey,80)  it,iz,-ix,pole(1,i)
   80             format ('multipole',1x,3i5,11x,f11.5)
               else
                  write (ikey,90)  it,iz,-ix,iy,pole(1,i)
   90             format ('multipole',1x,4i5,6x,f11.5)
               end if
            else if (polaxe(i) .eq. 'Z-Bisect') then
               write (ikey,100)  it,iz,-ix,-iy,pole(1,i)
  100          format ('multipole',1x,4i5,6x,f11.5)
            else if (polaxe(i) .eq. '3-Fold') then
               write (ikey,110)  it,-iz,-ix,-iy,pole(1,i)
  110          format ('multipole',1x,4i5,6x,f11.5)
            end if
            write (ikey,120)  pole(2,i),pole(3,i),pole(4,i)
  120       format (36x,3f11.5)
            write (ikey,130)  pole(5,i)
  130       format (36x,f11.5)
            write (ikey,140)  pole(8,i),pole(9,i)
  140       format (36x,2f11.5)
            write (ikey,150)  pole(11,i),pole(12,i),pole(13,i)
  150       format (36x,3f11.5)
         end if
      end do
      close (unit=ikey)
      return
      end