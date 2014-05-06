c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program xtalfit  --  fit parameters to crystal structures  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "xtalfit" computes an optimized set of potential energy
c     parameters for user specified van der Waals and electrostatic
c     interactions by fitting to crystal structure, lattice energy
c     and monomer dipole moment data
c
c
      program xtalfit
      use sizes
      use files
      use iounit
      use molcul
      use potent
      use xtals
      implicit none
      integer i,ixtal
      integer atom1,atom2
      integer nresid,prmtyp
      real*8 grdmin
      real*8 xx(maxlsq)
      real*8 xhi(maxlsq)
      real*8 xlo(maxlsq)
      real*8 f(maxrsd)
      real*8 g(maxlsq)
      real*8 jacobian(maxrsd,maxlsq)
      logical exist,query
      character*16 blank
      character*16 label(6)
      character*120 record
      character*120 string
      external xtalerr,xtalwrt
c
c
c     initialize some variables to be used during fitting
c
      call initial
      nvary = 0
      nresid = 0
      blank = '                '
      do i = 1, maxlsq
         vartyp(i) = blank
      end do
      do i = 1, maxrsd
         rsdtyp(i) = blank
      end do
c
c     print informational header about available parameters
c
      write (iout,10)
   10 format (/,' The Following Parameters can be Fit for',
     &           ' each Atom Type :',
     &        //,4x,'(1) Van der Waals Atomic Radius',
     &        /,4x,'(2) Van der Waals Well Depth',
     &        /,4x,'(3) Hydrogen Atom Reduction Factor',
     &        /,4x,'(4) Atomic Partial Charge',
     &        /,4x,'(5) Bond Dipole Moment Magnitude',
     &        /,4x,'(6) Bond Dipole Moment Position')
c
c     get types of potential parameters to be optimized
c
      query = .true.
      do while (query)
         prmtyp = -1
         atom1 = 0
         atom2 = 0
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  prmtyp
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=20,end=20)  atom2
   20    continue
         if (prmtyp .ne. 0) then
            prmtyp = 0
            write (iout,30)
   30       format (/,' Enter Parameter Type then Atom Class',
     &                 ' or Type(s) :  ',$)
            read (input,40)  record
   40       format (a120)
            read (record,*,err=50,end=50)  prmtyp,atom1,atom2
   50       continue
         end if
         if (prmtyp .eq. 0) then
            query = .false.
         else
            query = .true.
            nvary = nvary + 1
            ivary(nvary) = prmtyp
            if (prmtyp .lt. 5) then
               vary(1,nvary) = atom1
            else
               vary(1,nvary) = min(atom1,atom2)
               vary(2,nvary) = max(atom1,atom2)
            end if
         end if
      end do
c
c     get termination criterion as RMS gradient over parameters
c
      grdmin = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=60,end=60)  grdmin
   60 continue
      if (grdmin .le. 0.0d0) then
         write (iout,70)
   70    format (/,' Enter RMS Gradient Termination Criterion',
     &              ' [0.1] :  ',$)
         read (input,80)  grdmin
   80    format (f20.0)
      end if
      if (grdmin .le. 0.0d0)  grdmin = 0.1d0
c
c     get number of crystal structures to use in optimization
c
      nxtal = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=90,end=90)  nxtal
   90 continue
      if (nxtal .le. 0) then
         write (iout,100)
  100    format (/,' Enter Number of Crystals to be Used [1] :  ',$)
         read (input,110)  nxtal
  110    format (i10)
      end if
c
c     check for too few or too many crystal structures
c
      if (nxtal .eq. 0)  nxtal = 1
      if (nxtal .gt. maxref) then
         write (iout,120)
  120    format (/,' XTALFIT  --  Too many Crystals,',
     &              ' Increase the Value of MAXREF')
         call fatal
      end if
c
c     get the coordinates and parameters for current crystal
c
      do ixtal = 1, nxtal
         call getxyz
         call mechanic
c
c     get an ideal value for the crystal lattice energy
c
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=130,end=130)  e0_lattice
            query = .false.
         end if
  130    continue
         if (query) then
            write (iout,140)
  140       format (/,' Enter Lattice Energy Value [<CR> to omit]',
     &                 ' :  ',$)
            read (input,150)  e0_lattice
  150       format (f20.0)
         end if
         if (e0_lattice .gt. 0.0d0)  e0_lattice = -e0_lattice
c
c     get an ideal value for the monomer dipole moment
c
         query = .true.
         call nextarg (string,exist)
         if (exist) then
            read (string,*,err=160,end=160)  moment_0
            query = .false.
         end if
  160    continue
         if (query) then
            write (iout,170)
  170       format (/,' Enter Monomer Dipole Moment [<CR> to Omit]',
     &                 ' :   ',$)
            read (input,180)  moment_0
  180       format (f20.0)
         end if
c
c     set the types of residuals for use in optimization
c
         do i = 1, 6
            iresid(nresid+i) = ixtal
         end do
         rsdtyp(nresid+1) = 'Force a-Axis    '
         rsdtyp(nresid+2) = 'Force b-Axis    '
         rsdtyp(nresid+3) = 'Force c-Axis    '
         rsdtyp(nresid+4) = 'Force Alpha     '
         rsdtyp(nresid+5) = 'Force Beta      '
         rsdtyp(nresid+6) = 'Force Gamma     '
         nresid = nresid + 6
c
c     print molecules per unit cell, lattice energy and dipole
c
         write (iout,190)  ixtal,filename(1:35),nmol
  190    format (/,' File Name of Crystal Structure',i4,' :  ',5x,a35,
     &           /,' Number of Molecules per Unit Cell :',i13)
         if (e0_lattice .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            rsdtyp(nresid) = 'Lattice Energy  '
            write (iout,200)  e0_lattice
  200       format (' Value of Crystal Lattice Energy :  ',f13.2)
         end if
         if (moment_0 .ne. 0.0d0) then
            nresid = nresid + 1
            iresid(nresid) = ixtal
            rsdtyp(nresid) = 'Dipole Moment   '
            write (iout,210)  moment_0
  210       format (' Value of Molecular Dipole Moment : ',f13.2)
         end if
c
c     set the initial values of the parameters
c
         call xtalprm ('STORE',ixtal,xx)
      end do
c
c     turn off all local interactions and extra terms
c
      call potoff
      use_vdw = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
c
c     types of variables for use in optimization
c
      label(1) = 'Atomic Radius   '
      label(2) = 'Well Depth      '
      label(3) = 'H Reduction     '
      label(4) = 'Partial Charge  '
      label(5) = 'Dipole Magnitude'
      label(6) = 'Dipole Position '
      do i = 1, nvary
         vartyp(i) = label(ivary(i))
      end do
c
c     print the initial parameter values
c
      write (iout,220)
  220 format (/,' Initial Values of the Parameters :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,230)  i,vartyp(i),vary(1,i),xx(i)
  230       format (3x,'(',i2,')',2x,a16,4x,'Atom Class',i5,4x,f12.4)
         else if (ivary(i).eq.4 .or. ivary(i).eq.5) then
            write (iout,240)  i,vartyp(i),vary(1,i),xx(i)
  240       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,4x,f12.4)
         else if (ivary(i) .eq. 6) then
            write (iout,250)  i,vartyp(i),vary(1,i),vary(2,i),xx(i)
  250       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,f12.4)
         end if
      end do
c
c     set upper and lower bounds based on the parameter type
c
      do i = 1, nvary
         if (ivary(i) .eq. 1) then
            xlo(i) = 0.9d0 * xx(i)
            xhi(i) = 1.1d0 * xx(i)
         else if (ivary(i) .eq. 2) then
            xlo(i) = 0.8d0 * xx(i)
            xhi(i) = 1.2d0 * xx(i)
         else if (ivary(i) .eq. 3) then
            xlo(i) = 0.8d0 * xx(i)
            xhi(i) = 1.2d0 * xx(i)
            if (xhi(i) .gt. 1.0d0)  xhi(i) = 1.0d0
         else if (ivary(i) .eq. 4) then
            xlo(i) = min(xx(i)-0.2d0*abs(xx(i)),xx(i)-0.2d0)
            xhi(i) = max(xx(i)+0.2d0*abs(xx(i)),xx(i)+0.2d0)
         else if (ivary(i) .eq. 5) then
            xlo(i) = min(xx(i)-0.2d0*abs(xx(i)),xx(i)-0.2d0)
            xhi(i) = max(xx(i)+0.2d0*abs(xx(i)),xx(i)+0.2d0)
         else if (ivary(i) .eq. 6) then
            xlo(i) = 0.8d0 * xx(i)
            xhi(i) = 1.2d0 * xx(i)
            if (xlo(i) .lt. 0.2d0)  xhi(i) = 0.2d0
            if (xhi(i) .gt. 0.8d0)  xhi(i) = 0.8d0
         end if
      end do
c
c     use nonlinear least squares to refine the parameters
c
      call square (nresid,nvary,xlo,xhi,xx,f,g,jacobian,
     &                maxrsd,grdmin,xtalerr,xtalwrt)
c
c     print the final parameter values
c
      write (iout,260)
  260 format (/,' Final Values of the Parameters and Scaled',
     &            ' Derivatives :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,270)  i,vartyp(i),vary(1,i),xx(i),g(i)
  270       format (3x,'(',i2,')',2x,a16,4x,'Atom Class',i5,4x,2f12.4)
         else if (ivary(i).eq.4 .or. ivary(i).eq.5) then
            write (iout,280)  i,vartyp(i),vary(1,i),xx(i),g(i)
  280       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,4x,2f12.4)
         else if (ivary(i) .eq. 6) then
            write (iout,290)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),g(i)
  290       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,2f12.4)
         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalprm  --  energy/optimization conversion  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalprm" stores or retrieves a crystal structure; used
c     to make a previously stored structure the currently active
c     structure, or to store a structure for later use
c
c     the current version only provides for intermolecular potential
c     energy terms
c
c
      subroutine xtalprm (mode,ixtal,xx)
      use sizes
      use atoms
      use atomid
      use boxes
      use charge
      use dipole
      use files
      use fracs
      use kvdws
      use molcul
      use vdw
      use xtals
      implicit none
      integer i,j,k
      integer ixtal,prmtyp
      integer atom1,atom2
      real*8 xmid,ymid,zmid
      real*8 e0_lattices(maxref)
      real*8 moment_0s(maxref)
      real*8 xx(*)
      logical first
      character*5 mode
      save e0_lattices
      save moment_0s
      save first
      data first  / .true. /
c
c
c     save or restore the key values for the current crystal
c
      if (mode .eq. 'STORE') then
         call makeref (ixtal)
      else if (mode .eq. 'RESET') then
         call getref (ixtal)
         call basefile (filename)
         call mechanic
         call bounds
      end if
c
c     perform dynamic allocation of some global arrays
c
      if (mode .eq. 'RESET') then
         if (first) then
            first = .false.
            allocate (xfrac(nmol))
            allocate (yfrac(nmol))
            allocate (zfrac(nmol))
         end if
      end if
c
c     coordinates of molecular centers of mass
c
      if (mode .eq. 'RESET') then
         do i = 1, nmol
            xmid = 0.0d0
            ymid = 0.0d0
            zmid = 0.0d0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               xmid = xmid + x(k)*mass(k)
               ymid = ymid + y(k)*mass(k)
               zmid = zmid + z(k)*mass(k)
            end do
            zmid = zmid / gamma_term
            ymid = (ymid - zmid*beta_term) / gamma_sin
            xmid = xmid - ymid*gamma_cos - zmid*beta_cos
            xfrac(i) = xmid / (xbox * molmass(i))
            yfrac(i) = ymid / (ybox * molmass(i))
            zfrac(i) = zmid / (zbox * molmass(i))
         end do
      end if
c
c     values of ideal lattice energy and dipole moment
c
      if (mode .eq. 'STORE') then
         e0_lattices(ixtal) = e0_lattice
         moment_0s(ixtal) = moment_0
      else if (mode .eq. 'RESET') then
         e0_lattice = e0_lattices(ixtal)
         moment_0 = moment_0s(ixtal)
      end if
c
c     store or reset values of the optimization variables
c
      do j = 1, nvary
         prmtyp = ivary(j)
         atom1 = vary(1,j)
         if (prmtyp .eq. 1) then
            if (mode .eq. 'STORE') then
               xx(j) = rad(atom1)
            else if (mode .eq. 'RESET') then
               rad(atom1) = xx(j)
               do i = 1, maxclass
                  radmin(i,atom1) = rad(i) + rad(atom1)
                  radmin(atom1,i) = radmin(i,atom1)
               end do
            end if
         else if (prmtyp .eq. 2) then
            if (mode .eq. 'STORE') then
               xx(j) = eps(atom1)
            else if (mode .eq. 'RESET') then
               eps(atom1) = abs(xx(j))
               do i = 1, maxclass
                  epsilon(i,atom1) = sqrt(eps(i) * eps(atom1))
                  epsilon(atom1,i) = epsilon(i,atom1)
               end do
            end if
         else if (prmtyp .eq. 3) then
            if (mode .eq. 'STORE') then
               do i = 1, n
                  if (class(i) .eq. atom1) then
                     xx(j) = kred(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, n
                  if (class(i) .eq. atom1)  kred(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 4) then
            if (mode .eq. 'STORE') then
               do i = 1, nion
                  if (type(iion(i)) .eq. atom1) then
                     xx(j) = pchg(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, n
                  if (type(iion(i)) .eq. atom1)  pchg(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 5) then
            atom2 = vary(2,j)
            if (mode .eq. 'STORE') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2) then
                     xx(j) = bdpl(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2)  bdpl(i) = xx(j)
               end do
            end if
         else if (prmtyp .eq. 6) then
            atom2 = vary(2,j)
            if (mode .eq. 'STORE') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2) then
                     xx(j) = sdpl(i)
                     goto 10
                  end if
               end do
            else if (mode .eq. 'RESET') then
               do i = 1, ndipole
                  if (type(idpl(1,i)).eq.atom1 .and.
     &                type(idpl(2,i)).eq.atom2)  sdpl(i) = xx(j)
               end do
            end if
         end if
   10    continue
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine xtalerr  --  error function for xtalfit  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "xtalerr" computes an error function value derived from
c     derivatives with respect to lattice parameters, lattice
c     energy and monomer dipole moments
c
c
      subroutine xtalerr (nresid,nvaried,xx,resid)
      use sizes
      use atoms
      use boxes
      use bound
      use charge
      use dipole
      use math
      use molcul
      use vdw
      use xtals
      implicit none
      integer i,i1,i2,ixtal
      integer nresid,nvaried
      integer n_old,nvdw_old
      integer ndipole_old,nion_old
      real*8 energy,eps,temp
      real*8 e_xtal,e_monomer
      real*8 e_lattice,moment
      real*8 xi,yi,zi,fik
      real*8 x_dpm,y_dpm,z_dpm
      real*8 e_xbox,e_ybox,e_zbox
      real*8 e_alpha,e_beta,e_gamma
      real*8 g_xbox,g_ybox,g_zbox
      real*8 g_alpha,g_beta,g_gamma
      real*8 xx(*)
      real*8 resid(*)
c
c
c     zero out the number of residual functions
c
      nresid = 0
c
c     set the values of the potential energy parameters
c     and get the energy of the base structure
c
      do ixtal = 1, nxtal
         call xtalprm ('RESET',ixtal,xx)
         e_xtal = energy ()
c
c     get energy derivative with respect to lattice a-axis
c
         eps = 0.00001d0
         temp = xbox
         xbox = xbox + eps
         call xtalmove
         e_xbox = energy ()
         xbox = temp
         g_xbox = (e_xbox - e_xtal) / eps
c
c     get energy derivative with respect to lattice b-axis
c
         temp = ybox
         ybox = ybox + eps
         call xtalmove
         e_ybox = energy ()
         ybox = temp
         g_ybox = (e_ybox - e_xtal) / eps
c
c     get energy derivative with respect to lattice c-axis
c
         temp = zbox
         zbox = zbox + eps
         call xtalmove
         e_zbox = energy ()
         zbox = temp
         g_zbox = (e_zbox - e_xtal) / eps
c
c     get energy derivative with respect to lattice alpha
c
         temp = alpha
         alpha = alpha + radian*eps
         call xtalmove
         e_alpha = energy ()
         alpha = temp
         g_alpha = (e_alpha - e_xtal) / eps
c
c     get energy derivative with respect to lattice beta
c
         temp = beta
         beta = beta + radian*eps
         call xtalmove
         e_beta = energy ()
         beta = temp
         g_beta = (e_beta - e_xtal) / eps
c
c     get energy derivative with respect to lattice gamma
c
         temp = gamma
         gamma = gamma + radian*eps
         call xtalmove
         e_gamma = energy ()
         gamma = temp
         g_gamma = (e_gamma - e_xtal) / eps
         call xtalmove
c
c     setup to compute properties of monomer; assumes that
c     molecules are contiguous in the coordinates file
c
         use_bounds = .false.
         use_replica = .false.
         n_old = n
         nvdw_old = nvdw
         nion_old = nion
         ndipole_old = ndipole
         n = n / nmol
         nvdw = nvdw / nmol
         nion = nion / nmol
         ndipole = ndipole / nmol
c
c     compute the crystal lattice energy
c
         e_monomer = energy ()
         e_lattice = (e_xtal - nmol*e_monomer) / dble(nmol)
c
c     compute the dipole moment of the monomer
c
         x_dpm = 0.0d0
         y_dpm = 0.0d0
         z_dpm = 0.0d0
         do i = 1, ndipole
            i1 = idpl(1,i)
            i2 = idpl(2,i)
            xi = x(i2) - x(i1)
            yi = y(i2) - y(i1)
            zi = z(i2) - z(i1)
            fik = bdpl(i) / sqrt(xi*xi + yi*yi + zi*zi)
            x_dpm = x_dpm - xi*fik
            y_dpm = y_dpm - yi*fik
            z_dpm = z_dpm - zi*fik
         end do
         moment = sqrt(x_dpm**2 + y_dpm**2 + z_dpm**2)
c
c     convert back from monomer to full crystal
c
         n = n_old
         nvdw = nvdw_old
         nion = nion_old
         ndipole = ndipole_old
         use_bounds = .true.
         use_replica = .true.
c
c     compute the residual vector as a collection
c     of the forces on the crystal lattice, plus any
c     lattice energy and dipole moment deviations
c
         nresid = nresid + 1
         resid(nresid) = g_xbox
         nresid = nresid + 1
         resid(nresid) = g_ybox
         nresid = nresid + 1
         resid(nresid) = g_zbox
         nresid = nresid + 1
         resid(nresid) = g_alpha
         nresid = nresid + 1
         resid(nresid) = g_beta
         nresid = nresid + 1
         resid(nresid) = g_gamma
         if (e0_lattice .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = e_lattice - e0_lattice
         end if
         if (moment_0 .ne. 0.0d0) then
            nresid = nresid + 1
            resid(nresid) = 2.0d0 * (moment-moment_0)
         end if
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine xtalmove  --  translation of rigid molecules  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "xtalmove" converts fractional to Cartesian coordinates for
c     rigid molecules during fitting of force field parameters to
c     crystal structure data
c
c
      subroutine xtalmove
      use sizes
      use atoms
      use atomid
      use boxes
      use fracs
      use molcul
      implicit none
      integer i,j,k
      integer init,stop
      real*8 weigh
      real*8 xmid,ymid,zmid
      real*8, allocatable :: xoff(:)
      real*8, allocatable :: yoff(:)
      real*8, allocatable :: zoff(:)
c
c
c     get values for fractional coordinate interconversion
c
      call lattice
c
c     perform dynamic allocation of some local arrays
c
      allocate (xoff(n))
      allocate (yoff(n))
      allocate (zoff(n))
c
c     locate the center of mass of each molecule
c
      do i = 1, nmol
         init = imol(1,i)
         stop = imol(2,i)
         xmid = 0.0d0
         ymid = 0.0d0
         zmid = 0.0d0
         do j = init, stop
            k = kmol(j)
            weigh = mass(k)
            xmid = xmid + x(k)*weigh
            ymid = ymid + y(k)*weigh
            zmid = zmid + z(k)*weigh
         end do
         weigh = molmass(i)
         xmid = xmid / weigh
         ymid = ymid / weigh
         zmid = zmid / weigh
c
c     save atomic coordinates relative to center of mass
c
         do j = init, stop
            k = kmol(j)
            xoff(k) = x(k) - xmid
            yoff(k) = y(k) - ymid
            zoff(k) = z(k) - zmid
         end do
c
c     convert fractional center of mass to Cartesian coordinates
c
         xmid = xfrac(i)*xbox + yfrac(i)*ybox*gamma_cos
     &             + zfrac(i)*zbox*beta_cos
         ymid = yfrac(i)*ybox*gamma_sin + zfrac(i)*zbox*beta_term
         zmid = zfrac(i)*zbox*gamma_term
c
c     translate coordinates via offset from center of mass
c
         do j = init, stop
            k = kmol(j)
            x(k) = xoff(k) + xmid
            y(k) = yoff(k) + ymid
            z(k) = zoff(k) + zmid
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xoff)
      deallocate (yoff)
      deallocate (zoff)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine xtalwrt  --  output optimization parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "xtalwrt" is a utility that prints intermediate results
c     during fitting of force field parameters to crystal data
c
c
      subroutine xtalwrt (niter,xx,gs,nresid,f)
      use iounit
      use xtals
      implicit none
      integer i,niter
      integer nresid
      real*8 xx(*)
      real*8 gs(*)
      real*8 f(*)
c
c
c     write the values of parameters and scaled derivatives
c
      write (iout,10)  niter
   10 format (/,' Parameters and Scaled Derivatives at',
     &          ' Iteration',i4,' :',/)
      do i = 1, nvary
         if (ivary(i) .le. 3) then
            write (iout,20)  i,vartyp(i),vary(1,i),xx(i),gs(i)
   20       format (3x,'(',i2,')',2x,a16,4x,'Atom Class',i5,4x,2f12.4)
         else if (ivary(i).eq.4 .or. ivary(i).eq.5) then
            write (iout,30)  i,vartyp(i),vary(1,i),xx(i),gs(i)
   30       format (3x,'(',i2,')',2x,a16,4x,'Atom Type ',i5,4x,2f12.4)
         else if (ivary(i) .eq. 6) then
            write (iout,40)  i,vartyp(i),vary(1,i),vary(2,i),xx(i),gs(i)
   40       format (3x,'(',i2,')',2x,a16,4x,'Bond Type ',2i5,2f12.4)
         end if
      end do
c
c     write the values of the residual functions
c
      write (iout,50)  niter
   50 format (/,' Residual Error Function Values at Iteration',
     &           i4,' :',/)
      do i = 1, nresid
         write (iout,60)  i,rsdtyp(i),iresid(i),f(i)
   60    format (3x,'(',i2,')',2x,a16,4x,2x,'Crystal',i4,4x,f12.4)
      end do
      write (iout,70)
   70 format ()
      return
      end
