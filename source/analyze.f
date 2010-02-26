c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze  --  energy partitioning and analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential energy;
c     options are provided to display system and force field info,
c     partition the energy by atom or by potential function type,
c     show force field parameters by atom; output the large energy
c     interactions and find electrostatic and inertial properties
c
c
      program analyze
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,ixyz
      integer frame
      integer freeunit
      integer trimtext
      integer list(20)
      logical dosystem,doparam
      logical doenergy,doatom
      logical dolarge,dodetail
      logical doprops
      logical exist
      logical active(maxatm)
      character*1 letter
      character*120 record
      character*120 string
      character*120 xyzfile
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get the desired types of analysis to be performed
c
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' The TINKER Analysis Facility can Provide :',
     &           //,' General System and Force Field Information [G]',
     &           /,' Force Field Parameters for Interactions [P]',
     &           /,' Total Potential Energy and its Components [E]',
     &           /,' Energy Breakdown over Each of the Atoms [A]',
     &           /,' List of the Large Individual Interactions [L]',
     &           /,' Details for All Individual Interactions [D]',
     &           /,' Electrostatic, Inertial & Virial Properties [M]')
   20    continue
         write (iout,30)
   30    format (/,' Enter the Desired Analysis Types',
     &              ' [G,P,E,A,L,D,M] :  ',$)
         read (input,40,err=20)  string
   40    format (a120)
      end if
c
c     set option control flags based desired analysis types
c
      dosystem = .false.
      doparam = .false.
      doenergy = .false.
      doatom = .false.
      dolarge = .false.
      dodetail = .false.
      doprops = .false.
      call upcase (string)
      do i = 1, trimtext(string)
         letter = string(i:i)
         if (letter .eq. 'G')  dosystem = .true.
         if (letter .eq. 'P')  doparam = .true.
         if (letter .eq. 'E')  doenergy = .true.
         if (letter .eq. 'A')  doatom = .true.
         if (letter .eq. 'L')  dolarge = .true.
         if (letter .eq. 'D')  dodetail = .true.
         if (letter .eq. 'M')  doprops = .true.
      end do
c
c     get the list of atoms for which output is desired
c
      if (doatom .or. doparam) then
         do i = 1, 20
            list(i) = 0
         end do
         if (exist) then
            do i = 1, 20
               call nextarg (string,exist)
               if (.not. exist)  goto 50
               read (string,*,err=50,end=50)  list(i)
            end do
   50       continue
         else
            write (iout,60)
   60       format (/,' List Atoms for which Output is Desired',
     &                 ' [ALL] :  '/,'    >  ',$)
            read (input,70)  record
   70       format (a120)
            read (record,*,err=80,end=80)  (list(i),i=1,20)
   80       continue
         end if
         do i = 1, n
            active(i) = .true.
         end do
         i = 1
         dowhile (list(i) .ne. 0)
            if (i .eq. 1) then
               do j = 1, n
                  active(j) = .false.
               end do
            end if
            if (list(i) .gt. 0) then
               active(list(i)) = .true.
               i = i + 1
            else
               do j = abs(list(i)), abs(list(i+1))
                  active(j) = .true.
               end do
               i = i + 2
            end if
         end do
      end if
c
c     setup to write out the large individual energy terms
c
      if (dolarge) then
         verbose = .true.
      else
         verbose = .false.
      end if
c
c     setup to write out all of the individual energy terms
c
      if (dodetail) then
         doenergy = .true.
         debug = .true.
      else
         debug = .false.
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz')
      call version (xyzfile,'old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     get info on the molecular system and force field
c
      if (dosystem)  call systyze
c
c     get parameters used for molecular mechanics potentials
c
      if (doparam)  call paramyze (active)
c
c     perform analysis for each successive coordinate structure
c
      dowhile (.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            write (iout,90)  frame
   90       format (/,' Analysis for Archive Structure :',8x,i8)
         end if
c
c     make the call to compute the potential energy
c
         if (doenergy .or. doatom .or. dolarge)  call enrgyze
c
c     energy partitioning by potential energy components
c
         if (doenergy) then
            if (digits .ge. 8) then
               call analyz8
            else if (digits .ge. 6) then
               call analyz6
            else
               call analyz4
            end if
         end if
c
c     get various electrostatic and inertial properties
c
         if (doprops) then
            debug = .false.
            call propyze
            if (dodetail)  debug = .true.
         end if
c
c     energy partitioning over the individual atoms
c
         if (doatom)  call atomyze (active)
c
c     attempt to read next structure from the coordinate file
c
         call readxyz (ixyz)
      end do
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      if (dodetail)  debug = .false.
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine systyze  --  system & force field information  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "systyze" is an auxiliary routine for the analyze program
c     that prints general information about the molecular system
c     the the force field model
c
c
      subroutine systyze
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'ewald.i'
      include 'fields.i'
      include 'iounit.i'
      include 'molcul.i'
      include 'pme.i'
      include 'potent.i'
      include 'units.i'
      include 'vdwpot.i'
      integer i
      real*8 dens
      character*20 value
      character*20 label(5)
c
c
c     info on number of atoms, molecules and mass
c
      if (n .ne. 0) then
         write (iout,10)  n,nmol,totmass
   10    format (/,' Overall System Contents :',
     &           //,' Number of Atoms',25x,i8,
     &           /,' Number of Molecules',21x,i8,
     &           /,' Total System Mass',19x,f12.4)
         if (use_bounds) then
            dens = (1.0d24/volbox) * (totmass/avogadro)
            write (iout,20)  dens
   20       format (' System Density',22x,f12.4)
         end if
      end if
c
c     periodic box dimensions and crystal lattice type
c
      if (use_bounds) then
         value = 'ORTHOGONAL'
         if (monoclinic)  value = 'MONOCLINIC'
         if (triclinic)  value = 'TRICLINIC'
         if (octahedron)  value = 'TRUNCATED OCTAHEDRON'
         call justify (value)
         write (iout,30)  xbox,ybox,zbox,alpha,beta,gamma,volbox,value
   30    format (/,' Periodic Boundary Conditions :',
     &           //,' a-Axis Length',23x,f12.4,
     &           /,' b-Axis Length',23x,f12.4,
     &           /,' c-Axis Length',23x,f12.4,
     &           /,' Alpha Angle',25x,f12.4,
     &           /,' Beta Angle',26x,f12.4,
     &           /,' Gamma Angle',25x,f12.4,
     &           /,' Cell Volume',25x,f12.4,
     &           /,' Lattice Type',16x,a20)
         if (spacegrp .ne. '          ') then
            value = spacegrp
            call justify (value)
            write (iout,40)  value
   40       format (' Space Group',17x,a20)
         end if
      end if
c
c     info on force field potential energy function
c
      value = forcefield
      call justify (value)
      write (iout,50)  value
   50 format (/,' Force Field Name :',10x,a20)
c
c     details of vdw potential energy functional form
c
      if (use_vdw) then
         label(1) = vdwtyp
         label(2) = radtyp
         label(3) = radsiz
         label(4) = radrule
         label(5) = epsrule
         do i = 1, 5
            call justify (label(i))
         end do
         write (iout,60)  (label(i),i=1,5)
   60    format (/,' VDW Function',16x,a20,
     &           /,' Size Descriptor',13x,a20,
     &           /,' Size Unit Type',14x,a20,
     &           /,' Size Combining Rule',9x,a20,
     &           /,' Well Depth Rule',13x,a20)
         if (vdwcut .le. 1000.0d0) then
            write (iout,70)  vdwcut
   70       format (' VDW Cutoff',26x,f12.4)
         end if
      end if
c
c     details of electrostatic energy functional form
c
      if (use_charge .or. use_dipole .or. use_mpole .or. use_polar) then
         write (iout,80)
   80    format ()
      end if
      if (use_charge) then
         value = 'PARTIAL CHARGE'
         call justify (value)
         write (iout,90)  value
   90    format (' Electrostatics',14x,a20)
      end if
      if (use_dipole) then
         value = 'BOND DIPOLE'
         call justify (value)
         write (iout,100)  value
  100    format (' Electrostatics',14x,a20)
      end if
      if (use_mpole) then
         value = 'ATOMIC MULTIPOLE'
         call justify (value)
         write (iout,110)  value
  110    format (' Electrostatics',14x,a20)
      end if
      if (use_polar) then
         value = 'INDUCED DIPOLE'
         call justify (value)
         write (iout,120)  value
  120    format (' Electrostatics',14x,a20)
      end if
c
c     details of particle mesh Ewald calculation
c
      if (use_ewald) then
         value = boundary
         call justify (value)
         write (iout,130)  aewald,ewaldcut,nfft1,nfft2,nfft3,
     &                    bsorder,value
  130    format (/,' Particle Mesh Ewald :',
     &           //,' Ewald Coefficient',19x,f12.4,
     &           /,' Real-Space Cutoff',19x,f12.4,
     &           /,' Grid Dimensions',21x,3i4,
     &           /,' B-Spline Order',26x,i8,
     &           /,' Boundary Condition',10x,a20)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine paramyze  --  force field parameter analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "paramyze" prints the force field parameters used in the
c     computation of each of the potential energy terms
c
c
      subroutine paramyze (active)
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'angpot.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bitor.i'
      include 'bond.i'
      include 'charge.i'
      include 'dipole.i'
      include 'improp.i'
      include 'imptor.i'
      include 'iounit.i'
      include 'korbs.i'
      include 'ktrtor.i'
      include 'kvdws.i'
      include 'math.i'
      include 'mpole.i'
      include 'opbend.i'
      include 'opdist.i'
      include 'piorbs.i'
      include 'pistuf.i'
      include 'pitors.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'solute.i'
      include 'strbnd.i'
      include 'strtor.i'
      include 'tors.i'
      include 'tortor.i'
      include 'units.i'
      include 'urey.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k
      integer ia,ib,ic
      integer id,ie,ig
      integer ixaxe
      integer iyaxe
      integer izaxe
      integer fold(6)
      real*8 bla,blc
      real*8 radj,rad4j
      real*8 ampli(6)
      real*8 phase(6)
      real*8 mpl(13)
      logical header
      logical active(maxatm)
c
c
c     number of each type of interaction and site
c
      if (n .ne. 0) then
         write (iout,10)
   10    format (/,' Interactions and Sites :',/)
      end if
      if (nbond .ne. 0) then
         write (iout,20)  nbond
   20    format (' Bond Stretches',19x,i15)
      end if
      if (nangle .ne. 0) then
         write (iout,30)  nangle
   30    format (' Angle Bends',22x,i15)
      end if
      if (nstrbnd .ne. 0) then
         write (iout,40)  nstrbnd
   40    format (' Stretch-Bends',20x,i15)
      end if
      if (nurey .ne. 0) then
         write (iout,50)  nurey
   50    format (' Urey-Bradley',21x,i15)
      end if
      if (nangang .ne. 0) then
         write (iout,50)  nangang
   60    format (' Angle-Angles',21x,i15)
      end if
      if (nopbend .ne. 0) then
         write (iout,70)  nopbend
   70    format (' Out-of-Plane Bends',15x,i15)
      end if
      if (nopdist .ne. 0) then
         write (iout,80)  nopdist
   80    format (' Out-of-Plane Distances',11x,i15)
      end if
      if (niprop .ne. 0) then
         write (iout,90)  niprop
   90    format (' Improper Dihedrals',15x,i15)
      end if
      if (nitors .ne. 0) then
         write (iout,100)  nitors
  100    format (' Improper Torsions',16x,i15)
      end if
      if (ntors .ne. 0) then
         write (iout,110)  ntors
  110    format (' Torsional Angles',17x,i15)
      end if
      if (npitors .ne. 0) then
         write (iout,120)  npitors
  120    format (' Pi-Orbital Torsions',14x,i15)
      end if
      if (nstrtor .ne. 0) then
         write (iout,130)  nstrtor
  130    format (' Stretch-Torsions',17x,i15)
      end if
      if (ntortor .ne. 0) then
         write (iout,140)  ntortor
  140    format (' Torsion-Torsions',17x,i15)
      end if
      if (nvdw .ne. 0) then
         write (iout,150)  nvdw
  150    format (' Van der Waals Sites',14x,i15)
      end if
      if (nion .ne. 0) then
         write (iout,160)  nion
  160    format (' Atomic Partial Charges',11x,i15)
      end if
      if (ndipole .ne. 0) then
         write (iout,170)  ndipole
  170    format (' Bond Dipole Moments',14x,i15)
      end if
      if (npole .ne. 0) then
         write (iout,180)  npole
  180    format (' Atomic Multipoles',16x,i15)
      end if
      if (npolar .ne. 0) then
         write (iout,190)  npolar
  190    format (' Polarizable Sites',16x,i15)
      end if
      if (norbit .ne. 0) then
         write (iout,200)  norbit
  200    format (' Pisystem Atoms',19x,i15)
      end if
      if (nbpi .ne. 0) then
         write (iout,210)  nbpi
  210    format (' Conjugated Pi-Bonds',14x,i15)
      end if
c
c     parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,220)
  220          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',2x,'Symbol',2x,'Type',
     &                    2x,'Class',2x,'Atomic',3x,'Mass',
     &                    2x,'Valence',2x,'Description',/)
            end if
            write (iout,230)  i,name(i),type(i),class(i),atomic(i),
     &                        mass(i),valence(i),story(i)
  230       format (i6,5x,a3,2i7,i6,f10.3,i5,5x,a24)
         end if
      end do
c
c     parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,240)
  240             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',7x,'Size',
     &                       3x,'Epsilon',3x,'Size 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  radj = rad(j)
                  if (radsiz .eq. 'DIAMETER')  radj = 2.0d0 * radj
                  if (radtyp .eq. 'SIGMA')  radj = radj / twosix
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,250)  k,i,radj,eps(j)
  250                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,260)  k,i,radj,eps(j),reduct(j)
  260                format (i6,3x,i6,7x,2f10.4,22x,f10.4)
                  end if
               else
                  radj = rad(j)
                  rad4j = rad4(j)
                  if (radsiz .eq. 'DIAMETER') then
                     radj = 2.0d0 * radj
                     rad4j = 2.0d0 * rad4j
                  end if
                  if (radtyp .eq. 'SIGMA') then
                     radj = radj / twosix
                     rad4j = rad4j / twosix
                  end if
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,270)  k,i,radj,eps(j),rad4j,eps4(j)
  270                format (i6,3x,i6,7x,2f10.4,1x,2f10.4)
                  else
                     write (iout,280)  k,i,radj,eps(j),rad4j,
     &                                eps4(j),reduct(j)
  280                format (i6,3x,i6,7x,2f10.4,1x,2f10.4,1x,f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,290)
  290             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',7x,'Bond',/)
               end if
               write (iout,300)  i,ia,ib,bk(i),bl(i)
  300          format (i6,3x,2i6,19x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',3x,'Fold',4x,'Type',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,320)  i,ia,ib,ic,ak(i),anat(i)
  320             format (i6,3x,3i6,13x,2f10.3)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,330)  i,ia,ib,ic,ak(i),anat(i)
  330             format (i6,3x,3i6,13x,2f10.3,9x,'In-Plane')
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,340)  i,ia,ib,ic,ak(i),anat(i)
  340             format (i6,3x,3i6,13x,2f10.3,9x,'Linear')
               else if (angtyp(i) .eq. 'FOURIER ') then
                  write (iout,350)  i,ia,ib,ic,ak(i),anat(i),afld(i)
  350             format (i6,3x,3i6,13x,2f10.3,f7.1,2x,'Fourier')
               end if
            end if
         end do
      end if
c
c     parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,360)
  360             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',8x,'KSB 1',5x,'KSB 2',
     &                       6x,'Angle',3x,'Bond 1',3x,'Bond 2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               write (iout,370)  i,ia,ib,ic,sbk(1,i),sbk(2,i),
     &                           anat(k),bla,blc
  370          format (i6,3x,3i6,1x,2f10.3,2x,f9.3,2f9.4)
            end if
         end do
      end if
c
c     parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ib = iury(2,i)
            ic = iury(3,i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Urey-Bradley Parameters :',
     &                    //,13x,'Atom Numbers',21x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,390)  i,ia,ib,ic,uk(i),ul(i)
  390          format (i6,3x,3i6,13x,f10.3,f10.4)
            end if
         end do
      end if
c
c     parameters used for out-of-plane bend interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,400)
  400             format (/,' Out-of-Plane Bend Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               write (iout,410)  i,id,ib,ia,ic,opbk(i)
  410          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,420)
  420             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPD',/)
               end if
               write (iout,430)  i,ia,ib,ic,id,opdk(i)
  430          format (i6,3x,4i6,9x,f10.3)
            end if
         end do
      end if
c
c     parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,440)
  440             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,450)  i,ia,ib,ic,id,kprop(i),vprop(i)
  450          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,460)
  460             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,470)  i,ia,ib,ic,id
  470             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,480)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  480             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,490)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  490             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,500)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  500             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,510)
  510             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,520)  i,ia,ib,ic,id
  520             format (i6,3x,4i6)
               else
                  write (iout,530)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  530             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     parameters used for pi-orbital torsion interactions
c
      if (use_pitors) then
         header = .true.
         do i = 1, npitors
            ia = ipit(1,i)
            ib = ipit(2,i)
            ic = ipit(3,i)
            id = ipit(4,i)
            ie = ipit(5,i)
            ig = ipit(6,i)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &             active(id) .or. active(ie) .or. active(ig)) then
               if (header) then
                  header = .false.
                  write (iout,540)
  540             format (/,' Pi-Orbital Torsion Parameters :',
     &                    //,10x,'Atom Numbers',19x,'Amplitude',/)
               end if
               write (iout,550)  i,ic,id,kpit(i)
  550          format (i6,3x,2i6,19x,f10.4)
            end if
         end do
      end if
c
c     parameters used for stretch-torsion interactions
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,560)
  560             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               j = 0
               if (kst(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = kst(1,i)
                  phase(j) = tors1(2,k)
               end if
               if (kst(2,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = kst(2,i)
                  phase(j) = tors2(2,k)
               end if
               if (kst(3,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = kst(3,i)
                  phase(j) = tors3(2,k)
               end if
               write (iout,570)  i,ia,ib,ic,id,bl(ist(2,i)),
     &                           (ampli(k),nint(phase(k)),
     &                           fold(k),k=1,j)
  570          format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,'/',i1))
            end if
         end do
      end if
c
c     parameters used for torsion-torsion interactions
c
      if (use_tortor) then
         header = .true.
         do i = 1, ntortor
            k = itt(1,i)
            ia = ibitor(1,k)
            ib = ibitor(2,k)
            ic = ibitor(3,k)
            id = ibitor(4,k)
            ie = ibitor(5,k)
            if (active(ia) .or. active(ib) .or. active(ic) .or.
     &                active(id) .or. active(ie)) then
               if (header) then
                  header = .false.
                  write (iout,580)
  580             format (/,' Torsion-Torsion Parameters :',
     &                    //,20x,'Atom Numbers',18x,'Spline Grid',/)
               end if
               j = itt(2,i)
               write (iout,590)  i,ia,ib,ic,id,ie,tnx(j),tny(j)
  590          format (i6,3x,5i6,10x,2i6)
            end if
         end do
      end if
c
c     parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,600)
  600             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,610)  i,ia,pchg(i)
  610             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,620)  i,ia,pchg(i),ib,ic
  620             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,630)
  630             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,640)  i,ia,ib,bdpl(i),sdpl(i)
  640          format (i6,3x,2i6,19x,f10.4,f10.3)
            end if
         end do
      end if
c
c     parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,650)
  650             format (/,' Atomic Multipole Parameters :',
     &                    //,11x,'Atom',3x,'Z-Axis',1x,'X-Axis',
     &                       1x,'Y-Axis',2x,
     &                       'Frame',11x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (izaxe .gt. n)  izaxe = 0
               if (ixaxe .gt. n)  ixaxe = 0
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               mpl(1) = pole(1,i)
               do j = 2, 4
                  mpl(j) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  mpl(j) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (ixaxe .eq. 0) then
                  write (iout,660)  i,ia,izaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  660             format (i6,3x,i6,1x,i7,17x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else  if (iyaxe .eq. 0) then
                  write (iout,670)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  670             format (i6,3x,i6,1x,2i7,10x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               else
                  write (iout,680)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (mpl(j),j=1,5),mpl(8),mpl(9),
     &                              (mpl(j),j=11,13)
  680             format (i6,3x,i6,1x,3i7,3x,a8,2x,f9.5,/,50x,3f9.5,
     &                    /,50x,f9.5,/,50x,2f9.5,/,50x,3f9.5)
               end if
            end if
         end do
      end if
c
c     parameters used for dipole polarizability
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,690)
  690             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',5x,'Alpha',5x,'Damp',
     &                       6x,'Polarization Group',/)
               end if
               write (iout,700)  i,ia,polarity(i),thole(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  700          format (i6,3x,i6,6x,f10.4,f9.3,3x,20i6)
            end if
         end do
      end if
c
c     parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,710)
  710             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,720)  k,i,rsolv(i),asolv(i)
  720          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
c
c     parameters used for conjugated pisystem atoms
c
      if (use_orbit) then
         header = .true.
         do i = 1, norbit
            ia = iorbit(i)
            j = class(ia)
            if (header) then
               header = .false.
               write (iout,730)
  730          format (/,' Conjugated Pi-Atom Parameters :',
     &                 //,10x,'Atom Number',14x,'Nelect',
     &                    6x,'Ionize',4x,'Repulsion',/)
            end if
            write (iout,740)  i,ia,electron(j),ionize(j),repulse(j)
  740       format (i6,3x,i6,17x,f8.1,3x,f10.4,2x,f10.4)
         end do
      end if
c
c     parameters used for conjugated pibond interactions
c
      if (use_orbit) then
         header = .true.
         do i = 1, nbpi
            ia = ibpi(2,i)
            ib = ibpi(3,i)
            if (header) then
               header = .false.
               write (iout,750)
  750          format (/,' Conjugated Pi-Bond Parameters :',
     &                 //,10x,'Atom Numbers',21x,'K Slope',
     &                    3x,'L Slope',/)
            end if
            write (iout,760)  i,ia,ib,kslope(i),lslope(i)
  760       format (i6,3x,2i6,19x,2f10.4)
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enrgyze  --  compute & report energy analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "enrgyze" is an auxiliary routine for the analyze program
c     that performs the energy analysis and prints the total and
c     intermolecular energies
c
c
      subroutine enrgyze
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      real*8 energy
c
c
c     perform the energy analysis by atom and component
c
      call analysis (energy)
c
c     print out the total potential energy of the system
c
      if (digits .ge. 8) then
         if (abs(energy) .lt. 1.0d10) then
            write (iout,10)  energy
   10       format (/,' Total Potential Energy :',8x,f20.8,
     &                    ' Kcal/mole')
         else
            write (iout,20)  energy
   20       format (/,' Total Potential Energy :',8x,d20.8,
     &                    ' Kcal/mole')
         end if
      else if (digits .ge. 6) then
         if (abs(energy) .lt. 1.0d10) then
            write (iout,30)  energy
   30       format (/,' Total Potential Energy :',8x,f18.6,
     &                    ' Kcal/mole')
         else
            write (iout,40)  energy
   40       format (/,' Total Potential Energy :',8x,d18.6,
     &                    ' Kcal/mole')
         end if
      else
         if (abs(energy) .lt. 1.0d10) then
            write (iout,50)  energy
   50       format (/,' Total Potential Energy :',8x,f16.4,
     &                    ' Kcal/mole')
         else
            write (iout,60)  energy
   60       format (/,' Total Potential Energy :',8x,d16.4,
     &                    ' Kcal/mole')
         end if
      end if
c
c     intermolecular energy for systems with multiple molecules
c
      if (nmol.gt.1 .and. nmol.lt.n .and. .not.use_ewald) then
         if (digits .ge. 8) then
            if (abs(einter) .lt. 1.0d10) then
               write (iout,70)  einter
   70          format (/,' Intermolecular Energy :',9x,f20.8,
     &                       ' Kcal/mole')
            else
               write (iout,80)  einter
   80          format (/,' Intermolecular Energy :',9x,d20.8,
     &                       ' Kcal/mole')
            end if
         else if (digits .ge. 6) then
            if (abs(einter) .lt. 1.0d10) then
               write (iout,90)  einter
   90          format (/,' Intermolecular Energy :',9x,f18.6,
     &                       ' Kcal/mole')
            else
               write (iout,100)  einter
  100          format (/,' Intermolecular Energy :',9x,d18.6,
     &                       ' Kcal/mole')
            end if
         else
            if (abs(einter) .lt. 1.0d10) then
               write (iout,110)  einter
  110          format (/,' Intermolecular Energy :',9x,f16.4,
     &                       ' Kcal/mole')
            else
               write (iout,120)  einter
  120          format (/,' Intermolecular Energy :',9x,d16.4,
     &                       ' Kcal/mole')
            end if
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analyz4  --  low precision energy components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analyz4" prints the energy to 4 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz4
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',11x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',17x,f16.4,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',19x,f16.4,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',20x,f16.4,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',20x,f16.4,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',21x,f16.4,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',15x,f16.4,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',11x,f16.4,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',15x,f16.4,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',16x,f16.4,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',17x,f16.4,i15)
      end if
      if (use_pitors .and. nept.ne.0) then
         write (iout,120)  ept,nept
  120    format (' Pi-Orbital Torsion',14x,f16.4,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,130)  ebt,nebt
  130    format (' Stretch-Torsion',17x,f16.4,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,140)  ett,nett
  140    format (' Torsion-Torsion',17x,f16.4,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,150)  ev,nev
  150       format (' Van der Waals',19x,f16.4,i15)
         else
            write (iout,160)  ev,nev
  160       format (' Van der Waals',19x,d16.4,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         if (abs(ec) .lt. 1.0d10) then
            write (iout,170)  ec,nec
  170       format (' Charge-Charge',19x,f16.4,i15)
         else
            write (iout,180)  ec,nec
  180       format (' Charge-Charge',19x,d16.4,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',19x,f16.4,i15)
         else
            write (iout,200)  ecd,necd
  200       format (' Charge-Dipole',19x,d16.4,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',19x,f16.4,i15)
         else
            write (iout,220)  ed,ned
  220       format (' Dipole-Dipole',19x,d16.4,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         if (abs(em) .lt. 1.0d10) then
            write (iout,230)  em,nem
  230       format (' Atomic Multipoles',15x,f16.4,i15)
         else
            write (iout,240)  em,nem
  240       format (' Atomic Multipoles',15x,d16.4,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         if (abs(ep) .lt. 1.0d10) then
            write (iout,250)  ep,nep
  250       format (' Polarization',20x,f16.4,i15)
         else
            write (iout,260)  ep,nep
  260       format (' Polarization',20x,d16.4,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,270)  er,ner
  270    format (' Reaction Field',18x,f16.4,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,280)  es,nes
  280    format (' Continuum Solvation',13x,f16.4,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,290)  elf,nelf
  290    format (' Metal Ligand Field',14x,f16.4,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,300)  eg,neg
  300    format (' Geometric Restraints',12x,f16.4,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,310)  ex,nex
  310    format (' Extra Energy Terms',14x,f16.4,i15)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine analyz6  --  medium precision energy components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "analyz6" prints the energy to 6 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz6
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',11x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',15x,f18.6,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',17x,f18.6,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',18x,f18.6,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',18x,f18.6,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',19x,f18.6,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',13x,f18.6,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',9x,f18.6,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',13x,f18.6,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',14x,f18.6,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',15x,f18.6,i15)
      end if
      if (use_pitors .and. nept.ne.0) then
         write (iout,120)  ept,nept
  120    format (' Pi-Orbital Torsion',12x,f18.6,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,130)  ebt,nebt
  130    format (' Stretch-Torsion',15x,f18.6,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,140)  ett,nett
  140    format (' Torsion-Torsion',15x,f18.6,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,150)  ev,nev
  150       format (' Van der Waals',17x,f18.6,i15)
         else
            write (iout,160)  ev,nev
  160       format (' Van der Waals',17x,d18.6,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         if (abs(ec) .lt. 1.0d10) then
            write (iout,170)  ec,nec
  170       format (' Charge-Charge',17x,f18.6,i15)
         else
            write (iout,180)  ec,nec
  180       format (' Charge-Charge',17x,d18.6,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',17x,f18.6,i15)
         else
            write (iout,200)  ecd,necd
  200       format (' Charge-Dipole',17x,d18.6,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',17x,f18.6,i15)
         else
            write (iout,220)  ed,ned
  220       format (' Dipole-Dipole',17x,d18.6,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         if (abs(em) .lt. 1.0d10) then
            write (iout,230)  em,nem
  230       format (' Atomic Multipoles',13x,f18.6,i15)
         else
            write (iout,240)  em,nem
  240       format (' Atomic Multipoles',13x,d18.6,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         if (abs(ep) .lt. 1.0d10) then
            write (iout,250)  ep,nep
  250       format (' Polarization',18x,f18.6,i15)
         else
            write (iout,260)  ep,nep
  260       format (' Polarization',18x,d18.6,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,270)  er,ner
  270    format (' Reaction Field',16x,f18.6,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,280)  es,nes
  280    format (' Continuum Solvation',11x,f18.6,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,290)  elf,nelf
  290    format (' Metal Ligand Field',12x,f18.6,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,300)  eg,neg
  300    format (' Geometric Restraints',10x,f18.6,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,310)  ex,nex
  310    format (' Extra Energy Terms',12x,f18.6,i15)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine analyz8  --  high precision energy components  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "analyz8" prints the energy to 8 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz8
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',11x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',13x,f20.8,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',15x,f20.8,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',16x,f20.8,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',16x,f20.8,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',17x,f20.8,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',11x,f20.8,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',7x,f20.8,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',11x,f20.8,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',12x,f20.8,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',13x,f20.8,i15)
      end if
      if (use_pitors .and. nept.ne.0) then
         write (iout,120)  ept,nept
  120    format (' Pi-Orbital Torsion',10x,f20.8,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,130)  ebt,nebt
  130    format (' Stretch-Torsion',13x,f20.8,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,140)  ett,nett
  140    format (' Torsion-Torsion',13x,f20.8,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,150)  ev,nev
  150       format (' Van der Waals',15x,f20.8,i15)
         else
            write (iout,160)  ev,nev
  160       format (' Van der Waals',15x,d20.8,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         if (abs(ec) .lt. 1.0d10) then
            write (iout,170)  label,ec,nec
  170       format (' Charge-Charge',15x,f20.8,i15)
         else
            write (iout,180)  label,ec,nec
  180       format (' Charge-Charge',15x,d20.8,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',15x,f20.8,i15)
         else
            write (iout,200)  ecd,necd
  200       format (' Charge-Dipole',15x,d20.8,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',15x,f20.8,i15)
         else
            write (iout,220)  ed,ned
  220       format (' Dipole-Dipole',15x,d20.8,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         if (abs(em) .lt. 1.0d10) then
            write (iout,230)  em,nem
  230       format (' Atomic Multipoles',11x,f20.8,i15)
         else
            write (iout,240)  em,nem
  240       format (' Atomic Multipoles',11x,d20.8,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         if (abs(ep) .lt. 1.0d10) then
            write (iout,250)  ep,nep
  250       format (' Polarization',16x,f20.8,i15)
         else
            write (iout,260)  ep,nep
  260       format (' Polarization',16x,d20.8,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,270)  er,ner
  270    format (' Reaction Field',14x,f20.8,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,280)  es,nes
  280    format (' Continuum Solvation',9x,f20.8,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,290)  elf,nelf
  290    format (' Metal Ligand Field',10x,f20.8,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,300)  eg,neg
  300    format (' Geometric Restraints',8x,f20.8,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,310)  ex,nex
  310    format (' Extra Energy Terms',10x,f20.8,i15)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine propyze  --  electrostatic & inertial analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "propyze" finds and prints the total charge, dipole moment
c     components, radius of gyration and moments of inertia
c
c
      subroutine propyze
      include 'sizes.i'
      include 'chgpot.i'
      include 'iounit.i'
      include 'moment.i'
      include 'virial.i'
      real*8 rg,energy
      real*8 derivs(3,maxatm)
c
c
c     get the total charge, dipole and quadrupole moments
c
      call moments
      write (iout,10)  netchg
   10 format (/,' Total Electric Charge :',13x,f12.5,' Electrons')
      write (iout,20)  netdpl,xdpl,ydpl,zdpl
   20 format (/,' Dipole Moment Magnitude :',11x,f12.3,' Debyes',
     &        //,' Dipole X,Y,Z-Components :',11x,3f12.3)
      write (iout,30)  xxqdp,xyqdp,xzqdp,yxqdp,yyqdp,
     &                 yzqdp,zxqdp,zyqdp,zzqdp
   30 format (/,' Quadrupole Moment Tensor :',10x,3f12.3,
     &        /,6x,'(Buckinghams)',18x,3f12.3,
     &        /,37x,3f12.3)
      write (iout,40)  netqdp(1),netqdp(2),netqdp(3)
   40 format (/,' Principal Axes Quadrupole :',9x,3f12.3)
      if (dielec .ne. 1.0d0) then
         write (iout,50)  dielec
   50    format (/,' Dielectric Constant :',15x,f12.3)
         write (iout,60)  netchg/sqrt(dielec)
   60    format (' Effective Total Charge :',12x,f12.5,' Electrons')
         write (iout,70)  netdpl/sqrt(dielec)
   70    format (' Effective Dipole Moment :',11x,f12.3,' Debyes')
      end if
c
c     get the radius of gyration and moments of inertia
c
      call gyrate (rg)
      write (iout,80)  rg
   80 format (/,' Radius of Gyration :',16x,f12.3,' Angstroms')
      call inertia (1)
c
c     get the internal virial tensor via gradient calculation
c
      call gradient (energy,derivs)
      write (iout,90)  (vir(1,i),vir(2,i),vir(3,i),i=1,3)
   90 format (/,' Internal Virial Tensor :',12x,3f12.3,
     &        /,37x,3f12.3,/,37x,3f12.3)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine atomyze  --  individual atom energy analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "atomyze" prints the potential energy components broken
c     down by atom and to a choice of precision
c
c
      subroutine atomyze (active)
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      integer i
      logical active(maxatm)
c
c
c     energy partitioning over the individual atoms
c
      write (iout,10)
   10 format (/,' Potential Energy Breakdown over Atoms :')
      if (digits .ge. 8) then
         write (iout,20)
   20    format (/,'  Atom',9x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &           /,15x,'EAA',13x,'EOPB',12x,'EOPD',12x,'EID',
     &           /,15x,'EIT',13x,'ET',14x,'EPT',13x,'EBT',
     &           /,15x,'ETT',13x,'EV',14x,'EC',14x,'ECD',
     &           /,15x,'ED',14x,'EM',14x,'EP',14x,'ER',
     &           /,15x,'ES',14x,'ELF',13x,'EG',14x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,30)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aept(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   30          format (/,i6,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
     &                 /,6x,4f16.8,/,6x,4f16.8,/,6x,4f16.8)
            end if
         end do
      else if (digits .ge. 6) then
         write (iout,40)
   40    format (/,'  Atom',8x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &              11x,'EAA',
     &           /,14x,'EOPB',10x,'EOPD',10x,'EID',11x,'EIT',
     &              11x,'ET',
     &           /,14x,'EPT',11x,'EBT',11x,'ETT',11x,'EV',12x,'EC',
     &           /,14x,'ECD',11x,'ED',12x,'EM',12x,'EP',12x,'ER',
     &           /,14x,'ES',12x,'ELF',11x,'EG',12x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,50)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aept(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   50          format (/,i6,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
     &                 /,6x,5f14.6,/,6x,4f14.6)
            end if
         end do
      else
         write (iout,60)
   60    format (/,'  Atom',8x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &              9x,'EAA',9x,'EOPB',
     &           /,14x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',10x,'EPT',
     &              9x,'EBT',
     &           /,14x,'ETT',9x,'EV',10x,'EC',10x,'ECD',9x,'ED',
     &              10x,'EM',
     &           /,14x,'EP',10x,'ER',10x,'ES',10x,'ELF',9x,'EG',
     &              10x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,70)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aept(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   70          format (/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,
     &                 /,6x,6f12.4)
            end if
         end do
      end if
      return
      end
