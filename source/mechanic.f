c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mechanic  --  initialize molecular mechanics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mechanic" sets up needed parameters for the potential energy
c     calculation and reads in many of the user selectable options
c
c
      subroutine mechanic
      use inform
      use iounit
      implicit none
c
c
c     set the bonded connectivity lists and active atoms
c
      call attach
      call active
c
c     find bonds, angles, torsions, bitorsions and small rings
c
      call bonds
      call angles
      call torsions
      call bitors
      call rings
c
c     get the base force field from parameter file and keyfile
c
      call field
c
c     find unit cell type, lattice parameters and cutoff values
c
      call unitcell
      call lattice
      call polymer
      call cutoffs
c
c     setup needed for potential energy smoothing methods
c
      call flatten
c
c     assign atom types, classes and other atomic information
c
      call katom
c
c     assign atoms to molecules and set the atom groups
c
      call molecule
      call cluster
c
c     find any pisystem atoms, bonds and torsional angles
c
      call orbital
c
c     assign bond, angle and cross term potential parameters
c
      call kbond
      call kangle
      call kstrbnd
      call kurey
      call kangang
c
c     assign out-of-plane deformation potential parameters
c
      call kopbend
      call kopdist
      call kimprop
      call kimptor
c
c     assign torsion and torsion cross term potential parameters
c
      call ktors
      call kpitors
      call kstrtor
      call kangtor
      call ktortor
c
c     assign electrostatic interaction potential parameters
c
      call kcharge
      call kdipole
      call kmpole
      call kpolar
      call kchgtrn
      call kchgflx
c
c     assign van der Waals, repulsion and dispersion parameters
c
      call kvdw
      call krepel
      call kxrepel
      call kdisp
c
c     assign solvation, metal, pisystem and restraint parameters
c
      call ksolv
      call kmetal
      call korbit
      call kgeom
      call kextra
c
c     assign electrostatic and dispersion Ewald sum parameters
c
      call kewald
c
c     set any holonomic interatomic distance constraints
c
      call shakeup
c
c     set hybrid parameter values for free energy perturbation
c
      call mutate
c
c     quit if essential parameter information is missing
c
      if (abort) then
         write (iout,10)
   10    format (/,' MECHANIC  --  Some Required Potential Energy',
     &              ' Parameters are Undefined')
         call fatal
      end if
      return
      end
