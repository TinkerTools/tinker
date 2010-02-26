c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine initial  --  initial values and program setup  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "initial" sets up original values for some parameters
c     and variables that might not otherwise get initialized
c
c
      subroutine initial
      implicit none
      include 'sizes.i'
      include 'align.i'
      include 'argue.i'
      include 'atoms.i'
      include 'bath.i'
      include 'bound.i'
      include 'cell.i'
      include 'files.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'linmin.i'
      include 'minima.i'
      include 'molcul.i'
      include 'mutant.i'
      include 'neigh.i'
      include 'output.i'
      include 'params.i'
      include 'pdb.i'
      include 'precis.i'
      include 'rigid.i'
      include 'scales.i'
      include 'sequen.i'
      include 'socket.i'
      include 'warp.i'
      include 'zclose.i'
      real*8 precise
c
c
c     number of atoms used in superposition
c
      nfit = 0
c
c     number of command line arguments
c
      narg = 0
c
c     number of atoms in the system
c
      n = 0
c
c     flags for temperature and pressure baths
c
      isothermal = .false.
      isobaric = .false.
c
c     flags for periodic boundaries
c
      use_bounds = .false.
      use_replica = .false.
      use_polymer = .false.
c
c     number of unit cell replicates
c
      ncell = 0
c
c     flag for use of atom groups
c
      use_group = .false.
c
c     highest numbered previous cycle file
c
      nprior = 0
c
c     information levels within the program
c
      verbose = .false.
      debug = .false.
      abort = .false.
c
c     default input/output unit numbers
c
      input = 5
      iout = 6
c
c     number of lines in the keyfile
c
      nkey = 0
c
c     default parameters used by line search
c
      stpmin = 0.0d0
      stpmax = 0.0d0
      cappa = 0.0d0
      slpmax = 0.0d0
      angmax = 0.0d0
      intmax = 0
c
c     default parameters used by optimizations
c
      fctmin = 0.0d0
      maxiter = 0
      nextiter = 0
      iprint = -1
      iwrite = -1
c
c     number of molecules in the system
c
      nmol = 0
c
c     number of mutated atoms in the system
c
      nmut = 0
c
c     flags for rebuilding of neighbor lists
c
      dovlst = .true.
      doclst = .true.
      domlst = .true.
c
c     type of coordinates file
c
      coordtype = 'NONE'
c
c     number of lines in the parameter file
c
      nprm = 0
c
c     number of atoms in Protein Data Bank format
c
      npdb = 0
c
c     flag for use of rigid bodies
c
      use_rigid = .false.
c
c     flag to show setting of optimization scale factors
c
      set_scale = .false.
c
c     number of residues and chains in biopolymer sequence
c
      nseq = 0
      nchain = 0
c
c     flags for external Java socket communication
c
      skt_init = .false.
      use_socket = .false.
c
c     flags for potential energy smoothing
c
      use_smooth = .false.
      use_dem = .false.
      use_gda = .false.
      use_tophat = .false.
      use_stophat = .false.
c
c     number of bonds added or deleted from Z-matrix
c
      nadd = 0
      ndel = 0
c
c     display program info and copyright notice
c
      call promo
c
c     names of biopolymer residue types
c
      call initres
c
c     determine a set of machine precision values
c
      tiny = precise (1)
      small = precise (2)
      huge = precise (3)
c
c     get any command line arguments to the program
c
      call command
      return
      end
