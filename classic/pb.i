c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2006 by Michael Schnieders & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  pb.i  --  parameters for Poisson-Boltzmann solvation  ##
c     ##                                                        ##
c     ############################################################
c
c
c     pbe      Poisson-Boltzman permanent multipole solvation energy
c     apbe     Poisson-Boltzman permanent multipole energy over atoms
c     pbr      Poisson-Boltzman cavity radii for atom types
c     pbep     Poisson-Boltzman energies on permanent multipoles
c     pbfp     Poisson-Boltzman forces on permanent multipoles
c     pbtp     Poisson-Boltzman torques on permanent multipoles
c     pbeuind  Poisson-Boltzman field due to induced dipoles
c     pbeuinp  Poisson-Boltzman field due to non-local induced dipoles
c
c     APBS configuration parameters (see APBS documentation for more details)
c     In the column on the right are possible values for each variable, with
c     the default values given in brackets. Note that only a subset of APBS
c     options are supported and/or are appropriate for use with AMOEBA.
c
c     pbtyp                                   lpbe
c
c     At some point AMOEBA with the non-linear PBE could be supported, but
c     there is only have theory for energies (no gradients).
c
c     pbsoln                                  mg-auto, [mg-manual]
c
c     Currently there is only limited support for focusing calculations,
c     which is a powerful feature of APBS. The current requirement is
c     that energies and forces must all be calculated using the finest
c     solution.
c
c     bcfl     boundary conditions            zero, sdh, [mdh]
c     chgm     multipole discretization       spl4
c
c     other charge discretization methods are not appropriate for AMOEBA
c
c     srfm     surface method                 mol, smol, [spl4]
c
c     spl4 is required for forces calculations, although mol is useful for
c     comparison with generalized Kirkwood
c
c     dime     number of grid points          [65, 65, 65]
c     grid     grid spacing (mg-manual)       fxn of "dime"
c     cgrid    coarse grid spacing            fxn of "dime"
c     fgrid    fine grid spacing              cgrid / 2
c
c     stable results require grid spacing to be fine enough to keep
c     multipoles inside the dielectric boundary (2.5 * grid < PBR)
c
c     gcent    grid center  (mg-manual)       center of mass
c     cgcent   coarse grid center             center of mass
c     fgcent   fine grid center               center of mass
c     pdie     solute/homogeneous dieletric   [1.0]
c     sdie     solvent dieletric              [78.3]
c     ionn     number of ion species          [0]
c     ionc     ion concentration (M)          [0.0]
c     ionq     ion charge (electrons)         [1.0]
c     ionr     ion radius (A)                 [2.0]
c     srad     solvent probe radius (A)       [1.4]
c     swin     surface spline window width    [0.3]
c     sdens    density of surface points      [10.0]
c
c     additional parameter to facilitate default grid setup
c
c     smin     minimum distance between an    [10.0]
c              atomic center and the grid
c              boundary (A)
c
c
      integer maxion
      parameter (maxion=10)
      integer ionn,dime,ionq
      real*8 pbe,apbe,pbr
      real*8 pbep,pbfp,pbtp
      real*8 pbeuind,pbeuinp
      real*8 grid,gcent
      real*8 cgrid,cgcent
      real*8 fgrid,fgcent
      real*8 ionr,ionc
      real*8 pdie,sdie
      real*8 srad,swin
      real*8 sdens,smin
      character*20 pbtyp,pbsoln
      character*20 bcfl,srfm,chgm
      common /pb/ pbe,apbe(maxatm),pbr(maxatm),pbep(3,maxatm),
     &            pbfp(3,maxatm),pbtp(3,maxatm),pbeuind(3,maxatm),
     &            pbeuinp(3,maxatm),grid(3),gcent(3),cgrid(3),cgcent(3),
     &            fgrid(3),fgcent(3),ionr(maxion),ionc(maxion),pdie,
     &            sdie,srad,swin,sdens,smin,ionn,dime(3),ionq(maxion),
     &            pbtyp,pbsoln,bcfl,srfm,chgm
