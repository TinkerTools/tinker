c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  module polpot  --  polarization functional form details  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     politer      maximum number of induced dipole SCF iterations
c     poleps       induced dipole convergence criterion (rms Debye/atom)
c     p2scale      scale factor for 1-2 polarization energy interactions
c     p3scale      scale factor for 1-3 polarization energy interactions
c     p4scale      scale factor for 1-4 polarization energy interactions
c     p5scale      scale factor for 1-5 polarization energy interactions
c     p2iscale     scale factor for 1-2 intragroup polarization energy
c     p3iscale     scale factor for 1-3 intragroup polarization energy
c     p4iscale     scale factor for 1-4 intragroup polarization energy
c     p5iscale     scale factor for 1-5 intragroup polarization energy
c     d1scale      scale factor for intra-group direct induction
c     d2scale      scale factor for 1-2 group direct induction
c     d3scale      scale factor for 1-3 group direct induction
c     d4scale      scale factor for 1-4 group direct induction
c     u1scale      scale factor for intra-group mutual induction
c     u2scale      scale factor for 1-2 group mutual induction
c     u3scale      scale factor for 1-3 group mutual induction
c     u4scale      scale factor for 1-4 group mutual induction
c     w2scale      scale factor for 1-2 induced dipole interactions
c     w3scale      scale factor for 1-3 induced dipole interactions
c     w4scale      scale factor for 1-4 induced dipole interactions
c     w5scale      scale factor for 1-5 induced dipole interactions
c     uaccel       acceleration factor for induced dipole SCF iterations
c     polprt       flag to print summary of induced dipole iterations
c     dpequal      flag to set dscale values equal to pscale values
c     use_thole    flag to use Thole damped polarization interactions
c     use_tholed   flag to use alternate Thole for direct polarization
c     use_expol    flag to use damped exchange polarization correction
c     scrtyp       type of exchange polarization (S2U, S2 or G)
c     poltyp       type of polarization (MUTUAL, DIRECT, OPT or TCG)
c
c
      module polpot
      implicit none
      integer politer
      real*8 poleps
      real*8 p2scale,p3scale
      real*8 p4scale,p5scale
      real*8 p2iscale,p3iscale
      real*8 p4iscale,p5iscale
      real*8 d1scale,d2scale
      real*8 d3scale,d4scale
      real*8 u1scale,u2scale
      real*8 u3scale,u4scale
      real*8 w2scale,w3scale
      real*8 w4scale,w5scale
      real*8 uaccel
      logical polprt
      logical dpequal
      logical use_thole
      logical use_tholed
      logical use_expol
      character*3 scrtyp
      character*6 poltyp
      save
      end
