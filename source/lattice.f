c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine lattice  --  setup periodic boundary conditions  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "lattice" stores the periodic box dimensions and sets angle
c     values to be used in computing fractional coordinates
c
c
      subroutine lattice
      use bound
      use boxes
      use cell
      use inform
      use iounit
      use math
      implicit none
      real*8 boxmax
      real*8 ar1,ar2,ar3
      real*8 br1,br2,br3
      real*8 cr1,cr2,cr3
c
c
c     use periodic boundary conditions if a cell was defined
c
      boxmax = max(xbox,ybox,zbox)
      if (boxmax .ne. 0.0d0)  use_bounds = .true.
c
c     set unspecified periodic boundary box lengths and angles
c
      if (use_bounds) then
         if (xbox .eq. 0.0d0)  xbox = boxmax
         if (ybox .eq. 0.0d0)  ybox = boxmax
         if (zbox .eq. 0.0d0)  zbox = boxmax
         if (alpha .eq. 0.0d0)  alpha = 90.0d0
         if (beta .eq. 0.0d0)  beta = 90.0d0
         if (gamma .eq. 0.0d0)  gamma = 90.0d0
c
c     determine the general periodic boundary lattice type
c
         orthogonal = .false.
         monoclinic = .false.
         triclinic = .false.
         if (nosymm) then
            triclinic = .true.
         else if (alpha.eq.90.0d0 .and. beta.eq.90.0d0
     &               .and. gamma.eq.90.0d0) then
            orthogonal = .true.
         else if (alpha.eq.90.0d0 .and. gamma.eq.90.0d0) then
            monoclinic = .true.
         else
            triclinic = .true.
         end if
      end if
c
c     check for proper use of non-prism periodic boundaries
c
      if (octahedron .or. dodecadron) then
         orthogonal = .false.
         monoclinic = .false.
         triclinic = .false.
         nonprism = .true.
         ybox = xbox
         if (octahedron) then
            zbox = xbox
         else if (dodecadron) then
            zbox = xbox * root2
         end if         
         alpha = 90.0d0
         beta = 90.0d0
         gamma = 90.0d0
      end if
c
c     compute and store half box lengths and other lengths
c
      xbox2 = 0.5d0 * xbox
      ybox2 = 0.5d0 * ybox
      zbox2 = 0.5d0 * zbox
      if (octahedron)  box34 = 0.75d0 * xbox
c
c     set replicated cell dimensions equal to the unit cell
c
      xcell = xbox
      ycell = ybox
      zcell = zbox
      xcell2 = xbox2
      ycell2 = ybox2
      zcell2 = zbox2
c
c     get values needed for fractional coordinate computations
c
      if (triclinic) then
         alpha_sin = sin(alpha/radian)
         alpha_cos = cos(alpha/radian)
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = sin(gamma/radian)
         gamma_cos = cos(gamma/radian)
         beta_term = (alpha_cos - beta_cos*gamma_cos) / gamma_sin
         gamma_term = sqrt(beta_sin**2 - beta_term**2)
      else if (monoclinic) then
         alpha_sin = 1.0d0
         alpha_cos = 0.0d0
         beta_sin = sin(beta/radian)
         beta_cos = cos(beta/radian)
         gamma_sin = 1.0d0
         gamma_cos = 0.0d0
         beta_term = 0.0d0
         gamma_term = beta_sin
      else
         alpha_sin = 1.0d0
         alpha_cos = 0.0d0
         beta_sin = 1.0d0
         beta_cos = 0.0d0
         gamma_sin = 1.0d0
         gamma_cos = 0.0d0
         beta_term = 0.0d0
         gamma_term = 1.0d0
      end if
c
c     determine the volume of the parent periodic box
c
      volbox = 0.0d0
      if (triclinic) then
         volbox = (gamma_sin*gamma_term) * xbox * ybox * zbox
      else if (monoclinic) then
         volbox = beta_sin * xbox * ybox * zbox
      else
         volbox = xbox * ybox * zbox
      end if
c
c     compute and store real space lattice vectors as rows
c
      ar1 = xbox
      ar2 = 0.0d0
      ar3 = 0.0d0
      br1 = ybox * gamma_cos
      br2 = ybox * gamma_sin
      br3 = 0.0d0
      cr1 = zbox * beta_cos
      cr2 = zbox * beta_term
      cr3 = zbox * gamma_term
      lvec(1,1) = ar1
      lvec(1,2) = ar2
      lvec(1,3) = ar3
      lvec(2,1) = br1
      lvec(2,2) = br2
      lvec(2,3) = br3
      lvec(3,1) = cr1
      lvec(3,2) = cr2
      lvec(3,3) = cr3
c
c     compute and store reciprocal lattice vectors as columns
c
      if (volbox .ne. 0.0d0) then
         recip(1,1) = (br2*cr3 - cr2*br3) / volbox
         recip(2,1) = (br3*cr1 - cr3*br1) / volbox
         recip(3,1) = (br1*cr2 - cr1*br2) / volbox
         recip(1,2) = (cr2*ar3 - ar2*cr3) / volbox
         recip(2,2) = (cr3*ar1 - ar3*cr1) / volbox
         recip(3,2) = (cr1*ar2 - ar1*cr2) / volbox
         recip(1,3) = (ar2*br3 - br2*ar3) / volbox
         recip(2,3) = (ar3*br1 - br3*ar1) / volbox
         recip(3,3) = (ar1*br2 - br1*ar2) / volbox
      end if
c
c     correct volume of non-parallelepiped periodic cells
c
      if (nonprism)  volbox = 0.5d0 * volbox
      return
      end
