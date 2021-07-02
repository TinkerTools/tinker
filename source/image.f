c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine image  --  compute the minimum image distance  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "image" takes the components of pairwise distance between
c     two points in a periodic box and converts to the components
c     of the minimum image distance
c
c     literature reference:
c
c     U. K. Deiters, "Efficient Coding of the Minimum Image Convention",
c     Zeitschrift fur Physikalische Chemie, 227, 345-352 (2013)
c
c     note the "do while" clause below can be written using the "nint"
c     intrinsic, and the two forms give equivalent values:
c
c     do while (abs(xr) .gt. xbox2)
c        xr = xr - sign(xbox,xr)    vs.  xr = xr - xbox*nint(xr/xbox)
c     end do
c
c     which is faster depends on the specific machine and compiler
c     combination, and other implementations are also possible
c
c
      subroutine image (xr,yr,zr)
      use boxes
      use cell
      use math
      implicit none
      real*8 xr,yr,zr
      real*8 corr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
c
c     for monoclinic lattice, convert x and z to fractional,
c     find desired image, then translate back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert to fractional coordinates,
c     find image, then translate fractional back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
         corr = box23 * int(abs(xr/xbox)+abs(yr/ybox)+abs(zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)
c
c     for rhombic dodecahedron, align along the x- and y-axes
c
      else if (dodecadron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - root2*zbox*nint(zr/(zbox*root2))
         corr = xbox2 * int(abs(xr/xbox)+abs(yr/ybox)
     &                        +abs(root2*zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)*root2
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine imager  --  replicate minimum image distance  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "imager" takes the components of pairwise distance between
c     two points in the same or neighboring periodic boxes and
c     converts to the components of the minimum image distance
c
c
      subroutine imager (xr,yr,zr,i)
      use boxes
      use cell
      use math
      implicit none
      integer i
      real*8 xr,yr,zr
      real*8 xmove,ymove,zmove
      real*8 corr
c
c
c     set the distance to translate along each cell axis
c
      xmove = icell(1,i) * xbox
      ymove = icell(2,i) * ybox
      zmove = icell(3,i) * zbox
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr + xmove
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
c
c     for monoclinic lattice, convert x and z to fractional,
c     find desired image, then translate back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         xr = xr + xmove
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert to fractional coordinates,
c     find image, then translate fractional back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         xr = xr + xmove
         do while (abs(xr) .gt. xcell2)
            xr = xr - sign(xcell,xr)
         end do
         yr = yr + ymove
         do while (abs(yr) .gt. ycell2)
            yr = yr - sign(ycell,yr)
         end do
         zr = zr + zmove
         do while (abs(zr) .gt. zcell2)
            zr = zr - sign(zcell,zr)
         end do
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
         corr = box23 * int(abs(xr/xbox)+abs(yr/ybox)+abs(zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)
c
c     for rhombic dodecahedron, align along the x- and y-axes
c
      else if (dodecadron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - root2*zbox*nint(zr/(zbox*root2))
         corr = xbox2 * int(abs(xr/xbox)+abs(yr/ybox)
     &                        +abs(root2*zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)*root2
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine imagen  --  fast minimum image magnitude  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "imagen" takes the components of pairwise distance between
c     two points and converts to the components of the minimum
c     image distance
c
c     note this is a fast version for use in computing the 3D
c     distance during neighbor list construction
c
c
      subroutine imagen (xr,yr,zr)
      use boxes
      use math
      implicit none
      real*8 xr,yr,zr
      real*8 corr
c
c
c     for orthogonal lattice, find the desired image directly
c
      if (orthogonal) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
c
c     for monoclinic lattice, convert x and z to fractional,
c     find desired image, then translate back to Cartesian
c
      else if (monoclinic) then
         zr = zr / beta_sin
         xr = xr - zr*beta_cos
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
         xr = xr + zr*beta_cos
         zr = zr * beta_sin
c
c     for triclinic lattice, convert to fractional coordinates,
c     find image, then translate fractional back to Cartesian
c
      else if (triclinic) then
         zr = zr / gamma_term
         yr = (yr - zr*beta_term) / gamma_sin
         xr = xr - yr*gamma_cos - zr*beta_cos
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
         xr = xr + yr*gamma_cos + zr*beta_cos
         yr = yr*gamma_sin + zr*beta_term
         zr = zr * gamma_term
c
c     for truncated octahedron, remove the corner pieces
c
      else if (octahedron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - zbox*nint(zr/zbox)
         corr = box23 * int(abs(xr/xbox)+abs(yr/ybox)+abs(zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)
c
c     for rhombic dodecahedron, align along the x- and y-axes
c
      else if (dodecadron) then
         xr = xr - xbox*nint(xr/xbox)
         yr = yr - ybox*nint(yr/ybox)
         zr = zr - root2*zbox*nint(zr/(zbox*root2))
         corr = xbox2 * int(abs(xr/xbox)+abs(yr/ybox)
     &                        +abs(root2*zr/zbox))
         xr = xr - sign(corr,xr)
         yr = yr - sign(corr,yr)
         zr = zr - sign(corr,zr)*root2
      end if
      return
      end
