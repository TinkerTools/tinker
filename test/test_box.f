c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################
c     ##                                          ##
c     ##  subroutine test_box  --  box utilities  ##
c     ##                                          ##
c     ##############################################
c
c
c     "test_box" checks the volume and minimum-image cases exercised
c     by the tinker-gpu box.cpp tests, plus the rhombic dodecahedron
c
c
      subroutine test_box
      implicit none
c
c
      call test_box_orthogonal
      call test_box_monoclinic
      call test_box_triclinic
      call test_box_dodecahedron
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine test_box_orthogonal  --  cubic box  ##
c     ##                                                 ##
c     #####################################################
c
c
      subroutine test_box_orthogonal
      use boxes
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_box_orthogonal')
c
c
      if (skiptest(tname,'amoeba'))  return
      call initial
      xbox = 16.0d0
      ybox = 16.0d0
      zbox = 16.0d0
      alpha = 90.0d0
      beta = 90.0d0
      gamma = 90.0d0
      call lattice
      call assert_real (volbox,4096.0d0,1.0d-6,tname//' volume')
      call test_box_image (tname,-16.0d0,-12.0d0,-8.0d0,
     &                     0.0d0,4.0d0,-8.0d0)
      call test_box_image (tname,-4.0d0,0.0d0,4.0d0,
     &                     -4.0d0,0.0d0,4.0d0)
      call test_box_image (tname,8.0d0,12.0d0,16.0d0,
     &                     8.0d0,-4.0d0,0.0d0)
      call test_box_image (tname,51.0d0,-83.0d0,164.0d0,
     &                     3.0d0,-3.0d0,4.0d0)
      call final
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine test_box_monoclinic  --  monoclinic  ##
c     ##                                                  ##
c     ######################################################
c
c
      subroutine test_box_monoclinic
      use boxes
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_box_monoclinic')
c
c
      if (skiptest(tname,'amoeba'))  return
      call initial
      xbox = 32.0d0
      ybox = 24.0d0
      zbox = 20.0d0
      alpha = 90.0d0
      beta = 30.0d0
      gamma = 90.0d0
      call lattice
      call assert_real (volbox,7680.0d0,1.0d-6,tname//' volume')
      call assert_real (lvec(1,1),32.0d0,1.0d-6,tname//' lvec11')
      call assert_real (lvec(1,2),0.0d0,1.0d-6,tname//' lvec12')
      call assert_real (lvec(1,3),0.0d0,1.0d-6,tname//' lvec13')
      call assert_real (lvec(2,1),0.0d0,1.0d-6,tname//' lvec21')
      call assert_real (lvec(2,2),24.0d0,1.0d-6,tname//' lvec22')
      call assert_real (lvec(2,3),0.0d0,1.0d-6,tname//' lvec23')
      call assert_real (lvec(3,1),17.3205080757d0,1.0d-6,
     &                  tname//' lvec31')
      call assert_real (lvec(3,2),0.0d0,1.0d-6,tname//' lvec32')
      call assert_real (lvec(3,3),10.0d0,1.0d-6,tname//' lvec33')
      call test_box_image (tname,0.0d0,0.0d0,0.0d0,
     &                     0.0d0,0.0d0,0.0d0)
      call test_box_image (tname,-8.0d0,-6.0d0,0.0d0,
     &                     -8.0d0,-6.0d0,0.0d0)
      call test_box_image (tname,5.0d0,10.0d0,15.0d0,
     &                     2.3589838486d0,10.0d0,-5.0d0)
      call test_box_image (tname,-13.0d0,-30.0d0,20.0d0,
     &                     -15.6410161514d0,-6.0d0,0.0d0)
      call test_box_image (tname,-18.0d0,-40.0d0,5.0d0,
     &                     -3.3205080757d0,8.0d0,-5.0d0)
      call test_box_image (tname,-18.0d0,-16.0d0,5.0d0,
     &                     -3.3205080757d0,8.0d0,-5.0d0)
      call final
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine test_box_triclinic  --  triclinic  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine test_box_triclinic
      use boxes
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_box_triclinic')
c
c
      if (skiptest(tname,'amoeba'))  return
      call initial
      xbox = 32.0d0
      ybox = 24.0d0
      zbox = 20.0d0
      alpha = 75.0d0
      beta = 60.0d0
      gamma = 45.0d0
      call lattice
      call assert_real (volbox,9292.805125725014d0,1.0d-6,
     &                  tname//' volume')
      call test_box_image (tname,0.0d0,0.0d0,0.0d0,
     &                     0.0d0,0.0d0,0.0d0)
      call test_box_image (tname,5.0d0,10.0d0,15.0d0,
     &                     10.02943725d0,-4.29107082d0,
     &                     -2.11199354d0)
      call test_box_image (tname,-13.0d0,-30.0d0,20.0d0,
     &                     10.94112550d0,6.62061742d0,
     &                     2.88800646d0)
      call test_box_image (tname,-18.0d0,-40.0d0,5.0d0,
     &                     -16.05887450d0,-6.05887450d0,
     &                     5.0d0)
      call test_box_image (tname,0.91168825d0,10.91168824d0,
     &                     5.0d0,-16.05887450d0,
     &                     -6.05887450d0,5.0d0)
      call final
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_box_dodecahedron  --  rhombic dodec  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_box_dodecahedron" checks the rhombic dodecahedron with
c     edge length a, stored as the a x a x sqrt(2)a parent cell, by
c     hand-computed images and by a sweep against a brute-force
c     search over the face-centered cubic lattice vectors
c
c
      subroutine test_box_dodecahedron
      use boxes
      use math
      implicit none
      integer i,j,k
      integer n1,n2,n3
      integer nbad,ngrid
      parameter (ngrid=10)
      real*8 a,eps
      real*8 xr,yr,zr
      real*8 xn,yn,zn
      real*8 r2,rmin2
      real*8 c1,c2,c3
      real*8 x0,y0,z0
      real*8 lv(3,3)
      logical skiptest
      character*(*) tname
      parameter (tname='test_box_dodecahedron')
c
c
      if (skiptest(tname,'amoeba'))  return
      call initial
      a = 16.0d0
      xbox = a
      ybox = a
      zbox = a * root2
      alpha = 90.0d0
      beta = 90.0d0
      gamma = 90.0d0
      octahedron = .false.
      dodecadron = .true.
      call lattice
      call assert_real (volbox,2896.3093757400984d0,1.0d-6,
     &                  tname//' volume')
c
c     hand-computed images, the first two missed by a z wrap
c     with period 2a in place of sqrt(2)a
c
      call test_box_image (tname,1.0d0,1.0d0,15.0d0,
     &                     1.0d0,1.0d0,-7.627416997970d0)
      call test_box_image (tname,20.0d0,-3.0d0,30.0d0,
     &                     -4.0d0,5.0d0,-3.941125496954d0)
      call test_box_image (tname,7.0d0,7.0d0,7.0d0,
     &                     -1.0d0,-1.0d0,-4.313708498985d0)
      call test_box_image (tname,7.0d0,-7.0d0,-7.0d0,
     &                     -1.0d0,1.0d0,4.313708498985d0)
      call test_box_image (tname,-4.0d0,0.0d0,4.0d0,
     &                     -4.0d0,0.0d0,4.0d0)
c
c     face-centered cubic lattice vectors with nearest-neighbor
c     distance a, built here rather than taken from "lattice"
c
      do i = 1, 3
         do j = 1, 3
            lv(j,i) = 0.0d0
         end do
      end do
      lv(1,1) = a
      lv(2,2) = a
      lv(1,3) = 0.5d0 * a
      lv(2,3) = 0.5d0 * a
      lv(3,3) = 0.5d0 * a * root2
c
c     sweep a grid over [-1.5a,1.5a] with irrational offsets; each
c     image must be a lattice shift of the input, as short as the
c     brute-force minimum, and agree with "imagen" in magnitude
c
      eps = 1.0d-6
      nbad = 0
      do i = 1, ngrid
         do j = 1, ngrid
            do k = 1, ngrid
               x0 = 3.0d0*a*(dble(i)-0.381966d0)/dble(ngrid) - 1.5d0*a
               y0 = 3.0d0*a*(dble(j)-0.732051d0)/dble(ngrid) - 1.5d0*a
               z0 = 3.0d0*a*(dble(k)-0.414214d0)/dble(ngrid) - 1.5d0*a
               rmin2 = 1.0d10
               do n1 = -5, 5
                  do n2 = -5, 5
                     do n3 = -5, 5
                        xr = x0 + n1*lv(1,1) + n2*lv(1,2) + n3*lv(1,3)
                        yr = y0 + n1*lv(2,1) + n2*lv(2,2) + n3*lv(2,3)
                        zr = z0 + n1*lv(3,1) + n2*lv(3,2) + n3*lv(3,3)
                        r2 = xr*xr + yr*yr + zr*zr
                        rmin2 = min(r2,rmin2)
                     end do
                  end do
               end do
               xr = x0
               yr = y0
               zr = z0
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               c3 = (zr-z0) / lv(3,3)
               c1 = (xr-x0-c3*lv(1,3)) / a
               c2 = (yr-y0-c3*lv(2,3)) / a
               xn = x0
               yn = y0
               zn = z0
               call imagen (xn,yn,zn)
               if (abs(sqrt(r2)-sqrt(rmin2)) .gt. eps .or.
     &             abs(c1-dble(nint(c1))) .gt. eps .or.
     &             abs(c2-dble(nint(c2))) .gt. eps .or.
     &             abs(c3-dble(nint(c3))) .gt. eps .or.
     &             abs(abs(xn)-abs(xr)) .gt. eps .or.
     &             abs(abs(yn)-abs(yr)) .gt. eps .or.
     &             abs(abs(zn)-abs(zr)) .gt. eps) then
                  nbad = nbad + 1
               end if
            end do
         end do
      end do
      call assert_int (nbad,0,tname//' brute-force images')
      dodecadron = .false.
      nonprism = .false.
      call final
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine test_box_image  --  compare image  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine test_box_image (tname,xr0,yr0,zr0,xa,ya,za)
      implicit none
      real*8 xr0,yr0,zr0,xa,ya,za
      real*8 xr,yr,zr
      character*(*) tname
c
c
      xr = xr0
      yr = yr0
      zr = zr0
      call image (xr,yr,zr)
      call assert_real (xr,xa,1.0d-6,tname//' image x')
      call assert_real (yr,ya,1.0d-6,tname//' image y')
      call assert_real (zr,za,1.0d-6,tname//' image z')
      xr = xr0
      yr = yr0
      zr = zr0
      call imagen (xr,yr,zr)
      call assert_real (abs(xr),abs(xa),1.0d-6,tname//' imagen x')
      call assert_real (abs(yr),abs(ya),1.0d-6,tname//' imagen y')
      call assert_real (abs(zr),abs(za),1.0d-6,tname//' imagen z')
      return
      end
