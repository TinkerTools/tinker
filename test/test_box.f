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
c     by the tinker-gpu box.cpp tests
c
c
      subroutine test_box
      implicit none
c
c
      call test_box_orthogonal
      call test_box_monoclinic
      call test_box_triclinic
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
