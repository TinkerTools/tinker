c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  mpole.i  --  multipole components for current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c     dpole     derivative rotation matrix for each multipole
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     polsiz    number of mutipole components at each multipole site
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     polaxe    local axis type for each multipole site
c
c
      integer maxpole
      parameter (maxpole=13)
      integer npole,ipole,polsiz
      integer zaxis,xaxis
      real*8 pole,rpole,dpole
      character*8 polaxe
      common /mpole/ pole(maxpole,maxatm),rpole(maxpole,maxatm),
     &               dpole(maxpole,3,3,maxatm),npole,ipole(maxatm),
     &               polsiz(maxatm),zaxis(maxatm),xaxis(maxatm),
     &               polaxe(maxatm)
c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ###########################
c     ##                       ##
c     ##  subroutine drotpole  ##
c     ##                       ##
c     ###########################
c
c
c     "drotpole" computes the derivatives of the atomic multipoles
c     in the global coordinate frame with respect to motion of the
c     sites that define the local coordinate frame
c
c
      subroutine drotpole (imdq,a,d,p,q)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m,nn
      integer imdq,p,q
      real*8 a(3,3),d(3,3,3,3)
      real*8 m2(3,3),dm2(3,3)
c
c
c     derivative for rotation of the monopole term is zero
c
      dpole(1,p,q,imdq) = 0.0d0
c
c     derivative terms for rotation of the dipole components
c
      do i = 2, 4
         dpole(i,p,q,imdq) = 0.0d0
         do j = 2, 4
            dpole(i,p,q,imdq) = dpole(i,p,q,imdq)
     &                              + pole(j,imdq)*d(i-1,j-1,p,q)
         end do
      end do
c
c     derivative terms for rotation of the quadrupole components
c
      nn = 5
      do i = 1, 3
         do j = 1, 3
            m2(i,j) = pole(nn,imdq)
            dm2(i,j) = 0.0d0
            nn = nn + 1
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do m = 1, 3
                  dm2(i,j) = dm2(i,j) + m2(k,m) *
     &                         (d(i,k,p,q)*a(j,m)+a(i,k)*d(j,m,p,q))
               end do
            end do
          end do
      end do
      nn = 5
      do i = 1, 3
         do j = 1, 3
            dpole(nn,p,q,imdq) = dm2(i,j)
            nn = nn + 1
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine drotmat  --  multipole rotation matrix derivs  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "drotmat" finds the derivative rotation matrices that convert
c     multipoles from the local coordinate system to the global system
c
c
      subroutine drotmat (i,d)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i
      real*8 d(3,3,3,3)
c
c
      if (polaxe(i) .eq. 'Z-then-X') then
         call drotmat1 (i,d)
      else if (polaxe(i) .eq. 'Bisector') then
         call drotmat2 (i,d)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine drotmat1  --  Z-then-X local coordinate derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "drotmat1" finds the multipole rotation matrix derivatives
c     for local coordinates defined via the "Z-then-X" method
c
c
      subroutine drotmat1 (i,d)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k
      real*8 xi,yi,zi
      real*8 xz,yz,zz
      real*8 xx,yx,zx
      real*8 t1,t2,t3,t4,t5
      real*8 t6,t7,t8,t9,t10
      real*8 t12,t13,t14,t15
      real*8 t16,t17,t18,t19
      real*8 t20,t21,t22,t23
      real*8 t24,t25,t26,t27
      real*8 t28,t29,t30,t31
      real*8 t32,t33,t34,t35
      real*8 t36,t37,t38,t39
      real*8 t40,t41,t42,t43
      real*8 t44,t46,t47,t49
      real*8 t50,t52,t58,t59
      real*8 t60,t62,t64,t65
      real*8 t66,t68,t69,t70
      real*8 t72,t73,t74,t75
      real*8 t76,t77,t78,t80
      real*8 t81,t86,t91,t93
      real*8 t95,t97,t98,t101
      real*8 t105,t106,t111,t116
      real*8 t118,t120,t123,t124
      real*8 t126,t130,t131,t136
      real*8 t141,t143,t145,t146
      real*8 t148,t150,t152,t153
      real*8 t155,t157,t159,t163
      real*8 t167,t168,t173,t178
      real*8 t180,t182,t186,t190
      real*8 t191,t196,t201,t203
      real*8 t205,t209,t213,t214
      real*8 t219,t224,t226,t228
      real*8 t230,t232,t234,t236
      real*8 t238,t240,t242,t243
      real*8 t246,t249,t251,t253
      real*8 t254,t255,t256,t259
      real*8 t263,t265,t266,t271
      real*8 t273,t276,t278,t282
      real*8 t285,t287,t290,t292
      real*8 t295,t297,t300,t304
      real*8 t308,t310,t314,t317
      real*8 t324,t329,t332,t354
      real*8 t358,t362,t366,t370
      real*8 t374,t381,t388,t395
      real*8 t416,t418,t422,t425
      real*8 t428,t430,t432,t436
      real*8 t438,t441,t446
      real*8 d(3,3,3,3)
c
c
c     coordinates of main atom and those defining local axes
c
      xi = x(ipole(i))
      yi = y(ipole(i))
      zi = z(ipole(i))
      xz = x(zaxis(i))
      yz = y(zaxis(i))
      zz = z(zaxis(i))
      xx = x(xaxis(i))
      yx = y(xaxis(i))
      zx = z(xaxis(i))
c
c     temporary variables from Maple symbolic algebra derivation
c
      t1 = xz - xi
      t2 = t1 * t1
      t3 = yz - yi
      t4 = t3 * t3
      t5 = zz - zi
      t6 = t5 * t5
      t7 = t2 + t4 + t6
      t8 = sqrt(t7)
      t9 = 1.0d0 / t8
      t10 = t7 * t7
      t12 = t8 / t10
      t13 = t1 * t12
      t14 = 2.0d0 * (xi-xz)
      t15 = t13 * t14
      t16 = 2.0d0 * (yi-yz)
      t17 = t13 * t16
      t18 = 2.0d0 * (zi-zz)
      t19 = t13 * t18
      t20 = t3 * t12
      t21 = t20 * t14
      t22 = t20 * t16
      t23 = t20 * t18
      t24 = t5 * t12
      t25 = t24 * t14
      t26 = t24 * t16
      t27 = t24 * t18
      t28 = 2.0d0 * (xz-xi)
      t29 = t13 * t28
      t30 = 2.0d0 * (yz-yi)
      t31 = t13 * t30
      t32 = 2.0d0 * (zz-zi)
      t33 = t13 * t32
      t34 = t20 * t28
      t35 = t20 * t30
      t36 = t20 * t32
      t37 = t24 * t28
      t38 = t24 * t30
      t39 = t24 * t32
      t40 = t1 * t9
      t41 = xx - xi
      t42 = t41 * t9
      t43 = t41 * t1
      t44 = t12 * t14
      t46 = yx - yi
      t47 = t46 * t3
      t49 = zx - zi
      t50 = t49 * t5
      t52 = -t40 - t42 - 0.5d0*t44*(t43+t47+t50)
      t58 = t9 * (t43+t47+t50)
      t59 = t58 * t9
      t60 = t58 * t1
      t62 = -1.0d0 - t52*t1*t9 + t59 + 0.5d0*t60*t44
      t64 = xx - xi - t60*t9
      t65 = t64 * t64
      t66 = t58 * t3
      t68 = yx - yi - t66*t9
      t69 = t68 * t68
      t70 = t58 * t5
      t72 = zx - zi - t70*t9
      t73 = t72 * t72
      t74 = t65 + t69 + t73
      t75 = sqrt(t74)
      t76 = 1.0d0 / t75
      t77 = t62 * t76
      t78 = t74 * t74
      t80 = t75 / t78
      t81 = t64 * t80
      t86 = -t52*t3*t9 + 0.5d0*t66*t44
      t91 = -t52*t5*t9 + 0.5d0*t70*t44
      t93 = 2.0d0 * (t64*t62 + t68*t86 + t72*t91)
      t95 = t12 * t16
      t97 = t3 * t9
      t98 = t46 * t9
      t101 =  -t97 - t98 - 0.5d0*t95*(t43+t47+t50)
      t105 = -t101*t1*t9 + 0.5d0*t60*t95
      t106 = t105 * t76
      t111 = -1.0d0 - t101*t3*t9 + t59 + 0.5d0*t66*t95
      t116 = -t101*t5*t9 + 0.5d0*t70*t95
      t118 = 2.0d0 * (t64*t105 + t68*t111 + t72*t116)
      t120 = t12 * t18
      t123 = t5 * t9
      t124 = t49 * t9
      t126 = - t123 - t124 - 0.5d0*t120*(t43+t47+t50)
      t130 = -t126*t1*t9 + 0.5d0*t60*t120
      t131 = t130 * t76
      t136 = -t126*t3*t9 + 0.5d0*t66*t120
      t141 = -1.0d0 - t126*t5*t9 + t59 + 0.5d0*t70*t120
      t143 = 2.0d0 * (t64*t130 + t68*t136 + t72*t141)
      t145 = t86 * t76
      t146 = t68 * t80
      t148 = t111 * t76
      t150 = t136 * t76
      t152 = t91 * t76
      t153 = t72 * t80
      t155 = t116 * t76
      t157 = t141 * t76
      t159 = t12 * t28
      t163 = t42 - 0.5d0*t159*(t43+t47+t50)
      t167 = -t163*t1*t9 - t59 + 0.5d0*t60*t159
      t168 = t167 * t76
      t173 = -t163*t3*t9 + 0.5d0*t66*t159
      t178 = -t163*t5*t9 + 0.5d0*t70*t159
      t180 = 2.0d0 * (t64*t167 + t68*t173 + t72*t178)
      t182 = t12 * t30
      t186 = t98 - 0.5d0*t182*(t43+t47+t50)
      t190 = -t186*t1*t9 + 0.5d0*t60*t182
      t191 = t190 * t76
      t196 = -t186*t3*t9 - t59 + 0.5d0*t66*t182
      t201 = -t186*t5*t9 + 0.5d0*t70*t182
      t203 = 2.0d0 * (t64*t190 + t68*t196 + t72*t201)
      t205 = t12 * t32
      t209 = t124 - 0.5d0*t205*(t43+t47+t50)
      t213 = -t209*t1*t9 + 0.5d0*t60*t205
      t214 = t213 * t76
      t219 = -t209*t3*t9 + 0.5d0*t66*t205
      t224 = -t209*t5*t9 - t59 + 0.5d0*t70*t205
      t226 = 2.0d0 * (t64*t213 + t68*t219 + t72*t224)
      t228 = t173 * t76
      t230 = t196 * t76
      t232 = t219 * t76
      t234 = t178 * t76
      t236 = t201 * t76
      t238 = t224 * t76
      t240 = 1.0d0 / t7
      t242 = 1.0d0 - t2*t240
      t243 = t242 * t76
      t246 = t240 * t3
      t249 = t240 * t5
      t251 = 2.0d0 * (t64*t242 - t68*t1*t246 - t72*t1*t249)
      t253 = t1 * t240
      t254 = t3 * t76
      t255 = t253 * t254
      t256 = t64 * t1
      t259 = 1.0d0 - t4*t240
      t263 =  2.0d0 * (t68*t259 - t256*t246 - t72*t3*t249)
      t265 = t5 * t76
      t266 = t253 * t265
      t271 = 1.0d0 - t6*t240
      t273 =  2.0d0 * (t72*t271 - t256*t249 - t68*t3*t249)
      t276 = t259 * t76
      t278 = t246 * t265
      t282 = t271 * t76
      t285 = t97 * t93
      t287 = t72 * t76
      t290 = t123 * t93
      t292 = t68 * t76
      t295 = t97 * t118
      t297 = t287 * t9
      t300 = t123 * t118
      t304 = t97 * t143
      t308 = t123 * t143
      t310 = t292 * t9
      t314 = t64 * t76
      t317 = t40 * t93
      t324 = t40 * t118
      t329 = t314 * t9
      t332 = t40 * t143
      t354 = t97 * t180
      t358 = t123 * t180
      t362 = t97 * t203
      t366 = t123 * t203
      t370 = t97 * t226
      t374 = t123 * t226
      t381 = t40 * t180
      t388 = t40 * t203
      t395 = t40 * t226
      t416 = t97 * t251
      t418 = t123 * t251
      t422 = t97 * t263
      t425 = t123 * t263
      t428 = t97 * t273
      t430 = t6 * t76
      t432 = t123 * t273
      t436 = t2 * t12
      t438 = t40 * t251
      t441 = t40 * t263
      t446 = t40 * t273
c
c     set values for the derivatives of the rotation matrix
c
      d(1,3,1,1) = -t9 - 0.5d0*t15
      d(1,3,1,2) = -0.5d0 * t17
      d(1,3,1,3) = -0.5d0 * t19
      d(2,3,1,1) = -0.5d0 * t21
      d(2,3,1,2) = -t9 - 0.5d0*t22
      d(2,3,1,3) = -0.5d0 * t23
      d(3,3,1,1) = -0.5d0 * t25
      d(3,3,1,2) = -0.5d0 * t26
      d(3,3,1,3) = -t9 - 0.5d0*t27
      d(1,3,2,1) = t9 - 0.5d0*t29
      d(1,3,2,2) = -0.5d0 * t31
      d(1,3,2,3) = -0.5d0 * t33
      d(2,3,2,1) = -0.5d0 * t34
      d(2,3,2,2) = t9 - 0.5d0*t35
      d(2,3,2,3) = -0.5d0 * t36
      d(3,3,2,1) = -0.5d0 * t37
      d(3,3,2,2) = -0.5d0 * t38
      d(3,3,2,3) = t9 - 0.5d0*t39
      d(1,1,1,1) = t77 - 0.5d0*t81*t93
      d(1,1,1,2) = t106 - 0.5d0*t81*t118
      d(1,1,1,3) = t131 - 0.5d0*t81*t143
      d(2,1,1,1) = t145 - 0.5d0*t146*t93
      d(2,1,1,2) = t148 - 0.5d0*t146*t118
      d(2,1,1,3) = t150 - 0.5d0*t146*t143
      d(3,1,1,1) = t152 - 0.5d0*t153*t93
      d(3,1,1,2) = t155 - 0.5d0*t153*t118
      d(3,1,1,3) = t157 - 0.5d0*t153*t143
      d(1,1,2,1) = t168 - 0.5d0*t81*t180
      d(1,1,2,2) = t191 - 0.5d0*t81*t203
      d(1,1,2,3) = t214 - 0.5d0*t81*t226
      d(2,1,2,1) = t228 - 0.5d0*t146*t180
      d(2,1,2,2) = t230 - 0.5d0*t146*t203
      d(2,1,2,3) = t232 - 0.5d0*t146*t226
      d(3,1,2,1) = t234 - 0.5d0*t153*t180
      d(3,1,2,2) = t236 - 0.5d0*t153*t203
      d(3,1,2,3) = t238 - 0.5d0*t153*t226
      d(1,1,3,1) = t243 - 0.5d0*t81*t251
      d(1,1,3,2) = -t255 - 0.5d0*t81*t263
      d(1,1,3,3) = -t266 - 0.5d0*t81*t273
      d(2,1,3,1) = -t255 - 0.5d0*t146*t251
      d(2,1,3,2) = t276 - 0.5d0*t146*t263
      d(2,1,3,3) = -t278 - 0.5d0*t146*t273
      d(3,1,3,1) = -t266 - 0.5d0*t153*t251
      d(3,1,3,2) = -t278 - 0.5d0*t153*t263
      d(3,1,3,3) = t282 - 0.5d0*t153*t273
      d(1,2,1,1) = t152*t97 - 0.5d0*(t153*t285 + t287*t21)
     &                - t145*t123 + 0.5d0*(t146*t290 + t292*t25)
      d(1,2,1,2) = t155*t97 - t297 - 0.5d0*(t153*t295 + t287*t22)
     &                - t148*t123 + 0.5d0*(t146*t300 + t292*t26)
      d(1,2,1,3) = t157*t97 + t310 - 0.5d0*(t153*t304 + t287*t23)
     &                - t150*t123 + 0.5d0*(t146*t308 + t292*t27)
      d(2,2,1,1) = t77*t123 + t297 - 0.5d0*(t81*t290 + t314*t25)
     &                - t152*t40 + 0.5d0*(t153*t317 + t287*t15)
      d(2,2,1,2) = t106*t123 - 0.5d0*(t81*t300 + t314*t26)
     &                - t155*t40 + 0.5d0*(t153*t324 + t287*t17)
      d(2,2,1,3) = t131*t123 - t329 - 0.5d0*(t81*t308 + t314*t27)
     &                - t157*t40 + 0.5d0*(t153*t332 + t287*t19)
      d(3,2,1,1) = t145*t40 - t310 - 0.5d0*(t146*t317 + t292*t15)
     &                - t77*t97 + 0.5d0*(t81*t285 + t314*t21)
      d(3,2,1,2) = t148*t40 + t329 - 0.5d0*(t146*t324 + t292*t17)
     &                - t106*t97 + 0.5d0*(t81*t295 + t314*t22)
      d(3,2,1,3) = t150*t40 - 0.5d0*(t146*t332 + t292*t19)
     &                - t131*t97 + 0.5d0*(t81*t304 + t314*t23)
      d(2,2,1,3) = t131*t123 - t329 - 0.5d0*(t81*t308 + t314*t27)
     &                - t157*t40 + 0.5d0*(t153*t332 + t287*t19)
      d(3,2,1,1) = t145*t40 - t310 - 0.5d0*(t146*t317 + t292*t15)
     &                - t77*t97 + 0.5d0*(t81*t285 + t314*t21)
      d(3,2,1,2) = t148*t40 + t329 - 0.5d0*(t146*t324 + t292*t17)
     &                - t106*t97 + 0.5d0*(t81*t295 + t314*t22)
      d(3,2,1,3) = t150*t40 - 0.5d0*(t146*t332 + t292*t19)
     &                - t131*t97 + 0.5d0*(t81*t304 + t314*t23)
      d(1,2,2,1) = t234*t97 - 0.5d0*(t153*t354 + t287*t34)
     &                - t228*t123 + 0.5d0*(t146*t358 + t292*t37)
      d(1,2,2,2) = t236*t97 + t297 - 0.5d0*(t153*t362 + t287*t35)
     &                - t230*t123 + 0.5d0*(t146*t366 + t292*t38)
      d(1,2,2,3) = t238*t97 - t310 - 0.5d0*(t153*t370 + t287*t36)
     &                - t232*t123 + 0.5d0*(t146*t374 + t292*t39)
      d(2,2,2,1) = t168*t123 - t297 - 0.5d0*(t81*t358 + t314*t37)
     &                - t234*t40 + 0.5d0*(t153*t381 + t287*t29)
      d(2,2,2,2) = t191*t123 - 0.5d0*(t81*t366 + t314*t38)
     &                - t236*t40 + 0.5d0*(t153*t388 + t287*t31)
      d(2,2,2,3) = t214*t123 + t329 - 0.5d0*(t81*t374 + t314*t39)
     &                - t238*t40 + 0.5d0*(t153*t395 + t287*t33)
      d(3,2,2,1) = t228*t40 + t310 - 0.5d0*(t146*t381 + t292*t29)
     &                - t168*t97 + 0.5d0*(t81*t354 + t314*t34)
      d(3,2,2,2) = t230*t40 - t329 - 0.5d0*(t146*t388 + t292*t31)
     &                - t191*t97 + 0.5d0*(t81*t362 + t314*t35)
      d(3,2,2,3) = t232*t40 - 0.5d0*(t146*t395 + t292*t33)
     &                - t214*t97 + 0.5d0*(t81*t370 + t314*t36)
      d(1,2,3,1) = 0.5d0 * (t146*t418 - t153*t416)
      d(1,2,3,2) = -t4*t12*t265 - t276*t123
     &                + 0.5d0*(t146*t425 - t153*t422)
      d(1,2,3,3) = t282*t97 + t20*t430 + 0.5d0*(t146*t432-t153*t428)
      d(2,2,3,1) = t243*t123 + t436*t265 + 0.5d0*(t153*t438-t81*t418)
      d(2,2,3,2) = 0.5d0 * (t153*t441 - t81*t425)
      d(2,2,3,3) = -t13*t430 - t282*t40 + 0.5d0*(t153*t446-t81*t432)
      d(3,2,3,1) = -t436*t254 - t243*t97 + 0.5d0*(t81*t416-t146*t438)
      d(3,2,3,2) = t276*t40 + t13*t4*t76 + 0.5d0*(t81*t422-t146*t441)
      d(3,2,3,3) = 0.5d0 * (t81*t428 - t146*t446)
c
c     some of the derivative matrix values are always zero
c
      do j = 1, 3
         do k = 1, 3
            d(k,3,3,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine drotmat2  --  bisector local coordinate derivs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "drotmat2" finds the multipole rotation matrix derivatives
c     for local coordinates defined via the "bisector" method
c
c
      subroutine drotmat2 (i,d)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'mpole.i'
      integer i,j,k,m
      real*8 xi,yi,zi
      real*8 xz,yz,zz
      real*8 xx,yx,zx
      real*8 t1,t2,t3,t4,t5
      real*8 t6,t7,t8,t9,t10
      real*8 t12,t13,t14,t16
      real*8 t18,t19,t20,t21
      real*8 t22,t23,t24,t25
      real*8 t26,t27,t28,t29
      real*8 t31,t32,t33,t35
      real*8 t36,t37,t38,t39
      real*8 t40,t41,t42,t44
      real*8 t45,t47,t48,t50
      real*8 t52,t53,t54,t55
      real*8 t57,t58,t60,t62
      real*8 t65,t66,t67,t68
      real*8 t70,t73,t75,t77
      real*8 t78,t79,t81,t82
      real*8 t83,t84,t86,t88
      real*8 t89,t91,t92,t93
      real*8 t94,t95,t96,t97
      real*8 t98,t99,t100,t101
      real*8 t103,t104,t106,t108
      real*8 t109,t110,t111,t113
      real*8 t114,t116,t117,t118
      real*8 t121,t122,t123,t124
      real*8 t126,t129,t130,t131
      real*8 t133,t134,t135,t137
      real*8 t138,t139,t141,t143
      real*8 t145,t146,t147,t150
      real*8 t151,t154,t157,t160
      real*8 t162,t164,t166,t169
      real*8 t170,t172,t174,t175
      real*8 t176,t178,t179,t180
      real*8 t182,t183,t184,t185
      real*8 t186,t187,t188,t190
      real*8 t191,t195,t198,t202
      real*8 t205,t207,t209,t210
      real*8 t217,t220,t222,t224
      real*8 t225,t232,t238,t240
      real*8 t242,t249,t253,t255
      real*8 t256,t262,t269,t271
      real*8 t273,t274,t276,t278
      real*8 t280,t281,t283,t285
      real*8 t295,t300,t303,t308
      real*8 t310,t311,t316,t317
      real*8 t320,t325,t328,t330
      real*8 t335,t348,t351,t352
      real*8 t354,t356,t357,t364
      real*8 t371,t373,t390,t393
      real*8 t395,t397,t398,t405
      real*8 t412,t414,t416,t418
      real*8 t420,t422,t424,t426
      real*8 t429,t431,t432,t436
      real*8 t438,t439,t443,t448
      real*8 t453,t458,t464,t465
      real*8 t469,t478,t480,t488
      real*8 t498,t503,t519,t521
      real*8 t525,t527,t531,t536
      real*8 t541,t546,t552,t556
      real*8 t565,t567,t575
      real*8 t585,t590
      real*8 d(3,3,3,3)
c
c
c     coordinates of main atom and those defining local axes
c
      xi = x(ipole(i))
      yi = y(ipole(i))
      zi = z(ipole(i))
      xz = x(zaxis(i))
      yz = y(zaxis(i))
      zz = z(zaxis(i))
      xx = x(xaxis(i))
      yx = y(xaxis(i))
      zx = z(xaxis(i))
c
c     temporary variables from Maple symbolic algebra derivation
c
      t1 = xz - xi
      t2 = t1**2
      t3 = yz - yi
      t4 = t3**2
      t5 = zz - zi
      t6 = t5**2
      t7 = t2 + t4 + t6
      t8 = dsqrt(t7)
      t9 = 1.0d0 / t8
      t10 = t7**2
      t12 = t8 / t10
      t13 = t1 * t12
      t14 = 2.0d0 * (xz-xi)
      t16 = t9 - 0.5d0*t13*t14
      t18 = xx - xi
      t19 = t18**2
      t20 = yx - yi
      t21 = t20**2
      t22 = zx - zi
      t23 = t22**2
      t24 = t19 + t21 + t23
      t25 = dsqrt(t24)
      t26 = 1.0d0 / t25
      t27 = t18 * t26
      t28 = t1*t9 + t27
      t29 = t28**2
      t31 = t20 * t26
      t32 = t3*t9 + t31
      t33 = t32**2
      t35 = t22 * t26
      t36 = t5*t9 + t35
      t37 = t36**2
      t38 = t29 + t33 + t37
      t39 = dsqrt(t38)
      t40 = 1.0d0 / t39
      t41 = t16 * t40
      t42 = t38**2
      t44 = t39 / t42
      t45 = t28 * t44
      t47 = t32 * t3
      t48 = t12 * t14
      t50 = t36 * t5
      t52 = 2.0d0*t28*t16 - t47*t48 - t50*t48
      t53 = t45 * t52
      t54 = 2.0d0 * (yz-yi)
      t55 = t54 * t40
      t57 = t28 * t1
      t58 = t12 * t54
      t60 = t3 * t12
      t62 = t9 - 0.5d0*t60*t54
      t65 = -t57*t58 + 2.0d0*t32*t62 - t50*t58
      t66 = t45 * t65
      t67 = 2.0d0 * (zz-zi)
      t68 = t67 * t40
      t70 = t12 * t67
      t73 = t5 * t12
      t75 = t9 - 0.5d0*t73*t67
      t77 = -t57*t70 - t47*t70 + 2.0d0*t36*t75
      t78 = t45 * t77
      t79 = t14 * t40
      t81 = t32 * t44
      t82 = t81 * t52
      t83 = t62 * t40
      t84 = t81 * t65
      t86 = t81 * t77
      t88 = t36 * t44
      t89 = t88 * t52
      t91 = t88 * t65
      t92 = t75 * t40
      t93 = t88 * t77
      t94 = t24**2
      t95 = 1.0d0 / t94
      t96 = t25 * t95
      t97 = t18 * t96
      t98 = 2.0d0 * (xx-xi)
      t99 = t97 * t98
      t100 = t26 - 0.5d0*t99
      t101 = t100 * t40
      t103 = t32 * t20
      t104 = t96 * t98
      t106 = t36 * t22
      t108 = 2.0d0*t28*t100 - t103*t104 - t106*t104
      t109 = t45 * t108
      t110 = 2.0d0 * (yx-yi)
      t111 = t110 * t40
      t113 = t28 * t18
      t114 = t96 * t110
      t116 = t20 * t96
      t117 = t116 * t110
      t118 = t26 - 0.5d0*t117
      t121 = -t113*t114 + 2.0d0*t32*t118 - t106*t114
      t122 = t45 * t121
      t123 = 2.0d0 * (zx-zi)
      t124 = t123 * t40
      t126 = t96 * t123
      t129 = t22 * t96
      t130 = t129 * t123
      t131 = t26 - 0.5d0*t130
      t133 = -t113*t126 - t103*t126 + 2.0d0*t36*t131
      t134 = t45 * t133
      t135 = t98 * t40
      t137 = t81 * t108
      t138 = t118 * t40
      t139 = t81 * t121
      t141 = t81 * t133
      t143 = t88 * t108
      t145 = t88 * t121
      t146 = t131 * t40
      t147 = t88 * t133
      t150 = t31 * t3
      t151 = t48 * t40
      t154 = t35 * t5
      t157 = t27*t41 - 0.5d0*(t27*t53 + t150*t151 + t31*t82
     &                            + t154*t151 + t35*t89)
      t160 = t28 * t40
      t162 = t32 * t40
      t164 = t36 * t40
      t166 = t27*t160 + t31*t162 + t35*t164
      t169 = t166 * t28
      t170 = t44 * t52
      t172 = -t157*t28*t40 - t166*t16*t40 + 0.5d0*t169*t170
      t174 = t27 - t169*t40
      t175 = t174**2
      t176 = t166 * t32
      t178 = t31 - t176*t40
      t179 = t178**2
      t180 = t166 * t36
      t182 = t35 - t180*t40
      t183 = t182**2
      t184 = t175 + t179 + t183
      t185 = dsqrt(t184)
      t186 = 1.0d0 / t185
      t187 = t172 * t186
      t188 = t184**2
      t190 = t185 / t188
      t191 = t174 * t190
      t195 = t166 * t3
      t198 = -t157*t32*t40 + 0.5d0*(t195*t151 + t176*t170)
      t202 = t166 * t5
      t205 = -t157*t36*t40 + 0.5d0*(t202*t151 + t180*t170)
      t207 = 2.0d0 * (t174*t172 + t178*t198 + t182*t205)
      t209 = t27 * t1
      t210 = t58 * t40
      t217 = t31*t83 - 0.5d0*(t209*t210 + t27*t66 + t31*t84
     &                            + t154*t210 + t35*t91)
      t220 = t166 * t1
      t222 = t44 * t65
      t224 = -t217*t28*t40 + 0.5d0*(t220*t210 + t169*t222)
      t225 = t224 * t186
      t232 = -t217*t32*t40 - t166*t62*t40 + 0.5d0*t176*t222
      t238 = -t217*t36*t40 + 0.5d0*(t202*t210 + t180*t222)
      t240 = 2.0d0 * (t174*t224 + t178*t232 + t182*t238)
      t242 = t70 * t40
      t249 = t35*t92 - 0.5d0*(t209*t242 + t27*t78 + t150*t242
     &                              + t31*t86 + t35*t93)
      t253 = t44 * t77
      t255 = -t249*t28*t40 + 0.5d0*(t220*t242 + t169*t253)
      t256 = t255 * t186
      t262 = -t249*t32*t40 + 0.5d0*(t195*t242 + t176*t253)
      t269 = -t249*t36*t40 - t166*t75*t40 + 0.5d0*t180*t253
      t271 = 2.0d0 * (t174*t255 + t178*t262 + t182*t269)
      t273 = t198 * t186
      t274 = t178 * t190
      t276 = t232 * t186
      t278 = t262 * t186
      t280 = t205 * t186
      t281 = t182 * t190
      t283 = t238 * t186
      t285 = t269 * t186
      t295 = t21 * t95
      t300 = t23 * t95
      t303 = t26*t28*t40 + t27*t101
     &          - 0.5d0 * (t97*t160*t98 + t27*t109 + t116*t162*t98
     &                        + t295*t135 + t31*t137 + t300*t135
     &                        + t35*t143 + t129*t164*t98)
      t308 = t44 * t108
      t310 = t26 - t303*t28*t40 - t166*t100*t40
     &          + 0.5d0*(t169*t308 - t99)
      t311 = t310 * t186
      t316 = t166 * t20
      t317 = t104 * t40
      t320 = -t303*t32*t40 + 0.5d0*(t316*t317 + t176*t308 - t116*t98)
      t325 = t166 * t22
      t328 = -t303*t36*t40 + 0.5d0*(t325*t317 + t180*t308 - t129*t98)
      t330 = 2.0d0 * (t174*t310 + t178*t320 + t182*t328)
      t335 = t19 * t95
      t348 = t26*t32*t40 + t31*t138
     &          - 0.5d0 * (t97*t160*t110 + t335*t111 + t27*t122
     &                        + t116*t162*t110 + t129*t164*t110
     &                        + t300*t111 + t35*t145 + t31*t139)
      t351 = t166 * t18
      t352 = t114 * t40
      t354 = t44 * t121
      t356 = -t348*t28*t40 + 0.5d0*(t351*t352 + t169*t354 - t97*t110)
      t357 = t356 * t186
      t364 = t26 - t348*t32*t40 - t166*t118*t40
     &          + 0.5d0*(t176*t354 - t117)
      t371 = -t348*t36*t40 + 0.5d0*(t325*t352 + t180*t354 - t129*t110)
      t373 = 2.0d0 * (t174*t356 + t178*t364 + t182*t371)
      t390 = t26*t36*t40 + t35*t146
     &          - 0.5d0 * (t97*t160*t123 + t335*t124 + t27*t134
     &                        + t116*t162*t123 + t129*t164*t123
     &                        + t295*t124 + t31*t141 + t35*t147)
      t393 = t126 * t40
      t395 = t44 * t133
      t397 = -t390*t28*t40 + 0.5d0*(t351*t393 + t169*t395 - t97*t123)
      t398 = t397 * t186
      t405 = -t390*t32*t40 + 0.5d0*(t316*t393 + t176*t395 - t116*t123)
      t412 = t26 - t390*t36*t40 - t166*t131*t40
     &          + 0.5d0*(t180*t395 - t130)
      t414 = 2.0d0 * (t174*t397 + t178*t405 + t182*t412)
      t416 = t320 * t186
      t418 = t364 * t186
      t420 = t405 * t186
      t422 = t328 * t186
      t424 = t371 * t186
      t426 = t412 * t186
      t429 = t162 * t207
      t431 = t182 * t186
      t432 = t431 * t3
      t436 = t164 * t207
      t438 = t178 * t186
      t439 = t438 * t5
      t443 = t162 * t240
      t448 = t164 * t240
      t453 = t162 * t271
      t458 = t164 * t271
      t464 = t174 * t186
      t465 = t464 * t5
      t469 = t160 * t207
      t478 = t160 * t240
      t480 = t431 * t1
      t488 = t160 * t271
      t498 = t464 * t3
      t503 = t438 * t1
      t519 = t162 * t330
      t521 = t431 * t20
      t525 = t164 * t330
      t527 = t438 * t22
      t531 = t162 * t373
      t536 = t164 * t373
      t541 = t162 * t414
      t546 = t164 * t414
      t552 = t464 * t22
      t556 = t160 * t330
      t565 = t160 * t373
      t567 = t431 * t18
      t575 = t160 * t414
      t585 = t464 * t20
      t590 = t438 * t18
c
c     set values for the derivatives of the rotation matrix
c
      d(1,3,2,1) = t41 - 0.5d0*t53
      d(1,3,2,2) = -0.5d0 * (t13*t55 + t66)
      d(1,3,2,3) = -0.5d0 * (t13*t68 + t78)
      d(2,3,2,1) = -0.5d0 * (t60*t79 + t82)
      d(2,3,2,2) = t83 - 0.5d0*t84
      d(2,3,2,3) = -0.5d0 * (t60*t68 + t86)
      d(3,3,2,1) = -0.5d0 * (t73*t79 + t89)
      d(3,3,2,2) = -0.5d0 * (t73*t55 + t91)
      d(3,3,2,3) = t92 - 0.5d0*t93
      d(1,3,3,1) = t101 - 0.5d0*t109
      d(1,3,3,2) = -0.5d0 * (t97*t111 + t122)
      d(1,3,3,3) = -0.5d0 * (t97*t124 + t134)
      d(2,3,3,1) = -0.5d0 * (t116*t135 + t137)
      d(2,3,3,2) = t138 - 0.5d0*t139
      d(2,3,3,3) = -0.5d0 * (t116*t124 + t141)
      d(3,3,3,1) = -0.5d0 * (t129*t135 + t143)
      d(3,3,3,2) = -0.5d0 * (t129*t111 + t145)
      d(3,3,3,3) = t146 - 0.5d0*t147
      d(1,1,2,1) = t187 - 0.5d0*t191*t207
      d(1,1,2,2) = t225 - 0.5d0*t191*t240
      d(1,1,2,3) = t256 - 0.5d0*t191*t271
      d(2,1,2,1) = t273 - 0.5d0*t274*t207
      d(2,1,2,2) = t276 - 0.5d0*t274*t240
      d(2,1,2,3) = t278 - 0.5d0*t274*t271
      d(3,1,2,1) = t280 - 0.5d0*t281*t207
      d(3,1,2,2) = t283 - 0.5d0*t281*t240
      d(3,1,2,3) = t285 - 0.5d0*t281*t271
      d(1,1,3,1) = t311 - 0.5d0*t191*t330
      d(1,1,3,2) = t357 - 0.5d0*t191*t373
      d(1,1,3,3) = t398 - 0.5d0*t191*t414
      d(2,1,3,1) = t416 - 0.5d0*t274*t330
      d(2,1,3,2) = t418 - 0.5d0*t274*t373
      d(2,1,3,3) = t420 - 0.5d0*t274*t414
      d(3,1,3,1) = t422 - 0.5d0*t281*t330
      d(3,1,3,2) = t424 - 0.5d0*t281*t373
      d(3,1,3,3) = t426 - 0.5d0*t281*t414
      d(1,2,2,1) = t280*t162 - t273*t164
     &                - 0.5d0*(t281*t429 + t432*t151 + t431*t82)
     &                + 0.5d0*(t274*t436 + t439*t151 + t438*t89)
      d(1,2,2,2) = t283*t162 + t431*t83 - t276*t164
     &                - 0.5d0*(t281*t443 + t431*t84)
     &                + 0.5d0*(t274*t448 + t439*t210 + t438*t91)
      d(1,2,2,3) = t285*t162 - t278*t164 - t438*t92
     &                - 0.5d0*(t281*t453 + t432*t242 + t431*t86)
     &                + 0.5d0*(t274*t458 + t438*t93)
      d(2,2,2,1) = t187*t164 - t280*t160 - t431*t41
     &                - 0.5d0*(t191*t436 + t465*t151 + t464*t89)
     &                + 0.5d0*(t281*t469 + t431*t53)
      d(2,2,2,2) = t225*t164 - t283*t160
     &                - 0.5d0*(t191*t448 + t465*t210 + t464*t91)
     &                + 0.5d0*(t281*t478 + t480*t210 + t431*t66)
      d(2,2,2,3) = t256*t164 + t464*t92 - t285*t160
     &                - 0.5d0*(t191*t458 + t464*t93)
     &                + 0.5d0*(t281*t488 + t480*t242 + t431*t78)
      d(3,2,2,1) = t273*t160 + t438*t41 - t187*t162
     &                - 0.5d0*(t274*t469 + t438*t53)
     &                + 0.5d0*(t191*t429 + t498*t151 + t464*t82)
      d(3,2,2,2) = t276*t160 - t225*t162 - t464*t83
     &                - 0.5d0*(t274*t478 + t503*t210 + t438*t66)
     &                + 0.5d0*(t191*t443 + t464*t84)
      d(3,2,2,3) = t278*t160 - t256*t162
     &                - 0.5d0*(t274*t488 + t503*t242 + t438*t78)
     &                + 0.5d0*(t191*t453 + t498*t242 + t464*t86)
      d(1,2,3,1) = t422*t162 - t416*t164
     &                - 0.5d0*(t281*t519 + t521*t317 + t431*t137)
     &                + 0.5d0*(t274*t525 + t527*t317 + t438*t143)
      d(1,2,3,2) = t424*t162 + t431*t138 - t418*t164
     &                - 0.5d0*(t281*t531 + t431*t139)
     &                + 0.5d0*(t274*t536 + t527*t352 + t438*t145)
      d(1,2,3,3) = t426*t162 - t420*t164 - t438*t146
     &                - 0.5d0*(t281*t541 + t521*t393 + t431*t141)
     &                + 0.5d0*(t274*t546 + t438*t147)
      d(2,2,3,1) = t311*t164 - t422*t160 - t431*t101
     &                - 0.5d0*(t191*t525 + t552*t317 + t464*t143)
     &                + 0.5d0*(t281*t556 + t431*t109)
      d(2,2,3,2) = t357*t164 - t424*t160
     &                - 0.5d0*(t191*t536 + t552*t352 + t464*t145)
     &                + 0.5d0*(t281*t565 + t567*t352 + t431*t122)
      d(2,2,3,3) = t398*t164 + t464*t146 - t426*t160
     &                - 0.5d0*(t191*t546 + t464*t147)
     &                + 0.5d0*(t281*t575 + t567*t393 + t431*t134)
      d(3,2,3,1) = t416*t160 + t438*t101 - t311*t162
     &                - 0.5d0*(t274*t556 + t438*t109)
     &                + 0.5d0*(t191*t519 + t585*t317 + t464*t137)
      d(3,2,3,2) = t418*t160 - t357*t162 - t464*t138
     &                - 0.5d0*(t274*t565 + t590*t352 + t438*t122)
     &                + 0.5d0*(t191*t531 + t464*t139)
      d(3,2,3,3) = t420*t160 - t398*t162
     &                - 0.5d0*(t274*t575 + t590*t393 + t438*t134)
     &                + 0.5d0*(t191*t541 + t585*t393 + t464*t141)
c
c     some derivative matrix values are combinations of others
c
      do j = 1, 3
         do k = 1, 3
            do m = 1, 3
               d(m,k,1,j) = -d(m,k,2,j) - d(m,k,3,j)
            end do
         end do
      end do
      return
      end
