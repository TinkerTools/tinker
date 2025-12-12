c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine settle  --  SETTLE distance constraint method  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "settle" implements the SETTLE algorithm by correcting atomic
c     positions to maintain rigid three-site water models
c
c     literature reference:
c
c     S. Miyamoto and P. A. Kollman, "SETTLE: An Analytical Version
c     of the SHAKE and RATTLE Algorithm for Rigid Water Models",
c     Journal of Computational Chemistry, 13, 952-962 (1992)
c
c
      subroutine settle (dt,xold,yold,zold)
      use atomid
      use atoms
      use freeze
      use math
      use moldyn
      implicit none
      integer i,ia,ib,ic
      real*8 dt,mtot
      real*8 mia,mib,mic
      real*8 ra,rb,rc
      real*8 rc2,rhh2
      real*8 xcom,ycom,zcom
      real*8 xia0,yia0,zia0
      real*8 xib0,yib0,zib0
      real*8 xic0,yic0,zic0
      real*8 xb0,yb0,zb0
      real*8 xc0,yc0,zc0
      real*8 xa1,ya1,za1
      real*8 xb1,yb1,zb1
      real*8 xc1,yc1,zc1
      real*8 xaksx,yaksx,zaksx
      real*8 xaksy,yaksy,zaksy
      real*8 xaksz,yaksz,zaksz
      real*8 rot11,rot12,rot13
      real*8 rot21,rot22,rot23
      real*8 rot31,rot32,rot33
      real*8 norm,za1d
      real*8 xb0d,yb0d
      real*8 xc0d,yc0d
      real*8 xb1d,yb1d,zb1d
      real*8 xc1d,yc1d,zc1d
      real*8 sinphi,cosphi
      real*8 sinpsi,cospsi
      real*8 ya2d,xb2d,yb2d
      real*8 yc2d,xb2d2
      real*8 alpa,beta,gama
      real*8 al2be2
      real*8 sintheta
      real*8 costheta
      real*8 xa3d,ya3d,za3d
      real*8 xb3d,yb3d,zb3d
      real*8 xc3d,yc3d,zc3d
      real*8 xold(*)
      real*8 yold(*)
      real*8 zold(*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nwat,iwat,kwat,mass,
!$OMP& x,y,z,xold,yold,zold,dt)
!$OMP& shared(v)
!$OMP DO reduction(+:v)
c
c     get coordinates, masses and geometric parameters
c
      do i = 1, nwat
         ia = iwat(1,i)
         ib = iwat(2,i)
         ic = iwat(3,i)
         mia = mass(ia)
         mib = mass(ib)
         mic = mass(ic)
         mtot = mia + mib + mic
         ra = kwat(1,i) * cos(0.5d0*kwat(3,i)/radian)
         rb = ra * mia / mtot
         ra = ra - rb
         rc2 = kwat(2,i)
         rc = 0.5d0 * rc2
         rhh2 = rc2 * rc2
c
c     store current coordinates needed for velocity correction
c
         xia0 = x(ia)
         yia0 = y(ia)
         zia0 = z(ia)
         xib0 = x(ib)
         yib0 = y(ib)
         zib0 = z(ib)
         xic0 = x(ic)
         yic0 = y(ic)
         zic0 = z(ic)
c
c     find prior interatomic vectors and center of mass
c
         xb0 = xold(ib) - xold(ia)
         yb0 = yold(ib) - yold(ia)
         zb0 = zold(ib) - zold(ia)
         xc0 = xold(ic) - xold(ia)
         yc0 = yold(ic) - yold(ia)
         zc0 = zold(ic) - zold(ia)
         xcom = (mia*x(ia)+mib*x(ib)+mic*x(ic)) / mtot
         ycom = (mia*y(ia)+mib*y(ib)+mic*y(ic)) / mtot
         zcom = (mia*z(ia)+mib*z(ib)+mic*z(ic)) / mtot
         xa1 = x(ia) - xcom
         ya1 = y(ia) - ycom
         za1 = z(ia) - zcom
         xb1 = x(ib) - xcom
         yb1 = y(ib) - ycom
         zb1 = z(ib) - zcom
         xc1 = x(ic) - xcom
         yc1 = y(ic) - ycom
         zc1 = z(ic) - zcom
c
c     compute axes and rotation matrix for alternative frame
c
         xaksz = yb0*zc0 - zb0*yc0
         yaksz = zb0*xc0 - xb0*zc0
         zaksz = xb0*yc0 - yb0*xc0
         xaksx = ya1*zaksz - za1*yaksz
         yaksx = za1*xaksz - xa1*zaksz
         zaksx = xa1*yaksz - ya1*xaksz
         xaksy = yaksz*zaksx - zaksz*yaksx
         yaksy = zaksz*xaksx - xaksz*zaksx
         zaksy = xaksz*yaksx - yaksz*xaksx
         norm = sqrt(xaksx*xaksx + yaksx*yaksx + zaksx*zaksx)
         rot11 = xaksx / norm
         rot21 = yaksx / norm
         rot31 = zaksx / norm
         norm = sqrt(xaksy*xaksy + yaksy*yaksy + zaksy*zaksy)
         rot12 = xaksy / norm
         rot22 = yaksy / norm
         rot32 = zaksy / norm
         norm = sqrt(xaksz*xaksz + yaksz*yaksz + zaksz*zaksz)
         rot13 = xaksz / norm
         rot23 = yaksz / norm
         rot33 = zaksz / norm
         xb0d = rot11*xb0 + rot21*yb0 + rot31*zb0
         yb0d = rot12*xb0 + rot22*yb0 + rot32*zb0
         xc0d = rot11*xc0 + rot21*yc0 + rot31*zc0
         yc0d = rot12*xc0 + rot22*yc0 + rot32*zc0
         za1d = rot13*xa1 + rot23*ya1 + rot33*za1
         xb1d = rot11*xb1 + rot21*yb1 + rot31*zb1
         yb1d = rot12*xb1 + rot22*yb1 + rot32*zb1
         zb1d = rot13*xb1 + rot23*yb1 + rot33*zb1
         xc1d = rot11*xc1 + rot21*yc1 + rot31*zc1
         yc1d = rot12*xc1 + rot22*yc1 + rot32*zc1
         zc1d = rot13*xc1 + rot23*yc1 + rot33*zc1
c
c     transform via rotation about Y' (phi) and X' axis (psi)
c
         sinphi = za1d / ra
         cosphi = sqrt(1.0d0 - sinphi*sinphi)
         sinpsi = (zb1d-zc1d) / (rc2*cosphi)
         cospsi = sqrt(1.0d0 - sinpsi*sinpsi)
         ya2d = ra * cosphi
         xb2d = -rc * cospsi
         yb2d = -rb*cosphi - rc*sinpsi*sinphi
         yc2d = -rb*cosphi + rc*sinpsi*sinphi
         xb2d2 = xb2d * xb2d
         xb2d = -0.5d0 * sqrt(rhh2 - (yb2d-yc2d)*(yb2d-yc2d)
     &                          - (zb1d-zc1d)*(zb1d-zc1d))
c
c    transform via a rotation about the Z' axis (theta)
c
         alpa = (xb2d*(xb0d-xc0d) + yb0d*yb2d + yc0d*yc2d)
         beta = (xb2d*(yc0d-yb0d) + xb0d*yb2d + xc0d*yc2d)
         gama = xb0d*yb1d - xb1d*yb0d + xc0d*yc1d - xc1d*yc0d
         al2be2 = alpa*alpa + beta*beta
         sintheta = (alpa*gama - beta*sqrt(al2be2-gama*gama)) / al2be2
         costheta = sqrt(1.0d0 - sintheta*sintheta)
         xa3d = -ya2d * sintheta
         ya3d = ya2d * costheta
         za3d = za1d
         xb3d = xb2d*costheta - yb2d*sintheta
         yb3d = xb2d*sintheta + yb2d*costheta
         zb3d = zb1d
         xc3d = -xb2d*costheta - yc2d*sintheta
         yc3d = -xb2d*sintheta + yc2d*costheta
         zc3d = zc1d
c
c     translate and rotate back into original coordinate frame
c
         x(ia) = xcom + rot11*xa3d + rot12*ya3d + rot13*za3d
         y(ia) = ycom + rot21*xa3d + rot22*ya3d + rot23*za3d
         z(ia) = zcom + rot31*xa3d + rot32*ya3d + rot33*za3d
         x(ib) = xcom + rot11*xb3d + rot12*yb3d + rot13*zb3d
         y(ib) = ycom + rot21*xb3d + rot22*yb3d + rot23*zb3d
         z(ib) = zcom + rot31*xb3d + rot32*yb3d + rot33*zb3d
         x(ic) = xcom + rot11*xc3d + rot12*yc3d + rot13*zc3d
         y(ic) = ycom + rot21*xc3d + rot22*yc3d + rot23*zc3d
         z(ic) = zcom + rot31*xc3d + rot32*yc3d + rot33*zc3d
c
c     use velocity correction derived from position movement
c
         if (dt .ne. 0.0d0) then
            v(1,ia) = v(1,ia) + (x(ia)-xia0)/dt
            v(2,ia) = v(2,ia) + (y(ia)-yia0)/dt
            v(3,ia) = v(3,ia) + (z(ia)-zia0)/dt
            v(1,ib) = v(1,ib) + (x(ib)-xib0)/dt
            v(2,ib) = v(2,ib) + (y(ib)-yib0)/dt
            v(3,ib) = v(3,ib) + (z(ib)-zib0)/dt
            v(1,ic) = v(1,ic) + (x(ic)-xic0)/dt
            v(2,ic) = v(2,ic) + (y(ic)-yic0)/dt
            v(3,ic) = v(3,ic) + (z(ic)-zic0)/dt
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     set position and velocity for four-site water extra sites
c
      if (nwat4 .ne. 0) then
         call watfour (dt)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine settle2  --  SETTLE atom velocity corrections  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "settle2" implements the second portion of the SETTLE algorithm
c     by correcting the full-step velocities in order to maintain
c     rigid three-site water models
c
c
      subroutine settle2 (dt)
      use atomid
      use atoms
      use freeze
      use moldyn
      use units
      use virial
      implicit none
      integer i,ia,ib,ic
      real*8 dt,norm
      real*8 mo,mh,moh,mhh
      real*8 mohoh,mhmh
      real*8 momoh,mhmhh
      real*8 momhh,mohmoh
      real*8 xab,yab,zab
      real*8 xbc,ybc,zbc
      real*8 xca,yca,zca
      real*8 xvab,yvab,zvab
      real*8 xvbc,yvbc,zvbc
      real*8 xvca,yvca,zvca
      real*8 vabab,vbcbc,vcaca
      real*8 cosa,cosb,cosc
      real*8 abmc,bcma,camb
      real*8 tab,tbc,tca
      real*8 denom,vterm
      real*8 dvax,dvay,dvaz
      real*8 dvbx,dvby,dvbz
      real*8 dvcx,dvcy,dvcz
      real*8 tax,tay,taz
      real*8 tbx,tby,tbz
      real*8 tcx,tcy,tcz
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(n,nwat,iwat,mass,x,y,z,v,dt)
!$OMP& shared(vir)
!$OMP DO reduction(+:vir)
c
c     find atoms in the molecule and mass combinations
c
      do i = 1, nwat
         ia = iwat(1,i)
         ib = iwat(2,i)
         ic = iwat(3,i)
         mo = mass(ia)
         mh = 0.5d0 * (mass(ib)+mass(ic))
         moh = mo + mh
         mhh = mh + mh
         mohoh = moh + moh
         mhmh = mh * mh
         momoh = mo * moh
         mhmhh = mh * mhh
         momhh = mo * mhh
         mohmoh = moh * moh
c
c     determine the normalized interactomic vectors
c
         xab = x(ib) - x(ia)
         yab = y(ib) - y(ia)
         zab = z(ib) - z(ia)
         norm = sqrt(xab*xab + yab*yab + zab*zab)
         xab = xab / norm
         yab = yab / norm
         zab = zab / norm
         xbc = x(ic) - x(ib)
         ybc = y(ic) - y(ib)
         zbc = z(ic) - z(ib)
         norm = sqrt(xbc*xbc + ybc*ybc + zbc*zbc)
         xbc = xbc / norm
         ybc = ybc / norm
         zbc = zbc / norm
         xca = x(ia) - x(ic)
         yca = y(ia) - y(ic)
         zca = z(ia) - z(ic)
         norm = sqrt(xca*xca + yca*yca + zca*zca)
         xca = xca / norm
         yca = yca / norm
         zca = zca / norm
c
c     compute angle cosines between interatomic vectors
c
         cosa = -xab*xca - yab*yca - zab*zca
         cosb = -xbc*xab - ybc*yab - zbc*zab
         cosc = -xca*xbc - yca*ybc - zca*zbc
c
c     find velocity difference vectors between adjacent atoms
c
         xvab = v(1,ib) - v(1,ia)
         yvab = v(2,ib) - v(2,ia)
         zvab = v(3,ib) - v(3,ia)
         xvbc = v(1,ic) - v(1,ib)
         yvbc = v(2,ic) - v(2,ib)
         zvbc = v(3,ic) - v(3,ib)
         xvca = v(1,ia) - v(1,ic)
         yvca = v(2,ia) - v(2,ic)
         zvca = v(3,ia) - v(3,ic)
c
c     get dot product of interatomic and velocity vectors
c
         vabab = xvab*xab + yvab*yab + zvab*zab
         vbcbc = xvbc*xbc + yvbc*ybc + zvbc*zbc
         vcaca = xvca*xca + yvca*yca + zvca*zca
c
c    intermediates for velocity constraint corrections
c
         abmc = mh*cosa*cosb - moh*cosc
         bcma = mo*cosb*cosc - mhh*cosa
         camb = mh*cosc*cosa - moh*cosb
         tab = vabab*(mohoh-mo*cosc*cosc) + vbcbc*camb + vcaca*bcma
         tbc = vbcbc*(mohmoh-mhmh*cosa*cosa)
     &            + vcaca*abmc*mo + vabab*camb*mo
         tca = vcaca*(mohoh-mo*cosb*cosb) + vabab*bcma + vbcbc*abmc
         denom = 2.0d0*mohmoh + momhh*cosa*cosb*cosc
     &              - mhmhh*cosa*cosa - momoh*(cosb*cosb+cosc*cosc)
c
c     construct the velocity constraint correction components
c
         dvax = (xab*tab-xca*tca)*mh / denom
         dvay = (yab*tab-yca*tca)*mh / denom
         dvaz = (zab*tab-zca*tca)*mh / denom
         dvbx = (xbc*tbc-xab*tab*mo) / denom
         dvby = (ybc*tbc-yab*tab*mo) / denom
         dvbz = (zbc*tbc-zab*tab*mo) / denom
         dvcx = (xca*tca*mo-xbc*tbc) / denom
         dvcy = (yca*tca*mo-ybc*tbc) / denom
         dvcz = (zca*tca*mo-zbc*tbc) / denom
c
c     modify velocity components with constraint corrections
c
         v(1,ia) = v(1,ia) + dvax
         v(2,ia) = v(2,ia) + dvay
         v(3,ia) = v(3,ia) + dvaz
         v(1,ib) = v(1,ib) + dvbx
         v(2,ib) = v(2,ib) + dvby
         v(3,ib) = v(3,ib) + dvbz
         v(1,ic) = v(1,ic) + dvcx
         v(2,ic) = v(2,ic) + dvcy
         v(3,ic) = v(3,ic) + dvcz
c
c     increment the internal virial tensor components
c
         vterm = -2.0d0 / (dt*ekcal)
         tax = vterm * x(ia) * mo
         tay = vterm * y(ia) * mo
         taz = vterm * z(ia) * mo
         tbx = vterm * x(ib) * mh
         tby = vterm * y(ib) * mh
         tbz = vterm * z(ib) * mh
         tcx = vterm * x(ic) * mh
         tcy = vterm * y(ic) * mh
         tcz = vterm * z(ic) * mh
         vir(1,1) = vir(1,1) + tax*dvax + tbx*dvbx + tcx*dvcx
         vir(2,1) = vir(2,1) + tay*dvax + tby*dvbx + tcy*dvcx
         vir(3,1) = vir(3,1) + taz*dvax + tbz*dvbx + tcz*dvcx
         vir(1,2) = vir(1,2) + tax*dvay + tbx*dvby + tcx*dvcy
         vir(2,2) = vir(2,2) + tay*dvay + tby*dvby + tcy*dvcy
         vir(3,2) = vir(3,2) + taz*dvay + tbz*dvby + tcz*dvcy
         vir(1,3) = vir(1,3) + tax*dvaz + tbx*dvbz + tcx*dvcz
         vir(2,3) = vir(2,3) + tay*dvaz + tby*dvbz + tcy*dvcz
         vir(3,3) = vir(3,3) + taz*dvaz + tbz*dvbz + tcz*dvcz
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine settleg  --  SETTLE gradient vector correction  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "settleg" modifies the gradient to remove components along any
c     holonomic distance contraints for rigid three-site water using
c     a variant of the SETTLE algorithm
c
c
      subroutine settleg (derivs)
      use atoms
      use freeze
      implicit none
      integer i,j,ia,ib,ic
      real*8 r1x,r1y,r1z
      real*8 r2x,r2y,r2z
      real*8 r3x,r3y,r3z
      real*8 gxO,gyO,gzO
      real*8 gx1,gy1,gz1
      real*8 gx2,gy2,gz2
      real*8 c1o(3),c1h1(3),c1h2(3)
      real*8 c2o(3),c2h1(3),c2h2(3)
      real*8 c3o(3),c3h1(3),c3h2(3)
      real*8 a(3,3),b(3),p(3)
      real*8 derivs(3,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nwat,iwat,x,y,z,derivs)
!$OMP DO
c
c     determine the atom numbers and interactomic vectors
c
      do i = 1, nwat
         ia = iwat(1,i)
         ib = iwat(2,i)
         ic = iwat(3,i)
         r1x = x(ib) - x(ia)
         r1y = y(ib) - y(ia)
         r1z = z(ib) - z(ia)
         r2x = x(ic) - x(ia)
         r2y = y(ic) - y(ia)
         r2z = z(ic) - z(ia)
         r3x = x(ic) - x(ib)
         r3y = y(ic) - y(ib)
         r3z = z(ic) - z(ib)
c
c     construct the vectors for each of the constraints
c
         c1o(1) = -2.0d0 * r1x
         c1o(2) = -2.0d0 * r1y
         c1o(3) = -2.0d0 * r1z
         c1h1(1) = 2.0d0 * r1x
         c1h1(2) = 2.0d0 * r1y
         c1h1(3) = 2.0d0 * r1z
         c1h2(1) = 0.0d0
         c1h2(2) = 0.0d0
         c1h2(3) = 0.0d0
         c2o(1) = -2.0d0 * r2x
         c2o(2) = -2.0d0 * r2y
         c2o(3) = -2.0d0 * r2z
         c2h1(1) = 0.0d0
         c2h1(2) = 0.0d0
         c2h1(3) = 0.0d0
         c2h2(1) = 2.0d0 * r2x
         c2h2(2) = 2.0d0 * r2y
         c2h2(3) = 2.0d0 * r2z
         c3o(1) = 0.0d0
         c3o(2) = 0.0d0
         c3o(3) = 0.0d0
         c3h1(1) = -2.0d0 * r3x
         c3h1(2) = -2.0d0 * r3y
         c3h1(3) = -2.0d0 * r3z
         c3h2(1) = 2.0d0 * r3x
         c3h2(2) = 2.0d0 * r3y
         c3h2(3) = 2.0d0 * r3z
c
c     build the individual element of the constraint matrix
c
         a(1,1) = c1o(1)*c1o(1) + c1o(2)*c1o(2)
     &               + c1o(3)*c1o(3) + c1h1(1)*c1h1(1)
     &               + c1h1(2)*c1h1(2) + c1h1(3)*c1h1(3)
         a(1,2) = c1o(1)*c2o(1) + c1o(2)*c2o(2)
     &               + c1o(3)*c2o(3) + c1h1(1)*c2h1(1)
     &               + c1h1(2)*c2h1(2) + c1h1(3)*c2h1(3)
     &               + c1h2(1)*c2h2(1) + c1h2(2)*c2h2(2)
     &               + c1h2(3)*c2h2(3)
         a(1,3) = c1o(1)*c3o(1) + c1o(2)*c3o(2)
     &               + c1o(3)*c3o(3) + c1h1(1)*c3h1(1)
     &               + c1h1(2)*c3h1(2) + c1h1(3)*c3h1(3)
     &               + c1h2(1)*c3h2(1) + c1h2(2)*c3h2(2)
     &               + c1h2(3)*c3h2(3)
         a(2,1) = a(1,2)
         a(2,2) = c2o(1)*c2o(1) + c2o(2)*c2o(2)
     &               + c2o(3)*c2o(3) + c2h2(1)*c2h2(1)
     &               + c2h2(2)*c2h2(2) + c2h2(3)*c2h2(3)
         a(2,3) = c2o(1)*c3o(1) + c2o(2)*c3o(2)
     &               + c2o(3)*c3o(3) + c2h1(1)*c3h1(1)
     &               + c2h1(2)*c3h1(2) + c2h1(3)*c3h1(3)
     &               + c2h2(1)*c3h2(1) + c2h2(2)*c3h2(2)
     &               + c2h2(3)*c3h2(3)
         a(3,1) = a(1,3)
         a(3,2) = a(2,3)
         a(3,3) = c3h1(1)*c3h1(1) + c3h1(2)*c3h1(2)
     &               + c3h1(3)*c3h1(3) + c3h2(1)*c3h2(1)
     &               + c3h2(2)*c3h2(2) + c3h2(3)*c3h2(3)
c
c     copy the current gradient values into local variables
c
         gxO = -derivs(1,ia)
         gyO = -derivs(2,ia)
         gzO = -derivs(3,ia)
         gx1 = -derivs(1,ib)
         gy1 = -derivs(2,ib)
         gz1 = -derivs(3,ib)
         gx2 = -derivs(1,ic)
         gy2 = -derivs(2,ic)
         gz2 = -derivs(3,ic)
c
c     evaluate the constraint-gradient dot product sums
c
         b(1) = c1o(1)*gxO + c1o(2)*gyO + c1o(3)*gzO
     &             + c1h1(1)*gx1 + c1h1(2)*gy1 + c1h1(3)*gz1
     &             + c1h2(1)*gx2 + c1h2(2)*gy2 + c1h2(3)*gz2
         b(2) = c2o(1)*gxO + c2o(2)*gyO + c2o(3)*gzO
     &             + c2h1(1)*gx1 + c2h1(2)*gy1 + c2h1(3)*gz1
     &             + c2h2(1)*gx2 + c2h2(2)*gy2 + c2h2(3)*gz2
         b(3) = c3o(1)*gxO + c3o(2)*gyO + c3o(3)*gzO
     &             + c3h1(1)*gx1 + c3h1(2)*gy1 + c3h1(3)*gz1
     &             + c3h2(1)*gx2 + c3h2(2)*gy2 + c3h2(3)*gz2
c
c     use a 3x3 Gaussian elimination to solve A * p = b
c
         call solve3 (a,b,p)
c
c     gradient correction to remove projection on constraints
c
         derivs(1,ia) = derivs(1,ia) + p(1)*c1o(1)
     &                     + p(2)*c2o(1) + p(3)*c3o(1)
         derivs(2,ia) = derivs(2,ia) + p(1)*c1o(2)
     &                     + p(2)*c2o(2) + p(3)*c3o(2)
         derivs(3,ia) = derivs(3,ia) + p(1)*c1o(3)
     &                     + p(2)*c2o(3) + p(3)*c3o(3)
         derivs(1,ib) = derivs(1,ib) + p(1)*c1h1(1)
     &                     + p(2)*c2h1(1) + p(3)*c3h1(1)
         derivs(2,ib) = derivs(2,ib) + p(1)*c1h1(2)
     &                     + p(2)*c2h1(2) + p(3)*c3h1(2)
         derivs(3,ib) = derivs(3,ib) + p(1)*c1h1(3)
     &                     + p(2)*c2h1(3) + p(3)*c3h1(3)
         derivs(1,ic) = derivs(1,ic) + p(1)*c1h2(1)
     &                     + p(2)*c2h2(1) + p(3)*c3h2(1)
         derivs(2,ic) = derivs(2,ic) + p(1)*c1h2(2)
     &                     + p(2)*c2h2(2) + p(3)*c3h2(2)
         derivs(3,ic) = derivs(3,ic) + p(1)*c1h2(3)
     &                     + p(2)*c2h2(3) + p(3)*c3h2(3)

      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine solve3  --  3x3 Gaussian elimination with pivot  ## 
c     ##                                                              ##
c     ##################################################################
c
c
c     "solve3" uses 3x3 Gaussian elimination with pivoting to solve
c     A * p = b as a utility to find gradient corrections for rigid
c     three-site water models under the SETTLE algorithm
c
c     note the inputs "a" and "b" are overwritten, and must be saved
c     upon entering this routine if needed upon return; the solution
c     is returned in "p"
c
c
      subroutine solve3 (a,b,p)
      implicit none
      integer i,j,k
      integer pivot
      real*8 b(3)
      real*8 p(3)
      real*8 a(3,3)
      real*8 eps,amax
      real*8 swap,factor
c
c
c     use Gaussian elimination with partial pivoting
c
      eps = 0.00000001d0
      do k = 1, 2
         pivot = k
         amax = abs(a(k,k))
         do i = k+1, 3
            if (abs(a(i,k)) .gt. amax) then
               amax = abs(a(i,k))
               pivot = i
            end if
         end do
c
c     swap rows k and pivot in the "a" matrix and "b" vector
c
         if (pivot .ne. k) then
            do j = k, 3
               swap = a(k,j)
               a(k,j) = a(pivot,j)
               a(pivot,j) = swap
            end do
            swap = b(k)
            b(k) = b(pivot)
            b(pivot) = swap
         end if
c
c     if small pivot, then near singular and use zero solution
c
         if (abs(a(k,k)) .lt. eps) then
            p(1) = 0.0d0
            p(2) = 0.0d0
            p(3) = 0.0d0
            return
         end if
c
c     eliminate the row or rows below
c
         do i = k+1, 3
            factor = a(i,k) / a(k,k)
            a(i,k) = 0.0d0
            do j = k+1, 3
               a(i,j) = a(i,j) - factor*a(k,j)
            end do
            b(i) = b(i) - factor*b(k)
         end do
      end do
c
c     perform the back substitution phase
c
      if (abs(a(3,3)) .lt. eps) then
         p(1) = 0.0d0
         p(2) = 0.0d0
         p(3) = 0.0d0
      else
         p(3) = b(3) / a(3,3)
         p(2) = (b(2)-a(2,3)*p(3)) / a(2,2)
         p(1) = (b(1)-a(1,2)*p(2)-a(1,3)*p(3)) / a(1,1)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine watfour  --  set 4-site water extra position  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "watfour" sets the position and zeros the velocity of the
c     extra site of rigid planar four-site water molecules 
c
c
      subroutine watfour (dt)
      use atoms
      use freeze
      use moldyn
      implicit none
      integer i
      integer ia,ib,ic,id
      real*8 dt
      real*8 oterm,hterm
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nwat4,iwat4,kwat4,x,y,z,v,dt)
!$OMP DO
c
c     set the position of four-site water extra centers
c
      do i = 1, nwat4
         ia = iwat4(1,i)
         ib = iwat4(2,i)
         ic = iwat4(3,i)
         id = iwat4(4,i)
         oterm = kwat4(1,i)
         hterm = kwat4(2,i)
         x(id) = oterm*x(ia) + hterm*(x(ib)+x(ic))
         y(id) = oterm*y(ia) + hterm*(y(ib)+y(ic))
         z(id) = oterm*z(ia) + hterm*(z(ib)+z(ic))
         if (dt .ne. 0.0d0) then
            v(1,id) = 0.0d0
            v(2,id) = 0.0d0
            v(3,id) = 0.0d0
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine watfour2  --  distribute 4-site water gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "watfour2" transfers gradient components on the extra site of
c     rigid planar four-site water molecules to the H-O-H atoms
c
c
      subroutine watfour2 (derivs)
      use freeze
      implicit none
      integer i,j
      integer ia,ib,ic,id
      real*8 oterm,hterm
      real*8 derivs(3,*)
c
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nwat4,iwat4,kwat4,derivs)
!$OMP DO
c
c     distribute the gradient on four-site water extra centers
c
      do i = 1, nwat4
         ia = iwat4(1,i)
         ib = iwat4(2,i)
         ic = iwat4(3,i)
         id = iwat4(4,i)
         oterm = kwat4(1,i)
         hterm = kwat4(2,i)
         do j = 1, 3
            derivs(j,ia) = derivs(j,ia) + oterm*derivs(j,id)
            derivs(j,ib) = derivs(j,ib) + hterm*derivs(j,id)
            derivs(j,ic) = derivs(j,ic) + hterm*derivs(j,id)
            derivs(j,id) = 0.0d0
         end do
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
