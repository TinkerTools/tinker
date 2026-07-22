c
c
c     ######################################################
c     ##  COPYRIGHT (C) 2023 by Zhi Wang & Jay W. Ponder  ##
c     ##                All Rights Reserved               ##
c     ######################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine exfield  --  external electric field energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "exfield" calculates the electrostatic energy due to an
c     external electric field applied to the system
c
c
      subroutine exfield (mode,exf)
      use atoms
      use bound
      use charge
      use chgpot
      use energi
      use extfld
      use mpole
      use usage
      implicit none
      integer i,ii
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     calculate external field energy over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf)
!$OMP    DO reduction(+:exf)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate external field energy over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,f,rpole,texfld,exf)
!$OMP    DO reduction(+:exf)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exfield1  --  external field energy & gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exfield1" calculates the electrostatic energy, gradient and
c     virial due to an external electric field applied to the system
c
c
      subroutine exfield1 (mode,exf)
      use atoms
      use bound
      use charge
      use chgpot
      use deriv
      use energi
      use extfld
      use mpole
      use usage
      use virial
      implicit none
      integer i,ii
      integer ix,iy,iz
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 frx,fry,frz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 fix(3),fiy(3)
      real*8 fiz(3),tem(3)
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     calculate energy and derivatives over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf,dec,vir)
!$OMP    DO reduction(+:exf,dec,vir)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
c
c     gradient and virial components from charge interactions
c
               frx = -f * texfld(1) * ci
               fry = -f * texfld(2) * ci
               frz = -f * texfld(3) * ci
               dec(1,i) = dec(1,i) + frx
               dec(2,i) = dec(2,i) + fry
               dec(3,i) = dec(3,i) + frz
               vxx = xi * frx
               vyy = yi * fry
               vzz = zi * frz
               vxy = 0.5d0 * (yi*frx+xi*fry)
               vxz = 0.5d0 * (zi*frx+xi*frz)
               vyz = 0.5d0 * (zi*fry+yi*frz)
c
c     increment the total internal virial tensor components
c
               vir(1,1) = vir(1,1) + vxx
               vir(2,1) = vir(2,1) + vxy
               vir(3,1) = vir(3,1) + vxz
               vir(1,2) = vir(1,2) + vxy
               vir(2,2) = vir(2,2) + vyy
               vir(3,2) = vir(3,2) + vyz
               vir(1,3) = vir(1,3) + vxz
               vir(2,3) = vir(2,3) + vyz
               vir(3,3) = vir(3,3) + vzz
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate energy and derivatives over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,xaxis,yaxis,zaxis,f,rpole,texfld,exf,dem,emvir)
!$OMP    DO reduction(+:exf,dem,emvir)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
c
c     gradient and virial components from dipole interactions
c
               tem(1) = f * (diy*texfld(3)-diz*texfld(2))
               tem(2) = f * (diz*texfld(1)-dix*texfld(3))
               tem(3) = f * (dix*texfld(2)-diy*texfld(1))
               call torque (i,tem,fix,fiy,fiz,dem)
               iz = zaxis(i)
               ix = xaxis(i)
               iy = abs(yaxis(i))
               if (iz .eq. 0)  iz = i
               if (ix .eq. 0)  ix = i
               if (iy .eq. 0)  iy = i
               xiz = x(iz) - x(i)
               yiz = y(iz) - y(i)
               ziz = z(iz) - z(i)
               xix = x(ix) - x(i)
               yix = y(ix) - y(i)
               zix = z(ix) - z(i)
               xiy = x(iy) - x(i)
               yiy = y(iy) - y(i)
               ziy = z(iy) - z(i)
               vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
               vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                           + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
               vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                           + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
               vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
               vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                           + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
               vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
c
c     gradient and virial components from monopole interactions
c
               frx = -f * texfld(1) * ci
               fry = -f * texfld(2) * ci
               frz = -f * texfld(3) * ci
               dem(1,i) = dem(1,i) + frx
               dem(2,i) = dem(2,i) + fry
               dem(3,i) = dem(3,i) + frz
               vxx = vxx + xi*frx
               vyy = vyy + yi*fry
               vzz = vzz + zi*frz
               vxy = vxy + 0.5d0*(yi*frx+xi*fry)
               vxz = vxz + 0.5d0*(zi*frx+xi*frz)
               vyz = vyz + 0.5d0*(zi*fry+yi*frz)
c
c     increment the total internal virial tensor components
c
               emvir(1,1) = emvir(1,1) + vxx
               emvir(2,1) = emvir(2,1) + vxy
               emvir(3,1) = emvir(3,1) + vxz
               emvir(1,2) = emvir(1,2) + vxy
               emvir(2,2) = emvir(2,2) + vyy
               emvir(3,2) = emvir(3,2) + vyz
               emvir(1,3) = emvir(1,3) + vxz
               emvir(2,3) = emvir(2,3) + vyz
               emvir(3,3) = emvir(3,3) + vzz
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine exfield3  --  electric field energy & analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "exfield3" calculates the electrostatic energy and partitions
c     the energy among the atoms due to an external electric field
c     applied to the system
c
c
      subroutine exfield3 (mode,exf)
      use action
      use analyz
      use atoms
      use bound
      use charge
      use chgpot
      use energi
      use extfld
      use mpole
      use usage
      implicit none
      integer i,ii
      real*8 exf,e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      character*6 mode
c
c
c     zero out the external electric field energy
c
      exf = 0.0d0
      f = electric / dielec
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     calculate energy and partitioning over partial charges
c
      if (mode .eq. 'CHARGE') then
!$OMP    PARALLEL default(private) shared(nion,iion,use,
!$OMP&    x,y,z,f,pchg,texfld,exf,nec,aec)
!$OMP    DO reduction(+:exf,nec,aec)
         do ii = 1, nion
            i = iion(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = pchg(i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               e = -f * ci * phi
               exf = exf + e
               nec = nec + 1
               aec(i) = aec(i) + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     calculate energy and partitioning over atomic multipoles
c
      if (mode .eq. 'MPOLE') then
!$OMP    PARALLEL default(private) shared(npole,ipole,use,
!$OMP&    x,y,z,f,rpole,texfld,exf,nem,aem)
!$OMP    DO reduction(+:exf,nem,aem)
         do ii = 1, npole
            i = ipole(ii)
            if (use(i)) then
               xi = x(i)
               yi = y(i)
               zi = z(i)
               ci = rpole(1,i)
               phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               e = -f * (ci*phi + dix*texfld(1)
     &                      + diy*texfld(2) + diz*texfld(3))
               exf = exf + e
               nem = nem + 1
               aem(i) = aem(i) + e
            end if
         end do
!$OMP    END DO
!$OMP    END PARALLEL
      end if
c
c     save the external field energy for use by any caller that
c     combines several subsystem states
c
      exfe = exf
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine exfield4  --  external field lambda derivatives  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "exfield4" calculates the electrostatic energy, gradient, virial
c     and lambda derivatives due to an external electric field applied
c     to the system, for the single topology multipole routines
c
c
      subroutine exfield4
      use atoms
      use bound
      use chgpot
      use deriv
      use dlmda
      use energi
      use extfld
      use mpole
      use mutant
      use usage
      use virial
      implicit none
      integer i,ii,j
      integer ix,iy,iz
      real*8 e,f,phi
      real*8 xi,yi,zi
      real*8 ci,dix,diy,diz
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 frx,fry,frz
      real*8 vxx,vyy,vzz
      real*8 vxy,vxz,vyz
      real*8 scalelmda
      real*8 fix(3),fiy(3)
      real*8 fiz(3),tem(3)
      real*8 dlfix(3),dlfiy(3)
      real*8 dlfiz(3),stem(3)
      logical muti
c
c
      f = electric / dielec
c
c     maintain any periodic boundary conditions
c
      if (use_bounds)  call bounds
c
c     calculate energy, derivatives and lambda derivatives over
c     the atomic multipoles
c
!$OMP PARALLEL default(private) shared(npole,ipole,use,mut,elambda,
!$OMP& x,y,z,xaxis,yaxis,zaxis,f,rpole,texfld,em,dem,emvir,
!$OMP& demdl,dfmdl,demvirdl)
!$OMP DO reduction(+:em,dem,emvir,demdl,dfmdl,demvirdl)
      do ii = 1, npole
         i = ipole(ii)
         if (use(i)) then
            xi = x(i)
            yi = y(i)
            zi = z(i)
            ci = rpole(1,i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            phi = xi*texfld(1) + yi*texfld(2) + zi*texfld(3)
            e = -f * (ci*phi + dix*texfld(1)
     &                   + diy*texfld(2) + diz*texfld(3))
c
c     the mutated sites carry multipoles scaled by elambda
c
            muti = mut(i)
            scalelmda = 1.0d0
            if (muti)  scalelmda = elambda
c
c     unscaled force and torque from the external field
c
            frx = -f * texfld(1) * ci
            fry = -f * texfld(2) * ci
            frz = -f * texfld(3) * ci
            tem(1) = f * (diy*texfld(3)-diz*texfld(2))
            tem(2) = f * (diz*texfld(1)-dix*texfld(3))
            tem(3) = f * (dix*texfld(2)-diy*texfld(1))
c
c     get the local frame definition used by the virial
c
            iz = zaxis(i)
            ix = xaxis(i)
            iy = abs(yaxis(i))
            if (iz .eq. 0)  iz = i
            if (ix .eq. 0)  ix = i
            if (iy .eq. 0)  iy = i
            xiz = x(iz) - x(i)
            yiz = y(iz) - y(i)
            ziz = z(iz) - z(i)
            xix = x(ix) - x(i)
            yix = y(ix) - y(i)
            zix = z(ix) - z(i)
            xiy = x(iy) - x(i)
            yiy = y(iy) - y(i)
            ziy = z(iy) - z(i)
c
c     lambda derivative of the energy, gradient and virial, which
c     for a mutated site is the unscaled external field term
c
            if (muti) then
               demdl = demdl + e
               call torque (i,tem,dlfix,dlfiy,dlfiz,dfmdl)
               vxx = xix*dlfix(1) + xiy*dlfiy(1) + xiz*dlfiz(1)
               vxy = 0.5d0 * (yix*dlfix(1) + yiy*dlfiy(1) + yiz*dlfiz(1)
     &                  + xix*dlfix(2) + xiy*dlfiy(2) + xiz*dlfiz(2))
               vxz = 0.5d0 * (zix*dlfix(1) + ziy*dlfiy(1) + ziz*dlfiz(1)
     &                  + xix*dlfix(3) + xiy*dlfiy(3) + xiz*dlfiz(3))
               vyy = yix*dlfix(2) + yiy*dlfiy(2) + yiz*dlfiz(2)
               vyz = 0.5d0 * (zix*dlfix(2) + ziy*dlfiy(2) + ziz*dlfiz(2)
     &                  + yix*dlfix(3) + yiy*dlfiy(3) + yiz*dlfiz(3))
               vzz = zix*dlfix(3) + ziy*dlfiy(3) + ziz*dlfiz(3)
               dfmdl(1,i) = dfmdl(1,i) + frx
               dfmdl(2,i) = dfmdl(2,i) + fry
               dfmdl(3,i) = dfmdl(3,i) + frz
               vxx = vxx + xi*frx
               vyy = vyy + yi*fry
               vzz = vzz + zi*frz
               vxy = vxy + 0.5d0*(yi*frx+xi*fry)
               vxz = vxz + 0.5d0*(zi*frx+xi*frz)
               vyz = vyz + 0.5d0*(zi*fry+yi*frz)
               demvirdl(1,1) = demvirdl(1,1) + vxx
               demvirdl(2,1) = demvirdl(2,1) + vxy
               demvirdl(3,1) = demvirdl(3,1) + vxz
               demvirdl(1,2) = demvirdl(1,2) + vxy
               demvirdl(2,2) = demvirdl(2,2) + vyy
               demvirdl(3,2) = demvirdl(3,2) + vyz
               demvirdl(1,3) = demvirdl(1,3) + vxz
               demvirdl(2,3) = demvirdl(2,3) + vyz
               demvirdl(3,3) = demvirdl(3,3) + vzz
            end if
c
c     scale the energy, force and torque by lambda
c
            e = scalelmda * e
            frx = scalelmda * frx
            fry = scalelmda * fry
            frz = scalelmda * frz
            do j = 1, 3
               stem(j) = scalelmda * tem(j)
            end do
c
c     increment the energy, gradient and virial from the dipole
c     interactions with the external field
c
            em = em + e
            call torque (i,stem,fix,fiy,fiz,dem)
            vxx = xix*fix(1) + xiy*fiy(1) + xiz*fiz(1)
            vxy = 0.5d0 * (yix*fix(1) + yiy*fiy(1) + yiz*fiz(1)
     &                        + xix*fix(2) + xiy*fiy(2) + xiz*fiz(2))
            vxz = 0.5d0 * (zix*fix(1) + ziy*fiy(1) + ziz*fiz(1)
     &                        + xix*fix(3) + xiy*fiy(3) + xiz*fiz(3))
            vyy = yix*fix(2) + yiy*fiy(2) + yiz*fiz(2)
            vyz = 0.5d0 * (zix*fix(2) + ziy*fiy(2) + ziz*fiz(2)
     &                        + yix*fix(3) + yiy*fiy(3) + yiz*fiz(3))
            vzz = zix*fix(3) + ziy*fiy(3) + ziz*fiz(3)
c
c     increment the gradient and virial from the monopole
c     interactions with the external field
c
            dem(1,i) = dem(1,i) + frx
            dem(2,i) = dem(2,i) + fry
            dem(3,i) = dem(3,i) + frz
            vxx = vxx + xi*frx
            vyy = vyy + yi*fry
            vzz = vzz + zi*frz
            vxy = vxy + 0.5d0*(yi*frx+xi*fry)
            vxz = vxz + 0.5d0*(zi*frx+xi*frz)
            vyz = vyz + 0.5d0*(zi*fry+yi*frz)
c
c     increment the multipole virial tensor components
c
            emvir(1,1) = emvir(1,1) + vxx
            emvir(2,1) = emvir(2,1) + vxy
            emvir(3,1) = emvir(3,1) + vxz
            emvir(1,2) = emvir(1,2) + vxy
            emvir(2,2) = emvir(2,2) + vyy
            emvir(3,2) = emvir(3,2) + vyz
            emvir(1,3) = emvir(1,3) + vxz
            emvir(2,3) = emvir(2,3) + vyz
            emvir(3,3) = emvir(3,3) + vzz
         end if
      end do
!$OMP END DO
!$OMP END PARALLEL
      return
      end
