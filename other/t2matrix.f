c
c
c     #############################################################
c     ##                                                         ##
c     ##  sizes.i  --  parameter values to set array dimensions  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sizes.i" sets values for critical array dimensions used
c     throughout the software; these parameters will fix the size
c     of the largest systems that can be handled; values too large
c     for the computer's memory and/or swap space to accomodate
c     will result in failure to operate
c
c     parameter:      maximum allowed number of:
c
c     maxatm          atoms in the molecular system
c     maxval          atoms directly bonded to an atom
c     maxgrp          user-defined groups of atoms
c     maxtyp          force field atom type definitions
c     maxclass        force field atom class definitions
c     maxkey          lines in the keyword file
c     maxrot          bonds for torsional rotation
c     maxpolar        mutually dipole polarizable sites
c     maxhess         off-diagonal Hessian elements
c     maxopt          full-matrix optimization variables
c     maxvib          vibrational frequencies
c     maxgeo          distance geometry points
c     maxpi           atoms in conjugated pisystem
c     maxcell         unit cells in crystal
c     maxring         3-, 4-, or 5-membered rings
c     maxfix          geometric restraints
c     maxbio          biopolymer atom definitions
c     maxres          residues in the macromolecule
c     maxamino        defined amino acid types
c     maxvar          Cartesian degrees of freedom
c     maxbnd          covalent bonds in molecular system
c     maxang          bond angles in molecular system
c     maxtors         dihedral angles in molecular system
c     maxlight        sites for method of lights neighbors
c     maxpib          covalent bonds in pisystem
c     maxpit          dihedrals involving pisystem
c
c
      integer maxatm,maxval,maxgrp,maxtyp,maxclass,maxkey,maxrot
      integer maxpolar,maxhess,maxopt,maxvib,maxgeo,maxpi,maxcell
      integer maxring,maxfix,maxbio,maxres,maxamino,maxvar,maxbnd
      integer maxang,maxtors,maxlight,maxpib,maxpit
      parameter (maxatm=3000)
      parameter (maxval=4)
      parameter (maxgrp=8)
      parameter (maxtyp=500)
      parameter (maxclass=400)
      parameter (maxkey=5000)
      parameter (maxrot=500)
      parameter (maxpolar=1200)
      parameter (maxhess=1000000)
      parameter (maxopt=500)
      parameter (maxvib=500)
      parameter (maxgeo=2000)
      parameter (maxpi=50)
      parameter (maxcell=500)
      parameter (maxring=500)
      parameter (maxfix=500)
      parameter (maxbio=500)
      parameter (maxres=500)
      parameter (maxamino=26)
      parameter (maxvar=3*maxatm)
      parameter (maxbnd=6*maxatm/5)
      parameter (maxang=12*maxatm/5)
      parameter (maxtors=4*maxatm)
      parameter (maxlight=8*maxatm)
      parameter (maxpib=6*maxpi/5)
      parameter (maxpit=4*maxpi)
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  mpole.i  --  multipole components for current structure  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxpole   max components (monopole=1,dipole=4,quadrupole=13)
c     maxrank   maximum rank needed (energy=2,gradient=3,hessian=4)
c     maxrank2  square of the maximum required rank value
c
c     pole      multipole values for each site in the local frame
c     rpole     multipoles rotated to the global coordinate system
c     dpole     derivative rotation matrix for each multipole
c     mdqsiz    number of mutipole components per atom to be used
c     npole     total number of multipole sites in the system
c     ipole     number of the atom for each multipole site
c     zaxis     number of the z-axis defining atom for each site
c     xaxis     number of the x-axis defining atom for each site
c     c         coefficients for Applequist-Dijkstra T2 elements
c     polaxe    local axis type for each multipole site
c
c
      integer maxpole,maxrank,maxrank2
      parameter (maxpole=13)
      parameter (maxrank=4)
      parameter (maxrank2=maxrank*maxrank)
      integer mdqsiz,npole,ipole,zaxis,xaxis,c
      real*8 pole,rpole,dpole
      character*8 polaxe
      common /mpole/ pole(maxpole,maxatm),rpole(maxpole,maxatm),
     &               dpole(maxpole,3,3,maxatm),mdqsiz,npole,
     &               ipole(maxatm),zaxis(maxatm),xaxis(maxatm),
     &               c(2*maxrank2+1,0:maxrank2,0:maxrank2),
     &               polaxe(maxatm)
c
c
c     ################################################################
c     ##                                                            ##
c     ##  polar.i  --  polarizabilities and induced dipole moments  ##
c     ##                                                            ##
c     ################################################################
c
c
c     polarize  atomic dipole polarizability for each atom (Ang**3)
c     pdamp     value of polarizability damping factor for each atom
c     uind      induced dipole components for each atom in the system
c     t2save    T2 matrix elements for induced dipole-multipole pairs
c
c
      real*8 polarize,pdamp,uind,t2save
      common /polar/ polarize(maxatm),pdamp(maxatm),uind(3,maxatm),
     &               t2save(5,maxpolar*(maxpolar-1)/2)
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine t2coeff  --  Applequist T2 matrix coefficients  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "t2coeff" is a service routine that sets coefficients needed
c     for the T2 matrix elements during calculation of multipole
c     interactions; needed to be called once from routine "kmpole"
c
c
      subroutine t2coeff
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer i,j,k,nn,kk
c
c
c     get the number of T2 coefficients to be generated
c     note maxrank for energy=2, gradient=3, hessian=4
c
      kk = maxrank + 2
      nn = 2*kk + 1
c
c     set coefficients used in the Applequist T2 matrix elements
c
      do k = 1, nn, 2
         c(k,0,0) = 1
         do j = 1, kk
            c(k,j,j) = c(k,j-1,j-1) * (-k-(j-1)*2)
         end do
      end do
      do k = 1, nn, 2
         do i = 2, kk, 2
            c(k,i,0) = c(k,i-1,1)
            do j = 1, kk-i
               c(k,j+i,j) = c(k,j+i-1,j-1)*(-k-i-(j-1)*2)
     &                         + c(k,j+i-1,j+1)*(j+1)
            end do
         end do
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  function t2matrix  --  get individual T2 matrix elements  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "t2matrix" computes the value of a single T2 matrix element
c     used in a polarizable multipole treatment of electrostatics
c     via the Applequist-Dykstra polytensor formalism
c
c
      function t2matrix (di,dj,dk,rr,bx,by)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      integer k,di,dj,dk
      real*8 t2matrix,rr
      real*8 bx(0:5,0:5),by(11,0:5,0:1)
c
c
c     get a single element of the Applequist-Dykstra T2 matrix
c
      t2matrix = 0.0d0
      do k = dk, 0, -2
         t2matrix = t2matrix + bx(dk,k)*by(k+dk+1,dj,di)
      end do
      t2matrix = t2matrix / rr
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empik  --  multipole & polarization pair energy  ##
c     ##                                                              ##
c     ##################################################################
c
c
      subroutine empik (ii,kk,r,xr,yr,zr,rpi,rpk,indi,indk,eik,ei,ek)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'electr.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer i,j,k,ii,kk
      real*8 eik,ei,ek
      real*8 t2matrix,factor,damp
      real*8 r,r2,r3,r4,r5
      real*8 xr,yr,zr,zrr
      real*8 xrr,xrr2,xrr3,xrr4
      real*8 yrr,yrr2,yrr3,yrr4
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 t2(13,13),tt2(4,13)
      real*8 fieldi(3),fieldk(3),m2t2(13)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
c
c
c     zeroth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 1) then
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     zeroth order T2 matrix elements
c
         t2(1,1) = t2matrix (0,0,0,r,bx,by)
c
c     first order coefficients to speed the T2 calculation
c
         r2 = r * r
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
c
c     first order T2 matrix elements
c
         t2(2,1) = t2matrix (0,0,1,r2,bx,by)
         t2(3,1) = t2matrix (0,1,0,r2,bx,by)
         t2(4,1) = t2matrix (1,0,0,r2,bx,by)
         t2(1,2) = -t2(2,1)
         t2(1,3) = -t2(3,1)
         t2(1,4) = -t2(4,1)
      end if
c
c     second order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 4) then
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
c
c     second order T2 matrix elements
c
         t2(2,2) = -t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = -t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = -t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = -t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = -t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
         k = 5
         do i = 2, 4
            do j = 2, 4
               t2(1,k) = -t2(i,j)
               t2(k,1) = t2(1,k)
               k = k + 1
            end do
         end do
      end if
c
c     third order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 13) then
         r4 = r2 * r2
         bz(7,0) = c(7,0,0)
         bz(5,1) = c(5,1,1) * zrr
         yrr3 = yrr2 * yrr
         by(7,0,0) = c(7,0,0)*bz(7,0)
         by(5,0,1) = c(5,0,0)*bz(5,1)
         by(5,1,0) = c(5,1,1)*bz(7,0)*yrr
         by(3,1,1) = c(3,1,1)*bz(5,1)*yrr
         by(3,2,0) = c(3,2,0)*bz(5,0) + c(3,2,2)*bz(7,0)*yrr2
         by(1,2,1) = c(1,2,0)*bz(3,1) + c(1,2,2)*bz(5,1)*yrr2
         by(1,3,0) = c(1,3,1)*bz(5,0)*yrr + c(1,3,3)*bz(7,0)*yrr3
         xrr3 = xrr2 * xrr
         bx(3,1) = c(1,3,1) * xrr
         bx(3,3) = c(1,3,3) * xrr3
c
c     third order T2 matrix elements
c
         t2(2,5) = t2matrix (0,0,3,r4,bx,by)
         t2(3,5) = t2matrix (0,1,2,r4,bx,by)
         t2(4,5) = t2matrix (1,0,2,r4,bx,by)
         t2(3,6) = t2matrix (0,2,1,r4,bx,by)
         t2(4,6) = t2matrix (1,1,1,r4,bx,by)
         t2(4,7) = -t2(2,5) - t2(3,6)
         t2(2,6) = t2(3,5)
         t2(2,7) = t2(4,5)
         t2(3,7) = t2(4,6)
         t2(2,8) = t2(2,6)
         t2(2,9) = t2(3,6)
         t2(2,10) = t2(4,6)
         t2(3,9) = t2matrix (0,3,0,r4,bx,by)
         t2(3,10) = t2matrix (1,2,0,r4,bx,by)
         t2(4,10) = -t2(2,8) - t2(3,9)
         t2(3,8) = t2(2,9)
         t2(4,8) = t2(2,10)
         t2(4,9) = t2(3,10)
         t2(2,11) = t2(2,7)
         t2(2,12) = t2(3,7)
         t2(2,13) = t2(4,7)
         t2(3,12) = t2(3,10)
         t2(3,13) = t2(4,10)
         t2(4,13) = -t2(2,11) - t2(3,12)
         t2(3,11) = t2(2,12)
         t2(4,11) = t2(2,13)
         t2(4,12) = t2(3,13)
         do i = 5, 13
            do j = 2, 4
               t2(i,j) = -t2(j,i)
            end do
         end do
c
c     fourth order coefficients to speed the T2 calculation
c
         r5 = r2 * r3
         bz(9,0) = c(9,0,0)
         bz(7,1) = c(7,1,1) * zrr
         yrr4 = yrr2 * yrr2
         by(9,0,0) = c(9,0,0)*bz(9,0)
         by(7,0,1) = c(7,0,0)*bz(7,1)
         by(7,1,0) = c(7,1,1)*bz(9,0)*yrr
         by(5,1,1) = c(5,1,1)*bz(7,1)*yrr
         by(5,2,0) = c(5,2,0)*bz(7,0) + c(5,2,2)*bz(9,0)*yrr2
         by(3,2,1) = c(3,2,0)*bz(5,1) + c(3,2,2)*bz(7,1)*yrr2
         by(3,3,0) = c(3,3,1)*bz(7,0)*yrr + c(3,3,3)*bz(9,0)*yrr3
         by(1,3,1) = c(1,3,1)*bz(5,1)*yrr + c(1,3,3)*bz(7,1)*yrr3
         by(1,4,0) = c(1,4,0)*bz(5,0) + c(1,4,2)*bz(7,0)*yrr2
     &                  + c(1,4,4)*bz(9,0)*yrr4
         xrr4 = xrr2 * xrr2
         bx(4,0) = c(1,4,0)
         bx(4,2) = c(1,4,2) * xrr2
         bx(4,4) = c(1,4,4) * xrr4
c
c     fourth order T2 matrix elements
c
         t2(5,5) = t2matrix (0,0,4,r5,bx,by)
         t2(6,5) = t2matrix (0,1,3,r5,bx,by)
         t2(7,5) = t2matrix (1,0,3,r5,bx,by)
         t2(6,6) = t2matrix (0,2,2,r5,bx,by)
         t2(7,6) = t2matrix (1,1,2,r5,bx,by)
         t2(7,7) = -t2(5,5) - t2(6,6)
         t2(8,5) = t2(6,5)
         t2(9,5) = t2(6,6)
         t2(10,5) = t2(7,6)
         t2(9,6) = t2matrix (0,3,1,r5,bx,by)
         t2(10,6) = t2matrix (1,2,1,r5,bx,by)
         t2(10,7) = -t2(8,5) - t2(9,6)
         t2(8,6) = t2(9,5)
         t2(8,7) = t2(10,5)
         t2(9,7) = t2(10,6)
         t2(11,5) = t2(7,5)
         t2(12,5) = t2(7,6)
         t2(13,5) = t2(7,7)
         t2(12,6) = t2(10,6)
         t2(13,6) = t2(10,7)
         t2(13,7) = -t2(11,5) - t2(12,6)
         t2(11,6) = t2(12,5)
         t2(11,7) = t2(13,5)
         t2(12,7) = t2(13,6)
         t2(8,8) = t2(8,6)
         t2(9,8) = t2(9,6)
         t2(10,8) = t2(10,6)
         t2(9,9) = t2matrix (0,4,0,r5,bx,by)
         t2(10,9) = t2matrix (1,3,0,r5,bx,by)
         t2(10,10) = -t2(8,8) - t2(9,9)
         t2(11,8) = t2(11,6)
         t2(12,8) = t2(12,6)
         t2(13,8) = t2(13,6)
         t2(12,9) = t2(10,9)
         t2(13,9) = t2(10,10)
         t2(13,10) = -t2(11,8) - t2(12,9)
         t2(11,9) = t2(12,8)
         t2(11,10) = t2(13,8)
         t2(12,10) = t2(13,9)
         t2(11,11) = t2(11,7)
         t2(12,11) = t2(12,7)
         t2(13,11) = t2(13,7)
         t2(12,12) = t2(12,10)
         t2(13,12) = t2(13,10)
         t2(13,13) = -t2(11,11) - t2(12,12)
         do i = 5, 12
            do j = i+1, 13
               t2(i,j) = t2(j,i)
            end do
         end do
      end if
c
c     compute interaction energy between the two multipole sites
c
      factor = electric / dielec
      do i = 1, mdqsiz
         m2t2(i) = 0.0d0
         do j = 1, mdqsiz
            m2t2(i) = m2t2(i) + rpk(j)*t2(j,i)
         end do
      end do
      eik = 0.0d0
      do i = 1, mdqsiz
         eik = eik + m2t2(i)*rpi(i)
      end do
      eik = factor * eik
c
c     compute dipole polarization energy at both polarizable sites
c
      do i = 2, 4
         do j = 1, mdqsiz
            tt2(i,j) = t2(j,i)
         end do
         do j = 2, 4
            tt2(i,j) = -tt2(i,j)
         end do
      end do
      do i = 1, 3
         fieldi(i) = 0.0d0
         fieldk(i) = 0.0d0
         do j = 1, mdqsiz
            fieldi(i) = fieldi(i) - rpk(j)*t2(j,i+1)
            fieldk(i) = fieldk(i) + tt2(i+1,j)*rpi(j)
         end do
      end do
      ei = 0.0d0
      ek = 0.0d0
      do i = 1, 3
         ei = ei + indi(i)*fieldi(i)
         ek = ek + indk(i)*fieldk(i)
      end do
      ei = -0.5d0 * factor * ei
      ek = -0.5d0 * factor * ek
c
c     apply an exponential damping factor to polarization energy
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            ei = ei * damp
            ek = ek * damp
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empik1  --  mpole & polarization pair gradient  ##
c     ##                                                             ##
c     #################################################################
c
c
      subroutine empik1 (ii,kk,r,xr,yr,zr,rpi,rpk,indi,indk,eik,ei,ek,
     &                      gt,gmi,gmj,d1ik,d1ki,d2ik,d2ki,d3,utu)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'electr.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'units.i'
      integer i,j,k,m
      integer ii,kk,polsiz
      real*8 eik,ei,ek
      real*8 t2matrix,factor
      real*8 damp,ddamp,de,term
      real*8 r,r2,r3,r4,r5,r6
      real*8 xr,yr,zr,zrr
      real*8 xrr,xrr2,xrr3,xrr4,xrr5
      real*8 yrr,yrr2,yrr3,yrr4,yrr5
      real*8 rpi(13),rpk(13)
      real*8 indi(3),indk(3)
      real*8 gt(3),gmi(3,3),gmj(3,3)
      real*8 t2(40,13),tt2(4,13),t2m1(13),m2t2(13)
      real*8 m2t3x(13),m2t3y(13),m2t3z(13)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
      real*8 d1ik(3),d1ki(3),d2ik(3,3),d2ki(3,3),d3(3),utu
      real*8 p2t2(13),fieldi(3),fieldk(3)
      real*8 interx(3),intery(3),interz(3),m1t(13)
      real*8 interx2(3),intery2(3),interz2(3),m1t2(13)
      real*8 t3x(3,13),t3y(3,13),t3z(3,13)
      real*8 t3x2(3,13),t3y2(3,13),t3z2(3,13)
c
c
c     zeroth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 1) then
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     zeroth order T2 matrix elements
c
         t2(1,1) = t2matrix (0,0,0,r,bx,by)
c
c     first order coefficients to speed the T2 calculation
c
         r2 = r * r
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
c
c     first order T2 matrix elements
c
         t2(2,1) = t2matrix (0,0,1,r2,bx,by)
         t2(3,1) = t2matrix (0,1,0,r2,bx,by)
         t2(4,1) = t2matrix (1,0,0,r2,bx,by)
         t2(1,2) = -t2(2,1)
         t2(1,3) = -t2(3,1)
         t2(1,4) = -t2(4,1)
c
c     second order coefficients to speed the T2 calculation
c
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
c
c     second order T2 matrix elements
c
         t2(2,2) = -t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = -t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = -t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = -t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = -t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
         k = 5
         do i = 2, 4
            do j = 2, 4
               t2(1,k) = -t2(i,j)
               t2(k,1) = t2(1,k)
               k = k + 1
            end do
         end do
c
c     third order coefficients to speed the T2 calculation
c
         r4 = r2 * r2
         bz(7,0) = c(7,0,0)
         bz(5,1) = c(5,1,1) * zrr
         yrr3 = yrr2 * yrr
         by(7,0,0) = c(7,0,0)*bz(7,0)
         by(5,0,1) = c(5,0,0)*bz(5,1)
         by(5,1,0) = c(5,1,1)*bz(7,0)*yrr
         by(3,1,1) = c(3,1,1)*bz(5,1)*yrr
         by(3,2,0) = c(3,2,0)*bz(5,0) + c(3,2,2)*bz(7,0)*yrr2
         by(1,2,1) = c(1,2,0)*bz(3,1) + c(1,2,2)*bz(5,1)*yrr2
         by(1,3,0) = c(1,3,1)*bz(5,0)*yrr + c(1,3,3)*bz(7,0)*yrr3
         xrr3 = xrr2 * xrr
         bx(3,1) = c(1,3,1) * xrr
         bx(3,3) = c(1,3,3) * xrr3
c
c     third order T2 matrix elements
c
         t2(2,5) = t2matrix (0,0,3,r4,bx,by)
         t2(3,5) = t2matrix (0,1,2,r4,bx,by)
         t2(4,5) = t2matrix (1,0,2,r4,bx,by)
         t2(3,6) = t2matrix (0,2,1,r4,bx,by)
         t2(4,6) = t2matrix (1,1,1,r4,bx,by)
         t2(4,7) = -t2(2,5) - t2(3,6)
         t2(2,6) = t2(3,5)
         t2(2,7) = t2(4,5)
         t2(3,7) = t2(4,6)
         t2(2,8) = t2(2,6)
         t2(2,9) = t2(3,6)
         t2(2,10) = t2(4,6)
         t2(3,9) = t2matrix (0,3,0,r4,bx,by)
         t2(3,10) = t2matrix (1,2,0,r4,bx,by)
         t2(4,10) = -t2(2,8) - t2(3,9)
         t2(3,8) = t2(2,9)
         t2(4,8) = t2(2,10)
         t2(4,9) = t2(3,10)
         t2(2,11) = t2(2,7)
         t2(2,12) = t2(3,7)
         t2(2,13) = t2(4,7)
         t2(3,12) = t2(3,10)
         t2(3,13) = t2(4,10)
         t2(4,13) = -t2(2,11) - t2(3,12)
         t2(3,11) = t2(2,12)
         t2(4,11) = t2(2,13)
         t2(4,12) = t2(3,13)
         do i = 5, 13
            do j = 2, 4
               t2(i,j) = -t2(j,i)
            end do
         end do
      end if
c
c     fourth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 13) then
         r5 = r2 * r3
         bz(9,0) = c(9,0,0)
         bz(7,1) = c(7,1,1) * zrr
         yrr4 = yrr2 * yrr2
         by(9,0,0) = c(9,0,0)*bz(9,0)
         by(7,0,1) = c(7,0,0)*bz(7,1)
         by(7,1,0) = c(7,1,1)*bz(9,0)*yrr
         by(5,1,1) = c(5,1,1)*bz(7,1)*yrr
         by(5,2,0) = c(5,2,0)*bz(7,0) + c(5,2,2)*bz(9,0)*yrr2
         by(3,2,1) = c(3,2,0)*bz(5,1) + c(3,2,2)*bz(7,1)*yrr2
         by(3,3,0) = c(3,3,1)*bz(7,0)*yrr + c(3,3,3)*bz(9,0)*yrr3
         by(1,3,1) = c(1,3,1)*bz(5,1)*yrr + c(1,3,3)*bz(7,1)*yrr3
         by(1,4,0) = c(1,4,0)*bz(5,0) + c(1,4,2)*bz(7,0)*yrr2
     &                  + c(1,4,4)*bz(9,0)*yrr4
         xrr4 = xrr2 * xrr2
         bx(4,0) = c(1,4,0)
         bx(4,2) = c(1,4,2) * xrr2
         bx(4,4) = c(1,4,4) * xrr4
c
c     fourth order T2 matrix elements
c
         t2(5,5) = t2matrix (0,0,4,r5,bx,by)
         t2(6,5) = t2matrix (0,1,3,r5,bx,by)
         t2(7,5) = t2matrix (1,0,3,r5,bx,by)
         t2(6,6) = t2matrix (0,2,2,r5,bx,by)
         t2(7,6) = t2matrix (1,1,2,r5,bx,by)
         t2(7,7) = -t2(5,5) - t2(6,6)
         t2(8,5) = t2(6,5)
         t2(9,5) = t2(6,6)
         t2(10,5) = t2(7,6)
         t2(9,6) = t2matrix (0,3,1,r5,bx,by)
         t2(10,6) = t2matrix (1,2,1,r5,bx,by)
         t2(10,7) = -t2(8,5) - t2(9,6)
         t2(8,6) = t2(9,5)
         t2(8,7) = t2(10,5)
         t2(9,7) = t2(10,6)
         t2(11,5) = t2(7,5)
         t2(12,5) = t2(7,6)
         t2(13,5) = t2(7,7)
         t2(12,6) = t2(10,6)
         t2(13,6) = t2(10,7)
         t2(13,7) = -t2(11,5) - t2(12,6)
         t2(11,6) = t2(12,5)
         t2(11,7) = t2(13,5)
         t2(12,7) = t2(13,6)
         t2(8,8) = t2(8,6)
         t2(9,8) = t2(9,6)
         t2(10,8) = t2(10,6)
         t2(9,9) = t2matrix (0,4,0,r5,bx,by)
         t2(10,9) = t2matrix (1,3,0,r5,bx,by)
         t2(10,10) = -t2(8,8) - t2(9,9)
         t2(11,8) = t2(11,6)
         t2(12,8) = t2(12,6)
         t2(13,8) = t2(13,6)
         t2(12,9) = t2(10,9)
         t2(13,9) = t2(10,10)
         t2(13,10) = -t2(11,8) - t2(12,9)
         t2(11,9) = t2(12,8)
         t2(11,10) = t2(13,8)
         t2(12,10) = t2(13,9)
         t2(11,11) = t2(11,7)
         t2(12,11) = t2(12,7)
         t2(13,11) = t2(13,7)
         t2(12,12) = t2(12,10)
         t2(13,12) = t2(13,10)
         t2(13,13) = -t2(11,11) - t2(12,12)
         do i = 5, 12
            do j = i+1, 13
               t2(i,j) = t2(j,i)
            end do
         end do
c
c     fifth order coefficients to speed the T2 calculation
c
         r6 = r3 * r3
         bz(11,0) = c(11,0,0)
         bz(9,1) = c(9,1,1) * zrr
         yrr5 = yrr3 * yrr2
         by(11,0,0) = c(11,0,0)*bz(11,0)
         by(9,0,1) = c(9,0,0)*bz(9,1)
         by(9,1,0) = c(9,1,1)*bz(11,0)*yrr
         by(7,1,1) = c(7,1,1)*bz(9,1)*yrr
         by(7,2,0) = c(7,2,0)*bz(9,0) + c(7,2,2)*bz(11,0)*yrr2
         by(5,2,1) = c(5,2,0)*bz(7,1) + c(5,2,2)*bz(9,1)*yrr2
         by(5,3,0) = c(5,3,1)*bz(9,0)*yrr + c(5,3,3)*bz(11,0)*yrr3
         by(3,3,1) = c(3,3,1)*bz(7,1)*yrr + c(3,3,3)*bz(9,1)*yrr3
         by(3,4,0) = c(3,4,0)*bz(7,0) + c(3,4,2)*bz(9,0)*yrr2
     &                  + c(3,4,4)*bz(11,0)*yrr4
         by(1,4,1) = c(1,4,0)*bz(5,1) + c(1,4,2)*bz(7,1)*yrr2
     &                  + c(1,4,4)*bz(9,1)*yrr4
         by(1,5,0) = c(1,5,1)*bz(7,0)*yrr + c(1,5,3)*bz(9,0)*yrr3
     &                  + c(1,5,5)*bz(11,0)*yrr5
         xrr5 = xrr3 * xrr2
         bx(5,1) = c(1,5,1) * xrr
         bx(5,3) = c(1,5,3) * xrr3
         bx(5,5) = c(1,5,5) * xrr5
c
c     first block of fifth order T2 matrix elements
c
         k = 14
         do i = 2, 4
            do j = 5, 13
               t2(k,1) = t2(i,j)
               k = k + 1
            end do
         end do
         k = 14
         do i = 5, 13, 3
            do j = 5, 13
               t2(k,2) = -t2(i,j)
               t2(k,3) = -t2(i+1,j)
               t2(k,4) = -t2(i+2,j)
               k = k + 1
            end do
         end do
         t2(14,5) = t2matrix (0,0,5,r6,bx,by)
         t2(15,5) = t2matrix (0,1,4,r6,bx,by)
         t2(16,5) = t2matrix (1,0,4,r6,bx,by)
         t2(15,6) = t2matrix (0,2,3,r6,bx,by)
         t2(16,6) = t2matrix (1,1,3,r6,bx,by)
         t2(16,7) = -t2(14,5) - t2(15,6)
         t2(17,5) = t2(15,5)
         t2(18,5) = t2(15,6)
         t2(19,5) = t2(16,6)
         t2(18,6) = t2matrix (0,3,2,r6,bx,by)
         t2(19,6) = t2matrix (1,2,2,r6,bx,by)
         t2(19,7) = -t2(17,5) - t2(18,6)
         t2(17,6) = t2(18,5)
         t2(17,7) = t2(19,5)
         t2(18,7) = t2(19,6)
         t2(17,8) = t2(17,6)
         t2(18,8) = t2(18,6)
         t2(19,8) = t2(19,6)
         t2(18,9) = t2matrix (0,4,1,r6,bx,by)
         t2(19,9) = t2matrix (1,3,1,r6,bx,by)
         t2(19,10) = -t2(17,8) - t2(18,9)
         t2(20,5) = t2(16,5)
         t2(21,5) = t2(16,6)
         t2(22,5) = t2(16,7)
         t2(21,6) = t2(19,6)
         t2(22,6) = t2(19,7)
         t2(22,7) = -t2(20,5) - t2(21,6)
         t2(20,6) = t2(21,5)
         t2(20,7) = t2(22,5)
         t2(21,7) = t2(22,6)
         t2(20,8) = t2(17,7)
         t2(21,8) = t2(18,7)
         t2(22,8) = t2(19,7)
         t2(21,9) = t2(19,9)
         t2(22,9) = t2(19,10)
         t2(22,10) = -t2(20,8) - t2(21,9)
         t2(20,9) = t2(21,8)
         t2(20,10) = t2(22,8)
         t2(21,10) = t2(22,9)
         t2(20,11) = t2(20,7)
         t2(21,11) = t2(21,7)
         t2(22,11) = t2(22,7)
         t2(21,12) = t2(19,10)
         t2(22,12) = t2(22,10)
         t2(22,13) = -t2(20,11) - t2(21,12)
         do i = 14, 21
            do j = i-8, 13
               k = j + 9
               m = i - 9
               t2(i,j) = t2(k,m)
            end do
         end do
c
c     second block of fifth order T2 matrix elements
c
         t2(23,5) = t2(17,5)
         t2(24,5) = t2(18,5)
         t2(25,5) = t2(19,5)
         t2(24,6) = t2(18,6)
         t2(25,6) = t2(19,6)
         t2(25,7) = t2(19,7)
         t2(26,5) = t2(24,5)
         t2(27,5) = t2(24,6)
         t2(28,5) = t2(25,6)
         t2(27,6) = t2(18,9)
         t2(28,6) = t2(19,9)
         t2(28,7) = -t2(26,5) - t2(27,6)
         t2(26,6) = t2(27,5)
         t2(26,7) = t2(28,5)
         t2(27,7) = t2(28,6)
         t2(26,8) = t2(26,6)
         t2(27,8) = t2(27,6)
         t2(28,8) = t2(28,6)
         t2(27,9) = t2matrix (0,5,0,r6,bx,by)
         t2(28,9) = t2matrix (1,4,0,r6,bx,by)
         t2(28,10) = -t2(26,8) - t2(27,9)
         t2(29,5) = t2(25,5)
         t2(30,5) = t2(25,6)
         t2(31,5) = t2(25,7)
         t2(30,6) = t2(28,6)
         t2(31,6) = t2(28,7)
         t2(31,7) = -t2(29,5) - t2(30,6)
         t2(29,6) = t2(30,5)
         t2(29,7) = t2(31,5)
         t2(30,7) = t2(31,6)
         t2(29,8) = t2(26,7)
         t2(30,8) = t2(27,7)
         t2(31,8) = t2(28,7)
         t2(30,9) = t2(28,9)
         t2(31,9) = t2(28,10)
         t2(31,10) = -t2(29,8) - t2(30,9)
         t2(29,9) = t2(30,8)
         t2(29,10) = t2(31,8)
         t2(30,10) = t2(31,9)
         t2(29,11) = t2(29,7)
         t2(30,11) = t2(30,7)
         t2(31,11) = t2(31,7)
         t2(30,12) = t2(30,10)
         t2(31,12) = t2(31,10)
         t2(31,13) = -t2(29,11) - t2(30,12)
         do i = 23, 30
            do j = i-17, 13
               k = j + 18
               m = i - 18
               t2(i,j) = t2(k,m)
            end do
         end do
c
c     third block of fifth order T2 matrix elements
c
         do i = 32, 40
            do j = 5, 7
               k = i - 18
               m = j + 6
               t2(i,j) = t2(k,m)
            end do
         end do
         do i = 32, 40
            do j = 8, 10
               k = i - 9
               m = j + 3
               t2(i,j) = t2(k,m)
            end do
         end do
         do i = 32, 34
            do j = 11, 13
               k = i - 12
               t2(i,j) = t2(k,j)
            end do
         end do
         do i = 35, 37
            do j = 11, 13
               k = i - 6
               t2(i,j) = t2(k,j)
            end do
         end do
         t2(38,11) = t2(20,13)
         t2(39,11) = t2(21,13)
         t2(40,11) = t2(22,13)
         t2(39,12) = t2(39,10)
         t2(40,12) = t2(40,10)
         t2(40,13) = -t2(38,11) - t2(39,12)
         t2(38,12) = t2(39,11)
         t2(38,13) = t2(40,11)
         t2(39,13) = t2(40,12)
      end if
c
c     compute interaction energy between the two multipole sites
c
      factor = electric / dielec
      do i = 1, mdqsiz
         m2t2(i) = 0.0d0
         do j = 1, mdqsiz
            m2t2(i) = m2t2(i) + rpk(j)*t2(j,i)
         end do
      end do
      eik = 0.0d0
      do i = 1, mdqsiz
         eik = eik + m2t2(i)*rpi(i)
      end do
      eik = factor * eik
c
c     get the permanent mulitpole (dM2/dx)*T*M1 derivative terms
c
      do i = 1, mdqsiz
         t2m1(i) = 0.0d0
         do j = 1, mdqsiz
            t2m1(i) = t2m1(i) + t2(i,j)*rpi(j)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            gmi(j,i) = 0.0d0
            gmj(j,i) = 0.0d0
            do k = 1, mdqsiz
               gmi(j,i) = gmi(j,i) + m2t2(k)*dpole(k,j,i,ii)
               gmj(j,i) = gmj(j,i) + dpole(k,j,i,kk)*t2m1(k)
            end do
            gmi(j,i) = factor * gmi(j,i)
            gmj(j,i) = factor * gmj(j,i)
         end do
      end do
c
c     get the permanent mulitpole M2*(dT/dx)*M1 derivative terms
c
      do i = 1, mdqsiz
         m2t3x(i) = 0.0d0
         m2t3y(i) = 0.0d0
         m2t3z(i) = 0.0d0
         do j = 1, mdqsiz
            m2t3x(i) = m2t3x(i) + rpk(j)*t2(3*j-1,i)
            m2t3y(i) = m2t3y(i) + rpk(j)*t2(3*j,i)
            m2t3z(i) = m2t3z(i) + rpk(j)*t2(3*j+1,i)
         end do
      end do
      do i = 1, 3
         gt(i) = 0.0d0
      end do
      do i = 1, mdqsiz
         gt(1) = gt(1) + m2t3x(i)*rpi(i)
         gt(2) = gt(2) + m2t3y(i)*rpi(i)
         gt(3) = gt(3) + m2t3z(i)*rpi(i)
      end do
      do i = 1, 3
         gt(i) = factor * gt(i)
      end do
c
c     compute dipole polarization energy at both polarizable sites
c
      do i = 2, 4
         do j = 1, mdqsiz
            tt2(i,j) = t2(j,i)
         end do
         do j = 2, 4
            tt2(i,j) = -tt2(i,j)
         end do
      end do
      do i = 1, 3
         fieldi(i) = 0.0d0
         fieldk(i) = 0.0d0
         do j = 1, mdqsiz
            fieldi(i) = fieldi(i) - rpk(j)*t2(j,i+1)
            fieldk(i) = fieldk(i) + tt2(i+1,j)*rpi(j)
         end do
      end do
      ei = 0.0d0
      ek = 0.0d0
      do i = 1, 3
         ei = ei + indi(i)*fieldi(i)
         ek = ek + indk(i)*fieldk(i)
      end do
      ei = -0.5d0 * factor * ei
      ek = -0.5d0 * factor * ek
c
c     copy some additional T matrix elements to auxiliary areas
c
      polsiz = max(4,mdqsiz)
      k = 1
      m = 1
      do i = 5, 13
         do j = 1, polsiz
            if (k .eq. 1) then
               t3x(m,j) = t2(i,j)
               t3x2(m,j) = t2(j,i)
            end if
            if (k .eq. 2) then
               t3y(m,j) = t2(i,j)
               t3y2(m,j) = t2(j,i)
            end if
            if (k .eq. 3) then
               t3z(m,j) = t2(i,j)
               t3z2(m,j) = t2(j,i)
            end if
         end do
         k = k + 1
         if (k .eq. 4) then
            k = 1
            m = m + 1
         end if
      end do
c
c     compute the d1 components for the induced dipole gradient
c
      do i = 1, 3
         interx(i) = 0.0d0
         intery(i) = 0.0d0
         interz(i) = 0.0d0
         interx2(i) = 0.0d0
         intery2(i) = 0.0d0
         interz2(i) = 0.0d0
         do j = 1, mdqsiz
            interx(i) = interx(i) + t3x(i,j)*rpi(j)
            intery(i) = intery(i) + t3y(i,j)*rpi(j)
            interz(i) = interz(i) + t3z(i,j)*rpi(j)
            interx2(i) = interx2(i) + t3x2(i,j)*rpk(j)
            intery2(i) = intery2(i) + t3y2(i,j)*rpk(j)
            interz2(i) = interz2(i) + t3z2(i,j)*rpk(j)
         end do
      end do
      do i = 1, 3
         d1ik(i) = 0.0d0
         d1ki(i) = 0.0d0
      end do
      do i = 1, 3
         d1ik(1) = d1ik(1) + indk(i)*interx(i)
         d1ik(2) = d1ik(2) + indk(i)*intery(i)
         d1ik(3) = d1ik(3) + indk(i)*interz(i)
         d1ki(1) = d1ki(1) + indi(i)*interx2(i)
         d1ki(2) = d1ki(2) + indi(i)*intery2(i)
         d1ki(3) = d1ki(3) + indi(i)*interz2(i)
      end do
      do i = 1, 3
         d1ik(i) = factor * d1ik(i)
         d1ki(i) = factor * d1ki(i)
      end do
c
c     compute the d2 components for the induced dipole gradient
c
      do i = 1, mdqsiz
         m1t(i) = 0.0d0
         m1t2(i) = 0.0d0
         do j = 1, 3
            m1t(i) = m1t(i) + indk(j)*t2(j+1,i)
            m1t2(i) = m1t2(i) + indi(j)*t2(i,j+1)
         end do
      end do
      do i = 1, 3
         do j = 1, 3
            d2ik(i,j) = 0.0d0
            d2ki(i,j) = 0.0d0
            do k = 1, mdqsiz
               d2ik(i,j) = d2ik(i,j) + m1t(k)*dpole(k,i,j,ii)
               d2ki(i,j) = d2ki(i,j) + m1t2(k)*dpole(k,i,j,kk)
            end do
            d2ik(i,j) = factor * d2ik(i,j)
            d2ki(i,j) = factor * d2ki(i,j)
         end do
      end do
c
c     compute the mutual polarization induced dipole gradient terms
c
      do i = 1, 3
         d3(i) = 0.0d0
      end do
      utu = 0.0d0
      if (poltyp .eq. 'MUTUAL') then
         do i = 1, 3
            interx(i) = 0.0d0
            intery(i) = 0.0d0
            interz(i) = 0.0d0
            do j = 2, 4
               interx(i) = interx(i) + t3x(i,j)*indk(j-1)
               intery(i) = intery(i) + t3y(i,j)*indk(j-1)
               interz(i) = interz(i) + t3z(i,j)*indk(j-1)
            end do
         end do
         do i = 1, 3
            d3(1) = d3(1) + indi(i)*interx(i)
            d3(2) = d3(2) + indi(i)*intery(i)
            d3(3) = d3(3) + indi(i)*interz(i)
         end do
         do i = 1, 3
            d3(i) = factor * d3(i)
         end do
         do i = 2, 4
            p2t2(i) = 0.0d0
            do j = 2, 4
               p2t2(i) = p2t2(i) + indk(j-1)*t2(j,i)
            end do
         end do
         do i = 2, 4
            utu = utu + p2t2(i)*indi(i-1)
         end do
         utu = factor * utu
      end if
c
c     apply a damping factor to polarization energy and derivatives
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         term = -pgamma * (r/damp)**3
         if (term .gt. -50.0d0) then
            term = exp(term)
            ddamp = (3.0d0*pgamma*r2/damp**3) * term
            damp = 1.0d0 - term
            de = 2.0d0 * ei * ddamp
            ei = ei * damp
            d1ik(1) = d1ik(1)*damp + de*(xr/r)
            d1ik(2) = d1ik(2)*damp + de*(yr/r)
            d1ik(3) = d1ik(3)*damp + de*(zr/r)
            d2ik(1,1) = d2ik(1,1) * damp
            d2ik(1,2) = d2ik(1,2) * damp
            d2ik(1,3) = d2ik(1,3) * damp
            d2ik(2,1) = d2ik(2,1) * damp
            d2ik(2,2) = d2ik(2,2) * damp
            d2ik(2,3) = d2ik(2,3) * damp
            d2ik(3,1) = d2ik(3,1) * damp
            d2ik(3,2) = d2ik(3,2) * damp
            d2ik(3,3) = d2ik(3,3) * damp
            de = 2.0d0 * ek * ddamp
            ek = ek * damp
            d1ki(1) = d1ki(1)*damp - de*(xr/r)
            d1ki(2) = d1ki(2)*damp - de*(yr/r)
            d1ki(3) = d1ki(3)*damp - de*(zr/r)
            d2ki(1,1) = d2ki(1,1) * damp
            d2ki(1,2) = d2ki(1,2) * damp
            d2ki(1,3) = d2ki(1,3) * damp
            d2ki(2,1) = d2ki(2,1) * damp
            d2ki(2,2) = d2ki(2,2) * damp
            d2ki(2,3) = d2ki(2,3) * damp
            d2ki(3,1) = d2ki(3,1) * damp
            d2ki(3,2) = d2ki(3,2) * damp
            d2ki(3,3) = d2ki(3,3) * damp
            de = utu * ddamp
            d3(1) = d3(1)*damp + de*(xr/r)
            d3(2) = d3(2)*damp + de*(yr/r)
            d3(3) = d3(3)*damp + de*(zr/r)
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine t2direct  --  field for direct induced dipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "t2direct" evaluates the electric field at a polarizable atom
c     due to permanent atomic multipoles at a second atom, and vice
c     versa, for use in computation of direct induced dipole moments
c
c
      subroutine t2direct (ii,kk,xr,yr,zr,r,r2,rpi,rpk,
     &                         fieldi,fieldk,store)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,j,ii,kk
      real*8 r,r2,r3,r4,damp
      real*8 xr,yr,zr,zrr
      real*8 xrr,xrr2,xrr3
      real*8 yrr,yrr2,yrr3
      real*8 rpi(13),rpk(13)
      real*8 fieldi(3),fieldk(3)
      real*8 t2matrix,t2(13,4),tt2(4,13)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
      logical store
c
c
c     zeroth order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 1) then
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     first order coefficients to speed the T2 calculation
c
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
c
c     first order T2 matrix elements
c
         t2(1,2) = t2matrix (0,0,1,r2,bx,by)
         t2(1,3) = t2matrix (0,1,0,r2,bx,by)
         t2(1,4) = t2matrix (1,0,0,r2,bx,by)
      end if
c
c     second order coefficients to speed the T2 calculation
c
      if (mdqsiz.ge.4 .or. store) then
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
c
c     second order T2 matrix elements
c
         t2(2,2) = t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
      end if
c
c     third order coefficients to speed the T2 calculation
c
      if (mdqsiz .ge. 13) then
         r4 = r2 * r2
         bz(7,0) = c(7,0,0)
         bz(5,1) = c(5,1,1) * zrr
         yrr3 = yrr2 * yrr
         by(7,0,0) = c(7,0,0)*bz(7,0)
         by(5,0,1) = c(5,0,0)*bz(5,1)
         by(5,1,0) = c(5,1,1)*bz(7,0)*yrr
         by(3,1,1) = c(3,1,1)*bz(5,1)*yrr
         by(3,2,0) = c(3,2,0)*bz(5,0)
     &                  + c(3,2,2)*bz(7,0)*yrr2
         by(1,2,1) = c(1,2,0)*bz(3,1)
     &                  + c(1,2,2)*bz(5,1)*yrr2
         by(1,3,0) = c(1,3,1)*bz(5,0)*yrr
     &                  + c(1,3,3)*bz(7,0)*yrr3
         xrr3 = xrr2 * xrr
         bx(3,1) = c(1,3,1) * xrr
         bx(3,3) = c(1,3,3) * xrr3
c
c     third order T2 matrix elements
c
         t2(5,2) = t2matrix (0,0,3,r4,bx,by)
         t2(5,3) = t2matrix (0,1,2,r4,bx,by)
         t2(5,4) = t2matrix (1,0,2,r4,bx,by)
         t2(6,3) = t2matrix (0,2,1,r4,bx,by)
         t2(6,4) = t2matrix (1,1,1,r4,bx,by)
         t2(9,3) = t2matrix (0,3,0,r4,bx,by)
         t2(10,3) = t2matrix (1,2,0,r4,bx,by)
         t2(7,4) = -t2(5,2) - t2(6,3)
         t2(6,2) = t2(5,3)
         t2(7,2) = t2(5,4)
         t2(7,3) = t2(6,4)
         t2(8,2) = t2(6,2)
         t2(9,2) = t2(6,3)
         t2(10,2) = t2(6,4)
         t2(10,4) = -t2(8,2) - t2(9,3)
         t2(8,3) = t2(9,2)
         t2(8,4) = t2(10,2)
         t2(9,4) = t2(10,3)
         t2(11,2) = t2(7,2)
         t2(12,2) = t2(7,3)
         t2(13,2) = t2(7,4)
         t2(12,3) = t2(10,3)
         t2(13,3) = t2(10,4)
         t2(13,4) = -t2(11,2) - t2(12,3)
         t2(11,3) = t2(12,2)
         t2(11,4) = t2(13,2)
         t2(12,4) = t2(13,3)
      end if
c
c     store T2 elements for mutual induced dipole iterations
c
      if (store) then
         j = npole*(ii-1) - ii*(ii+1)/2 + kk
         t2save(1,j) = t2(2,2)
         t2save(2,j) = t2(3,2)
         t2save(3,j) = t2(4,2)
         t2save(4,j) = t2(3,3)
         t2save(5,j) = t2(4,3)
      end if
c
c     set the transpose elements for the inverse interaction
c
      do i = 2, 4
         do j = 1, mdqsiz
            tt2(i,j) = t2(j,i)
         end do
         do j = 2, 4
            tt2(i,j) = -tt2(i,j)
         end do
      end do
c
c     compute the field at atom i due to atom k and vice versa
c
      do i = 1, 3
         fieldi(i) = 0.0d0
         fieldk(i) = 0.0d0
         do j = 1, mdqsiz
            fieldi(i) = fieldi(i) + rpk(j)*t2(j,i+1)
            fieldk(i) = fieldk(i) - tt2(i+1,j)*rpi(j)
         end do
      end do
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine t2mutual  --  field for mutual induced dipoles  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "t2mutual" evaluates the electric field at a polarizable atom
c     due to the induced atomic dipoles at a second atom, and vice
c     versa, for use in computation of mutual induced dipole moments
c
c
      subroutine t2mutual (ii,kk,xr,yr,zr,r,r2,rpi,rpk,
     &                         fieldi,fieldk,store)
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      integer i,j,ii,kk
      real*8 r,r2,r3,damp
      real*8 xr,yr,zr,zrr
      real*8 xrr,xrr2,yrr,yrr2
      real*8 rpi(13),rpk(13)
      real*8 fieldi(3),fieldk(3)
      real*8 t2matrix,t2(4,4),tt2(4,4)
      real*8 bx(0:5,0:5),by(11,0:5,0:1),bz(11,0:1)
      logical store
c
c
c     set the T2 elements to previously stored values
c
      if (store) then
         j = npole*(ii-1) - ii*(ii+1)/2 + kk
         t2(2,2) = t2save(1,j)
         t2(3,2) = t2save(2,j)
         t2(4,2) = t2save(3,j)
         t2(3,3) = t2save(4,j)
         t2(4,3) = t2save(5,j)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
c
c     zeroth order coefficients to speed the T2 calculation
c
      else
         bz(1,0) = c(1,0,0)
         by(1,0,0) = c(1,0,0) * bz(1,0)
         bx(0,0) = c(1,0,0)
c
c     first order coefficients to speed the T2 calculation
c
         zrr = zr / r
         bz(3,0) = c(3,0,0)
         bz(1,1) = c(1,1,1) * zrr
         yrr = yr / r
         by(3,0,0) = c(3,0,0) * bz(3,0)
         by(1,0,1) = c(1,0,0) * bz(1,1)
         by(1,1,0) = c(1,1,1) * bz(3,0) * yrr
         xrr = xr / r
         bx(1,1) = c(1,1,1) * xrr
c
c     second order coefficients to speed the T2 calculation
c
         r3 = r2 * r
         bz(5,0) = c(5,0,0)
         bz(3,1) = c(3,1,1) * zrr
         yrr2 = yrr * yrr
         by(5,0,0) = c(5,0,0) * bz(5,0)
         by(3,0,1) = c(3,0,0) * bz(3,1)
         by(3,1,0) = c(3,1,1) * bz(5,0) * yrr
         by(1,1,1) = c(1,1,1) * bz(3,1) * yrr
         by(1,2,0) = c(1,2,0)*bz(3,0) + c(1,2,2)*bz(5,0)*yrr2
         xrr2 = xrr * xrr
         bx(2,0) = c(1,2,0)
         bx(2,2) = c(1,2,2) * xrr2
c
c     second order T2 matrix elements
c
         t2(2,2) = t2matrix (0,0,2,r3,bx,by)
         t2(3,2) = t2matrix (0,1,1,r3,bx,by)
         t2(4,2) = t2matrix (1,0,1,r3,bx,by)
         t2(3,3) = t2matrix (0,2,0,r3,bx,by)
         t2(4,3) = t2matrix (1,1,0,r3,bx,by)
         t2(4,4) = -t2(2,2) - t2(3,3)
         t2(2,3) = t2(3,2)
         t2(2,4) = t2(4,2)
         t2(3,4) = t2(4,3)
      end if
c
c     set the transpose elements for the inverse interaction
c
      do i = 2, 4
         do j = 2, 4
            tt2(i,j) = t2(j,i)
         end do
         do j = 2, 4
            tt2(i,j) = -tt2(i,j)
         end do
      end do
c
c     compute the field at atom i due to atom k and vice versa
c
      do i = 1, 3
         fieldi(i) = 0.0d0
         fieldk(i) = 0.0d0
         do j = 2, 4
            fieldi(i) = fieldi(i) + rpk(j)*t2(j,i+1)
            fieldk(i) = fieldk(i) - tt2(i+1,j)*rpi(j)
         end do
      end do
c
c     apply a damping factor to reduce the field at short range
c
      damp = pdamp(ii) * pdamp(kk)
      if (damp .ne. 0.0d0) then
         damp = -pgamma * (r/damp)**3
         if (damp .gt. -50.0d0) then
            damp = 1.0d0 - exp(damp)
            do i = 1, 3
               fieldi(i) = fieldi(i) * damp
               fieldk(i) = fieldk(i) * damp
            end do
         end if
      end if
      return
      end
