c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2010 by Chuanjie Wu and Jay William Ponder  ##
c     ##                    All Rights Reserved                     ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine evcorr  --  long range vdw energy correction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "evcorr" computes the long range van der Waals correction
c     to the energy via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
      subroutine evcorr (elrc)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'math.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,it,kt
      integer nstep,ndelta
      real*8 elrc,etot
      real*8 range,rdelta
      real*8 termi,termik
      real*8 e,eps
      real*8 offset,taper
      real*8 rv,rv2,rv6,rv7
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 p,p6,p12
      real*8 rho,tau,tau7
      real*8 expterm
      character*6 mode
c
c
c     zero out the long range van der Waals correction
c
      elrc = 0.0d0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 20
      range = 200.0d0
      ndelta = int(dble(nstep)*(range-cut))
      rdelta = (range-cut) / dble(ndelta)
      offset = cut - 0.5d0*rdelta
c
c     find the van der Waals energy via double loop search
c
      elrc = 0.0d0
      do i = 1, nvt
         it = ivt(i)
         termi = 2.0d0 * pi * dble(jvt(i))
         do k = 1, nvt
            kt = ivt(k)
            termik = termi * dble(jvt(k))
	    rv = radmin(kt,it) 
	    eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = 0.0d0
            do j = 1, ndelta
               r = offset + dble(j)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = 0.0d0
               if (vdwtyp .eq. 'LENNARD-JONES') then     
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
	          rho = r7 + ghal*rv7
	          tau = (dhal+1.0d0) / (r+dhal*rv)
                  tau7 = tau**7
	          e = eps * rv7 * tau7 * ((ghal+1.0d0)*rv7/rho-2.0d0)
               else if (vdwtyp.eq.'BUCKINGHAM' .or.
     &                  vdwtyp.eq.'MM3-HBOND') then
                  p = sqrt(rv2/r2)
                  p6 = rv6 / r6
                  expterm = abuck * exp(-bbuck/p)
                  e = eps * (expterm - cbuck*p6)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  e = e * (1.0d0-taper)
               end if
               etot = etot + e*rdelta*r2
            end do
            elrc = elrc + termik*etot
         end do
      end do
      elrc = elrc / volbox
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine evcorr1  --  long range vdw energy & virial  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "evcorr1" computes the long range van der Waals correction
c     to the energy and virial via numerical integration
c
c     literature reference:
c
c     M. P. Allen and D. J. Tildesley, "Computer Simulation of
c     Liquids", Oxford University Press, 1987, Section 2.8
c
c
      subroutine evcorr1 (elrc,vlrc)
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cutoff.i'
      include 'math.i'
      include 'shunt.i'
      include 'vdw.i'
      include 'vdwpot.i'
      integer i,j,k,it,kt
      integer nstep,ndelta
      real*8 elrc,vlrc
      real*8 etot,vtot
      real*8 range,rdelta
      real*8 termi,termik
      real*8 e,de,eps
      real*8 offset
      real*8 taper,dtaper
      real*8 rv,rv2,rv6,rv7
      real*8 r,r2,r3,r4
      real*8 r5,r6,r7
      real*8 p,p6,p12
      real*8 rho,tau,tau7
      real*8 dtau,gtau
      real*8 rvterm,expterm
      character*6 mode
c
c
c     zero out the long range van der Waals corrections
c
      elrc = 0.0d0
      vlrc = 0.0d0
c
c     only applicable if periodic boundaries are in use
c
      if (.not. use_bounds)  return
c
c     set the coefficients for the switching function
c
      mode = 'VDW'
      call switch (mode)
c
c     set number of steps and range for numerical integration
c
      nstep = 20
      range = 200.0d0
      ndelta = int(dble(nstep)*(range-cut))
      rdelta = (range-cut) / dble(ndelta)
      offset = cut - 0.5d0*rdelta
c
c     find the van der Waals energy via double loop search
c
      do i = 1, nvt
         it = ivt(i)
         termi = 2.0d0 * pi * dble(jvt(i))
         do k = 1, nvt
            kt = ivt(k)
            termik = termi * dble(jvt(k))
	    rv = radmin(kt,it) 
	    eps = epsilon(kt,it)
            rv2 = rv * rv
            rv6 = rv2 * rv2 * rv2
            rv7 = rv6 * rv
            etot = 0.0d0
            vtot = 0.0d0
            do j = 1, ndelta
               r = offset + dble(j)*rdelta
               r2 = r * r
               r3 = r2 * r
               r6 = r3 * r3
               r7 = r6 * r
               e = 0.0d0
               de = 0.0d0
               if (vdwtyp .eq. 'LENNARD-JONES') then     
                  p6 = rv6 / r6
                  p12 = p6 * p6
                  e = eps * (p12 - 2.0d0*p6)
                  de = eps * (p12-p6) * (-12.0d0/r)
               else if (vdwtyp .eq. 'BUFFERED-14-7') then
	          rho = r7 + ghal*rv7
	          tau = (dhal+1.0d0) / (r+dhal*rv)
                  tau7 = tau**7
                  dtau = tau / (dhal+1.0d0)
                  gtau = eps*tau7*r6*(ghal+1.0d0)*(rv7/rho)**2
	          e = eps * rv7 * tau7 * ((ghal+1.0d0)*rv7/rho-2.0d0)
                  de = -7.0d0 * (dtau*e+gtau)
               else if (vdwtyp.eq.'BUCKINGHAM' .or.
     &                  vdwtyp.eq.'MM3-HBOND') then
                  p = sqrt(rv2/r2)
                  p6 = rv6 / r6
                  rvterm = -bbuck / rv
                  expterm = abuck * exp(-bbuck/p)
                  e = eps * (expterm - cbuck*p6)
                  de = eps * (rvterm*expterm+6.0d0*cbuck*p6/r)
               end if
               if (r .lt. off) then
                  r4 = r2 * r2
                  r5 = r2 * r3
                  taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
                  dtaper = 5.0d0*c5*r4 + 4.0d0*c4*r3
     &                        + 3.0d0*c3*r2 + 2.0d0*c2*r + c1
                  de = de*(1.0d0-taper) - e*dtaper
                  e = e*(1.0d0-taper)
               end if
               etot = etot + e*rdelta*r2
               vtot = vtot + de*rdelta*r3
            end do
            elrc = elrc + termik*etot
            vlrc = vlrc + termik*vtot
         end do
      end do
      elrc = elrc / volbox
      vlrc = vlrc / (3.0d0*volbox)
      return
      end
