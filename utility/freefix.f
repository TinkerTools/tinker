c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2016  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program freefix  --  free energy restraint thermodynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "freefix" computes via a volume integral the free energy, enthalpy
c     and entropy correction for use of a flat-bottom harmonic restraint
c     to keep a guest close to a host during decoupling simulations
c
c     literature reference:
c
c     D. Hamelberg and J. A. McCammon, Standard Free Energy of Releasing
c     a Localized Water Molecule from the Binding Pockets of Proteins:
c     Double-Decoupling Method, Journal of the American Chemical Society,
c     126, 7683-7689 (2004)  [see equations 15-19]
c
c
      program freefix
      use iounit
      use math
      use units
      implicit none
      integer i,j,k
      integer i2,j2,k2
      integer maxgrid
      real*8 rt,force,flat
      real*8 temper,stdcon
      real*8 spacing,cube
      real*8 dist,dist2
      real*8 term,expterm
      real*8 dexpterm
      real*8 vol,dvol
      real*8 dg,ds,dh
      logical exist
      character*240 string
c
c
c     get force constant, restraint radius and temperature values
c
      call initial
      force = 1.0d0
      flat = 0.0d0
      temper = 298.0d0
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  force
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  flat
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  temper
   10    continue
      else
         write (iout,20)
   20    format (/,' Force Constant, Flat-Bottom Radius &',
     &              ' Temperature :  ',$)
         read (input,30)  string
   30    format (a240)
         read (string,*,err=40,end=40)  force,flat,temper
   40    continue
      end if
c
c     print the restraint parameter values and temperature
c
      write (iout,50)  force
   50 format (/,' Restraint Force Constant :',3x,f14.4,
     &           ' Kcal/mole/Ang^2')
      write (iout,60)  flat
   60 format (' Flat-Bottom Radius Value :',3x,f14.4,' Ang')
      write (iout,70)  temper
   70 format (' System Temperature Value :',3x,f14.4,' Kelvin')
c
c     find RT and Ang^3 per molecule at 1 mole/L concentration
c
      rt = temper * gasconst
      stdcon = 1.0d27 / avogadro
c
c     set values for the numerical integration grid spacing
c
      spacing = 0.03d0
      cube = spacing**3
      maxgrid = int((10.0d0+flat)/spacing)
c
c     numerically estimate the restraint volume integral
c
      vol = 0.0d0
      dvol = 0.0d0
      do i = -maxgrid, maxgrid
         i2 = i * i
         do j = -maxgrid, maxgrid
            j2 = j * j
            do k = -maxgrid, maxgrid
               k2 = k * k
               dist = spacing * sqrt(dble(i2+j2+k2))
               dist = max(0.0d0,dist-flat)
               dist2 = dist * dist
               term = -force * dist2 / rt
               expterm = 0.0d0
               dexpterm = 0.0d0
               if (term .ge. -20.0d0) then
                  expterm = cube * exp(term)
                  dexpterm = expterm * (force*dist2)/(rt*temper)
               end if
               vol = vol + expterm
               dvol = dvol + dexpterm
            end do
         end do
      end do
c
c     ouput the value of the restraint volume integral
c
      write (iout,80)  spacing
   80 format (/,' Numerical Grid Spacing :',5x,f14.6,' Ang')
      write (iout,90)  vol
   90 format (' Numerical Volume Integral :',2x,f14.6,' Ang^3')
c
c     use analytical solution if restraint is not flat-bottom
c
      if (flat .eq. 0.0d0) then
         vol = (pi*rt/force)**1.5d0
         dvol = 1.5d0 * vol / temper
         write (iout,100)  vol
  100    format (' Analytical Volume Integral :',1x,f14.6,' Ang^3')
      end if
c
c     calculate and print the constraint thermodynamic values
c
      dg = rt * log(vol/stdcon)
      ds = -dg/temper - rt*dvol/vol
      dh = dg + temper*ds
      write (iout,110)  dg
  110 format (/,' Restraint Free Energy :',6x,f14.4,' Kcal/mole')
      write (iout,120)  ds
  120 format (' Restraint Entropy Value :',4x,f14.4,' Kcal/mole/K')
      write (iout,130)  dh
  130 format (' Restraint Enthalpy Value :',3x,f14.4,' Kcal/mole')
c
c     perform any final tasks before program exit
c
      call final
      end
