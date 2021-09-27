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
c     "freefix" computes a restraint correction term for a single
c     flat-bottomed harmonic or six translation-rotation restraints
c     used to keep a ligand bound during free energy simulations
c
c
      program freefix
      use iounit
      implicit none
      integer next
      logical exist
      character*1 answer
      character*8 method
      character*240 string
c
c
c     get type of ligand binding restraint to be analyzed
c
      call initial
      method = 'HARMONIC'
      call nextarg (answer,exist)
      if (.not. exist) then
         answer = 'H'
         write (iout,10)  answer
   10    format (/,' Choose Harmonic Restraint or Boresch Restraint',
     &              ' [',a1,'] :  ',$)
         read (input,20)  string
   20    format (a240)
         next = 1
         call gettext (string,answer,next)
      end if
      call upcase (answer)
      if (answer .eq. 'B')  method = 'BORESCH '
c
c     compute the values of the free energy correction
c
      if (method .eq. 'HARMONIC')  call hfix
      if (method .eq. 'BORESCH ')  call bfix
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine hfix  --  find harmonic restraint correction  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "hfix" computes via a volume integral the free energy, enthalpy
c     and entropy correction for a flat-bottom harmonic restraint
c
c     literature reference:
c
c     D. Hamelberg and J. A. McCammon, "Standard Free Energy of
c     Releasing a Localized Water Molecule from the Binding Pockets
c     of Proteins: Double-Decoupling Method", Journal of the American
c     Chemical Society, 126, 7683-7689 (2004)  [equations 15-19]
c
c     enthalpy and entropy values added by Aaron Gordon, February 2017;
c     analytical volume integral derived by Zhi Wang, December 2017
c
c
      subroutine hfix
      use iounit
      use math
      use units
      implicit none
      integer i,j,k
      integer i2,j2,k2
      integer maxgrid
      real*8 kt,stdcon
      real*8 temp,force
      real*8 ri,ro,fi,fo
      real*8 spacing,cube
      real*8 dist,dist2
      real*8 term,expterm
      real*8 dexpterm
      real*8 v1,v2,v3
      real*8 dv1,dv2,dv3
      real*8 vol,dvol
      real*8 dg,ds,dh
      real*8 erf
      logical exist,donumer
      character*240 string
      external erf
c
c
c     get force constant, restraint radius and temperature values
c
      fi = 0.0d0
      ri = 0.0d0
      fo = 1.0d0
      ro = 0.0d0
      temp = 298.0d0
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  ri
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  fi
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ro
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  fo
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  temp
   10    continue
      else
         write (iout,20)  ri,fi
   20    format (/,' Enter Inner Radius & Force Constant [',
     &              f4.2,',',f4.1'] :  ',$)
         read (input,30)  string
   30    format (a240)
         read (string,*,err=40,end=40)  ri,fi
   40    continue
         ro = ri
         if (fi .ne. 0.0d0)  fo = fi
         write (iout,50)  ro,fo
   50    format (/,' Enter Outer Radius & Force Constant [',
     &              f4.2,',',f4.1'] :  ',$)
         read (input,60)  string
   60    format (a240)
         read (string,*,err=70,end=70)  ro,fo
   70    continue
         write (iout,80)  temp
   80    format (/,' Enter System Temperature Value [',f5.1,'K] :  ',$)
         read (input,90)  string
   90    format (a240)
         read (string,*,err=100,end=100)  temp
  100    continue
      end if
c
c     print the restraint parameter values and temperature
c
      write (iout,110)  ri
  110 format (/,' Inner Flat-Bottom Radius :',5x,f12.4,' Ang')
      write (iout,120)  fi
  120 format (' Inner Force Constant :',9x,f12.4,
     &           ' Kcal/mole/Ang^2')
      write (iout,130)  ro
  130 format (' Outer Flat-Bottom Radius :',5x,f12.4,' Ang')
      write (iout,140)  fo
  140 format (' Outer Force Constant :',9x,f12.4,
     &           ' Kcal/mole/Ang^2')
      write (iout,150)  temp
  150 format (' System Temperature Value :',5x,f12.4,' Kelvin')
c
c     zero values for force constants are not allowed
c
      if (fi .eq. 0.0d0) then
         fi = 1.0d0
         ri = 0.0d0
      end if
      if (fo .eq. 0.0d0)  fo = 1.0d0
c
c     find RT and Ang^3 per molecule at 1 mole/L concentration
c
      kt = temp * gasconst
      stdcon = 1.0d27 / avogadro
c
c     numerical estimation of the restraint volume integral
c
      donumer = .false.
      if (donumer) then
         spacing = 0.03d0
         cube = spacing**3
         maxgrid = int((10.0d0+ro)/spacing)
         vol = 0.0d0
         dvol = 0.0d0
         do i = -maxgrid, maxgrid
            i2 = i * i
            do j = -maxgrid, maxgrid
               j2 = j * j
               do k = -maxgrid, maxgrid
                  k2 = k * k
                  dist = spacing * sqrt(dble(i2+j2+k2))
                  if (dist .gt. ro) then
                     dist = dist - ro
                     force = fo
                  else if (dist .lt. ri) then
                     dist = ri - dist
                     force = fi
                  else
                     dist = 0.0d0
                     force = 0.0d0
                  end if
                  dist2 = dist * dist
                  term = -force * dist2 / kt
                  expterm = 0.0d0
                  dexpterm = 0.0d0
                  if (term .ge. -20.0d0) then
                     expterm = cube * exp(term)
                     dexpterm = expterm * (force*dist2)/(kt*temp)
                  end if
                  vol = vol + expterm
                  dvol = dvol + dexpterm
               end do
            end do
         end do
         write (iout,160)  spacing
  160    format (/,' Numerical Grid Spacing :',7x,f12.4,' Ang')
         write (iout,170)  vol
  170    format (' Numerical Volume Integral :',4x,f12.4,' Ang^3')
         write (iout,180)  dvol
  180    format (' Numerical dVol/dT Value :',6x,f12.4,' Ang^3/K')
      end if
c
c     analytical evaluation of the restraint volume integral
c
      v1 = 2.0d0*pi*ri*(-2.0d0+exp(-ri**2*fi/kt))*kt/fi
     &        + sqrt(kt*(pi/fi)**3)*(2.0d0*fi*ri*ri+kt)
     &             *erf(ri*sqrt(fi/kt))
      v2 = (4.0d0*pi/3.0d0) * (ro**3-ri**3)
      v3 = sqrt(kt*(pi/fo)**3)
     &        * (2.0d0*fo*ro*ro+kt+4.0d0*ro*sqrt(kt*fo/pi))
      vol = v1 + v2 + v3
      dv1 = 2.0d0*pi*ri**3*exp(-ri**2*fi/kt)/temp
     &         + 2.0d0*pi*ri*(-2.0d0+exp(-ri**2*fi/kt))*kt/(fi*temp)
     &         + 0.5d0*sqrt((pi/fi)**3)*sqrt(kt)*(2.0d0*ri**2*fi+kt)
     &             *erf(ri*sqrt(fi/kt))/temp
     &         - pi*ri*exp(-ri**2*fi/kt)*(2.0d0*ri**2*fi+kt)/(fi*temp)
     &         + sqrt((kt*pi/fi)**3)*erf(ri*sqrt(fi/kt))/temp
      dv2 = 0.0d0
      dv3 = sqrt(kt*(pi/fo)**3)*fo*ro*ro/temp
     &         + 4.0d0*kt*(pi/fo)*ro/temp
     &         + 1.5d0*sqrt((kt*pi/fo)**3)/temp
      dvol = dv1 + dv2 + dv3
      write (iout,190)  vol
  190 format (/,' Analytical Volume Integral :',3x,f12.4,' Ang^3')
      write (iout,200)  dvol
  200 format (' Analytical dVol/dT Value :',5x,f12.4,' Ang^3/K')
c
c     calculate and print the restraint thermodynamic values
c
      dg = -kt * log(vol/stdcon)
      ds = -dg/temp + kt*dvol/vol
      dh = dg + temp*ds
      write (iout,210)  dg
  210 format (/,' Restraint Free Energy :',8x,f12.4,' Kcal/mole')
      write (iout,220)  ds
  220 format (' Restraint Entropy Value :',6x,f12.4,' Kcal/mole/K')
      write (iout,230)  dh
  230 format (' Restraint Enthalpy Value :',5x,f12.4,' Kcal/mole')
      write (iout,240)  -temp*ds
  240 format (' Restraint -T deltaS Value :',4x,f12.4,' Kcal/mole')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine bfix  --  find Boresch restraint correction  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "bfix" computes via a volume integral the free energy, enthalpy
c     and entropy correction for a set of six harmonic translation-
c     rotation "Boresch" restraints
c
c     literature reference:
c
c     S. Boresch, F. Tettinger, M. Leitgeb and M. Karplus, "Absolute
c     Binding Free Energies: A Quantitative Approach for Their
c     Calculation", Journal of Physical Chemistry B, 107, 9535-9551
c     (2003)  [equation 14]
c
c
      subroutine bfix
      use iounit
      use math
      use units
      implicit none
      real*8 dist,temp
      real*8 ang1,ang2
      real*8 fd,fa1,fa2
      real*8 ft1,ft2,ft3
      real*8 kt,stdcon
      real*8 sine1,sine2
      real*8 term1,term2
      real*8 term3,term4
      real*8 dg,ds,dh
      logical exist
      character*240 string
c
c
c     get distance, angle, force constant and temperature values
c
      dist = 0.0d0
      fd = 0.0d0
      ang1 = 0.0d0
      fa1 = 0.0d0
      ang2 = 0.0d0
      fa2 = 0.0d0
      ft1 = 0.0d0
      ft2 = 0.0d0
      ft3 = 0.0d0
      temp = 298.0d0
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  dist
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  fd
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ang1
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  fa1
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ang2
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  fa2
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ft1
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ft2
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  ft3
         call nextarg (string,exist)
         read (string,*,err=10,end=10)  temp
   10    continue
      else
         write (iout,20)
   20    format (/,' Enter Distance Value & Force Constant :  ',$)
         read (input,30)  string
   30    format (a240)
         read (string,*,err=40,end=40)  dist,fd
   40    continue
         write (iout,50)
   50    format (/,' Enter 1st Angle Value & Force Constant :  ',$)
         read (input,60)  string
   60    format (a240)
         read (string,*,err=70,end=70)  ang1,fa1
   70    continue
         write (iout,80)
   80    format (/,' Enter 2nd Angle Value & Force Constant :  ',$)
         read (input,90)  string
   90    format (a240)
         read (string,*,err=100,end=100)  ang2,fa2
  100    continue
         write (iout,110)
  110    format (/,' Enter Three Torsional Force Constants :  ',$)
         read (input,120)  string
  120    format (a240)
         read (string,*,err=130,end=130)  ft1,ft2,ft3
  130    continue
         write (iout,140)  temp
  140    format (/,' Enter System Temperature Value [',f5.1,'K] :  ',$)
         read (input,150)  string
  150    format (a240)
         read (string,*,err=160,end=160)  temp
  160    continue
      end if
c
c     print the restraint parameter values and temperature
c
      write (iout,170)  dist
  170 format (/,' Distance Reference Value :',5x,f12.4,' Ang')
      write (iout,180)  fd
  180 format (' Distance Force Constant :',6x,f12.4,' Kcal/mole/Ang^2')
      write (iout,190)  ang1
  190 format (' Angle 1 Reference Value :',6x,f12.4,' Deg')
      write (iout,200)  fa1
  200 format (' Angle 1 Force Constant :',7x,f12.4,' Kcal/mole/Rad^2')
      write (iout,210)  ang2
  210 format (' Angle 2 Reference Value :',6x,f12.4,' Deg')
      write (iout,220)  fa2
  220 format (' Angle 2 Force Constant :',7x,f12.4,' Kcal/mole/Rad^2')
      write (iout,230)  ft1
  230 format (' Torsion 1 Force Constant :',5x,f12.4,' Kcal/mole/Rad^2')
      write (iout,240)  ft2
  240 format (' Torsion 2 Force Constant :',5x,f12.4,' Kcal/mole/Rad^2')
      write (iout,250)  ft3
  250 format (' Torsion 3 Force Constant :',5x,f12.4,' Kcal/mole/Rad^2')
      write (iout,260)  temp
  260 format (' System Temperature Value :',5x,f12.4,' Kelvin')
c
c     find RT and Ang^3 per molecule at 1 mole/L concentration
c
      kt = temp * gasconst
      stdcon = 1.0d27 / avogadro
c
c     compute the free energy correction due to Boresch restraints
c
      sine1 = sin(ang1/radian)
      sine2 = sin(ang2/radian)
      term1 = 8.0d0 * pi * pi * stdcon
      term2 = sqrt(fd*fa1*fa2*ft1*ft2*ft3)
      term3 = dist * dist * sine1 * sine2
      term4 = (2.0d0 * pi * kt)**3
c
c     calculate and print the restraint thermodynamic values
c
      dg = kt * log((term1*term2)/(term3*term4))
      ds = dg/temp - 3.0d0*gasconst
      dh = dg + temp*ds
      write (iout,270)  dg
  270 format (/,' Restraint Free Energy :',8x,f12.4,' Kcal/mole')
      write (iout,280)  ds
  280 format (' Restraint Entropy Value :',6x,f12.4,' Kcal/mole/K')
      write (iout,290)  dh
  290 format (' Restraint Enthalpy Value :',5x,f12.4,' Kcal/mole')
      write (iout,300)  -temp*ds
  300 format (' Restraint -T deltaS Value :',4x,f12.4,' Kcal/mole')
      return
      end
