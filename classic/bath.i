c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  bath.i  --  temperature and pressure control parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     maxnose     maximum length of the Nose-Hoover chain
c
c     kelvin0     target value for the system temperature (K)
c     kelvin      variable target temperature for thermostat (K)
c     atmsph      target value for the system pressure (atm)
c     tautemp     time constant for Berendsen thermostat (psec)
c     taupres     time constant for Berendsen barostat (psec)
c     compress    isothermal compressibility of medium (atm-1)
c     collide     collision frequency for Andersen thermostat
c     xnh         position of each chained Nose-Hoover thermostat
c     vnh         velocity of each chained Nose-Hoover thermostat
c     qnh         mass for each chained Nose-Hoover thermostat
c     gnh         coupling between chained Nose-Hoover thermostats
c     volmove     maximum volume move for Monte Carlo barostat (Ang**3)
c     voltrial    mean number of steps between Monte Carlo moves
c     isothermal  logical flag governing use of temperature control
c     isobaric    logical flag governing use of pressure control
c     anisotrop   logical flag governing use of anisotropic pressure
c     thermostat  choice of temperature control method to be used
c     barostat    choice of pressure control method to be used
c     volscale    choice of scaling method for Monte Carlo barostat
c
c
      integer maxnose
      parameter (maxnose=2)
      integer voltrial
      real*8 kelvin0,kelvin
      real*8 atmsph
      real*8 tautemp,taupres
      real*8 compress,collide
      real*8 xnh,vnh,qnh,gnh
      real*8 volmove
      logical isothermal
      logical isobaric
      logical anisotrop
      character*9 volscale
      character*10 barostat
      character*11 thermostat
      common /bath/ kelvin0,kelvin,atmsph,tautemp,taupres,compress,
     &              collide,xnh(maxnose),vnh(maxnose),qnh(maxnose),
     &              gnh(maxnose),volmove,voltrial,isothermal,isobaric,
     &              anisotrop,thermostat,barostat,volscale
