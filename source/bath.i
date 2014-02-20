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
c     maxnose     maximum length of Nose-Hoover thermostat chain
c
c     kelvin      target value for the system temperature (K)
c     atmsph      target value for the system pressure (atm)
c     tautemp     time constant for Berendsen thermostat (psec)
c     taupres     time constant for Berendsen barostat (psec)
c     compress    isothermal compressibility of medium (atm-1)
c     collide     collision frequency for Andersen thermostat
c     vnh         velocity of each chained Nose-Hoover thermostat
c     qnh         mass for each chained Nose-Hoover thermostat
c     gnh         force for each chained Nose-Hoover thermostat
c     vbar        velocity of log volume for Nose-Hoover barostat
c     qbar        mass of the volume for Nose-Hoover barostat
c     gbar        force for the volume for Nose-Hoover barostat
c     eta         velocity value for Bussi-Parrinello barostat
c     volmove     maximum volume move for Monte Carlo barostat (Ang**3)
c     voltrial    mean number of steps between Monte Carlo moves
c     isothermal  logical flag governing use of temperature control
c     isobaric    logical flag governing use of pressure control
c     anisotrop   logical flag governing use of anisotropic pressure
c     thermostat  choice of temperature control method to be used
c     barostat    choice of pressure control method to be used
c     volscale    choice of scaling method for Monte Carlo barostat
c     q_iso1/2    masses of isokinetic thermostats
c     v_iso1/2    velocities of isokinetic thermostats
c     f_iso1/2    forces acting on isokinetic thermostats
c     len_nhc     In isokinetic dynamics, this is the Nose-Hoover chain length M. 
c                 In stochastic isokinetic dynamics, this is the chain length L
c
      real*8 q_iso1,q_iso2,v_iso1,v_iso2
      real*8 vnh,qnh,gnh
      real*8 kelvin,atmsph
      real*8 tautemp,taupres
      real*8 compress,collide
      real*8 vbar,qbar,gbar
      real*8 eta,volmove
      integer maxnose
      parameter (maxnose=4)
      integer len_nhc,isok_L,isok_M
      integer isok_chain
      parameter (isok_chain=4)
      integer voltrial
      logical isothermal
      logical isobaric
      logical anisotrop
      character*9 volscale
      character*11 barostat
      character*11 thermostat
      common /bath/ q_iso1(isok_chain,3,maxatm),
     &              q_iso2(isok_chain,3,maxatm),
     &              v_iso1(isok_chain,3,maxatm),
     &              v_iso2(isok_chain,3,maxatm),
     &              f_iso1(isok_chain,3,maxatm),
     &              kelvin,atmsph,tautemp,taupres,compress,collide,
     &              vnh(maxnose),qnh(maxnose),gnh(maxnose),vbar,qbar,
     &              gbar,eta,volmove,voltrial,len_nhc,isothermal,
     &              isobaric,anisotrop,thermostat,
     &              barostat,volscale            
