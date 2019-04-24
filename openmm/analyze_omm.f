c
c
c     ###########################################################
c     ##  COPYRIGHT (C) 2019 by Zhi Wang & Jay William Ponder  ##
c     ##                  All Rights Reserved                  ##
c     ###########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program analyze_omm  --  force field energy via OpenMM API  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "analyze_omm" computes the force field energy via an interface to
c     the OpenMM GPU code
c
c     modifications by Zhi Wang, Ponder Lab, St. Louis, April 2019
c
c
      program analyze_omm
      use files
      use inform
      use iounit
      use mdstuf
      use openmm
      implicit none
      integer iarc
      integer freeunit
      real*8 energy,dt
      character*240 arcfile
c
c
c     setup the molecular mechanics calculation
c
      call initial
      call getarc (iarc)
      close (unit=iarc)
      call mechanic
c
c     set needed variables for OpenMM energy calculations
c
      call kopenmm
      call mdinit
      dt = 0.001d0
      integrate = 'VERLET'
c
c     reopen the dynamics trajectory file
c
      iarc = freeunit ()
      arcfile = filename
      call suffix (arcfile,'arc','old')
      open (unit=iarc,file=arcfile,status='old')
c
c     find potential energies for the trajectory
c
      call readxyz (iarc)
      call ommdata
      call openmm_init (ommHandle,dt)
      do while (.not. abort)
         call openmm_bar_energy (ommHandle,energy)
         write (iout,10) energy
   10    format(/,' Total Potential Energy :',8x,f20.8,' Kcal/mole')
         call readxyz (iarc)
      end do
c
c     perform any final tasks before program exit
c
      call openmm_cleanup (ommHandle)
      close (unit=iarc)
      call final
      end
