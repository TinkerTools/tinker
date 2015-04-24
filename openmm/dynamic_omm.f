c
c
c     ############################################################
c     ##                  COPYRIGHT (C) 2015                    ##
c     ##     by Mark Friedrichs, Lee-Ping Wang & Jay Ponder     ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  program dynamic_omm  --  molecular dynamics via OpenMM API  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "dynamic_omm" computes a molecular or stochastic dynamics
c     trajectory via an interface to the OpenMM GPU code for the
c     computation of forces and dynamics integration steps
c
c
      program dynamic_omm
      use sizes
      use atoms
      use bath
      use bndstr
      use bound
      use inform
      use iounit
      use keys
      use mdstuf
      use openmm
      use openmp
      use potent
      use solute
      use stodyn
      use usage
      implicit none
      integer i,istep,nstep
      integer mode,next
      real*8 e,dt,dtdump
      real*8, allocatable :: derivs(:,:)
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
c
c     additional variables used by the OpenMM interface
c
      integer nextStep
      integer nextUpdate
      integer updateCalls
      integer callMdStat
      integer callMdSave
      real*8 elapsed,nsPerDay,cpu
      logical oneTimeStepPerUpdate
c
c
c     set up the structure and molecular mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,160)
  160          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     initialize any holonomic constraints and setup dynamics
c
      call shakeup
      call mdinit
c
c     get TINKER energy/gradient values for initial structure
c
      allocate (derivs(3,n))
      call gradient (e,derivs)
      deallocate (derivs)
c
c     map TINKER data structures to OpenMM wrapper structures
c
      call map_tinker_to_openmm ()
c
c     check required potentials forces are available in OpenMM
c
      call openmm_init (ommHandle,dt)
c
c     compare the TINKER and OpenMM energy/gradient values
c
      call openmm_test ()
c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,330)
  330    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,340)
  340    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,350)
  350    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,360)
  360    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,370)
  370    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,380)
  380    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if
c
c     initialize some counters used during the MD steps
c
      istep = 0
      nextStep = 1
      updateCalls = 0
c
c     if oneTimeStepPerUpdate=true, then take one step on the GPU
c     and retrieve data (position/velocities/energies) from GPU;
c     if oneTimeStepPerUpdate=false, then take multiple time steps
c     on the GPU before retrieving data to the CPU; the number of
c     steps taken per iteration is determined by "iwrite"
c
c     oneTimestepPerUpdate = .true.
      oneTimestepPerUpdate = .false.
c
c     integrate equations of motion to take the MD time steps
c
      call settime
      if (oneTimeStepPerUpdate) then
         callMdStat = 1
         callMdSave = 1
         do while (istep .lt. nstep)
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      else
         nextUpdate = iwrite 
         callMdStat = 0
         callMdSave = 1
         do while (istep .lt. nstep)
            nextStep = nextUpdate - istep
            nextUpdate = nextUpdate + iwrite
            if (nextStep+istep .gt. nstep) then
               nextStep = nstep - istep
            end if
            call openmm_take_steps (ommHandle,nextStep)
            istep = istep + nextStep
            updateCalls = updateCalls + 1
            call openmm_update (ommHandle,dt,istep,
     &                          callMdStat,callMdSave)
         end do
      end if
      call gettime (elapsed,cpu)
c
c     print performance and timing information
c
      nsPerDay = 86.4d0 * nstep * dt / elapsed
      write (iout,410)  nsPerDay,elapsed,nstep,updateCalls,
     &                  1000.0d0*dt,n,nthread
  410 format (/,' Performance:  ns/day',9x,f12.4,
     &        /,15x,'Wall Time',6x,f12.4,
     &        /,15x,'Steps',14x,i8,
     &        /,15x,'Updates',12x,i8,
     &        /,15x,'Time Step',6x,f12.4,
     &        /,15x,'Atoms',14x,i8,
     &        /,15x,'Threads',12x,i8)
c
c     perform any final tasks before program exit
c
      call openmm_cleanup (ommHandle)
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine map_tinker_to_openmm  --  transfer structures  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "map_tinker_to_openmm" uses calls to the OpenMM C++ interface
c     to transfer TINKER data to the corresponding OpenMM structures
c
c
      subroutine map_tinker_to_openmm ()
      use sizes
      use angbnd
      use angpot
      use atomid
      use atoms
      use bath
      use bitor
      use bndpot
      use bndstr
      use boxes
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use freeze
      use inform
      use ktrtor
      use kvdws
      use limits
      use mdstuf
      use moldyn
      use mplpot
      use mpole
      use nonpol
      use opbend
      use pitors
      use pme
      use polar
      use polgrp
      use polpot
      use potent
      use solute
      use stodyn
      use strbnd
      use torpot
      use tors
      use tortor
      use urey
      use urypot
      use usage
      use vdw
      use vdwpot
      implicit none
c
c
c     use C++ interface calls to map TINKER variables to OpenMM
c
      call set_parameters (maxval)
      call set_angbnd_data (ak,anat,afld,nangle,iang)
      call set_angpot_data (angunit,stbnunit,aaunit,opbunit,opdunit,
     &                      cang,qang,pang,sang,copb,qopb,popb,sopb,
     &                      copd,qopd,popd,sopd,angtyp,opbtyp)
      call set_atomid_data (mass,tag,class,atomic,valence,name,story)
      call set_atoms_data (x,y,z,n,type)
      call set_bath_data (kelvin,atmsph,tautemp,taupres,compress,
     &                    collide,vnh,qnh,gnh,volmove,voltrial,
     &                    isothermal,isobaric,anisotrop,thermostat,
     &                    barostat,volscale)
      call set_bitor_data (ibitor,nbitor)
      call set_bndpot_data (cbnd,qbnd,bndunit,bndtyp)
      call set_bndstr_data (bk,bl,nbond,ibnd)
      call set_boxes_data (xbox,ybox,zbox,alpha,beta,gamma,xbox2,ybox2,
     &                     zbox2,box34,lvec,recip,volbox,beta_sin,
     &                     beta_cos,gamma_sin,gamma_cos,beta_term,
     &                     gamma_term,orthogonal,monoclinic,triclinic,
     &                     octahedron,spacegrp)
      call set_chgpot_data (electric,dielec,ebuffer,c2scale,c3scale,
     &                      c4scale,c5scale,neutnbr,neutcut)
      call set_couple_data (n12,i12,n13,i13,n14,i14,n15,i15,maxval)
      call set_deriv_data (desum,deb,dea,deba,deub,deaa,deopb,deopd,
     &                     deid,deit,det,dept,debt,deat,dett,dev,dec,
     &                     decd,ded,dem,dep,der,des,delf,deg,dex)
      call set_energi_data (esum,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,
     &                      et,ept,ebt,eat,ett,ev,ec,ecd,ed,em,ep,er,
     &                      es,elf,eg,ex)
      call set_ewald_data (aewald,boundary)
      call set_freeze_data (krat,nrat,nratx,irat,iratx,kratx,
     &                      ratimage,use_rattle)
      call set_inform_data (digits,iprint,iwrite,isend,verbose,debug,
     &                      holdup,abort)
      call set_ktrtor_data (ttx,tty,tbf,tbx,tby,tbxy,tnx,tny,
     &                      ktt,maxntt,maxtgrd)
      call set_kvdws_data (rad,eps,rad4,eps4,reduct)
      call set_limits_data (vdwcut,chgcut,dplcut,mpolecut,vdwtaper,
     &                      chgtaper,dpltaper,mpoletaper,ewaldcut,
     &                      use_ewald,use_lights,use_list,use_vlist,
     &                      use_clist,use_mlist)
      call set_mdstuf_data (nfree,irest,velsave,frcsave,uindsave,
     &                      integrate)
      call set_moldyn_data (v,a,aalt)
      call set_mplpot_data (m2scale,m3scale,m4scale,m5scale)
      call set_mpole_data (pole,rpole,npole,ipole,polsiz,pollist,
     &                     polaxe,zaxis,xaxis,yaxis,maxpole)
      call set_nonpol_data (solvprs,surften,spcut,spoff,stcut,stoff,
     &                      rcav,rdisp,cdisp)
      call set_opbend_data (opbk,iopb,nopbend)
      call set_pitors_data (kpit,ipit,npitors)
      call set_pme_data (bsmod1,bsmod2,bsmod3,thetai1,thetai2,thetai3,
     &                   qgrid,qfac,nfft1,nfft2,nfft3,bsorder,igrid)
      call set_polar_data (polarity,thole,pdamp,uind,uinp,
     &                     uinds,uinps,npolar)
      call set_polgrp_data (np11,ip11,np12,ip12,np13,ip13,np14,ip14,
     &                      maxp11,maxp12,maxp13,maxp14)
      call set_polpot_data (poleps,p2scale,p3scale,p4scale,p41scale,
     &                      p5scale,d1scale,d2scale,d3scale,d4scale,
     &                      u1scale,u2scale,u3scale,u4scale,poltyp)
      call set_potent_data (use_bond,use_angle,use_strbnd,use_urey,
     &                      use_angang,use_opbend,use_opdist,use_improp,
     &                      use_imptor,use_tors,use_pitors,use_strtor,
     &                      use_angtor,use_tortor,use_vdw,use_charge,
     &                      use_chgdpl,use_dipole,use_mpole,use_polar,
     &                      use_rxnfld,use_solv,use_metal,use_geom,
     &                      use_extra,use_born,use_orbit)
      call set_solute_data (rsolv,asolv,rborn,drb,drbp,drobc,doffset,
     &                      p1,p2,p3,p4,p5,gpol,shct,aobc,bobc,gobc,
     &                      vsolv,wace,s2ace,uace,solvtyp,borntyp)
      call set_stodyn_data (friction,fgamma,use_sdarea)
      call set_strbnd_data (sbk,isb,nstrbnd)
      call set_torpot_data (idihunit,itorunit,torsunit,ptorunit,
     &                      storunit,atorunit,ttorunit)
      call set_tors_data (tors1,tors2,tors3,tors4,tors5,tors6,
     &                    ntors,itors)
      call set_tortor_data (itt,ntortor)
      call set_urey_data (uk,ul,iury,nurey)
      call set_urypot_data (cury,qury,ureyunit)
      call set_usage_data (nuse,iuse,use)
      call set_vdw_data (radmin,epsilon,radmin4,epsilon4,radhbnd,
     &                   epshbnd,kred,ired,nvdw,ivdw,jvdw)
      call set_vdwpot_data (abuck,bbuck,cbuck,ghal,dhal,v2scale,
     &                      v3scale,v4scale,v5scale,igauss,ngauss,
     &                      use_vcorr,vdwindex,vdwtyp,radtyp,radsiz,
     &                      radrule,epsrule,gausstyp)
      return 
      end
