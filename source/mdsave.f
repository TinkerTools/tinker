c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mdsave  --  save trajectory and restart files  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mdsave" writes molecular dynamics trajectory snapshots and
c     auxiliary files with velocity, force or induced dipole data;
c     also checks for user requested termination of a simulation
c
c
      subroutine mdsave (istep,dt,epot,eksum)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use group
      use inform
      use iounit
      use mdstuf
      use mpole
      use output
      use ost
      use polar
      use potent
      use rgddyn
      use socket
      use titles
      implicit none
      integer i,j,lext
      integer istep
      integer ixyz,iind
      integer ivel,ifrc
      integer iend,isave
      integer freeunit
      integer trimtext
      integer modsave
      real*8 dt,pico
      real*8 epot,eksum
      logical exist,first
      character*7 ext
      character*240 endfile
      character*240 xyzfile
      character*240 velfile
      character*240 frcfile
      character*240 indfile
c
c
c     send data via external socket communication if desired
c
      if (.not.sktstart .or. use_socket)  call sktdyn (istep,dt,epot)
c
c     check number of steps between trajectory file saves
c
      modsave = mod(istep,iwrite)
      if (modsave .ne. 0)  return
c
c     get the sequence number of the current trajectory frame
c
      isave = nprior + istep/iwrite
      lext = 3
      call numeral (isave,ext,lext)
c
c     print header for the instantaneous values at current step
c
      pico = dble(istep) * dt
      write (iout,10)  istep
   10 format (/,' Instantaneous Values for Frame Saved at',
     &           i10,' Dynamics Steps')
c
c     print the current time, potential and kinetic energies
c
      if (digits .ge. 8) then
         write (iout,20)  pico
   20    format (/,' Current Time',8x,f19.8,' Picosecond')
         write (iout,30)  epot
   30    format (' Current Potential',3x,f19.8,' Kcal/mole')
         write (iout,40)  eksum
   40    format (' Current Kinetic',5x,f19.8,' Kcal/mole')
      else if (digits .ge. 6) then
         write (iout,50)  pico
   50    format (/,' Current Time',8x,f17.6,' Picosecond')
         write (iout,60)  epot
   60    format (' Current Potential',3x,f17.6,' Kcal/mole')
         write (iout,70)  eksum
   70    format (' Current Kinetic',5x,f17.6,' Kcal/mole')
      else
         write (iout,80)  pico
   80    format (/,' Current Time',8x,f15.4,' Picosecond')
         write (iout,90)  epot
   90    format (' Current Potential',3x,f15.4,' Kcal/mole')
         write (iout,100)  eksum
  100    format (' Current Kinetic',5x,f15.4,' Kcal/mole')
      end if
c
c     print adaptive lambda bias information if present
c
      if (use_ostdyn) then
         if (digits .ge. 8) then
            write (iout,110)  ostlambda
  110       format (' Current Lambda',6x,f19.8)
            write (iout,120)  ostdedl
  120       format (' Current dU/dLambda',2x,f19.8,' Kcal/mole')
            write (iout,130)  ostdgdl
  130       format (' Current dg/dLambda',2x,f19.8,' Kcal/mole')
            write (iout,140)  -ostddgdl
  140       format (' Current -ddG/dLambda',f19.8,' Kcal/mole')
            write (iout,150)  eosttot
  150       format (' Estimated Delta G',3x,f19.8,' Kcal/mole')
         else if (digits .ge. 6) then
            write (iout,160)  ostlambda
  160       format (' Current Lambda',6x,f17.6)
            write (iout,170)  ostdedl
  170       format (' Current dU/dLambda',2x,f17.6,' Kcal/mole')
            write (iout,180)  ostdgdl
  180       format (' Current dg/dLambda',2x,f17.6,' Kcal/mole')
            write (iout,190)  -ostddgdl
  190       format (' Current -ddG/dLambda',f17.6,' Kcal/mole')
            write (iout,200)  eosttot
  200       format (' Estimated Delta G',3x,f17.6,' Kcal/mole')
         else
            write (iout,210)  ostlambda
  210       format (' Current Lambda',6x,f15.4)
            write (iout,220)  ostdedl
  220       format (' Current dU/dLambda',2x,f15.4,' Kcal/mole')
            write (iout,230)  ostdgdl
  230       format (' Current dg/dLambda',2x,f15.4,' Kcal/mole')
            write (iout,240)  -ostddgdl
  240       format (' Current -ddG/dLambda',f15.4,' Kcal/mole')
            write (iout,250)  eosttot
  250       format (' Estimated Delta G',3x,f15.4,' Kcal/mole')
         end if
      else if (use_metadyn) then
         if (digits .ge. 8) then
            write (iout,110)  ostlambda
            write (iout,150)  eosttot
         else if (digits .ge. 6) then
            write (iout,160)  ostlambda
            write (iout,200)  eosttot
         else
            write (iout,210)  ostlambda
            write (iout,250)  eosttot
         end if
      end if
c
c     print the values of the lattice lengths and angles
c
      if (use_bounds) then
         if (digits .le. 6) then
            write (iout,260)  xbox,ybox,zbox
  260       format (' Lattice Lengths',6x,3f14.6)
            write (iout,270)  alpha,beta,gamma
  270       format (' Lattice Angles',7x,3f14.6)
         else if (digits .le. 8) then
            write (iout,280)  xbox,ybox,zbox
  280       format (' Lattice Lengths',6x,3f16.8)
            write (iout,290)  alpha,beta,gamma
  290       format (' Lattice Angles',7x,3f16.8)
         else
            write (iout,300)  xbox,ybox,zbox
  300       format (' Lattice Lengths',6x,3f18.10)
            write (iout,310)  alpha,beta,gamma
  310       format (' Lattice Angles',7x,3f18.10)
         end if
      end if
c
c     move stray molecules into periodic box if desired
c
      if (use_wrap)  call bounds
c
c     save coordinates to archive or numbered structure file
c
      ixyz = freeunit ()
      if (cyclesave) then
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
         call prtxyz (ixyz)
      else if (dcdsave) then
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'dcd','old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            first = .false.
            open (unit=ixyz,file=xyzfile,form='unformatted',
     &               status='old',position='append')
         else
            first = .true.
            open (unit=ixyz,file=xyzfile,form='unformatted',
     &               status='new')
         end if
         call prtdcd (ixyz,first)
      else
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc','old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            call openend (ixyz,xyzfile)
         else
            open (unit=ixyz,file=xyzfile,status='new')
         end if
         call prtxyz (ixyz)
      end if
      close (unit=ixyz)
      write (iout,320)  isave
  320 format (' Frame Number',13x,i10)
      write (iout,330)  xyzfile(1:trimtext(xyzfile))
  330 format (' Coordinate File',13x,a)
c
c     update the information needed to restart the trajectory
c
      call prtdyn
      if (use_ost)  call saveost
      if (use_meta)  call savemeta
c
c     save the velocity vector components at the current step
c
      if (velsave) then
         ivel = freeunit ()
         if (cyclesave) then
            velfile = filename(1:leng)//'.'//ext(1:lext)//'v'
            call version (velfile,'new')
            open (unit=ivel,file=velfile,status='new')
         else if (dcdsave) then
            velfile = filename(1:leng)
            call suffix (velfile,'dcdv','old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=ivel,file=velfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=ivel,file=velfile,form='unformatted',
     &                  status='new')
            end if
         else
            velfile = filename(1:leng)
            call suffix (velfile,'vel','old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               call openend (ivel,velfile)
            else
               open (unit=ivel,file=velfile,status='new')
            end if
         end if
         if (integrate .eq. 'RIGIDBODY') then
            write (ivel,340)  ngrp,title(1:ltitle)
  340       format (i6,2x,a)
            do i = 1, ngrp
               write (ivel,350)  i,(vcm(j,i),j=1,3)
  350          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
               write (ivel,360)  i,(wcm(j,i),j=1,3)
  360          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         else if (dcdsave) then
            call prtdcdv (ivel,first)
         else
            call prtvel (ivel)
         end if
         close (unit=ivel)
         write (iout,370)  velfile(1:trimtext(velfile))
  370    format (' Velocity File',15x,a)
      end if
c
c     save the force vector components for the current step; not
c     available for rigid body or multiple time step integrators
c
      if (frcsave .and. integrate.ne.'RIGIDBODY'
     &       .and. integrate.ne.'VRESPA'
     &       .and. integrate.ne.'BRESPA'
     &       .and. integrate.ne.'SRESPA') then
         ifrc = freeunit ()
         if (cyclesave) then
            frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
            call version (frcfile,'new')
            open (unit=ifrc,file=frcfile,status='new')
         else if (dcdsave) then
            frcfile = filename(1:leng)
            call suffix (frcfile,'dcdf','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=ifrc,file=frcfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=ifrc,file=frcfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcdf (ifrc,first)
         else
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
            call prtfrc (ifrc)
         end if
         close (unit=ifrc)
         write (iout,380)  frcfile(1:trimtext(frcfile))
  380    format (' Force Vector File',11x,a)
      end if
c
c     save the induced dipole components for the current step
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (cyclesave) then
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
         else if (dcdsave) then
            indfile = filename(1:leng)
            call suffix (indfile,'dcdu','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=iind,file=indfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=iind,file=indfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcdu (iind,first)
         else
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
            call prtuind (iind)
         end if
         close (unit=iind)
         write (iout,390)  indfile(1:trimtext(indfile))
  390    format (' Induced Dipole File',9x,a)
      end if
c
c     test for requested termination of the dynamics calculation
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,400)
  400    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     skip an extra line to keep the output formating neat
c
      modsave = mod(istep,iprint)
      if (verbose .and. modsave.ne.0) then
         write (iout,410)
  410    format ()
      end if
      return
      end
