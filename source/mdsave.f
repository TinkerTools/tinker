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
      use dlmda
      use extfld
      use files
      use group
      use inform
      use iounit
      use mdstuf
      use moment
      use mpole
      use output
      use ost
      use polar
      use potent
      use rgddyn
      use socket
      use titles
      use uatom
      use units
      implicit none
      integer i,j,lext
      integer istep
      integer ixyz,iind
      integer ivel,ifrc,istc
      integer iend,isave
      integer freeunit
      integer trimtext
      integer modsave
      real*8 dt,pico
      real*8 epot,eksum
      real*8 xustc,yustc,zustc
      real*8 xuind,yuind,zuind
      real*8 xuchg,yuchg,zuchg
      real*8 xf,yf,zf
      real*8 xm,ym,zm
      logical exist,first
      character*7 ext
      character*240 endfile
      character*240 xyzfile
      character*240 velfile
      character*240 frcfile
      character*240 indfile
      character*240 stcfile
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
c     print the external electric field values
c
      if (use_exfld) then
         xf = texfld(1) * elefield
         yf = texfld(2) * elefield
         zf = texfld(3) * elefield
         if (digits .le. 6) then
            write (iout,320)  xf,yf,zf
  320       format (' External Field',7x,3f14.6)
         else if (digits .le. 8) then
            write (iout,330)  xf,yf,zf
  330       format (' External Field',7x,3f16.8)
         else
            write (iout,340)  xf,yf,zf
  340       format (' External Field',7x,3f18.10)
         end if
      end if
c
c     move stray molecules into periodic box if desired
c
      if (use_wrap)  call bounds
c
c     compute center of mass if saving dipole moment
c
      if (usyssave .or. uchgsave)  call compcent(xm,ym,zm)
c
c     compute dipole moment of system
c
      if (usyssave) then
         call dmoments (xm,ym,zm,xustc,yustc,zustc,xuind,yuind,zuind,
     &                           xuchg,yuchg,zuchg)
         if (digits .le. 6) then
            write (iout,350)  xuchg,yuchg,zuchg
  350       format (' System Charge Dipole',1x,3f20.6)
            write (iout,360)  xustc,yustc,zustc
  360       format (' System Static Dipole',1x,3f20.6)
         else if (digits .le. 8) then
            write (iout,370)  xuchg,yuchg,zuchg
  370       format (' System Charge Dipole',1x,3f22.8)
            write (iout,380)  xustc,yustc,zustc
  380       format (' System Static Dipole',1x,3f22.8)
         else
            write (iout,390)  xuchg,yuchg,zuchg
  390       format (' System Charge Dipole',1x,3f24.10)
            write (iout,400)  xustc,yustc,zustc
  400       format (' System Static Dipole',1x,3f24.10)
         end if
         if (use_polar) then
            if (digits .le. 6) then
               write (iout,410)  xuind,yuind,zuind
  410          format (' System Induced Dipole',3f20.6)
            else if (digits .le. 8) then
               write (iout,420)  xuind,yuind,zuind
  420          format (' System Induced Dipole',3f22.8)
            else
               write (iout,430)  xuind,yuind,zuind
  430          format (' System Induced Dipole',3f24.10)
            end if
         end if
         write (iout,440)
  440    format (' Charge Dipole by Atom Type:',
     &        /,' Type',11x,'X-UCharge',11x,'Y-UCharge',11x,'Z-UCharge')
         do i = 1, nunique
            write (iout,450)  utype(i),utv3(1,i),utv3(2,i),utv3(3,i)
  450       format (i5,3f20.6)
         end do
         write (iout,460)
  460    format (' Static Dipole by Atom Type:',
     &        /,' Type',11x,'X-UStatic',11x,'Y-UStatic',11x,'Z-UStatic')
         do i = 1, nunique
            write (iout,470)  utype(i),utv1(1,i),utv1(2,i),utv1(3,i)
  470       format (i5,3f20.6)
         end do
         if (use_polar) then
            write (iout,480)
  480       format (' Induced Dipole by Atom Type:',
     &        /,' Type',11x,'X-UInduce',11x,'Y-UInduce',11x,'Z-UInduce')
            do i = 1, nunique
               write (iout,490)  utype(i),utv2(1,i),utv2(2,i),utv2(3,i)
  490          format (i5,3f20.6)
            end do
         end if
      end if
c
c     compute velocity of unique atom types in the system
c
      if (vsyssave) then
         call velunique
         write (iout,500)
  500    format (' Velocity by Atom Type:',
     &     /,' Type',10x,'X-Velocity',10x,'Y-Velocity',10x,'Z-Velocity')
         do i = 1, nunique
            write (iout,510)  utype(i),utv1(1,i),utv1(2,i),utv1(3,i)
  510       format (i5,3f20.6)
         end do
      end if
c
c     save coordinates to archive or numbered structure file
c
      write (iout,520)  isave
  520 format (' Frame Number',13x,i10)
      if (coordsave) then
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
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=ixyz,file=xyzfile,form='unformatted',
     &                  status='new')
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
         write (iout,530)  xyzfile(1:trimtext(xyzfile))
  530    format (' Coordinate File',13x,a)
      end if
c
c     update the information needed to restart the trajectory
c
      if (use_ostdyn)  call saveost
      if (use_metadyn)  call savemeta
      if (dynsave)  call prtdyn
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
            write (ivel,540)  ngrp,title(1:ltitle)
  540       format (i6,2x,a)
            do i = 1, ngrp
               write (ivel,550)  i,(vcm(j,i),j=1,3)
  550          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
               write (ivel,560)  i,(wcm(j,i),j=1,3)
  560          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         else if (dcdsave) then
            call prtdcdv3 (ivel,first,'VEL')
         else
            call prtvec3 (ivel,'VEL')
         end if
         close (unit=ivel)
         write (iout,570)  velfile(1:trimtext(velfile))
  570    format (' Velocity File',15x,a)
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
            call prtdcdv3 (ifrc,first,'FRC')
         else
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
            call prtvec3 (ifrc,'FRC')
         end if
         close (unit=ifrc)
         write (iout,580)  frcfile(1:trimtext(frcfile))
  580    format (' Force Vector File',11x,a)
      end if
c
c     save the induced dipole components for the current step
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (cyclesave) then
            indfile = filename(1:leng)//'.'//ext(1:lext)//'ui'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
            call prtvec3 (iind,'UIN')
         else if (dcdsave) then
            indfile = filename(1:leng)
            call suffix (indfile,'dcdui','old')
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
            call prtdcdv3 (iind,first,'UIN')
         else
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
            call prtvec3 (iind,'UIN')
         end if
         close (unit=iind)
         write (iout,590)  indfile(1:trimtext(indfile))
  590    format (' Induced Dipole File',9x,a)
      end if
c
c     save the static dipole components for the current step
c
      if (ustcsave) then
         istc = freeunit ()
         if (cyclesave) then
            stcfile = filename(1:leng)//'.'//ext(1:lext)//'us'
            call version (stcfile,'new')
            open (unit=istc,file=stcfile,status='new')
            call prtvec3 (istc,'UST')
         else if (dcdsave) then
            stcfile = filename(1:leng)
            call suffix (stcfile,'dcdus','old')
            inquire (file=stcfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=istc,file=stcfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=istc,file=stcfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcdv3 (istc,first,'UST')
         else
            stcfile = filename(1:leng)
            call suffix (stcfile,'ustc','old')
            inquire (file=stcfile,exist=exist)
            if (exist) then
               call openend (istc,stcfile)
            else
               open (unit=istc,file=stcfile,status='new')
            end if
            call prtvec3 (istc,'UST')
         end if
         close (unit=istc)
         write (iout,600)  stcfile(1:trimtext(stcfile))
  600    format (' Static Dipole File',10x,a)
      end if
c
c     save the charge dipole components for the current step
c
      if (uchgsave) then
         istc = freeunit ()
         if (cyclesave) then
            stcfile = filename(1:leng)//'.'//ext(1:lext)//'uc'
            call version (stcfile,'new')
            open (unit=istc,file=stcfile,status='new')
            call prtvec3 (istc,'UCH')
         else if (dcdsave) then
            stcfile = filename(1:leng)
            call suffix (stcfile,'dcduc','old')
            inquire (file=stcfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=istc,file=stcfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=istc,file=stcfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcdv3 (istc,first,'UCH')
         else
            stcfile = filename(1:leng)
            call suffix (stcfile,'uchg','old')
            inquire (file=stcfile,exist=exist)
            if (exist) then
               call openend (istc,stcfile)
            else
               open (unit=istc,file=stcfile,status='new')
            end if
            call prtvec3 (istc,'UCH')
         end if
         close (unit=istc)
         write (iout,610)  stcfile(1:trimtext(stcfile))
  610    format (' Charge Dipole File',10x,a)
      end if
c
c     save the direct induced dipole components for the current step
c
      if (udirsave .and. use_polar) then
         iind = freeunit ()
         if (cyclesave) then
            indfile = filename(1:leng)//'.'//ext(1:lext)//'ud'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
            call prtvec3 (iind,'UDR')
         else if (dcdsave) then
            indfile = filename(1:leng)
            call suffix (indfile,'dcdud','old')
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
            call prtdcdv3 (iind,first,'UDR')
         else
            indfile = filename(1:leng)
            call suffix (indfile,'udir','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
            call prtvec3 (iind,'UDR')
         end if
         close (unit=iind)
         write (iout,620)  indfile(1:trimtext(indfile))
  620    format (' Direct Induced Dipole File',2x,a)
      end if
c
c     save the direct atomic electric field for the current step
c
      if (defsave .and. use_polar) then
         iind = freeunit ()
         if (cyclesave) then
            indfile = filename(1:leng)//'.'//ext(1:lext)//'de'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
            call prtvec3 (iind,'DEF')
         else if (dcdsave) then
            indfile = filename(1:leng)
            call suffix (indfile,'dcdde','old')
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
            call prtdcdv3 (iind,first,'DEF')
         else
            indfile = filename(1:leng)
            call suffix (indfile,'def','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
            call prtvec3 (iind,'DEF')
         end if
         close (unit=iind)
         write (iout,630)  indfile(1:trimtext(indfile))
  630    format (' Direct Electric Field File',2x,a)
      end if
c
c     save the total atomic electric field for the current step
c
      if (tefsave .and. use_polar) then
         iind = freeunit ()
         if (cyclesave) then
            indfile = filename(1:leng)//'.'//ext(1:lext)//'te'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
            call prtvec3 (iind,'TEF')
         else if (dcdsave) then
            indfile = filename(1:leng)
            call suffix (indfile,'dcdte','old')
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
            call prtdcdv3 (iind,first,'TEF')
         else
            indfile = filename(1:leng)
            call suffix (indfile,'tef','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
            call prtvec3 (iind,'TEF')
         end if
         close (unit=iind)
         write (iout,640)  indfile(1:trimtext(indfile))
  640    format (' Total Electric Field File',3x,a)
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
         write (iout,650)
  650    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     skip an extra line to keep the output formating neat
c
      modsave = mod(istep,iprint)
      if (verbose .and. modsave.ne.0) then
         write (iout,660)
  660    format ()
      end if
      return
      end
