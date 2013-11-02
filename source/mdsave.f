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
      subroutine mdsave (istep,dt,epot)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'files.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'mpole.i'
      include 'output.i'
      include 'polar.i'
      include 'potent.i'
      include 'rgddyn.i'
      include 'socket.i'
      include 'titles.i'
      include 'units.i'
      integer i,j,k,istep
      integer ixyz,iind
      integer ivel,ifrc
      integer iend,idump,lext
      integer freeunit,trimtext
      integer moddump
      real*8 dt,epot,pico,wt
      logical exist
      character*7 ext
      character*120 endfile
      character*120 xyzfile
      character*120 velfile
      character*120 frcfile
      character*120 indfile
c     LPW
      integer ibox
      character*120 boxfile
c
c
c     send data via external socket communication if desired
c
      if (.not.skt_init .or. use_socket)  call sktdyn (istep,dt,epot)
c
c     check number of steps between trajectory file dumps
c
      moddump = mod(istep,iwrite)
      if (moddump .ne. 0)  return
c
c     get the sequence number of the current trajectory frame
c
      idump = nprior + istep/iwrite
      lext = 3
      call numeral (idump,ext,lext)
c
c     print header for the instantaneous values at current step
c
      pico = dble(istep) * dt
      write (iout,10)  istep
   10 format (/,' Instantaneous Values for Frame saved at',
     &           i10,' Dynamics Steps')
      if (digits .ge. 10) then
         write (iout,20)  pico
 20      format (/,' Current Time',6x,f21.10,' Picosecond')
         write (iout,30)  epot
 30      format (' Current Potential',1x,f21.10,' Kcal/mole')
      else if (digits .ge. 8) then
         write (iout,21)  pico
 21      format (/,' Current Time',6x,f19.8,' Picosecond')
         write (iout,31)  epot
 31      format (' Current Potential',1x,f19.8,' Kcal/mole')
      else if (digits .ge. 6) then
         write (iout,22)  pico
 22      format (/,' Current Time',6x,f17.6,' Picosecond')
         write (iout,32)  epot
 32      format (' Current Potential',1x,f17.6,' Kcal/mole')
      else 
         write (iout,23)  pico
 23      format (/,' Current Time',8x,f15.4,' Picosecond')
         write (iout,33)  epot
 33      format (' Current Potential',3x,f15.4,' Kcal/mole')
      end if
      if (use_bounds) then
         if (digits .ge. 10) then
            write (iout,40)  xbox,ybox,zbox
 40         format (' Lattice Lengths',6x,3f18.10)
            write (iout,50)  alpha,beta,gamma
 50         format (' Lattice Angles',7x,3f18.10)
            write (iout,60)  idump
 60         format (' Frame Number',17x,i10)
         else if (digits .ge. 8) then
            write (iout,41)  xbox,ybox,zbox
 41         format (' Lattice Lengths',6x,3f16.8)
            write (iout,51)  alpha,beta,gamma
 51         format (' Lattice Angles',7x,3f16.8)
            write (iout,61)  idump
 61         format (' Frame Number',15x,i10)
         else 
            write (iout,42)  xbox,ybox,zbox
 42         format (' Lattice Lengths',6x,3f14.6)
            write (iout,52)  alpha,beta,gamma
 52         format (' Lattice Angles',7x,3f14.6)
            write (iout,62)  idump
 62         format (' Frame Number',13x,i10)
         end if
      end if
c
c     update the information needed to restart the trajectory
c
      call prtdyn
c
c     save coordinates to an archive or numbered structure file
c
      ixyz = freeunit ()
      if (archive) then
         xyzfile = filename(1:leng)
         call suffix (xyzfile,'arc','old')
         inquire (file=xyzfile,exist=exist)
         if (exist) then
            call openend (ixyz,xyzfile)
         else
            open (unit=ixyz,file=xyzfile,status='new')
         end if
      else
         xyzfile = filename(1:leng)//'.'//ext(1:lext)
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
      end if
      call prtxyz (ixyz)
      close (unit=ixyz)
      write (iout,70)  xyzfile(1:trimtext(xyzfile))
   70 format (' Coordinate File',12x,a)
c
c     LPW save the box vectors to a .box file
c
      if (boxsave .and. use_bounds) then
         ibox = freeunit ()
         if (archive) then
            boxfile = filename(1:leng)
            call suffix (boxfile,'box','old')
            inquire (file=boxfile,exist=exist)
            if (exist) then
               call openend (ibox,boxfile)
            else
               open (unit=ibox,file=boxfile,status='new')
            end if
         else
            boxfile = filename(1:leng)//'.'//ext(1:lext)//'b'
            call version (boxfile,'new')
            open (unit=ibox,file=boxfile,status='new')
         end if
         write (ibox,183)  pico,xbox,ybox,zbox,alpha,beta,gamma
 183     format (f16.4,6f16.10)
         close (unit=ibox)
         write (iout,193)  boxfile(1:trimtext(boxfile))
 193     format (' Periodic Box File',10x,a)
      end if
c     end LPW modifications
c
c     save the velocity vector components at the current step
c
      if (velsave) then
         ivel = freeunit ()
         if (archive) then
            velfile = filename(1:leng)
            call suffix (velfile,'vel','old')
            inquire (file=velfile,exist=exist)
            if (exist) then
               call openend (ivel,velfile)
            else
               open (unit=ivel,file=velfile,status='new')
            end if
         else
            velfile = filename(1:leng)//'.'//ext(1:lext)//'v'
            call version (velfile,'new')
            open (unit=ivel,file=velfile,status='new')
         end if
         if (integrate .eq. 'RIGIDBODY') then
            write (ivel,80)  ngrp,title(1:ltitle)
   80       format (i6,2x,a)
            do i = 1, ngrp
               write (ivel,90)  i,(vcm(j,i),j=1,3)
   90          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
               write (ivel,100)  i,(wcm(j,i),j=1,3)
  100          format (i6,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         else
            write (ivel,110)  n,title(1:ltitle)
  110       format (i6,2x,a)
            do i = 1, n
               write (ivel,120)  i,name(i),(v(j,i),j=1,3)
  120          format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
            end do
         end if
         close (unit=ivel)
         write (iout,130)  velfile(1:trimtext(velfile))
  130    format (' Velocity File',15x,a)
      end if
c
c     save the force vector components for the current step
c
      if (frcsave .and. integrate.ne.'RIGIDBODY') then
         ifrc = freeunit ()
         if (archive) then
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               call openend (ifrc,frcfile)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
         else
            frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
            call version (frcfile,'new')
            open (unit=ifrc,file=frcfile,status='new')
         end if
         write (ifrc,140)  n,title(1:ltitle)
  140    format (i6,2x,a)
         do i = 1, n
            wt = mass(i) / convert
            write (ifrc,150)  i,name(i),(wt*a(j,i),j=1,3)
  150       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ifrc)
         write (iout,160)  frcfile(1:trimtext(frcfile))
  160    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar) then
         iind = freeunit ()
         if (archive) then
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               call openend (iind,indfile)
            else
               open (unit=iind,file=indfile,status='new')
            end if
         else
            indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
            call version (indfile,'new')
            open (unit=iind,file=indfile,status='new')
         end if
         write (iind,170)  n,title(1:ltitle)
  170    format (i6,2x,a)
         do i = 1, npole
            if (polarity(i) .ne. 0.0d0) then
               k = ipole(i)
               write (iind,180)  k,name(k),(debye*uind(j,i),j=1,3)
  180          format (i6,2x,a3,3f12.6)
            end if
         end do
         close (unit=iind)
         write (iout,190)  indfile(1:trimtext(indfile))
  190    format (' Induced Dipole File',10x,a)
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
         write (iout,200)
  200    format (/,' MDSAVE  --  Dynamics Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
c
c     skip an extra line to keep the output formating neat
c
      moddump = mod(istep,iprint)
      if (verbose .and. moddump.ne.0) then
         write (iout,210)
  210    format ()
      end if
      return
      end
