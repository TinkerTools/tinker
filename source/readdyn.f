c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine readdyn  --  input of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "readdyn" get the positions, velocities and accelerations
c     for a molecular dynamics restart from an external disk file
c
c
      subroutine readdyn (idyn,idynaux,idynauxp)!ALBAUGH
      use sizes
      use atoms
      use boxes
      use files
      use group
      use iounit
      use mdstuf
      use moldyn
      use rgddyn
      use ielscf!ALBAUGH
      use mpole!ALBAUGH
      use polar!ALBAUGH
      implicit none
      integer i,idyn,ndyn
      integer idynaux,idynauxp!ALBAUGH
      logical exist,opened,quit
      character*240 dynfile
      character*240 record
c
c
c     open the input file if it has not already been done
c
      inquire (unit=idyn,opened=opened)
      if (.not. opened) then
         dynfile = filename(1:leng)//'.dyn'
         call version (dynfile,'old')
         inquire (file=dynfile,exist=exist)
         if (exist) then
            open (unit=idyn,file=dynfile,status='old')
            rewind (unit=idyn)
         else
            write (iout,10)
   10       format (/,' READDYN  --  Unable to Find the Dynamics',
     &                 ' Restart File')
            call fatal
         end if
      end if
c
c     initialize error handling during reading of the file
c
      i = 0
      quit = .true.
c
c     get the number of atoms and check for consistency
c
      read (idyn,20)
   20 format ()
      read (idyn,30)  record
   30 format (a240)
      read (record,*,err=240,end=240)  ndyn
      if (ndyn .ne. n) then
         write (iout,40)
   40    format (/,' READDYN  --  Restart File has Incorrect',
     &              ' Number of Atoms')
         call fatal
      end if
c
c     get the periodic box edge lengths and angles
c
      read (idyn,50)
   50 format ()
      read (idyn,60)  record
   60 format (a240)
      read (record,*,err=240,end=240)  xbox,ybox,zbox
      read (idyn,70)  record
   70 format (a240)
      read (record,*,err=240,end=240)  alpha,beta,gamma
c
c     set the box volume and additional periodic box values
c
      call lattice
c
c     get rigid body positions, translational and angular velocities
c
      if (integrate .eq. 'RIGIDBODY') then
         read (idyn,80)
   80    format ()
         do i = 1, n
            read (idyn,90)  record
   90       format (a240)
            read (record,*,err=240,end=240)  x(i),y(i),z(i)
         end do
         read (idyn,100)
  100    format ()
         do i = 1, ngrp
            read (idyn,110)  record
  110       format (a240)
            read (record,*,err=240,end=240)  vcm(1,i),vcm(2,i),vcm(3,i)
         end do
         read (idyn,120)
  120    format ()
         do i = 1, ngrp
            read (idyn,130)  record
  130       format (a240)
            read (record,*,err=240,end=240)  wcm(1,i),wcm(2,i),wcm(3,i)
         end do
         read (idyn,140)
  140    format ()
         do i = 1, ngrp
            read (idyn,150)  record
  150       format (a240)
            read (record,*,err=240,end=240)  lm(1,i),lm(2,i),lm(3,i)
         end do
c
c     get the atomic positions, velocities and accelerations
c
      else
         read (idyn,160)
  160    format ()
         do i = 1, n
            read (idyn,170)  record
  170       format (a240)
            read (record,*,err=240,end=240)  x(i),y(i),z(i)
         end do
         read (idyn,180)
  180    format ()
         do i = 1, n
            read (idyn,190)  record
  190       format (a240)
            read (record,*,err=240,end=240)  v(1,i),v(2,i),v(3,i)
         end do
         read (idyn,200)
  200    format ()
         do i = 1, n
            read (idyn,210)  record
  210       format (a240)
            read (record,*,err=240,end=240)  a(1,i),a(2,i),a(3,i)
         end do
         read (idyn,220)
  220    format ()
         do i = 1, n
            read (idyn,230)  record
  230       format (a240)
            read (record,*,err=240,end=240)  aalt(1,i),aalt(2,i),
     &                                       aalt(3,i)
         end do
      end if
      quit = .false.
  240 continue
      if (.not. opened)  close (unit=idyn)
c
c     report any error in reading the dynamics restart file
c
      if (quit) then
         write (iout,250)  i
  250    format (/,' READDYN  --  Error in Dynamics Restart',
     &              ' File at Atom',i6)
         call fatal
      end if
      
      if (use_ielscf .or. use_iel0scf) then
c
c
c     open the input file if it has not already been done
c
         inquire (unit=idynaux,opened=opened)
         if (.not. opened) then
             dynfile = filename(1:leng)//'.auxdyn'
           call version (dynfile,'old')
            inquire (file=dynfile,exist=exist)
            if (exist) then
               open (unit=idynaux,file=dynfile,status='old')
               rewind (unit=idynaux)
            else
               write (iout,260)
  260          format (/,' READDYN  --  Unable to Find the Auxiliary',
     &                 ' Restart File')
               call fatal
            end if
         end if
c
c     initialize error handling during reading of the file
c
         i = 0
         quit = .true.
c
c     get the number of atoms and check for consistency
c
         read (idynaux,270)
  270    format ()
         read (idynaux,280)  record
  280    format (a240)
         read (record,*,err=410,end=410)  ndyn
         if (ndyn .ne. npole) then
            write (iout,290)
  290       format (/,' READDYN  --  Restart File has Incorrect',
     &              ' Number of Auxiliaries')
             call fatal
         end if
c
c     get the periodic box edge lengths and angles
c
         read (idynaux,300)
  300    format ()
         read (idynaux,310)  record
  310    format (a240)
         read (record,*,err=410,end=410)
         read (idynaux,320)  record
  320    format (a240)
         read (record,*,err=410,end=410)

         read (idynaux,330)
  330    format ()
         do i = 1, npole
            read (idynaux,340)  record
  340       format (a240)
            read (record,*,err=410,end=410)  uaux(1,i),uaux(2,i),
     &                                        uaux(3,i)
         end do
         read (idynaux,350)
  350    format ()
         do i = 1, npole
            read (idynaux,360)  record
  360       format (a240)
            read (record,*,err=410,end=410)  vaux(1,i),vaux(2,i),
     &                                        vaux(3,i)
         end do
         read (idynaux,370)
  370    format ()
         do i = 1, n
            read (idynaux,380)  record
  380       format (a240)
            read (record,*,err=410,end=410)  aaux(1,i),aaux(2,i),
     &                                        aaux(3,i)
         end do
         read (idynaux,390)
  390    format ()
         do i = 1, n
            read (idynaux,400)  record
  400       format (a240)
            read (record,*,err=410,end=410)  uind(1,i),uind(2,i),
     &                                         uind(3,i)
         end do
         
         quit = .false.
  410    continue
         if (.not. opened)  close (unit=idynaux)
c
c     report any error in reading the dynamics restart file
c
         if (quit) then
            write (iout,420)  i
  420       format (/,' READDYN  --  Error in Auxiliary Restart',
     &              ' File at Atom',i6)
            call fatal
         end if

c
c
c     open the input file if it has not already been done
c
         inquire (unit=idynauxp,opened=opened)
         if (.not. opened) then
             dynfile = filename(1:leng)//'.auxpdyn'
           call version (dynfile,'old')
            inquire (file=dynfile,exist=exist)
            if (exist) then
               open (unit=idynauxp,file=dynfile,status='old')
               rewind (unit=idynauxp)
            else
               write (iout,430)
  430          format (/,' READDYN  --  Unable to Find the P-Auxiliary',
     &                 ' Restart File')
               call fatal
            end if
         end if
c
c     initialize error handling during reading of the file
c
         i = 0
         quit = .true.
c
c     get the number of atoms and check for consistency
c
         read (idynauxp,440)
  440    format ()
         read (idynauxp,450)  record
  450    format (a240)
         read (record,*,err=580,end=580)  ndyn
         if (ndyn .ne. npole) then
            write (iout,460)
  460       format (/,' READDYN  --  Restart File has Incorrect',
     &              ' Number of P-Auxiliaries')
             call fatal
         end if
c
c     get the periodic box edge lengths and angles
c
         read (idynauxp,470)
  470    format ()
         read (idynauxp,480)  record
  480    format (a240)
         read (record,*,err=580,end=580)
         read (idynauxp,490)  record
  490    format (a240)
         read (record,*,err=580,end=580)

         read (idynauxp,500)
  500    format ()
         do i = 1, npole
            read (idynauxp,510)  record
  510       format (a240)
            read (record,*,err=580,end=580)  upaux(1,i),upaux(2,i),
     &                                        upaux(3,i)
         end do
         read (idynauxp,520)
  520    format ()
         do i = 1, npole
            read (idynauxp,530)  record
  530       format (a240)
            read (record,*,err=580,end=580)  vpaux(1,i),vpaux(2,i),
     &                                        vpaux(3,i)
         end do
         read (idynauxp,540)
  540    format ()
         do i = 1, n
            read (idynauxp,550)  record
  550       format (a240)
            read (record,*,err=580,end=580)  apaux(1,i),apaux(2,i),
     &                                        apaux(3,i)
         end do
         read (idynauxp,560)
  560    format ()
         do i = 1, n
            read (idynauxp,570)  record
  570       format (a240)
            read (record,*,err=580,end=580)  uinp(1,i),uinp(2,i),
     &                                         uinp(3,i)
         end do
         
         quit = .false.
  580    continue
         if (.not. opened)  close (unit=idynauxp)
c
c     report any error in reading the dynamics restart file
c
         if (quit) then
            write (iout,590)  i
  590       format (/,' READDYN  --  Error in P-Auxiliary Restart',
     &              ' File at Atom',i6)
            call fatal
         end if  
      end if
      
      return
      end
