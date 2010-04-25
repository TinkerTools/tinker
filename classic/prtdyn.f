c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdyn  --  output of MD restart information  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdyn" writes out the information needed to restart a
c     molecular dynamics trajectory to an external disk file
c
c
      subroutine prtdyn
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'files.i'
      include 'group.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'rgddyn.i'
      include 'titles.i'
      integer i,idyn
      integer freeunit
      logical exist
      character*120 dynfile
c
c
c     update an existing restart file or open a new one
c
      idyn = freeunit ()
      dynfile = filename(1:leng)//'.dyn'
      inquire (file=dynfile,exist=exist)
      if (exist) then
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
      else
         open (unit=idyn,file=dynfile,status='new')
      end if
c
c     save the number of atoms and the title string
c
      write (idyn,10)
   10 format (' Number of Atoms and Title :')
      if (ltitle .eq. 0) then
         write (idyn,20)  n
   20    format (i6)
      else
         write (idyn,30)  n,title(1:ltitle)
   30    format (i6,2x,a)
      end if
c
c     save the periodic box edge lengths and angles
c
      write (idyn,40)
   40 format (' Periodic Box Dimensions :')
      write (idyn,50)  xbox,ybox,zbox
   50 format (3d26.16)
      write (idyn,60)  alpha,beta,gamma
   60 format (3d26.16)
c
c     save rigid body positions, translational and angular velocities
c
      if (integrate .eq. 'RIGIDBODY') then
         write (idyn,70)
   70    format (' Current Atomic Positions :')
         do i = 1, n
            write (idyn,80)  x(i),y(i),z(i)
   80       format (3d26.16)
         end do
         write (idyn,90)
   90    format (' Current Translational Velocities :')
         do i = 1, ngrp
            write (idyn,100)  vcm(1,i),vcm(2,i),vcm(3,i)
  100       format (3d26.16)
         end do
         write (idyn,110)
  110    format (' Current Angular Velocities :')
         do i = 1, ngrp
            write (idyn,120)  wcm(1,i),wcm(2,i),wcm(3,i)
  120       format (3d26.16)
         end do
         write (idyn,130)
  130    format (' Current Angular Momenta :')
         do i = 1, ngrp
            write (idyn,140)  lm(1,i),lm(2,i),lm(3,i)
  140       format (3d26.16)
         end do
c
c     save the atomic positions, velocities and accelerations
c
      else
         write (idyn,150)
  150    format (' Current Atomic Positions :')
         do i = 1, n
            write (idyn,160)  x(i),y(i),z(i)
  160       format (3d26.16)
         end do
         write (idyn,170)
  170    format (' Current Atomic Velocities :')
         do i = 1, n
            write (idyn,180)  v(1,i),v(2,i),v(3,i)
  180       format (3d26.16)
         end do
         write (idyn,190)
  190    format (' Current Atomic Accelerations :')
         do i = 1, n
            write (idyn,200)  a(1,i),a(2,i),a(3,i)
  200       format (3d26.16)
         end do
         write (idyn,210)
  210    format (' Previous Atomic Accelerations :')
         do i = 1, n
            write (idyn,220)  aold(1,i),aold(2,i),aold(3,i)
  220       format (3d26.16)
         end do
      end if
c
c     close the dynamics trajectory restart file
c
      close (unit=idyn)
      return
      end
