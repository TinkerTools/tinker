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
      use sizes
      use atoms
      use boxes
      use files
      use group
      use mdstuf
      use moldyn
      use rgddyn
      use titles
      use ielscf!ALBAUGH
      use mpole!ALBAUGH
      use polar!ALBAUGH
      implicit none
      integer i,idyn
      integer idynaux,idynauxp!ALBAUGH
      integer freeunit
      logical exist
      character*2 atmc
      character*50 fstr
      character*240 dynfile
      character*240 dynfileaux,dynfileauxp!ALBAUGH
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
      if (use_ielscf .or. use_iel0scf) then!ALBAUGH
         idynaux = freeunit ()
         dynfileaux = filename(1:leng)//'.auxdyn'
         inquire (file=dynfileaux,exist=exist)
         if (exist) then
            open (unit=idynaux,file=dynfileaux,status='old')
            rewind (unit=idynaux)
         else
            open (unit=idynaux,file=dynfileaux,status='new')
         end if
         
         idynauxp = freeunit ()
         dynfileauxp = filename(1:leng)//'.auxpdyn'
         inquire (file=dynfileauxp,exist=exist)
         if (exist) then
            open (unit=idynauxp,file=dynfileauxp,status='old')
            rewind (unit=idynauxp)
         else
            open (unit=idynauxp,file=dynfileauxp,status='new')
         end if
      end if
c
c     save the number of atoms and the title string
c
      fstr = '('' Number of Atoms and Title :'')'
      write (idyn,fstr(1:32))
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (idyn,fstr(1:4))  n
      else
         fstr = '('//atmc//',2x,a)'
         write (idyn,fstr(1:9))  n,title(1:ltitle)
      end if
      
      if (use_ielscf .or. use_iel0scf) then!ALBAUGH
         fstr = '('' Number of Atoms and Title :'')'
         write (idynaux,fstr(1:32))
         atmc = 'i6'
         if (n .ge. 100000)  atmc = 'i7'
         if (n .ge. 1000000)  atmc = 'i8'
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (idynaux,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (idynaux,fstr(1:9))  n,title(1:ltitle)
         end if
         
         fstr = '('' Number of Atoms and Title :'')'
         write (idynauxp,fstr(1:32))
         atmc = 'i6'
         if (n .ge. 100000)  atmc = 'i7'
         if (n .ge. 1000000)  atmc = 'i8'
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (idynauxp,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (idynauxp,fstr(1:9))  n,title(1:ltitle)
         end if
      end if
c
c     save the periodic box edge lengths and angles
c
      fstr = '('' Periodic Box Dimensions :'')'
      write (idyn,fstr(1:30))
      fstr = '(3d26.16)'
      write (idyn,fstr(1:9))  xbox,ybox,zbox
      write (idyn,fstr(1:9))  alpha,beta,gamma
      
      if (use_ielscf .or. use_iel0scf) then!ALBAUGH
         fstr = '('' Periodic Box Dimensions :'')'
         write (idynaux,fstr(1:30))
         fstr = '(3d26.16)'
         write (idynaux,fstr(1:9))  xbox,ybox,zbox
         write (idynaux,fstr(1:9))  alpha,beta,gamma
         
         
         fstr = '('' Periodic Box Dimensions :'')'
         write (idynauxp,fstr(1:30))
         fstr = '(3d26.16)'
         write (idynauxp,fstr(1:9))  xbox,ybox,zbox
         write (idynauxp,fstr(1:9))  alpha,beta,gamma
      end if
c
c     save rigid body positions, translational and angular velocities
c
      if (integrate .eq. 'RIGIDBODY') then
         fstr = '('' Current Atomic Positions :'')'
         write (idyn,fstr(1:31))
         fstr = '(3d26.16)'
         do i = 1, n
            write (idyn,fstr(1:9))  x(i),y(i),z(i)
         end do
         fstr = '('' Current Translational Velocities :'')'
         write (idyn,fstr(1:39))
         fstr = '(3d26.16)'
         do i = 1, ngrp
            write (idyn,fstr(1:9))  vcm(1,i),vcm(2,i),vcm(3,i)
         end do
         fstr = '('' Current Angular Velocities :'')'
         write (idyn,fstr(1:33))
         fstr = '(3d26.16)'
         do i = 1, ngrp
            write (idyn,fstr(1:9))  wcm(1,i),wcm(2,i),wcm(3,i)
         end do
         fstr = '('' Current Angular Momenta :'')'
         write (idyn,fstr(1:30))
         fstr = '(3d26.16)'
         do i = 1, ngrp
            write (idyn,fstr(1:9))  lm(1,i),lm(2,i),lm(3,i)
         end do
c
c     save the atomic positions, velocities and accelerations
c
      else
         fstr = '('' Current Atomic Positions :'')'
         write (idyn,fstr(1:31))
         fstr = '(3d26.16)'
         do i = 1, n
            write (idyn,fstr(1:9))  x(i),y(i),z(i)
         end do
         fstr = '('' Current Atomic Velocities :'')'
         write (idyn,fstr(1:32))
         fstr = '(3d26.16)'
         do i = 1, n
            write (idyn,fstr(1:9))  v(1,i),v(2,i),v(3,i)
         end do
         fstr =  '('' Current Atomic Accelerations :'')'
         write (idyn,fstr(1:36))
         fstr = '(3d26.16)'
         do i = 1, n
            write (idyn,fstr(1:9))  a(1,i),a(2,i),a(3,i)
         end do
         fstr =  '('' Alternate Atomic Accelerations :'')'
         write (idyn,fstr(1:38))
         fstr = '(3d26.16)'
         if (integrate .eq. 'VERLET') then
            do i = 1, n
               write (idyn,fstr(1:9))  a(1,i),a(2,i),a(3,i)
            end do
         else
            do i = 1, n
               write (idyn,fstr(1:9))  aalt(1,i),aalt(2,i),aalt(3,i)
            end do
         end if
      end if
      
      if (use_iel0scf .or. use_ielscf) then!ALBAUGH
         fstr = '('' Current Auxiliary Dipole Positions :'')'
         write (idynaux,fstr(1:41))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynaux,fstr(1:9))  uaux(1,i),uaux(2,i),uaux(3,i)
         end do
         fstr = '('' Current Auxiliary Dipole Velocities :'')'
         write (idynaux,fstr(1:42))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynaux,fstr(1:9))  vaux(1,i),vaux(2,i),vaux(3,i)
         end do
         fstr =  '('' Current Auxiliary Dipole Accelerations :'')'
         write (idynaux,fstr(1:45))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynaux,fstr(1:9))  aaux(1,i),aaux(2,i),aaux(3,i)
         end do
         fstr =  '('' Current Real Dipole Accelerations :'')'
         write (idynaux,fstr(1:40))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynaux,fstr(1:9))  uind(1,i),uind(2,i),uind(3,i)
         end do
         
         fstr = '('' Current P-Auxiliary Dipole Positions :'')'
         write (idynauxp,fstr(1:43))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynauxp,fstr(1:9)) upaux(1,i),upaux(2,i),upaux(3,i)
         end do
         fstr = '('' Current P-Auxiliary Dipole Velocities :'')'
         write (idynauxp,fstr(1:44))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynauxp,fstr(1:9)) vpaux(1,i),vpaux(2,i),vpaux(3,i)
         end do
         fstr =  '('' Current P-Auxiliary Dipole Accelerations :'')'
         write (idynauxp,fstr(1:47))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynauxp,fstr(1:9)) apaux(1,i),apaux(2,i),apaux(3,i)
         end do
         fstr =  '('' Current Real P-Dipole Accelerations :'')'
         write (idynauxp,fstr(1:42))
         fstr = '(3d26.16)'
         do i = 1, npole
            write (idynauxp,fstr(1:9)) uinp(1,i),uinp(2,i),uinp(3,i)
         end do
      end if
c
c     close the dynamics trajectory restart file
c
      close (unit=idyn)
      if(use_ielscf .or. use_iel0scf) then!ALBAUGH
         close (unit=idynaux)
         close (unit=idynauxp)
      end if
      return
      end
