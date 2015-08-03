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
      subroutine readdynaux (idyn)
      use sizes
      use atoms
      use boxes
      use files
      use group
      use iounit
      use mdstuf
      use moldyn
      use rgddyn
      use iELSCF
      use polar
      implicit none
      integer i
      integer idyn,ndyn
      logical exist,opened,quit
      character*120 dynfile
      character*120 record
      real*8 temp
c
c
c     open the input file if it has not already been done
c
      inquire (unit=idyn,opened=opened)
      if (.not. opened) then
         dynfile = filename(1:leng)//'.auxdyn'
         call version (dynfile,'old')
         inquire (file=dynfile,exist=exist)
         if (exist) then
            open (unit=idyn,file=dynfile,status='old')
            rewind (unit=idyn)
         else
            write (iout,250)
  250       format (/,' READDYNAUX  --  Unable to Find the Auxiliary',
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
      read (idyn,260)
  260 format ()
      read (idyn,270)  record
  270 format (a120)
      read (record,*,err=400,end=400)  ndyn
      if (ndyn .ne. n) then
         write (iout,280)
  280    format (/,' READDYN  --  Auxiliary File has Incorrect',
     &              ' Number of Atoms')
         call fatal
      end if
      
      read (idyn,290)
  290 format ()
c      read (idynnik,300)  record
c  300 format (a120)
c      read (record,*,err=400,end=400)
c      read (idynnik,310)  record
c  310 format (a120)
c      read (record,*,err=400,end=400)
c      read (idynnik,320)
c  320 format ()
c      quit = .true.
c
c     get the auxiliary dipole positions, velocities and accelerations and real dipoles
c
      do i = 1, n
       read (idyn,330)  record
  330  format (a120)
       read (record,*,err=400,end=400) uind_aux(1,i),uind_aux(2,i),
     &                                 uind_aux(3,i)
      end do
      read (idyn,340)
  340 format ()
      do i = 1, n
         read (idyn,350)  record
  350    format (a120)
         read (record,*,err=400,end=400)  v_aux(1,i),v_aux(2,i),
     &                                    v_aux(3,i)
      end do
      read (idyn,360)
  360 format ()
      do i = 1, n
         read (idyn,370)  record
  370    format (a120)
         read (record,*,err=400,end=400)  a_aux(1,i),a_aux(2,i),
     &                                    a_aux(3,i)
      end do
      read (idyn,380)
  380 format ()
      do i = 1, n
         read (idyn,390)  record
  390    format (a120)
         read (record,*,err=400,end=400) uind(1,i),uind(2,i),uind(3,i)
      end do
      quit = .false.
  400 continue
      if (.not. opened)  close (unit=idyn)
c
c     report any error in reading the dynamics restart file
c
      if (quit) then
         write (iout,410)  i
  410    format (/,' READDYN  --  Error in Auxiliary Restart',
     &              ' File at Atom',i6)
         call fatal
      end if
      
      return
      end
