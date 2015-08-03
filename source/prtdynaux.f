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
      subroutine prtdynaux
      use sizes
      use atoms
      use boxes
      use files
      use group
      use mdstuf
      use moldyn
      use rgddyn
      use titles
      use iELSCF
      use polar
      implicit none
      integer i,idyn
      integer freeunit
      logical exist
      character*2 atmc
      character*50 fstr
      character*120 dynfile
      
      idyn = freeunit ()
      dynfile = filename(1:leng)//'.auxdyn'
      inquire (file=dynfile,exist=exist)
      if (exist) then
         open (unit=idyn,file=dynfile,status='old')
         rewind (unit=idyn)
      else
         open (unit=idyn,file=dynfile,status='new')
      end if
      
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
      
      fstr = '('' Current Auxiliary Dipole Positions :'')'
      write (idyn,fstr(1:41))
      fstr = '(3d26.16)'
      do i = 1, n
        write (idyn,fstr(1:9)) uind_aux(1,i),uind_aux(2,i),uind_aux(3,i)
      end do
      fstr = '('' Current Auxiliary Dipole Velocities :'')'
      write (idyn,fstr(1:42))
      fstr = '(3d26.16)'
      do i = 1, n
        write (idyn,fstr(1:9)) v_aux(1,i),v_aux(2,i),v_aux(3,i)
      end do
      fstr =  '('' Current Auxiliary Dipole Accelerations :'')'
      write (idyn,fstr(1:45))
      fstr = '(3d26.16)'
      do i = 1, n
        write (idyn,fstr(1:9)) a_aux(1,i),a_aux(2,i),a_aux(3,i)
      end do
      fstr =  '('' Current Real Dipole Positions :'')'
      write (idyn,fstr(1:36))
      fstr = '(3d26.16)'
      do i = 1, n
        write (idyn,fstr(1:9)) uind(1,i),uind(2,i),uind(3,i)
      end do
c
c     close the dynamics trajectory restart file
c
      close (unit=idyn)
      return
      end
