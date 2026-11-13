c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_verlet  --  short MD trajectory tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_verlet" runs twenty steps of "dynamic" on the argon box
c     with the velocity Verlet integrator at constant energy and with
c     the Nose-Hoover integrator at constant pressure, as in the
c     tinker-gpu verlet.cpp and nh.cpp tests, and checks the restart
c     file written at the end against one made by "dynamic"
c
c
      subroutine test_verlet
      implicit none
c
c
      call test_verlet_arbox
      call test_verlet_nosehoover
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_verlet_arbox  --  NVE velocity Verlet  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_verlet_arbox
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_verlet_arbox')
c
c
      if (skiptest(tname,'arbox verlet nve',4))  return
      call mdprep ('kinetic','verlet_arbox','arbox')
      call pushdir ('file/verlet_arbox')
      call tnist_append ('arbox.key','integrator verlet')
      call run_prog ('dynamic','arbox 20 1 0.02 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,tname//' status')
         call refpath ('verlet','arbox_verlet.dyn',rpath)
         call assert_dynfile (rpath,'arbox.dyn',1.0d-6,1.0d-6,tname)
      end if
      call popdir
      call tnist_clean ('verlet_arbox')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_verlet_nosehoover  --  NPT Nose-Hoover  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_verlet_nosehoover
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_verlet_nosehoover')
c
c
      if (skiptest(tname,'arbox nosehoover npt',4))  return
      call mdprep ('kinetic','verlet_nosehoover','arbox')
      call pushdir ('file/verlet_nosehoover')
      call tnist_append ('arbox.key','integrator nose-hoover')
      call tnist_append ('arbox.key','thermostat nose-hoover')
      call tnist_append ('arbox.key','barostat nose-hoover')
      call run_prog ('dynamic','arbox 20 1 0.02 4 298 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,tname//' status')
         call refpath ('verlet','arbox_nosehoover.dyn',rpath)
         call assert_dynfile (rpath,'arbox.dyn',1.0d-6,1.0d-6,tname)
      end if
      call popdir
      call tnist_clean ('verlet_nosehoover')
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine mdprep  --  scratch dynamics fixture  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "mdprep" copies the coordinates, keyfile and restart file of a
c     fixture into a fresh scratch directory under the file tree
c
c
      subroutine mdprep (src,work,base)
      implicit none
      character*(*) src,work,base
      character*512 cmd
c
c
      call pushdir ('file/'//src)
      cmd = 'rm -rf ../'//work//' ; mkdir -p ../'//work//' ; cp '//
     &      base//'.xyz '//base//'.key '//base//'.dyn ../'//work//'/'
      call execute_command_line (cmd)
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine assert_dynfile  --  compare MD restart files  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "assert_dynfile" compares the box, positions and velocities of
c     a dynamics restart file against a reference restart file
c
c
      subroutine assert_dynfile (reffile,outfile,epsq,epsv,label)
      implicit none
      integer nr,no
      real*8 epsq,epsv,d
      character*(*) reffile,outfile,label
c
c
      call dynatoms (reffile,nr)
      call dynatoms (outfile,no)
      call assert_int (no,nr,label//' atoms')
      if (no .ne. nr)  return
      call dyndiff (reffile,outfile,'Periodic Box Dimensions',2,d)
      call assert_real (d,0.0d0,epsq,label//' box')
      call dyndiff (reffile,outfile,'Current Atomic Positions',nr,d)
      call assert_real (d,0.0d0,epsq,label//' positions')
      call dyndiff (reffile,outfile,'Current Atomic Velocities',nr,d)
      call assert_real (d,0.0d0,epsv,label//' velocities')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine dynatoms  --  atom count of restart file  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine dynatoms (fname,nat)
      implicit none
      integer nat,iunit,ios
      integer freeunit
      character*(*) fname
      character*240 line
c
c
      nat = -1
      iunit = freeunit ()
      open (unit=iunit,file=fname,status='old',iostat=ios)
      if (ios .ne. 0)  return
      read (iunit,'(a)',iostat=ios)  line
      read (iunit,*,iostat=ios)  nat
      if (ios .ne. 0)  nat = -1
      close (unit=iunit)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine dyndiff  --  largest restart section diff  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "dyndiff" finds the largest absolute difference between the
c     rows that follow a section label in two restart files; a
c     missing section gives a huge value
c
c
      subroutine dyndiff (reffile,outfile,title,nrow,dmax)
      implicit none
      integer i,j,nrow
      integer ir,io,iosr,ioso
      integer freeunit
      real*8 dmax
      real*8 rv(3),ov(3)
      character*(*) reffile,outfile,title
      character*240 lr,lo
c
c
      dmax = 1.0d30
      ir = freeunit ()
      open (unit=ir,file=reffile,status='old',iostat=iosr)
      io = freeunit ()
      open (unit=io,file=outfile,status='old',iostat=ioso)
      if (iosr.ne.0 .or. ioso.ne.0)  goto 90
      do while (.true.)
         read (ir,'(a)',iostat=iosr)  lr
         if (iosr .ne. 0)  goto 90
         if (index(lr,title) .gt. 0)  goto 10
      end do
   10 continue
      do while (.true.)
         read (io,'(a)',iostat=ioso)  lo
         if (ioso .ne. 0)  goto 90
         if (index(lo,title) .gt. 0)  goto 20
      end do
   20 continue
      dmax = 0.0d0
      do i = 1, nrow
         read (ir,*,iostat=iosr)  (rv(j),j=1,3)
         read (io,*,iostat=ioso)  (ov(j),j=1,3)
         if (iosr.ne.0 .or. ioso.ne.0) then
            dmax = 1.0d30
            goto 90
         end if
         do j = 1, 3
            dmax = max(dmax,abs(ov(j)-rv(j)))
         end do
      end do
   90 continue
      close (unit=ir)
      close (unit=io)
      return
      end
