c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program diffuse  --  find liquid self-diffusion constant  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "diffuse" finds the self-diffusion constant for a homogeneous
c     liquid via the Einstein relation from a set of stored molecular
c     dynamics frames; molecular centers of mass are unfolded and mean
c     squared displacements are computed versus time separation
c
c     the estimate for the self-diffusion constant in 10-5 cm**2/sec
c     is printed in the far right column of output and can be checked
c     by plotting mean square displacements as a function of the time
c     separation
c
c     diffusion values for very large time separation are inaccurate
c     due to the small amount of data; the current version requires
c     an orthogonal unit cell
c
c
      program diffuse
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'molcul.i'
      integer maxframe
      parameter (maxframe=2000)
      integer i,j,k,m
      integer nframe,iframe
      integer iarc,freeunit
      integer start,stop
      integer step,skip
      integer ntime(maxframe)
      real*8 xmid,ymid,zmid
      real*8 xold,yold,zold
      real*8 xdiff,ydiff,zdiff
      real*8 dx,dy,dz,weigh
      real*8 tstep,dunits,delta
      real*8 xvalue,yvalue,zvalue
      real*8 rvalue,dvalue,counts
      real*8 xmsd(maxframe)
      real*8 ymsd(maxframe)
      real*8 zmsd(maxframe)
      real*8, allocatable :: xcm(:,:)
      real*8, allocatable :: ycm(:,:)
      real*8, allocatable :: zcm(:,:)
      logical exist,query
      character*120 arcfile
      character*120 record
      character*120 string
c
c
c     perform the standard initialization functions
c
      call initial
c
c     try to get a filename from the command line arguments
c
      call nextarg (arcfile,exist)
      if (exist) then
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end if
c
c     ask for the user specified input structure filename
c
      do while (.not. exist)
         write (iout,10)
   10    format (/,' Enter Coordinate Archive File Name :  ',$)
         read (input,20)  arcfile
   20    format (a120)
         call basefile (arcfile)
         call suffix (arcfile,'arc','old')
         inquire (file=arcfile,exist=exist)
      end do
c
c     read the first coordinate set in the archive
c
      iarc = freeunit ()
      open (unit=iarc,file=arcfile,status='old')
      call readxyz (iarc)
      rewind (unit=iarc)
c
c     get numbers of the coordinate frames to be processed
c
      start = 0
      stop = 0
      step = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=30,end=30)  start
         query = .false.
      end if
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  stop
      call nextarg (string,exist)
      if (exist)  read (string,*,err=30,end=30)  step
   30 continue
      if (query) then
         write (iout,40)
   40    format (/,' Numbers of First & Last Frame and Step',
     &              ' Increment :  ',$)
         read (input,50)  record
   50    format (a120)
         read (record,*,err=60,end=60)  start,stop,step
   60    continue
      end if
      if (stop .eq. 0)  stop = start
      if (step .eq. 0)  step = 1
c
c     get the time increment between frames in picoseconds
c
      tstep = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  tstep
   70 continue
      if (tstep .le. 0.0d0) then
         write (iout,80)
   80    format (/,' Enter the Time Increment in Picoseconds',
     &              ' [1.0] :  ',$)
         read (input,90)  tstep
   90    format (f20.0)
      end if
      if (tstep .le. 0.0d0)  tstep = 1.0d0
c
c     try to find the unit cell axis lengths in the keyfile
c
      call unitcell
c
c     get cell axis lengths from command line or interactively
c
      do while (xbox .eq. 0.0d0)
         write (iout,100)
  100    format (/,' Enter Unit Cell Axis Lengths :  ',$)
         read (input,110)  record
  110    format (a120)
         read (record,*,err=120,end=120)  xbox,ybox,zbox
  120    continue
         if (ybox .eq. 0.0d0)  ybox = xbox
         if (zbox .eq. 0.0d0)  zbox = xbox
      end do
c
c     check for possible non-orthogonal system unit cell
c
      if (.not. orthogonal) then
         write (iout,130)
  130    format (/,' DIFFUSE  --  Non-Orthogonal Unit Cell is',
     &              ' Not Supported')
         call fatal
      end if
c
c     set the half width values for the periodic box
c
      xbox2 = 0.5d0 * xbox
      ybox2 = 0.5d0 * ybox
      zbox2 = 0.5d0 * zbox
c
c     assign the atom parameters and count the molecules
c
      call field
      call katom
      call molecule
c
c     perform dynamic allocation of some local arrays
c
      allocate (xcm(nmol,maxframe))
      allocate (ycm(nmol,maxframe))
      allocate (zcm(nmol,maxframe))
c
c     get the archived coordinates for each frame in turn
c
      write (iout,140)
  140 format (/,' Reading the Coordinates Archive File :',/)
      nframe = 0
      iframe = start
      skip = start
      do while (iframe.ge.start .and. iframe.le.stop)
         skip = (skip-1) * (n+1)
         do j = 1, skip
            read (iarc,150,err=160,end=160)
  150       format ()
         end do
  160    continue
         iframe = iframe + step
         skip = step
         call readxyz (iarc)
         if (n .eq. 0)  goto 190
         nframe = nframe + 1
         if (mod(nframe,100) .eq. 0) then
            write (iout,170)  nframe
  170       format (4x,'Processing Coordinate Frame',i13)
         end if
         if (nframe .gt. maxframe) then
            write (iout,180)  maxframe
  180       format (/,' DIFFUSE  --  The Maximum of',i6,
     &                 ' Frames has been Exceeded')
            call fatal
         end if
c
c     unfold each molecule to get its corrected center of mass
c
         do i = 1, nmol
            xmid = 0.0d0
            ymid = 0.0d0
            zmid = 0.0d0
            do j = imol(1,i), imol(2,i)
               k = kmol(j)
               weigh = mass(k)
               xmid = xmid + x(k)*weigh
               ymid = ymid + y(k)*weigh
               zmid = zmid + z(k)*weigh
            end do
            weigh = molmass(i)
            xmid = xmid / weigh
            ymid = ymid / weigh
            zmid = zmid / weigh
            if (nframe .eq. 1) then
               xold = xmid
               yold = ymid
               zold = zmid
            else
               xold = xcm(i,nframe-1)
               yold = ycm(i,nframe-1)
               zold = zcm(i,nframe-1)
            end if
            dx = xmid - xold
            dy = ymid - yold
            dz = zmid - zold
            do while (dx .gt. xbox2)
               dx = dx - xbox
            end do
            do while (dx .lt. -xbox2)
               dx = dx + xbox
            end do
            do while (dy .gt. ybox2)
               dy = dy - ybox
            end do
            do while (dy .lt. -ybox2)
               dy = dy + ybox
            end do
            do while (dz .gt. zbox2)
               dz = dz - zbox
            end do
            do while (dz .lt. -zbox2)
               dz = dz + zbox
            end do
            xcm(i,nframe) = xold + dx
            ycm(i,nframe) = yold + dy
            zcm(i,nframe) = zold + dz
         end do
      end do
  190 continue
      close (unit=iarc)
      write (iout,200)  nframe
  200 format (/,' Total Number of Coordinate Frames :',i8)
c
c     increment the squared displacements for each frame pair
c
      do i = 1, nframe
         ntime(i) = 0
         xmsd(i) = 0.0d0
         ymsd(i) = 0.0d0
         zmsd(i) = 0.0d0
      end do
      do i = 1, nframe-1
         do j = i+1, nframe
            m = j - i
            ntime(m) = ntime(m) + 1
            do k = 1, nmol
               xdiff = xcm(k,j) - xcm(k,i)
               ydiff = ycm(k,j) - ycm(k,i)
               zdiff = zcm(k,j) - zcm(k,i)
               xmsd(m) = xmsd(m) + xdiff*xdiff
               ymsd(m) = ymsd(m) + ydiff*ydiff
               zmsd(m) = zmsd(m) + zdiff*zdiff
            end do
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xcm)
      deallocate (ycm)
      deallocate (zcm)
c
c     get mean squared displacements and convert units;
c     conversion is from sq. Ang/ps to 10-5 sq. cm/sec
c
      dunits = 10.0d0
      do i = 1, nframe-1
         counts = dble(nmol) * dble(ntime(i))
         xmsd(i) = xmsd(i) * (dunits/counts)
         ymsd(i) = ymsd(i) * (dunits/counts)
         zmsd(i) = zmsd(i) * (dunits/counts)
      end do
c
c     estimate the diffusion constant via the Einstein relation
c
      write (iout,210)
  210 format (/,' Mean Squared Diffusion Distance and Self-Diffusion',
     &           ' Constant :',
     &        //,5x,'Time Step',5x,'X MSD',7x,'Y MSD',7x,'Z MSD',
     &           7x,'R MSD',4x,'Diff Const',/)
      do i = 1, nframe-1
         delta = tstep * dble(i)
         xvalue = xmsd(i) / 2.0d0
         yvalue = ymsd(i) / 2.0d0
         zvalue = zmsd(i) / 2.0d0
         rvalue = (xmsd(i) + ymsd(i) + zmsd(i)) / 6.0d0
         dvalue = rvalue / delta
         write (iout,220)  delta,xvalue,yvalue,zvalue,rvalue,dvalue
  220    format (f12.2,4f12.2,f12.4)
      end do
c
c     perform any final tasks before program exit
c
      call final
      end
