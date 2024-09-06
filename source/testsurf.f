c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2023  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine testsurf  --  find & compare area and volume  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "testsurf" finds the accessible surface area, excluded volume
c     and their derivatives for a molecular system via the methods of
c     Tim Richmond, Michael Connolly, Craig Kundrot and Patrice Koehl
c
c
      program testsurf
      use atomid
      use atoms
      use files
      use iounit
      use kvdws
      use nonpol
      use ptable
      use vdwpot
      implicit none
      integer i,icrd
      integer nsize,nfudge
      integer freeunit
      real*8 surf,vol
      real*8 probe,rmax
      real*8 reentrant
      real*8 wall,cpu
      real*8, allocatable :: rsolv(:)
      real*8, allocatable :: weight(:)
      real*8, allocatable :: asurf(:)
      real*8, allocatable :: avol(:)
      real*8, allocatable :: dsurf(:,:)
      real*8, allocatable :: dvol(:,:)
      logical exist,query
      logical docrd
      logical doderiv,dovol
      character*240 crdfile
      character*240 string
c
c
c     set up the structure and values for the computation;
c     solute radii can be changed via the keyword mechanism
c
      call initial
      call getxyz
      call active
      call field
      call katom
      call kvdw
c
c     get probe radius for accessible area/excluded volume
c
      probe = 0.0d0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  probe
         query = .false.
      end if
   10 continue
      if (query) then
         probe = -1.0d0
         write (iout,20)
   20    format (/,' Enter a Value for the Probe Radius',
     &              ' [1.4 Ang] :  ',$)
         read (input,30)  string
   30    format (a240)
         read (string,*,err=40,end=40)  probe
         goto 50
   40    continue
         probe = 1.4d0
   50    continue
      end if
c
c     print out the total number of atoms
c
      write (iout,60)
   60 format (/,' Alternative Surface Area & Volume Methods')
      write (iout,70)  n,probe
   70 format (/,' Number of Atoms :',15x,i8,
     &        /,' Probe Size :',16x,f12.4)
c
c     perform dynamic allocation of some local arrays
c
      nfudge = 10
      nsize = n + nfudge
      allocate (rsolv(nsize))
      allocate (weight(nsize))
      allocate (asurf(nsize))
      allocate (avol(nsize))
      allocate (dsurf(3,nsize))
      allocate (dvol(3,nsize))
c
c     if all radii are zero then switch to generic vdw radii
c
      rmax = 0.0d0
      do i = 1, n
         rmax = rad(i)
         if (rmax .gt. 0.0d0)  goto 80
      end do
   80 continue
      if (rmax .eq. 0.0d0) then
         write (iout,90)
   90    format (/,' Atomic Radii not Set, Using Generic VDW Values')
         do i = 1, n
            rad(i) = vdwrad(atomic(i))
         end do
      end if
c
c     set radii to use for surface area and volume calculation
c
      do i = 1, n
         if (vdwindex .eq. 'CLASS') then
            rsolv(i) = rad(class(i))
         else
            rsolv(i) = rad(type(i))
         end if
         weight(i) = 1.0d0
      end do
c
c     initialize variables for Richmond and Connolly routines
c
      surf = 0.0d0
      vol = 0.0d0
      reentrant = 0.0d0
      do i = 1, n
         asurf(i) = 0.0d0
         dsurf(1,i) = 0.0d0
         dsurf(2,i) = 0.0d0
         dsurf(3,i) = 0.0d0
      end do
c
c     compute accessible surface area via Richmond method
c
      surf = 0.0d0
      write (iout,100)
  100 format (/,' Timothy Richmond Accessible Surface Area Method :')
      call settime
      call richmond (n,x,y,z,rsolv,weight,probe,surf,asurf)
      call gettime (wall,cpu)
      write (iout,110)  cpu,wall
  110 format (/,' CPU and Wall Times :',8x,2f12.4) 
      write (iout,120)  surf
  120 format (/,' Total Surface Area :',8x,f12.4)
c
c     compute accessible surface and derivatives via Richmond
c
      surf = 0.0d0
      write (iout,130)
  130 format (/,' Timothy Richmond Surface Area Derivative Method :')
      call settime
      call richmond1 (n,x,y,z,rsolv,weight,probe,surf,asurf,dsurf)
      call gettime (wall,cpu)
      write (iout,140)  cpu,wall
  140 format (/,' CPU and Wall Times :',8x,2f12.4)
      write (iout,150)  surf
  150 format (/,' Total Surface Area :',8x,f12.4)
      write (iout,160)
  160 format (/,' Surface Area Derivatives :  (First Ten Atoms)',
     &        //,5x,'Atom',11x,'dAx',7x,'dAy',7x'dAz',/)
      do i = 1, min(10,n)
         write (iout,170)  i,dsurf(1,i),dsurf(2,i),dsurf(3,i)
  170    format (i8,6x,3f10.4)
      end do
c
c     compute surface area and excluded volume via Connolly
c
      surf = 0.0d0
      vol = 0.0d0
      write (iout,180)
  180 format (/,' Michael Connolly Molecular Area-Volume Method :')
      call settime
      call connolly (n,x,y,z,rsolv,probe,reentrant,surf,vol)
      call gettime (wall,cpu)
      write (iout,190)  cpu,wall
  190 format (/,' CPU and Wall Times :',8x,2f12.4)
      write (iout,200)  surf
  200 format (/,' Total Surface Area :',8x,f12.4)
      write (iout,210)  vol
  210 format (/,' Total Excluded Volume :',5x,f12.4)
c
c     compute excluded volume derivatives via Kundrot method
c
      do i = 1, n
         dvol(1,i) = 0.0d0
         dvol(2,i) = 0.0d0
         dvol(3,i) = 0.0d0
      end do
      write (iout,220)
  220 format (/,' Craig Kundrot Excluded Volume Derivative Method :')
      call settime
      call kundrot1 (n,x,y,z,rsolv,probe,dvol)
      call gettime (wall,cpu)
      write (iout,230)  cpu,wall
  230 format (/,' CPU and Wall Times :',8x,2f12.4)
      write (iout,240)
  240 format (/,' Excluded Volume Derivatives :  (First Ten Atoms)',
     &        //,5x,'Atom',11x,'dVx',7x,'dVy',7x'dVz',/)
      do i = 1, min(10,n)
         write (iout,250)  i,dvol(1,i),dvol(2,i),dvol(3,i)
  250    format (i8,6x,3f10.4)
      end do

c
c     initialize variables for Koehl UnionBall routines
c
      doderiv = .true.
      dovol = .true.
      surf = 0.0d0
      vol = 0.0d0
      do i = 1, n
         asurf(i) = 0.0d0
         avol(i) = 0.0d0
         dsurf(1,i) = 0.0d0
         dsurf(2,i) = 0.0d0
         dsurf(3,i) = 0.0d0
         dvol(1,i) = 0.0d0
         dvol(2,i) = 0.0d0
         dvol(3,i) = 0.0d0
      end do
c
c     print out structure in UnionBall coordinate format
c
      docrd = .false.
      if (docrd) then
         icrd = freeunit ()
         crdfile = filename(1:leng)//'.crd'
         call version (crdfile,'new')
         open (unit=icrd,file=crdfile,status='new')
         write (icrd,260)  n
  260    format (i8,/)
         do i = 1, n
            write (icrd,270)  i,x(i),y(i),z(i),rsolv(i)
  270       format (i8,3f14.6,f12.4)
         end do
         close (unit=icrd)
      end if
c
c     compute area, volume and derivatives via UnionBall
c
      write (iout,280)
  280 format (/,' Patrice Koehl UnionBall Alpha Shape Method :')
      call settime
      call unionball (n,x,y,z,rsolv,weight,probe,doderiv,dovol,
     &                   surf,vol,asurf,avol,dsurf,dvol)
      call gettime (wall,cpu)
      write (iout,290)  cpu,wall
  290 format (/,' CPU and Wall Times :',8x,2f12.4)
      write (iout,300)  surf
  300 format (/,' Total Surface Area :',8x,f12.4)
      write (iout,310)  vol
  310 format (/,' Total Excluded Volume :',5x,f12.4)
      write (iout,320)
  320 format (/,' Surface Area & Volume Derivatives :',
     &           '  (First Ten Atoms)',
     &        //,5x,'Atom',11x,'dAx',7x,'dAy',7x'dAz',
     &           7x,'dVx',7x,'dVy',7x,'dVz',/)
      do i = 1, min(10,n)
         write (iout,330)  i,dsurf(1,i),dsurf(2,i),dsurf(3,i),
     &                     dvol(1,i),dvol(2,i),dvol(3,i)
  330    format (i8,6x,6f10.4)
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (rsolv)
      deallocate (weight)
      deallocate (asurf)
      deallocate (avol)
      deallocate (dsurf)
      deallocate (dvol)
      end
