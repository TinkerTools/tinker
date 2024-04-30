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
   70 format (/,' Number of Atoms :',12x,i8,
     &        /,' Probe Size :',11x,f14.4)
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
         if (rmax .gt. 0.0d0)  goto 99
      end do
   99 continue
      if (rmax .eq. 0.0d0) then
         write (iout,98)
   98    format (/,' Atomic Radii not Set, Using Generic VDW Values')
         do i = 1, n
            rad(i) = vdwrad(atomic(i))
         end do
      end if
c
c     set radii to use for surface area and volume calculation
c
      do i = 1, n
         if (vdwindex .eq. 'CLASS') then
            rsolv(i) = rad(class(i)) + cavoff
         else  if (vdwindex .eq. 'TYPE ') then
            rsolv(i) = rad(type(i)) + cavoff
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
      write (iout,80)
   80 format (/,' Timothy Richmond Accessible Surface Area Method :')
      call settime
      call richmond (n,x,y,z,rsolv,weight,probe,surf,asurf)
      call gettime (wall,cpu)
      write (iout,90)  cpu,wall
   90 format (/,' CPU and Wall Times :',5x,2f12.4) 
      write (iout,100)  surf
  100 format (/,' Total Surface Area :',5x,f12.4)
c
c     compute accessible surface and derivatives via Richmond
c
      surf = 0.0d0
      write (iout,110)
  110 format (/,' Timothy Richmond Surface Area Derivative Method :')
      call settime
      call richmond1 (n,x,y,z,rsolv,weight,probe,surf,asurf,dsurf)
      call gettime (wall,cpu)
      write (iout,120)  cpu,wall
  120 format (/,' CPU and Wall Times :',5x,2f12.4)
      write (iout,130)  surf
  130 format (/,' Total Surface Area :',5x,f12.4)
      write (iout,140)
  140 format (/,' Surface Area Derivatives :  (First Ten Atoms)',/)
      do i = 1, min(10,n)
         write (iout,150)  i,dsurf(1,i),dsurf(2,i),dsurf(3,i)
  150    format (i8,6x,3f10.4)
      end do
c
c     compute surface area and excluded volume via Connolly
c
      surf = 0.0d0
      vol = 0.0d0
      write (iout,160)
  160 format (/,' Michael Connolly Molecular Area-Volume Method :')
      call settime
      call connolly (n,x,y,z,rsolv,probe,reentrant,surf,vol)
      call gettime (wall,cpu)
      write (iout,170)  cpu,wall
  170 format (/,' CPU and Wall Times :',5x,2f12.4)
      write (iout,180)  surf
  180 format (/,' Total Surface Area :',5x,f12.4)
      write (iout,190)  vol
  190 format (/,' Total Excluded Vol :',5x,f12.4)
c
c     compute excluded volume derivatives via Kundrot method
c
      do i = 1, n
         dvol(1,i) = 0.0d0
         dvol(2,i) = 0.0d0
         dvol(3,i) = 0.0d0
      end do
      write (iout,200)
  200 format (/,' Craig Kundrot Excluded Volume Derivative Method :')
      call settime
      call kundrot1 (n,x,y,z,rsolv,probe,dvol)
      call gettime (wall,cpu)
      write (iout,210)  cpu,wall
  210 format (/,' CPU and Wall Times :',5x,2f12.4)
      write (iout,220)
  220 format (/,' Excluded Volume Derivatives :  (First Ten Atoms)',/)
      do i = 1, min(10,n)
         write (iout,230)  i,dvol(1,i),dvol(2,i),dvol(3,i)
  230    format (i8,6x,3f10.4)
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
         write (icrd,240)  n
  240    format (i8,/)
         do i = 1, n
            write (icrd,250)  i,x(i),y(i),z(i),rsolv(i)
  250       format (i8,3f14.6,f12.4)
         end do
         close (unit=icrd)
      end if
c
c     compute area, volume and derivatives via UnionBall
c
      write (iout,260)
  260 format (/,' Patrice Koehl UnionBall Alpha Shape Method :')
      call settime
      call unionball (n,x,y,z,rsolv,weight,probe,doderiv,dovol,
     &                   surf,vol,asurf,avol,dsurf,dvol)
      call gettime (wall,cpu)
      write (iout,270)  cpu,wall
  270 format (/,' CPU and Wall Times :',5x,2f12.4)
      write (iout,280)  surf
  280 format (/,' Total Surface Area :',5x,f12.4)
      write (iout,290)  vol
  290 format (/,' Total Excluded Vol :',5x,f12.4)
      write (iout,300)
  300 format (/,' Surface Area & Volume Derivatives :',
     &           '  (First Ten Atoms)',/)
      do i = 1, min(10,n)
         write (iout,310)  i,dsurf(1,i),dsurf(2,i),dsurf(3,i),
     &                     dvol(1,i),dvol(2,i),dvol(3,i)
  310    format (i8,6x,6f10.4)
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
