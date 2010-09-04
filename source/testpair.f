c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1993  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program testpair  --  time various neighbor pair schemes  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "testpair" performs a set of timing tests to compare the
c     evaluation of potential energy and energy/gradient using
c     different methods for finding pairwise neighbors
c
c
      program testpair
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inform.i'
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
      include 'potent.i'
      include 'vdwpot.i'
      integer i,j,k,m
      integer kgy,kgz
      integer start,stop
      integer ncalls,lmax
      integer npair,nterm
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
      real*8 wall,cpu,delta
      real*8 vrms,erms
      real*8 off,off2
      real*8 eloop,elight,elist
      real*8, allocatable :: xsort(:)
      real*8, allocatable :: ysort(:)
      real*8, allocatable :: zsort(:)
      real*8, allocatable :: dloop(:,:)
      real*8, allocatable :: dlight(:,:)
      real*8, allocatable :: dlist(:,:)
      logical exist,query
      logical header,match
      logical repeat
      character*1 axis(3)
      character*120 string
      data axis  / 'X','Y','Z' /
c
c
c     read the molecular system and setup molecular mechanics
c
      call initial
      call getxyz
      call mechanic
c
c     set difference threshhold via the energy precision
c
      delta = 1.0d-4
      if (digits .ge. 6)  delta = 1.0d-6
      if (digits .ge. 8)  delta = 1.0d-8
c
c     get the number of calculation cycles to be performed
c
      ncalls = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  ncalls
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter Desired Number of Repetitions [1] :  ',$)
         read (input,30,err=10)  ncalls
   30    format (i10)
      end if
      if (ncalls .eq. 0)  ncalls = 1
c
c     initialize number of pairs and generic cutoff distance
c
      npair = 0
      nterm = 0
      if (use_vdw)  nterm = nterm + 1
      if (use_charge)  nterm = nterm + 1
      if (use_mpole .or. use_polar)  nterm = nterm + 1
      nterm = nterm * ncalls
      off = 5.0d0
      off2 = off * off
c
c     perform dynamic allocation of some local arrays
c
      lmax = 8 * n
      allocate (xsort(lmax))
      allocate (ysort(lmax))
      allocate (zsort(lmax))
      allocate (dloop(3,n))
      allocate (dlight(3,n))
      allocate (dlist(3,n))
c
c     get the timing for setup of double nested loop
c
      call settime
      do m = 1, nterm
         do i = 1, n-1
            xi = x(i)
            yi = y(i)
            zi = z(i)
            do j = i+1, n
               xr = x(j) - xi
               yr = y(j) - yi
               zr = z(j) - zi
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .lt. off2)  npair = npair + 1
            end do
         end do
      end do
      call gettime (wall,cpu)
      write (iout,40)  ncalls
   40 format (/,' Total Wall Clock and CPU Time in Seconds for',
     &           i6,' Evaluations :')
      write (iout,50)
   50 format (/,' Computation Overhead',8x,'Wall',8x,'CPU')
      write (iout,60)  wall,cpu
   60 format (/,' Double Nested Loop',3x,2f11.3)
c
c     get the timing for setup of method of lights
c
      call settime
      do m = 1, nterm
         do i = 1, n
            xsort(i) = x(i)
            ysort(i) = y(i)
            zsort(i) = z(i)
         end do
         call lights (off,n,xsort,ysort,zsort)
         do i = 1, n
            xi = xsort(rgx(i))
            yi = ysort(rgy(i))
            zi = zsort(rgz(i))
            if (kbx(i) .le. kex(i)) then
               repeat = .false.
               start = kbx(i) + 1
               stop = kex(i)
            else
               repeat = .true.
               start = 1
               stop = kex(i)
            end if
   70       continue
            do j = start, stop
               k = locx(j)
               kgy = rgy(k)
               if (kby(i) .le. key(i)) then
                  if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 80
               else
                  if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 80
               end if
               kgz = rgz(k)
               if (kbz(i) .le. kez(i)) then
                  if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 80
               else
                  if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 80
               end if
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .lt. off2)  npair = npair + 1
   80          continue
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 70
            end if
         end do
      end do
      call gettime (wall,cpu)
      write (iout,90)  wall,cpu
   90 format (' Method of Lights',5x,2f11.3)
      if (npair .lt. 0)  call fatal
c
c     get the timing for setup of pair neighbor list
c
      call settime
      do m = 1, ncalls
         dovlst = .true.
         doclst = .true.
         domlst = .true.
         call nblist
      end do
      call gettime (wall,cpu)
      write (iout,100)  wall,cpu
  100 format (' Pair Neighbor List',3x,2f11.3)
c
c     get the timing for energy terms via double nested loop
c
      use_lights = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss
         end if
         if (use_charge)  call echarge
         if (use_mpole .or. use_polar)  call empole
      end do
      call gettime (wall,cpu)
      write (iout,110)
  110 format (/,' Potential Energy Only',7x,'Wall',8x,'CPU',
     &           13x,'Evdw',11x,'Eelect')
      eloop = ev + ec + em + ep
      if (digits .ge. 8) then
         write (iout,120)  wall,cpu,ev,ec+em+ep
  120    format (/,' Double Nested Loop',3x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,130)  wall,cpu,ev,ec+em+ep
  130    format (/,' Double Nested Loop',3x,2f11.3,2f17.6)
      else
         write (iout,140)  wall,cpu,ev,ec+em+ep
  140    format (/,' Double Nested Loop',3x,2f11.3,2f17.4)
      end if
c
c     get the timing for energy terms via method of lights
c
      use_lights = .true.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss
         end if
         if (use_charge)  call echarge
         if (use_mpole .or. use_polar)  call empole
      end do
      call gettime (wall,cpu)
      elight = ev + ec + em + ep
      if (digits .ge. 8) then
         write (iout,150)  wall,cpu,ev,ec+em+ep
  150    format (' Method of Lights',5x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,160)  wall,cpu,ev,ec+em+ep
  160    format (' Method of Lights',5x,2f11.3,2f17.6)
      else
         write (iout,170)  wall,cpu,ev,ec+em+ep
  170    format (' Method of Lights',5x,2f11.3,2f17.4)
      end if
c
c     get the timing for energy terms via pair neighbor list
c
      use_lights = .false.
      use_vlist = .true.
      use_clist = .true.
      use_mlist = .true.
      call nblist
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss
         end if
         if (use_charge)  call echarge
         if (use_mpole .or. use_polar)  call empole
      end do
      call gettime (wall,cpu)
      elist = ev + ec + em + ep
      if (digits .ge. 8) then
         write (iout,180)  wall,cpu,ev,ec+em+ep
  180    format (' Pair Neighbor List',3x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,190)  wall,cpu,ev,ec+em+ep
  190    format (' Pair Neighbor List',3x,2f11.3,2f17.6)
      else
         write (iout,200)  wall,cpu,ev,ec+em+ep
  200    format (' Pair Neighbor List',3x,2f11.3,2f17.4)
      end if
c
c     compare the nonbond energies from the various methods
c
      match = .true.
      if (abs(elight-eloop).gt.delta .or. abs(elist-eloop).gt.delta)
     &   match = .false.
      if (match) then
         write (iout,210)
  210    format (/,' Energies Computed via all Neighbor Methods',
     &              ' are Identical')
      end if
c
c     get the timing for gradient via double nested loop
c
      use_lights = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
         if (use_charge)  call echarge1
         if (use_mpole .or. use_polar)  call empole1
      end do
      call gettime (wall,cpu)
c
c     store the double loop gradient and get rms values
c
      vrms = 0.0d0
      erms = 0.0d0
      do i = 1, n
         do j = 1, 3
            dloop(j,i) = dev(j,i) + dec(j,i) + dem(j,i) + dep(j,i)
            vrms = vrms + dev(j,i)**2
            erms = erms + dec(j,i)**2 + dem(j,i)**2 + dep(j,i)**2
         end do
      end do
      vrms = sqrt(vrms/dble(n))
      erms = sqrt(erms/dble(n))
      write (iout,220)
  220 format (/,' Energy and Gradient',9x,'Wall',8x,'CPU',
     &           13x,'Dvdw',11x,'Delect')
      if (digits .ge. 8) then
         write (iout,230)  wall,cpu,vrms,erms
  230    format (/,' Double Nested Loop',3x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,240)  wall,cpu,vrms,erms
  240    format (/,' Double Nested Loop',3x,2f11.3,2f17.6)
      else
         write (iout,250)  wall,cpu,vrms,erms
  250    format (/,' Double Nested Loop',3x,2f11.3,2f17.4)
      end if
c
c     get the timing for gradient via method of lights
c
      use_lights = .true.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
         if (use_charge)  call echarge1
         if (use_mpole .or. use_polar)  call empole1
      end do
      call gettime (wall,cpu)
c
c     store the method of lights gradient and get rms values
c
      vrms = 0.0d0
      erms = 0.0d0
      do i = 1, n
         do j = 1, 3
            dlight(j,i) = dev(j,i) + dec(j,i) + dem(j,i) + dep(j,i)
            vrms = vrms + dev(j,i)**2
            erms = erms + dec(j,i)**2 + dem(j,i)**2 + dep(j,i)**2
         end do
      end do
      vrms = sqrt(vrms/dble(n))
      erms = sqrt(erms/dble(n))
      if (digits .ge. 8) then
         write (iout,260)  wall,cpu,vrms,erms
  260    format (' Method of Lights',5x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,270)  wall,cpu,vrms,erms
  270    format (' Method of Lights',5x,2f11.3,2f17.6)
      else
         write (iout,280)  wall,cpu,vrms,erms
  280    format (' Method of Lights',5x,2f11.3,2f17.4)
      end if
c
c     get the timing for gradient via pair neighbor list
c
      use_lights = .false.
      use_vlist = .true.
      use_clist = .true.
      use_mlist = .true.
      call settime
      do k = 1, ncalls
         if (use_vdw) then
            if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
            if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
            if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
            if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
            if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
         if (use_charge)  call echarge1
         if (use_mpole .or. use_polar)  call empole1
      end do
      call gettime (wall,cpu)
c
c     get the pair neighbor list gradient rms values
c
      vrms = 0.0d0
      erms = 0.0d0
      do i = 1, n
         do j = 1, 3
            dlist(j,i) = dev(j,i) + dec(j,i) + dem(j,i) + dep(j,i)
            vrms = vrms + dev(j,i)**2
            erms = erms + dec(j,i)**2 + dem(j,i)**2 + dep(j,i)**2
         end do
      end do
      vrms = sqrt(vrms/dble(n))
      erms = sqrt(erms/dble(n))
      if (digits .ge. 8) then
         write (iout,290)  wall,cpu,vrms,erms
  290    format (' Pair Neighbor List',3x,2f11.3,2f17.8)
      else if (digits .ge. 6) then
         write (iout,300)  wall,cpu,vrms,erms
  300    format (' Pair Neighbor List',3x,2f11.3,2f17.6)
      else
         write (iout,310)  wall,cpu,vrms,erms
  310    format (' Pair Neighbor List',3x,2f11.3,2f17.4)
      end if
c
c     compare the nonbond gradients from the various methods
c
      match = .true.
      header = .true.
      do i = 1, n
         do j = 1, 3
            if (abs(dlight(j,i)-dloop(j,i)).gt.delta .or.
     &          abs(dlist(j,i)-dloop(j,i)).gt.delta) then
               if (header) then
                  match = .false.
                  header = .false.
                  write (iout,320)
  320             format (/,' Comparison of Nonbond Gradients from',
     &                       ' the Neighbor Methods :',
     &                    //,11x,'Component',14x,'Loop',12x,'Lights',
     &                       14x,'List',/)
               end if
               if (digits .ge. 8) then
                  write (iout,330)  i,axis(j),dloop(j,i),dlight(j,i),
     &                              dlist(j,i)
  330             format (10x,i6,' (',a1,')',3f18.8)
               else if (digits .ge. 6) then
                  write (iout,340)  i,axis(j),dloop(j,i),dlight(j,i),
     &                              dlist(j,i)
  340             format (10x,i6,' (',a1,')',3f18.6)
               else
                  write (iout,350)  i,axis(j),dloop(j,i),dlight(j,i),
     &                              dlist(j,i)
  350             format (10x,i6,' (',a1,')',3f18.4)
               end if
            end if
         end do
      end do
      if (match) then
         write (iout,360)
  360    format (/,' Gradients Computed via all Neighbor Methods',
     &              ' are Identical')
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xsort)
      deallocate (ysort)
      deallocate (zsort)
      deallocate (dloop)
      deallocate (dlight)
      deallocate (dlist)
c
c     perform any final tasks before program exit
c
      call final
      end
