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
      include 'iounit.i'
      include 'light.i'
      include 'neigh.i'
      include 'potent.i'
      include 'vdwpot.i'
      integer i,j,k,m
      integer kgy,kgz
      integer start,stop
      integer ncalls
      integer npair,nterm
      real*8 xi,yi,zi
      real*8 xr,yr,zr,r2
      real*8 elapsed,eps
      real*8 vrms,erms
      real*8 off,off2
      real*8 xsort(maxlight)
      real*8 ysort(maxlight)
      real*8 zsort(maxlight)
      real*8 dloop(3,maxatm)
      real*8 dlight(3,maxatm)
      real*8 dlist(3,maxatm)
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
c     get the number calculation cycles to be performed
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
c     get the timing for setup of double nested loop
c
      call setime
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
      call getime (elapsed)
      write (iout,40)  ncalls
   40 format (/,' Computation Overhead :',11x,'Time',
     &           7x,i5,' Evaluations')
      write (iout,50)  elapsed
   50 format (/,' Double Nested Loop',9x,f10.3)
c
c     get the timing for setup of method of lights
c
      call setime
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
   60       continue
            do j = start, stop
               k = locx(j)
               kgy = rgy(k)
               if (kby(i) .le. key(i)) then
                  if (kgy.lt.kby(i) .or. kgy.gt.key(i))  goto 70
               else
                  if (kgy.lt.kby(i) .and. kgy.gt.key(i))  goto 70
               end if
               kgz = rgz(k)
               if (kbz(i) .le. kez(i)) then
                  if (kgz.lt.kbz(i) .or. kgz.gt.kez(i))  goto 70
               else
                  if (kgz.lt.kbz(i) .and. kgz.gt.kez(i))  goto 70
               end if
               xr = xi - xsort(j)
               yr = yi - ysort(kgy)
               zr = zi - zsort(kgz)
               call image (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .lt. off2)  npair = npair + 1
   70          continue
            end do
            if (repeat) then
               repeat = .false.
               start = kbx(i) + 1
               stop = nlight
               goto 60
            end if
         end do
      end do
      call getime (elapsed)
      write (iout,80)  elapsed
   80 format (' Method of Lights',11x,f10.3)
      if (npair .lt. 0)  call fatal
c
c     get the timing for setup of pair neighbor list
c
      call setime
      do m = 1, ncalls
         dovlst = .true.
         doclst = .true.
         domlst = .true.
         call nblist
      end do
      call getime (elapsed)
      write (iout,90)  elapsed
   90 format (' Pair Neighbor List',9x,f10.3)
c
c     get the timing for energy terms via double nested loop
c
      use_lights = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call setime
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
      call getime (elapsed)
      write (iout,100)
  100 format (/,' Potential Energy Only :',10x,'Time',14x,'Evdw',
     &           12x,'Eelect')
      write (iout,110)  elapsed,ev,ec+em+ep
  110 format (/,' Double Nested Loop',9x,f10.3,3x,f15.4,3x,f15.4)
c
c     get the timing for energy terms via method of lights
c
      use_lights = .true.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call setime
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
      call getime (elapsed)
      write (iout,120)  elapsed,ev,ec+em+ep
  120 format (' Method of Lights',11x,f10.3,3x,f15.4,3x,f15.4)
c
c     get the timing for energy terms via pair neighbor list
c
      use_lights = .false.
      use_vlist = .true.
      use_clist = .true.
      use_mlist = .true.
      call nblist
      call setime
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
      call getime (elapsed)
      write (iout,130)  elapsed,ev,ec+em+ep
  130 format (' Pair Neighbor List',9x,f10.3,3x,f15.4,3x,f15.4)
c
c     get the timing for gradient via double nested loop
c
      use_lights = .false.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call setime
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
      call getime (elapsed)
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
      write (iout,140)
  140 format (/,' Energy and Gradient :',12x,'Time',14x,'Dvdw',
     &           12x,'Delect')
      write (iout,150)  elapsed,vrms,erms
  150 format (/,' Double Nested Loop',9x,f10.3,3x,f15.4,3x,f15.4)
c
c     get the timing for gradient via method of lights
c
      use_lights = .true.
      use_vlist = .false.
      use_clist = .false.
      use_mlist = .false.
      call setime
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
      call getime (elapsed)
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
      write (iout,160)  elapsed,vrms,erms
  160 format (' Method of Lights',11x,f10.3,3x,f15.4,3x,f15.4)
c
c     get the timing for gradient via pair neighbor list
c
      use_lights = .false.
      use_vlist = .true.
      use_clist = .true.
      use_mlist = .true.
      call setime
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
      call getime (elapsed)
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
      write (iout,170)  elapsed,vrms,erms
  170 format (' Pair Neighbor List',9x,f10.3,3x,f15.4,3x,f15.4)
c
c     compare the nonbond gradients from the various methods
c
      eps = 1.0d-4
      match = .true.
      header = .true.
      do i = 1, n
         do j = 1, 3
            if (abs(dlight(j,i)-dloop(j,i)).gt.eps .or.
     &          abs(dlist(j,i)-dloop(j,i)).gt.eps) then
               if (header) then
                  match = .false.
                  header = .false.
                  write (iout,180)
  180             format (/,' Comparison of Nonbond Gradients from',
     &                       ' the Neighbor Methods:',
     &                    //,11x,'Component',16x,'Loop',10x,'Lights',
     &                       12x,'List',/)
               end if
               write (iout,190)  i,axis(j),dloop(j,i),dlight(j,i),
     &                           dlist(j,i)
  190          format (10x,i6,' (',a1,')',4x,3f16.4)
            end if
         end do
      end do
c
c     success if all the neighbor methods give the same gradient
c
      if (match) then
         write (iout,200)
  200    format (/,' Gradients Computed via all the Neighbor Methods',
     &              ' are Identical')
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
