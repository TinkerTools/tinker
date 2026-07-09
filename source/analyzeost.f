c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program analyzeost  --  analyze an OST restart file     ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "analyzeost" reads an orthogonal space tempering restart file
c     and prints saved history, final free energy, or the g kernel
c
c
      program analyzeost
      use bath
      use files
      use iounit
      use keys
      use ost
      implicit none
      integer next
      integer trimtext
      logical exist
      character*12 mode
      character*240 ostfile
      character*240 string
c
c
c     get the name of the OST restart file
c
      call initial
      call nextarg (ostfile,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter OST Restart File Name :  ',$)
         read (input,20)  ostfile
   20    format (a240)
      end if
c
c     set the base filename used by the OST restart reader
c
      call ostbasename (ostfile)
      call getkey
c
c     use room temperature if no simulation temperature is available
c
      if (kelvin .le. 0.0d0)  kelvin = 298.0d0
c
c     get the requested analysis mode
c
      mode = 'FREEENERGY'
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,30)
   30    format (/,' Choose SERIES, FREEENERGY or GKERNEL',
     &              ' [FREEENERGY] :  ',$)
         read (input,20)  string
      end if
      next = 1
      call gettext (string,mode,next)
      call upcase (mode)
      if (mode .eq. '    ')  mode = 'FREEENERGY'
c
c     read the OST file and rebuild its kernels
c
      call rdost
      if (.not. allocated(osthist)) then
         write (iout,40)  filename(1:leng)//'.ost'
   40    format (/,' ANALYZEOST  --  Unable to Read OST File :  ',a)
         call fatal
      end if
      call ostgridkey
c
c     perform the requested analysis
c
      if (mode(1:1) .eq. 'S' .or. mode .eq. 'HISTORY') then
         call ostseries
      else if (mode(1:1) .eq. 'F' .or. mode .eq. 'DG') then
         call ostfreeenergy
      else if (mode(1:1) .eq. 'G' .or. mode .eq. 'GRID') then
         call ostgkernel
      else
         write (iout,50)  mode(1:trimtext(mode))
   50    format (/,' ANALYZEOST  --  Unknown Analysis Mode :  ',a)
         call fatal
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostgridkey  --  override ost analysis grid    ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostgridkey" checks the keyfile for requested analysis grid
c     settings and rebuilds the ost kernels on the requested grid
c
c
      subroutine ostgridkey
      use keys
      use ost
      implicit none
      integer i,next
      integer nlmda1
      real*8 wflmda1
      logical setnlmda,setwflmda
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get requested analysis grid settings from the keyfile
c
      setnlmda = .false.
      setwflmda = .false.
      nlmda1 = nlmda
      wflmda1 = wflmda
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            read (string,*,err=10,end=10)  nlmda1
            setnlmda = .true.
         else if (keyword(1:14) .eq. 'FLAMBDA-WIDTH ') then
            read (string,*,err=10,end=10)  wflmda1
            setwflmda = .true.
         end if
   10    continue
      end do
c
c     normalize the requested lambda grid to the standard convention
c
      if (setnlmda) then
         if (nlmda1 .lt. 2)  nlmda1 = nlmda
         if (mod(nlmda1,2) .eq. 0)  nlmda1 = nlmda1 + 1
      end if
      if (setwflmda) then
         if (wflmda1 .le. 0.0d0)  wflmda1 = wflmda
      end if
c
c     rebuild kernels if the requested grid differs from restart grid
c
      if ((setnlmda .and. nlmda1.ne.nlmda) .or.
     &    (setwflmda .and. wflmda1.ne.wflmda)) then
         call remeshost (nlmda1,wflmda1)
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine remeshost  --  rebuild ost on a new grid      ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "remeshost" changes the analysis grid while preserving the
c     saved gaussian centers, widths and heights from the restart
c
c
      subroutine remeshost (nlmda1,wflmda1)
      use ost
      implicit none
      integer i,j
      integer lowbin,highbin
      integer nlmda1
      integer nflmda0,fli00
      real*8 wflmda1
      real*8 wflmda0
      real*8 flmin,flmax
      real*8 flow,fhigh
      real*8 etotfkernel
c
c
c     save old flambda grid range before changing the grid spacing
c
      nflmda0 = nflmda
      fli00 = fli0
      wflmda0 = wflmda
      flmin = dble(1-fli00) * wflmda0
      flmax = dble(nflmda0-fli00) * wflmda0
c
c     set the requested lambda and flambda grid spacings
c
      nlmda = nlmda1
      wlmda = 1.0d0 / dble(nlmda-1)
      wlmda2 = 0.5d0 * wlmda
      wflmda = wflmda1
      wflmda2 = 0.5d0 * wflmda
c
c     preserve old flambda range and all saved gaussian cutoffs
c
      do i = 1, nosthist
         flow = ostfhist(i) - oststdev*ostwfhist(i) - wflmda2
         fhigh = ostfhist(i) + oststdev*ostwfhist(i) + wflmda2
         flmin = min(flmin,flow)
         flmax = max(flmax,fhigh)
      end do
c
c     derive flambda bin limits that cover the requested range
c
      lowbin = int(flmin/wflmda)
      if (dble(lowbin)*wflmda .gt. flmin)  lowbin = lowbin - 1
      highbin = int(flmax/wflmda)
      if (dble(highbin)*wflmda .lt. flmax)  highbin = highbin + 1
      fli0 = 1 - lowbin
      nflmda = highbin - lowbin + 1
      if (nflmda .lt. 1)  nflmda = 1
c
c     reallocate grid-dependent arrays for the new analysis grid
c
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(fkernel))  deallocate (fkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      allocate (osthead(nlmda,nflmda))
      allocate (fkernel(nlmda))
      allocate (gkernel(nlmda,nflmda))
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         do j = 1, nflmda
            gkernel(i,j) = 0.0d0
            osthead(i,j) = 0
         end do
      end do
c
c     rebuild lookup table and kernels from saved gaussian centers
c
      call buildostindex
      call buildgkernel
      call buildfkernel
      eosttot = etotfkernel()
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostbasename  --  set base name for OST file   ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostbasename" strips a trailing ".ost" suffix, if present, and
c     sets the global base filename used by the OST restart reader
c
c
      subroutine ostbasename (ostfile)
      use files
      implicit none
      integer trimtext
      character*240 ostfile
      character*240 suffix
c
c
c     remove a trailing ".ost" suffix, accepting either case
c
      filename = ostfile
      leng = trimtext (filename)
      if (leng .ge. 4) then
         suffix = filename(leng-3:leng)
         call upcase (suffix)
         if (suffix(1:4) .eq. '.OST') then
            filename(leng-3:leng) = '    '
            leng = leng - 4
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostseries  --  print OST time series data     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostseries" prints saved OST history entries and recomputes
c     the cumulative free energy estimate after each deposited bias
c
c
      subroutine ostseries
      use bath
      use iounit
      use ost
      implicit none
      integer i,j
      integer ihist
      integer nsave
      integer step
      real*8 etotfkernel
      real*8 freeeng
c
c
c     write a column header for the history table
c
      nsave = nosthist
      write (iout,10)
   10 format (/,' OST Time Series :',
     &        //,3x,'Hist',8x,'Step',9x,'Lambda',11x,'dU/dLambda',
     &           10x,'Free Energy',10x,'Height',
     &           8x,'Width-Lambda',7x,'Width-dU/dL',/)
      write (iout,15)  kelvin
   15 format (3x,'Temperature Used',6x,d20.10,' K',/)
c
c     rebuild the kernels cumulatively over the saved history
c
      freeeng = 0.0d0
      nosthist = 0
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = 0.0d0
         end do
      end do
      do ihist = 1, nsave
         nosthist = ihist
         call updategkernel
         call buildfkernel
         freeeng = etotfkernel()
         step = ihist * iosthist
         write (iout,20)  ihist,step,ostlhist(ihist),ostfhist(ihist),
     &      freeeng,osthhist(ihist),ostwlhist(ihist),ostwfhist(ihist)
   20    format (i7,i12,6d20.10)
      end do
c
c     restore the full saved history free energy
c
      eosttot = freeeng
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostfreeenergy  --  print final OST free energy##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostfreeenergy" prints the final free energy estimate from the
c     f kernel rebuilt from the full saved OST history
c
c
      subroutine ostfreeenergy
      use bath
      use iounit
      use ost
      implicit none
      real*8 etotfkernel
c
c
c     recompute and print the total free energy estimate
c
      call buildgkernel
      call buildfkernel
      eosttot = etotfkernel()
      write (iout,10)  nosthist,eosttot
   10 format (/,' OST Free Energy Estimate :',
     &        //,4x,'Number of Gaussians',8x,i12,
     &         /,4x,'Delta G',20x,d20.10)
      write (iout,20)  kelvin
   20 format (4x,'Temperature Used',11x,d20.10,' K')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine ostgkernel  --  print OST g kernel grid       ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "ostgkernel" prints the OST g kernel on its lambda/flambda grid
c     with header lines describing the grid origin and spacing
c
c
      subroutine ostgkernel
      use iounit
      use ost
      implicit none
      integer iflmda
      integer ilmda
      real*8 flstart
      real*8 flmda
      real*8 lambda
c
c
c     rebuild the full g kernel before printing
c
      call buildgkernel
      flstart = dble(1-fli0) * wflmda
      write (iout,10)  nlmda,nflmda
      write (iout,20)  0.0d0,wlmda
      write (iout,30)  flstart,wflmda
      write (iout,40)
   10 format ('# OST gkernel grid',/,
     &        '# nlmda ',i12,' nflmda ',i12)
   20 format ('# lambda_start ',d20.10,' lambda_width ',d20.10)
   30 format ('# flambda_start',d20.10,' flambda_width',d20.10)
   40 format ('# ilmda iflmda lambda flambda gkernel')
c
c     print all grid values as one row per lambda/flambda point
c
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         do iflmda = 1, nflmda
            flmda = dble(iflmda-fli0) * wflmda
            write (iout,50)  ilmda,iflmda,lambda,flmda,
     &                       gkernel(ilmda,iflmda)
   50       format (2i8,3d22.12)
         end do
      end do
      return
      end
