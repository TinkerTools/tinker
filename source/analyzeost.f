c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  program analyzeost  --  analyze an OST restart file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "analyzeost" reads an orthogonal space tempering or adaptive
c     biasing force history file and prints the saved history, the
c     final free energy, or the ost g kernel
c
c
      program analyzeost
      use bath
      use dlmda
      use files
      use iounit
      use keys
      use ost
      implicit none
      integer next
      integer ihis
      integer freeunit
      integer trimtext
      logical exist
      character*12 mode
      character*240 ostfile
      character*240 string
      character*240 title
c
c
c     get the name of the OST history file
c
      call initial
      call nextarg (ostfile,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter OST History File Name :  ',$)
         read (input,20)  ostfile
   20    format (a240)
      end if
c
c     set the base filename and read the keyfile, then read
c     exactly the file that was named on the command line
c
      call basefile (ostfile)
      lmdasavefile = ostfile
c
c     set default temperature
c
      kelvin = 298.0d0
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
c     read the history header once, keeping the ost header fields
c     until the title shows an abf history, then read the rows of the
c     method that wrote the file and rebuild its kernels
c
      inquire (file=lmdasavefile,exist=exist)
      if (.not. exist) then
         write (iout,40)  lmdasavefile(1:trimtext(lmdasavefile))
   40    format (/,' ANALYZEOST  --  History File Not Found :  ',a)
         call fatal
      end if
      ihis = freeunit ()
      open (unit=ihis,file=lmdasavefile,status='old')
      rewind (unit=ihis)
      use_ost = .true.
      call rdbiashead (ihis,lmdasavefile,title)
      if (index(title,abftitle(1:trimtext(abftitle))) .gt. 0) then
         use_ost = .false.
         use_abf = .true.
         call rdabfhist (ihis,lmdasavefile)
      else if (index(title,osttitle(1:trimtext(osttitle))) .gt. 0) then
         call rdosthist (ihis,lmdasavefile)
      else
         close (unit=ihis)
         write (iout,42)  lmdasavefile(1:trimtext(lmdasavefile))
   42    format (/,' ANALYZEOST  --  Unknown History File :  ',a)
         call fatal
      end if
      if (use_abf) then
         call abfkey
      else
         call ostkey
      end if
c
c     perform the requested analysis
c
      if (mode(1:1) .eq. 'S' .or. mode .eq. 'HISTORY') then
         if (use_abf) then
            call abfseries
         else
            call ostseries
         end if
      else if (mode(1:1) .eq. 'F' .or. mode .eq. 'DG') then
         if (use_abf) then
            call abffreeenergy
         else
            call ostfreeenergy
         end if
      else if (mode(1:1) .eq. 'G' .or. mode .eq. 'GRID') then
         if (use_abf) then
            write (iout,45)
   45       format (/,' ANALYZEOST  --  GKERNEL is not available for',
     &                 ' an ABF history')
            call fatal
         end if
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
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ostkey  --  apply the ost analysis keys  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ostkey" checks the keyfile for requested analysis grid settings,
c     since this program never runs the setup that reads these keywords
c     for dynamics, and rebuilds the ost kernels on the requested grid
c     when it differs from the grid of the restart
c
c
      subroutine ostkey
      use dlmda
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
c     get the requested settings from the keyfile, keeping the grid
c     the restart was read with as the default
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
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine remeshost  --  rebuild ost on a new grid  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "remeshost" changes the analysis grid while preserving the
c     saved gaussian centers, widths and heights from the restart
c
c
      subroutine remeshost (nlmda1,wflmda1)
      use dlmda
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
      real*8 efreetot
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
      do i = 1, nlmdahist
         flow = lmdafhist(i) - oststdev*ostwfhist(i) - wflmda2
         fhigh = lmdafhist(i) + oststdev*ostwfhist(i) + wflmda2
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
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      if (allocated(vkernelmax))  deallocate (vkernelmax)
      allocate (osthead(nlmda,nflmda))
      allocate (lmdafmean(nlmda))
      allocate (lmdafsum(nlmda))
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      allocate (lmdafwt(nlmda))
      allocate (vkernelmax(nlmda))
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            gfkernel(i,j) = 0.0d0
            gkernel(i,j) = 0.0d0
            glfkernel(i,j) = 0.0d0
            glkernel(i,j) = 0.0d0
            osthead(i,j) = 0
         end do
      end do
c
c     rebuild lookup table and kernels from saved gaussian centers
c
      call buildostindex
      if (fastkernel) then
         call buildkernels
      else
         call buildgkernel
         call buildfkernel
      end if
      lmdadeltag = efreetot()
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine ostseries  --  print OST time series data  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ostseries" prints saved OST history entries and recomputes
c     the cumulative free energy estimate after each deposited bias
c
c
      subroutine ostseries
      use bath
      use dlmda
      use iounit
      use ost
      implicit none
      integer i,j
      integer ihist
      integer nsave
      integer step
      real*8 efreetot
      real*8 freeeng
c
c
c     write a column header for the history table
c
      nsave = nlmdahist
      write (iout,10)  kelvin
   10 format (/,' OST Time Series :',
     &        //,3x,'Temperature Used',6x,1p,d20.10,' K')
      write (iout,20)
   20 format (/,3x,'Hist',8x,'Step',4x,'Lambda',14x,'dU/dLambda',
     &           10x,'Free Energy',9x,'Height',
     &           14x,'Width-Lambda',8x,'Width-dU/dL',/)
c
c     rebuild the kernels cumulatively over the saved history
c
      freeeng = 0.0d0
      nlmdahist = 0
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            gfkernel(i,j) = 0.0d0
            gkernel(i,j) = 0.0d0
            glfkernel(i,j) = 0.0d0
            glkernel(i,j) = 0.0d0
         end do
      end do
      do ihist = 1, nsave
         nlmdahist = ihist
         if (fastkernel) then
            call updatekernels
         else
            call updategkernel
            call buildfkernel
         end if
         freeeng = efreetot()
         step = lmdaihist(ihist)
         write (iout,30)  ihist,step,lmdalhist(ihist),lmdafhist(ihist),
     &      freeeng,osthhist(ihist),ostwlhist(ihist),ostwfhist(ihist)
   30    format (i7,i12,1p,6d20.10)
      end do
c
c     restore the full saved history free energy
c
      lmdadeltag = freeeng
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ostfreeenergy  --  print final OST free energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ostfreeenergy" prints the final free energy estimate from the
c     f kernel rebuilt from the full saved OST history, followed by
c     the f kernel mean force at each lambda bin
c
c
      subroutine ostfreeenergy
      use bath
      use dlmda
      use iounit
      use ost
      implicit none
      integer ilmda
      real*8 lambda
      real*8 efreetot
c
c
c     recompute and print the total free energy estimate
c
      if (fastkernel) then
         call buildkernels
      else
         call buildgkernel
         call buildfkernel
      end if
      lmdadeltag = efreetot()
      write (iout,10)  nlmdahist,lmdadeltag,kelvin
   10 format (/,' OST Free Energy Estimate :',
     &        //,1x,'Number of Gaussians',i16,
     &         /,1x,'Delta G',20x,1p,d20.10,
     &         /,1x,'Temperature Used',11x,d20.10,' K')
c
c     print the f kernel mean force at each lambda bin, whose
c     trapezoid integral is the free energy estimate
c
      write (iout,20)
   20 format (/,' OST Mean Force dU/dL(L) :',
     &        //,1x,'Lambda',10x,'dU/dLambda',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         write (iout,30)  lambda,lmdafmean(ilmda)
   30    format (1x,1p,d10.4,d16.4)
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ostgkernel  --  print OST g kernel grid  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ostgkernel" prints the OST g kernel on its lambda/flambda grid
c     with header lines describing the grid origin and spacing
c
c
      subroutine ostgkernel
      use dlmda
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
      if (fastkernel) then
         call buildkernels
      else
         call buildgkernel
      end if
      flstart = dble(1-fli0) * wflmda
c
c     print the grid dimensions, spacings and flambda origin
c
      write (iout,10)
   10 format (/,' OST g Kernel Setting :',
     &        //,1x,'nLambda',3x,'nFLambda',2x,'wLambda',
     &           15x,'wFLambda',14x,'sFLambda',/)
      write (iout,20)  nlmda,nflmda,wlmda,wflmda,flstart
   20 format (i8,i11,1p,d20.12,2d22.12)
c
c     print all grid values as one row per lambda/flambda point
c
      write (iout,30)
   30 format (/,' OST g Kernel Grid :',
     &        //,1x,'iLambda',7x,'iFLambda',6x,'Lambda',
     &           16x,'FLambda',15x,'gKernel',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         do iflmda = 1, nflmda
            flmda = dble(iflmda-fli0) * wflmda
            write (iout,40)  ilmda,iflmda,lambda,flmda,
     &                       gkernel(ilmda,iflmda)
   40       format (i8,i15,1p,d24.12,2d22.12)
         end do
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine abfkey  --  apply the abf analysis keys  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "abfkey" checks the keyfile for a requested lambda grid, since
c     this program never runs the setup that reads the keyword for
c     dynamics, and rebuilds the abf mean force on that grid when it
c     differs from the grid of the history file
c
c
      subroutine abfkey
      use dlmda
      use keys
      implicit none
      integer i,next
      integer nlmda1
      real*8 efreetot
      logical setnlmda
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get the requested lambda grid, keeping that of the history file
c     as the default
c
      setnlmda = .false.
      nlmda1 = nlmda
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            read (string,*,err=10,end=10)  nlmda1
            setnlmda = .true.
         end if
   10    continue
      end do
      if (setnlmda) then
         if (nlmda1 .lt. 2)  nlmda1 = nlmda
         if (mod(nlmda1,2) .eq. 0)  nlmda1 = nlmda1 + 1
      end if
c
c     rebuild the lambda bins on the requested grid
c
      if (setnlmda .and. nlmda1.ne.nlmda) then
         nlmda = nlmda1
         wlmda = 1.0d0 / dble(nlmda-1)
         wlmda2 = 0.5d0 * wlmda
         if (allocated(lmdafmean))  deallocate (lmdafmean)
         if (allocated(lmdafsum))  deallocate (lmdafsum)
         if (allocated(lmdafwt))  deallocate (lmdafwt)
         allocate (lmdafmean(nlmda))
         allocate (lmdafsum(nlmda))
         allocate (lmdafwt(nlmda))
         call buildabfkernel
         lmdadeltag = efreetot()
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine abfseries  --  print ABF time series data  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "abfseries" prints the saved ABF interval samples and recomputes
c     the cumulative free energy estimate after each sample
c
c
      subroutine abfseries
      use dlmda
      use iounit
      implicit none
      integer i
      integer ihist
      real*8 efreetot
      real*8 freeeng
c
c
c     write a column header for the sample table
c
      write (iout,10)
   10 format (/,' ABF Time Series :')
      write (iout,20)
   20 format (/,3x,'Hist',8x,'Step',4x,'Lambda',14x,'dU/dLambda',
     &           10x,'Free Energy',/)
c
c     rebuild the mean force cumulatively over the saved samples
c
      freeeng = 0.0d0
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
      end do
      do ihist = 1, nlmdahist
         call addabfhist (ihist)
         freeeng = efreetot()
         write (iout,30)  ihist,lmdaihist(ihist),lmdalhist(ihist),
     &                    lmdafhist(ihist),freeeng
   30    format (i7,i12,1p,3d20.10)
      end do
c
c     keep the free energy of the full saved history
c
      lmdadeltag = freeeng
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine abffreeenergy  --  print final ABF free energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "abffreeenergy" prints the final free energy estimate from the
c     mean force rebuilt from the full saved ABF history, followed by
c     the mean force at each lambda bin
c
c
      subroutine abffreeenergy
      use dlmda
      use iounit
      implicit none
      integer ilmda
      real*8 lambda
      real*8 efreetot
c
c
c     recompute and print the total free energy estimate
c
      call buildabfkernel
      lmdadeltag = efreetot()
      write (iout,10)  nlmdahist,lmdadeltag
   10 format (/,' ABF Free Energy Estimate :',
     &        //,1x,'Number of Samples',i18,
     &         /,1x,'Delta G',20x,1p,d20.10)
c
c     print the mean force at each lambda bin, whose trapezoid
c     integral is the free energy estimate
c
      write (iout,20)
   20 format (/,' ABF Mean Force dU/dL(L) :',
     &        //,1x,'Lambda',10x,'dU/dLambda',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         write (iout,30)  lambda,lmdafmean(ilmda)
   30    format (1x,1p,d10.4,d16.4)
      end do
      return
      end
