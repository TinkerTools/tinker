c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses Chung, Pengyu Ren, Jay Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  program analyzelmda  --  lambda free energy analysis  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "analyzelmda" finds the free energy along lambda from the file
c     written by dynamics for the method chosen by the user; for an
c     orthogonal space tempering, metadynamics or adaptive biasing
c     force history, it writes the saved history and the final free
c     energy, along with the g kernel for orthogonal space tempering,
c     to separate files; for sequential free energy perturbation or
c     thermodynamic integration over lambda windows, it writes the
c     free energy of each window and the totals to a file
c
c
      program analyzelmda
      use dlmda
      use iounit
      use ost
      implicit none
      integer mode
      integer trimtext
      logical exist,query
      character*4 ext
      character*12 method
      character*24 prompt
      character*40 histtitle
      character*240 histfile
      character*240 string
c
c
c     find the lambda free energy method to analyze
c
      call initial
      mode = 0
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  mode
         if (mode.ge.1 .and. mode.le.5)  query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' The Tinker Lambda Free Energy Analysis Utility',
     &              ' Can :',
     &           //,4x,'(1) Analyze Orthogonal Space Tempering (OST)',
     &              ' History',
     &           /,4x,'(2) Analyze Metadynamics History',
     &           /,4x,'(3) Analyze Adaptive Biasing Force (ABF)',
     &              ' History',
     &           /,4x,'(4) Analyze Sequential Free Energy',
     &              ' Perturbation (FEP)',
     &           /,4x,'(5) Analyze Sequential Thermodynamic',
     &              ' Integration (TI)')
         do while (mode.lt.1 .or. mode.gt.5)
            mode = 0
            write (iout,30)
   30       format (/,' Enter the Number of the Desired Choice :  ',$)
            read (input,40,err=50,end=50)  mode
   40       format (i10)
   50       continue
         end do
      end if
c
c     set the file extension and the name of the chosen method
c
      histtitle = ' '
      if (mode .eq. 1) then
         ext = 'ost'
         method = 'OST'
         prompt = 'OST History'
         histtitle = osttitle
         use_ost = .true.
      else if (mode .eq. 2) then
         ext = 'meta'
         method = 'Metadynamics'
         prompt = 'Metadynamics History'
         histtitle = metatitle
         use_meta = .true.
      else if (mode .eq. 3) then
         ext = 'abf'
         method = 'ABF'
         prompt = 'ABF History'
         histtitle = abftitle
         use_abf = .true.
      else if (mode .eq. 4) then
         ext = 'fep'
         method = 'FEP'
         prompt = 'FEP Trial Lambda Energy'
      else
         ext = 'ti'
         method = 'TI'
         prompt = 'TI Block Average'
      end if
c
c     get the name of the input file, then set the base filename
c     and read the keyfile
c
      call nextarg (histfile,exist)
      if (exist) then
         call basefile (histfile)
         call suffix (histfile,ext(1:trimtext(ext)),'old')
         inquire (file=histfile,exist=exist)
      end if
      do while (.not. exist)
         write (iout,60)  prompt(1:trimtext(prompt))
   60    format (/,' Enter ',a,' File Name :  ',$)
         read (input,70)  histfile
   70    format (a240)
         call basefile (histfile)
         call suffix (histfile,ext(1:trimtext(ext)),'old')
         inquire (file=histfile,exist=exist)
      end do
c
c     analyze the file with the chosen method
c
      if (mode .le. 3) then
         call biascalc (histfile,histtitle,method)
      else if (mode .eq. 4) then
         call fepcalc (histfile)
      else
         call ticalc (histfile)
      end if
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine lmdaopen  --  open a new analysis output file  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "lmdaopen" opens a new file named from the base filename and
c     the given extension, using the next unused version so earlier
c     output is kept
c
c
      subroutine lmdaopen (ext,iunit,anafile)
      use files
      implicit none
      integer iunit
      integer freeunit
      character*(*) ext
      character*240 anafile
c
c
c     find the next unused version of the file and open it
c
      anafile = filename(1:leng)//'.'//ext
      call version (anafile,'new')
      iunit = freeunit ()
      open (unit=iunit,file=anafile,status='new')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine biascalc  --  free energy from a bias history  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "biascalc" reads an orthogonal space tempering, metadynamics or
c     adaptive biasing force history, checks that the chosen method
c     wrote it, and writes the time series and the final free energy,
c     along with the g kernel for orthogonal space tempering
c
c
      subroutine biascalc (histfile,histtitle,method)
      use bath
      use dlmda
      use iounit
      use ost
      implicit none
      integer ihis
      integer iser,ifree,igk
      integer freeunit
      integer trimtext
      logical found
      character*240 serfile
      character*240 freefile
      character*240 gkfile
      character*240 title
      character*(*) histfile
      character*(*) histtitle
      character*(*) method
c
c
c     set default temperature
c
      kelvin = 298.0d0
c
c     read the history header, check that the chosen method wrote
c     the file, then read its rows and rebuild its kernels
c
      ihis = freeunit ()
      open (unit=ihis,file=histfile,status='old')
      rewind (unit=ihis)
      call rdbiashead (ihis,histfile,title)
      found = (index(title,histtitle(1:trimtext(histtitle))) .gt. 0)
      if (.not. found) then
         close (unit=ihis)
         write (iout,10)  method(1:trimtext(method)),
     &                    histfile(1:trimtext(histfile))
   10    format (/,' BIASCALC  --  ',a,' History not Found in File :  ',
     &              a)
         call fatal
      end if
c
c     an ost history uses the fused kernel updates of dynamics from
c     the first build made as its rows are read
c
      if (use_ost) then
         fastkernel = .true.
         call rdosthist (ihis,histfile)
         call ostkey
      else if (use_meta) then
         call rdmetahist (ihis,histfile)
         call metakey
      else
         call rdabfhist (ihis,histfile)
         call abfkey
      end if
c
c     write the time series and the final free energy
c
      call lmdaopen ('series',iser,serfile)
      if (use_ost) then
         call ostseries (iser)
      else if (use_meta) then
         call metaseries (iser)
      else
         call abfseries (iser)
      end if
      close (unit=iser)
      call lmdaopen ('freeenergy',ifree,freefile)
      if (use_ost) then
         call ostfreeenergy (ifree)
      else if (use_meta) then
         call metafreeenergy (ifree)
      else
         call abffreeenergy (ifree)
      end if
      close (unit=ifree)
      write (iout,20)  serfile(1:trimtext(serfile))
   20 format (/,' Time Series Written To :  ',a)
      write (iout,30)  freefile(1:trimtext(freefile))
   30 format (' Free Energy Written To :  ',a)
c
c     write the g kernel of an orthogonal space tempering history
c
      if (use_ost) then
         call lmdaopen ('gkernel',igk,gkfile)
         call ostgkernel (igk)
         close (unit=igk)
         write (iout,40)  gkfile(1:trimtext(gkfile))
   40    format (' g Kernel Written To :  ',a)
      end if
      return
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
      integer lmdanbin
      real*8 wflmda1
      logical setnlmda,setwflmda
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get the requested settings from the keyfile, keeping the grid
c     the restart was read with when none is given
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
      if (setnlmda)  nlmda1 = lmdanbin(nlmda1)
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
c     ##  subroutine ostseries  --  write OST time series data  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "ostseries" writes saved OST history entries and recomputes
c     the cumulative free energy estimate after each deposited bias
c
c
      subroutine ostseries (iunit)
      use bath
      use dlmda
      use ost
      implicit none
      integer iunit
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
      write (iunit,10)  kelvin
   10 format (/,' OST Time Series :',
     &        //,3x,'Temperature Used',6x,1p,d20.10,' K')
      write (iunit,20)
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
         write (iunit,30)  ihist,step,lmdalhist(ihist),lmdafhist(ihist),
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
c     ##  subroutine ostfreeenergy  --  write final OST free energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ostfreeenergy" writes the final free energy estimate from the
c     f kernel of the full saved OST history, followed by the f kernel
c     mean force at each lambda bin
c
c     note the kernels are read as "ostseries" leaves them, built
c     from the full history
c
c
      subroutine ostfreeenergy (iunit)
      use bath
      use dlmda
      use ost
      implicit none
      integer iunit
      integer ilmda
      real*8 lambda
      real*8 efreetot
c
c
c     recompute and print the total free energy estimate
c
      lmdadeltag = efreetot()
      write (iunit,10)  nlmdahist,lmdadeltag,kelvin
   10 format (/,' OST Free Energy Estimate :',
     &        //,1x,'Number of Gaussians',i16,
     &         /,1x,'Delta G',20x,1p,d20.10,
     &         /,1x,'Temperature Used',11x,d20.10,' K')
c
c     print the f kernel mean force at each lambda bin, whose
c     trapezoid integral is the free energy estimate
c
      write (iunit,20)
   20 format (/,' OST Mean Force dU/dL(L) :',
     &        //,1x,'Lambda',10x,'dU/dLambda',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         write (iunit,30)  lambda,lmdafmean(ilmda)
   30    format (1x,1p,d10.4,d16.4)
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine ostgkernel  --  write OST g kernel grid  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "ostgkernel" writes the OST g kernel on its lambda/flambda grid
c     with header lines describing the grid origin and spacing
c
c     note the g kernel is read as "ostseries" leaves it, built from
c     the full history
c
c
      subroutine ostgkernel (iunit)
      use dlmda
      use ost
      implicit none
      integer iunit
      integer iflmda
      integer ilmda
      real*8 flstart
      real*8 flmda
      real*8 lambda
c
c
c     get the flambda origin of the grid
c
      flstart = dble(1-fli0) * wflmda
c
c     print the grid dimensions, spacings and flambda origin
c
      write (iunit,10)
   10 format (/,' OST g Kernel Setting :',
     &        //,1x,'nLambda',3x,'nFLambda',2x,'wLambda',
     &           15x,'wFLambda',14x,'sFLambda',/)
      write (iunit,20)  nlmda,nflmda,wlmda,wflmda,flstart
   20 format (i8,i11,1p,d20.12,2d22.12)
c
c     print all grid values as one row per lambda/flambda point
c
      write (iunit,30)
   30 format (/,' OST g Kernel Grid :',
     &        //,1x,'iLambda',7x,'iFLambda',6x,'Lambda',
     &           16x,'FLambda',15x,'gKernel',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         do iflmda = 1, nflmda
            flmda = dble(iflmda-fli0) * wflmda
            write (iunit,40)  ilmda,iflmda,lambda,flmda,
     &                       gkernel(ilmda,iflmda)
   40       format (i8,i15,1p,d24.12,2d22.12)
         end do
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine metakey  --  apply the metadynamics keys  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "metakey" checks the keyfile for a requested lambda grid, for
c     interpolation of the bias and for the theta map of the lambda
c     particle, since this program never runs the setup that reads
c     these keywords for dynamics, and rebuilds the bias grid when it
c     differs from the grid of the history file
c
c     note the theta map is not saved in the history, so the keyfile
c     has to name the same map and sharpness the dynamics used
c
c
      subroutine metakey
      use dlmda
      use keys
      use ost
      implicit none
      integer i,next
      integer ihist
      integer nlmda1
      integer lmdanbin
      real*8 metadeltag
      logical setnlmda
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get the requested lambda grid, keeping that of the history file
c     when none is given, and whether to interpolate the bias
c
      setnlmda = .false.
      nlmda1 = nlmda
      ostinterpol = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            read (string,*,err=10,end=10)  nlmda1
            setnlmda = .true.
         else if (keyword(1:16) .eq. 'OST-INTERPOLATE ') then
            ostinterpol = .true.
         end if
   10    continue
      end do
      if (setnlmda)  nlmda1 = lmdanbin(nlmda1)
c
c     get the theta map with the defaults of dynamics
c
      call lmdathkey
c
c     rebuild the bin lookup and the bias grid on the requested
c     lambda grid
c
      if (setnlmda .and. nlmda1.ne.nlmda) then
         nlmda = nlmda1
         wlmda = 1.0d0 / dble(nlmda-1)
         wlmda2 = 0.5d0 * wlmda
         if (allocated(metahead))  deallocate (metahead)
         if (allocated(vmetagrid))  deallocate (vmetagrid)
         if (allocated(dvmetagrid))  deallocate (dvmetagrid)
         allocate (metahead(nlmda))
         allocate (vmetagrid(nlmda))
         allocate (dvmetagrid(nlmda))
         do i = 1, nlmda
            vmetagrid(i) = 0.0d0
            dvmetagrid(i) = 0.0d0
         end do
         call buildmetaindex
         do ihist = 1, nmetahist
            call addmetagrid (ihist)
         end do
      end if
      lmdadeltag = metadeltag()
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine metaseries  --  write metadynamics time series  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "metaseries" writes the saved metadynamics gaussians and
c     recomputes the cumulative free energy estimate after each
c     deposited bias
c
c     note the free energy is read from the end points of the bias
c     grid, which holds every gaussian added so far
c
c
      subroutine metaseries (iunit)
      use bath
      use dlmda
      use ost
      implicit none
      integer iunit
      integer i
      integer ihist
      integer nsave
      real*8 metadeltag
      real*8 freeeng
c
c
c     write a column header for the history table
c
      nsave = nmetahist
      write (iunit,10)  kelvin
   10 format (/,' Metadynamics Time Series :',
     &        //,3x,'Temperature Used',6x,1p,d20.10,' K')
      write (iunit,20)
   20 format (/,3x,'Hist',8x,'Step',4x,'Lambda',14x,'Free Energy',
     &           9x,'Height',14x,'Width-Lambda',/)
c
c     rebuild the bias grid cumulatively over the saved history
c
      freeeng = 0.0d0
      do i = 1, nlmda
         vmetagrid(i) = 0.0d0
         dvmetagrid(i) = 0.0d0
      end do
      do ihist = 1, nsave
         nmetahist = ihist
         call addmetagrid (ihist)
         freeeng = metadeltag()
         write (iunit,30)  ihist,metaihist(ihist),metalhist(ihist),
     &                     freeeng,metahhist(ihist),metawhist(ihist)
   30    format (i7,i12,1p,4d20.10)
      end do
c
c     restore the full saved history free energy
c
      nmetahist = nsave
      lmdadeltag = freeeng
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine metafreeenergy  --  write final free energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "metafreeenergy" writes the final free energy estimate from the
c     bias deposited over the full saved metadynamics history, followed
c     by the free energy profile at each lambda bin,
c
c        G(L) = Vbias(0) - Vbias(L) + kT ln(J(L)/J(0))
c
c     where J is the factor dtheta/dlambda by which the theta map of
c     the lambda particle weights its visits to each lambda, smoothed
c     by the deposited gaussian as found by "lmdathjac"
c
c     note the factor is the same at both ends for the symmetric
c     theta maps, so the total free energy needs no correction
c
c
      subroutine metafreeenergy (iunit)
      use bath
      use dlmda
      use ost
      use units
      implicit none
      integer iunit
      integer ihist
      integer ilmda
      real*8 lambda
      real*8 vbias0,vbias
      real*8 dvdl,rt
      real*8 sig,hsum
      real*8 jac0,jac
      real*8 metadeltag
      real*8 lmdathjac
c
c
c     recompute and write the total free energy estimate
c
      lmdadeltag = metadeltag()
      write (iunit,10)  nmetahist,lmdadeltag,kelvin
   10 format (/,' Metadynamics Free Energy Estimate :',
     &        //,1x,'Number of Gaussians',i16,
     &         /,1x,'Delta G',20x,1p,d20.10,
     &         /,1x,'Temperature Used',11x,d20.10,' K')
c
c     smooth the theta map factor with the height weighted mean
c     width of the deposited gaussians
c
      sig = 0.0d0
      hsum = 0.0d0
      do ihist = 1, nmetahist
         sig = sig + metahhist(ihist)*metawhist(ihist)
         hsum = hsum + metahhist(ihist)
      end do
      if (hsum .gt. 0.0d0)  sig = sig / hsum
      rt = gasconst * kelvin
      jac0 = lmdathjac (0.0d0,sig,oststdev*sig)
c
c     write the free energy at each lambda bin as the negative of
c     the deposited bias with the theta map factor, shifted to zero
c     at lambda zero
c
      write (iunit,20)
   20 format (/,' Metadynamics Free Energy Profile G(L) :',
     &        //,1x,'Lambda',14x,'G(Lambda)',/)
      call emetabias (0.0d0,vbias0,dvdl)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         call emetabias (lambda,vbias,dvdl)
         jac = lmdathjac (lambda,sig,oststdev*sig)
         if (jac0.gt.0.0d0 .and. jac.gt.0.0d0) then
            vbias = vbias - rt*log(jac/jac0)
         end if
         write (iunit,30)  lambda,vbias0-vbias
   30    format (1x,1p,d10.4,d20.10)
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
      integer lmdanbin
      real*8 efreetot
      logical setnlmda
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     get the requested lambda grid, keeping that of the history file
c     when none is given
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
      if (setnlmda)  nlmda1 = lmdanbin(nlmda1)
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
c     ##  subroutine abfseries  --  write ABF time series data  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "abfseries" writes the saved ABF interval samples and recomputes
c     the cumulative free energy estimate after each sample
c
c
      subroutine abfseries (iunit)
      use dlmda
      implicit none
      integer iunit
      integer i
      integer ihist
      real*8 efreetot
      real*8 freeeng
c
c
c     write a column header for the sample table
c
      write (iunit,10)
   10 format (/,' ABF Time Series :')
      write (iunit,20)
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
         write (iunit,30)  ihist,lmdaihist(ihist),lmdalhist(ihist),
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
c     ##  subroutine abffreeenergy  --  write final ABF free energy  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "abffreeenergy" writes the final free energy estimate from the
c     mean force rebuilt from the full saved ABF history, followed by
c     the mean force at each lambda bin
c
c
      subroutine abffreeenergy (iunit)
      use dlmda
      implicit none
      integer iunit
      integer ilmda
      real*8 lambda
      real*8 efreetot
c
c
c     recompute and print the total free energy estimate
c
      call buildabfkernel
      lmdadeltag = efreetot()
      write (iunit,10)  nlmdahist,lmdadeltag
   10 format (/,' ABF Free Energy Estimate :',
     &        //,1x,'Number of Samples',i18,
     &         /,1x,'Delta G',20x,1p,d20.10)
c
c     print the mean force at each lambda bin, whose trapezoid
c     integral is the free energy estimate
c
      write (iunit,20)
   20 format (/,' ABF Mean Force dU/dL(L) :',
     &        //,1x,'Lambda',10x,'dU/dLambda',/)
      do ilmda = 1, nlmda
         lambda = dble(ilmda-1) * wlmda
         write (iunit,30)  lambda,lmdafmean(ilmda)
   30    format (1x,1p,d10.4,d16.4)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine fepcalc  --  free energy from Tinker FEP file  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "fepcalc" reads the trial lambda energies of a Tinker FEP file
c     written by dynamics over a schedule of lambda windows, finds the
c     free energy change between each pair of adjacent windows via the
c     FEP and BAR methods, sums them over the whole schedule, and
c     writes the values to a file
c
c     note the windows are taken from the lambda values of the data
c     records, so a window without samples of its own is still found
c     from the trial energies of its neighbors; the windows are put
c     in order of increasing lambda for either direction of the
c     simulation, so the values are for lambda 0 --> 1, and all of
c     the samples of each window are used; the native and trial
c     energies of a sample share the same periodic box, whose volume
c     is in a fifth column of newer files, and the volumes are used
c     for the PV enthalpy terms only if the file header also gives
c     the pressure of an isobaric simulation; the errors come from
c     bootstrap trials that resample the samples of every window at
c     once, so the errors of the totals include the correlation of
c     pairs that share a window
c
c
      subroutine fepcalc (fepfile)
      use iounit
      use units
      implicit none
      integer i,j,k,m,n
      integer ifep,nrow
      integer nsamp,nwin,nwhead
      integer npair,norphan
      integer nvol,nzero,izero
      integer nbst,nbstin,nsmax
      integer na,nb,iw,istat
      integer kt,maxiter,ifree
      integer ntota,ntotb
      integer npboot,nwork
      integer ntf,ntc,nlist
      integer barnbst
      integer freeunit
      integer trimtext
      integer nboot(2)
      integer, allocatable :: rstep(:)
      integer, allocatable :: sfirst(:)
      integer, allocatable :: slast(:)
      integer, allocatable :: swin(:)
      integer, allocatable :: wsel(:)
      integer, allocatable :: woff(:)
      integer, allocatable :: wsamp(:)
      integer, allocatable :: draw(:)
      integer, allocatable :: mapa(:)
      integer, allocatable :: mapb(:)
      integer, allocatable :: ia(:)
      integer, allocatable :: ib(:)
      integer, allocatable :: pstat(:)
      integer, allocatable :: pna(:)
      integer, allocatable :: pnb(:)
      integer, allocatable :: paoff(:)
      integer, allocatable :: pboff(:)
      integer, allocatable :: pskip(:)
      integer, allocatable :: pfail(:)
      integer, allocatable :: pconv(:)
      integer, allocatable :: ptry(:)
      integer, allocatable :: plist(:)
      real*8 temp,tol
      real*8 patm,pres,rt
      real*8 u0,u1,v0,v1
      real*8 lfirst,lsecond
      real*8 cumul,random
      real*8 barsd,hdir,tds
      real*8 cfore,cback
      real*8 hfore,hback
      real*8 ctrial,htrial
      real*8 sum(4),val(5)
      real*8 tsum(5),tsum2(5)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8, allocatable :: rsim(:)
      real*8, allocatable :: rtrl(:)
      real*8, allocatable :: rene(:)
      real*8, allocatable :: rvol(:)
      real*8, allocatable :: slmda(:)
      real*8, allocatable :: wlmda(:)
      real*8, allocatable :: pa0(:)
      real*8, allocatable :: pa1(:)
      real*8, allocatable :: pva(:)
      real*8, allocatable :: pb0(:)
      real*8, allocatable :: pb1(:)
      real*8, allocatable :: pvb(:)
      real*8, allocatable :: wa(:)
      real*8, allocatable :: hwa(:)
      real*8, allocatable :: ea(:)
      real*8, allocatable :: wb(:)
      real*8, allocatable :: hwb(:)
      real*8, allocatable :: eb(:)
      real*8, allocatable :: pxa(:)
      real*8, allocatable :: pxb(:)
      real*8, allocatable :: pcref(:)
      real*8, allocatable :: pfep(:,:)
      real*8, allocatable :: pbar(:,:)
      real*8, allocatable :: psdu(:,:)
      real*8, allocatable :: phs(:,:)
      real*8, allocatable :: psum(:,:)
      logical exist,header
      logical found,totals
      logical hasvol,usevol
      logical allf,allc,tconv
      character*240 record
      character*240 string
      character*240 freefile
      character*(*) fepfile
      character*80 line
      external random
c
c
c     get any requested number of bootstrap trials for each pair
c
      nbstin = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=70,end=70)  nbstin
   70 continue
      if (nbstin .lt. 0)  nbstin = 0
      if (nbstin .eq. 1)  nbstin = 2
c
c     read the header values and count the data records
c
      ifep = freeunit ()
      open (unit=ifep,file=fepfile,status='old')
      rewind (unit=ifep)
      temp = 0.0d0
      patm = 0.0d0
      nwhead = 0
      nrow = 0
      hasvol = .false.
      header = .true.
   80 continue
      read (ifep,90,err=100,end=100)  record
   90 format (a240)
      if (header) then
         k = index(record,'Number of Lambda Windows')
         if (k .ne. 0)  read (record(k+24:240),*,err=80,end=80)  nwhead
         k = index(record,'Simulation Temperature')
         if (k .ne. 0)  read (record(k+22:240),*,err=80,end=80)  temp
         k = index(record,'Simulation Pressure')
         if (k .ne. 0)  read (record(k+19:240),*,err=80,end=80)  patm
         if (index(record,'Sim Lambda') .ne. 0) then
            header = .false.
            hasvol = (index(record,'Volume') .ne. 0)
         end if
      else if (record .ne. ' ') then
         nrow = nrow + 1
      end if
      goto 80
  100 continue
      if (header) then
         write (iout,110)
  110    format (/,' FEPCALC  --  No Trial Lambda Energies Found',
     &              ' in the FEP File')
         call fatal
      end if
      if (temp .le. 0.0d0) then
         write (iout,120)
  120    format (/,' FEPCALC  --  Simulation Temperature not Found',
     &              ' in the FEP File')
         call fatal
      end if
      usevol = (hasvol .and. patm.gt.0.0d0)
      if (.not. usevol)  patm = 1.0d0
      pres = 0.0d0
      if (usevol)  pres = patm / prescon
c
c     perform dynamic allocation of some local arrays
c
      n = max(1,nrow)
      allocate (rstep(n))
      allocate (rsim(n))
      allocate (rtrl(n))
      allocate (rene(n))
      allocate (rvol(n))
      allocate (sfirst(n))
      allocate (slast(n))
      allocate (swin(n))
      allocate (wlmda(2*n))
c
c     read the data records, ignoring a final record that was cut
c     short; a whole record holds 77 columns as written by "savefep",
c     or 59 columns in older files without the volume column
c
      rewind (unit=ifep)
  130 continue
      read (ifep,140,err=180,end=180)  record
  140 format (a240)
      if (index(record,'Sim Lambda') .eq. 0)  goto 130
      m = 0
      n = 0
  150 continue
      read (ifep,160,err=180,end=180)  record
  160 format (a240)
      if (record .eq. ' ')  goto 150
      m = m + 1
      rvol(m) = 0.0d0
      if (hasvol) then
         if (m.eq.nrow .and. trimtext(record).lt.77)  goto 180
         read (record,*,err=170,end=170)  rstep(m),rsim(m),rtrl(m),
     &                                    rene(m),rvol(m)
      else
         if (m.eq.nrow .and. trimtext(record).lt.59)  goto 180
         read (record,*,err=170,end=170)  rstep(m),rsim(m),rtrl(m),
     &                                    rene(m)
      end if
      n = m
      goto 150
  170 continue
      if (m .lt. nrow) then
         write (iout,175)  record(1:trimtext(record))
  175    format (/,' FEPCALC  --  Unreadable Record in the FEP File :',
     &           //,' ',a)
         call fatal
      end if
  180 continue
      close (unit=ifep)
      nrow = n
c
c     a file with a pressure but only zero box volumes gets no PV
c     terms, while volumes missing from only some samples are an error
c
      tol = 1.0d-6
      if (usevol) then
         nvol = 0
         nzero = 0
         izero = 0
         do i = 1, nrow
            if (abs(rtrl(i)-rsim(i)) .lt. tol) then
               if (rvol(i) .ne. 0.0d0) then
                  nvol = nvol + 1
               else
                  nzero = nzero + 1
                  if (izero .eq. 0)  izero = i
               end if
            end if
         end do
         if (nvol .eq. 0) then
            usevol = .false.
            patm = 1.0d0
            pres = 0.0d0
            write (iout,182)
  182       format (/,' FEPCALC  --  File gives a Pressure but No Box',
     &                 ' Volumes;',
     &              /,' PV Terms are not Included')
         else if (nzero .gt. 0) then
            write (iout,184)  rstep(izero)
  184       format (/,' FEPCALC  --  Box Volume Missing at Step',i11,
     &              /,' in a File with Volumes for Other Samples')
            call fatal
         end if
      end if
c
c     group the records into samples, each with its native record
c     first, followed by those at the trial lambda values
c
      nsamp = 0
      norphan = 0
      do i = 1, nrow
         if (abs(rtrl(i)-rsim(i)) .lt. tol) then
            nsamp = nsamp + 1
            sfirst(nsamp) = i
            slast(nsamp) = i
         else if (nsamp .eq. 0) then
            norphan = norphan + 1
         else
            k = sfirst(nsamp)
            if (rstep(i).eq.rstep(k) .and. slast(nsamp).eq.i-1
     &             .and. abs(rsim(i)-rsim(k)).lt.tol) then
               slast(nsamp) = i
            else
               norphan = norphan + 1
            end if
         end if
      end do
c
c     at least two of the windows must have samples of their own
c
      found = .false.
      if (nsamp .gt. 0) then
         lfirst = rsim(sfirst(1))
         do i = 2, nsamp
            lsecond = rsim(sfirst(i))
            if (abs(lsecond-lfirst) .ge. tol) then
               found = .true.
               goto 190
            end if
         end do
      end if
  190 continue
      if (.not. found) then
         write (iout,200)
  200    format (/,' FEPCALC  --  Fewer than Two Lambda Windows',
     &              ' have Samples')
         return
      end if
c
c     find the windows from the simulation and trial lambda values,
c     in order of increasing lambda
c
      nwin = 0
      do i = 1, nrow
         call addlmdawin (rsim(i),tol,nwin,wlmda)
         call addlmdawin (rtrl(i),tol,nwin,wlmda)
      end do
      call sortlmdawin (nwin,wlmda)
c
c     assign each sample to its window, count those of each and
c     list the samples of each window in file order
c
      npair = nwin - 1
      allocate (slmda(max(1,nsamp)))
      allocate (wsel(nwin))
      allocate (woff(nwin+1))
      allocate (wsamp(max(1,nsamp)))
      do i = 1, nsamp
         slmda(i) = rsim(sfirst(i))
      end do
      call grouplmdawin (nsamp,slmda,tol,nwin,wlmda,swin,wsel,
     &                   woff,wsamp)
      deallocate (slmda)
      nsmax = 0
      do i = 1, nwin
         nsmax = max(nsmax,wsel(i))
      end do
c
c     perform dynamic allocation of some local arrays
c
      n = max(1,woff(nwin+1))
      allocate (draw(n))
      allocate (mapa(n))
      allocate (mapb(n))
      allocate (pa0(n))
      allocate (pa1(n))
      allocate (pva(n))
      allocate (pb0(n))
      allocate (pb1(n))
      allocate (pvb(n))
      allocate (wa(n))
      allocate (hwa(n))
      allocate (ea(n))
      allocate (wb(n))
      allocate (hwb(n))
      allocate (eb(n))
      n = max(1,nsmax)
      allocate (ia(n))
      allocate (ib(n))
      allocate (pstat(npair))
      allocate (pna(npair))
      allocate (pnb(npair))
      allocate (paoff(npair))
      allocate (pboff(npair))
      allocate (pskip(npair))
      allocate (pfail(npair))
      allocate (pconv(npair))
      allocate (ptry(npair))
      allocate (plist(npair))
      allocate (pfep(4,npair))
      allocate (pbar(4,npair))
      allocate (psdu(2,npair))
      allocate (phs(3,npair))
      allocate (pxa(npair))
      allocate (pxb(npair))
      allocate (pcref(npair))
      allocate (psum(8,npair))
c
c     gather the frames of each pair of adjacent windows, where a
c     sample of window "iw" gives a state 0 frame if it has the energy
c     at window "iw+1", and a sample of window "iw+1" gives a state 1
c     frame if it has the energy at window "iw"; then find the free
c     energy values from all of the frames of each pair
c
      rt = gasconst * temp
      npboot = 0
      nwork = 0
      ntota = 0
      ntotb = 0
      do iw = 1, npair
         paoff(iw) = ntota
         pboff(iw) = ntotb
         na = 0
         nb = 0
         pskip(iw) = 0
         do m = 1, wsel(iw)
            i = wsamp(woff(iw)+m)
            mapa(woff(iw)+m) = 0
            j = 0
            do k = sfirst(i)+1, slast(i)
               if (j.eq.0 .and. abs(rtrl(k)-wlmda(iw+1)).lt.tol)  j = k
            end do
            if (j .ne. 0) then
               na = na + 1
               pa0(ntota+na) = rene(sfirst(i))
               pa1(ntota+na) = rene(j)
               pva(ntota+na) = 0.0d0
               if (usevol)  pva(ntota+na) = rvol(sfirst(i))
               mapa(woff(iw)+m) = na
            else
               pskip(iw) = pskip(iw) + 1
            end if
         end do
         do m = 1, wsel(iw+1)
            i = wsamp(woff(iw+1)+m)
            mapb(woff(iw+1)+m) = 0
            j = 0
            do k = sfirst(i)+1, slast(i)
               if (j.eq.0 .and. abs(rtrl(k)-wlmda(iw)).lt.tol)  j = k
            end do
            if (j .ne. 0) then
               nb = nb + 1
               pb1(ntotb+nb) = rene(sfirst(i))
               pb0(ntotb+nb) = rene(j)
               pvb(ntotb+nb) = 0.0d0
               if (usevol)  pvb(ntotb+nb) = rvol(sfirst(i))
               mapb(woff(iw+1)+m) = nb
            else
               pskip(iw) = pskip(iw) + 1
            end if
         end do
         ntota = ntota + na
         ntotb = ntotb + nb
         pna(iw) = na
         pnb(iw) = nb
         if (na.eq.0 .or. nb.eq.0) then
            pstat(iw) = -1
         else
            k = paoff(iw) + 1
            j = pboff(iw) + 1
            call barpair (na,nb,pa0(k),pa1(k),pb0(j),pb1(j),pva(k),
     &                    pvb(j),temp,temp,patm,0,.false.,dgfep,dgbar,
     &                    nboot,sdu,uave,hbar,istat)
            pstat(iw) = istat
            pfep(1,iw) = dgfep(1)
            pfep(3,iw) = dgfep(3)
            pbar(1,iw) = dgbar(1)
            pbar(2,iw) = dgbar(2)
            psdu(1,iw) = sdu(1)
            psdu(2,iw) = sdu(2)
            phs(1,iw) = hbar(1)
            phs(3,iw) = hbar(3)
            pcref(iw) = 0.0d0
            if (istat .eq. 0)  pcref(iw) = dgbar(1)
            call bootprep (na,nb,pa0(k),pa1(k),pb0(j),pb1(j),pva(k),
     &                     pvb(j),rt,rt,pres,pcref(iw),pxa(iw),pxb(iw),
     &                     wa(k),hwa(k),wb(j),hwb(j),ea(k),eb(j))
            if (istat .eq. 0) then
               npboot = npboot + 1
               nwork = nwork + na + nb
            end if
         end if
      end do
c
c     set the number of bootstrap trials, shared by all of the pairs,
c     from the frames of the pairs used in each trial, with no trials
c     if no pair has a BAR free energy
c
      nbst = nbstin
      if (nbst .eq. 0)  nbst = barnbst (nwork)
      if (npboot .eq. 0)  nbst = 0
      maxiter = 100
c
c     apply the estimators of every pair to bootstrap trials that
c     resample the samples of each window once, so that the totals
c     over the pairs include the correlation between pairs sharing a
c     window; each BAR trial starts from the BAR value of its pair
c
      do iw = 1, npair
         do k = 1, 8
            psum(k,iw) = 0.0d0
         end do
         pconv(iw) = 0
         ptry(iw) = 0
      end do
      do k = 1, 5
         tsum(k) = 0.0d0
         tsum2(k) = 0.0d0
      end do
      ntf = 0
      ntc = 0
      do kt = 1, nbst
         do iw = 1, nwin
            do m = 1, wsel(iw)
               draw(woff(iw)+m) = int(dble(wsel(iw))*random()) + 1
            end do
         end do
         allf = .true.
         allc = .true.
         do k = 1, 5
            val(k) = 0.0d0
         end do
         do iw = 1, npair
            if (pstat(iw) .ne. 0)  goto 230
            na = 0
            do m = 1, wsel(iw)
               j = mapa(woff(iw)+draw(woff(iw)+m))
               if (j .ne. 0) then
                  na = na + 1
                  ia(na) = j
               end if
            end do
            nb = 0
            do m = 1, wsel(iw+1)
               j = mapb(woff(iw+1)+draw(woff(iw+1)+m))
               if (j .ne. 0) then
                  nb = nb + 1
                  ib(nb) = j
               end if
            end do
            if (na.eq.0 .or. nb.eq.0) then
               allf = .false.
               allc = .false.
               goto 230
            end if
            k = paoff(iw) + 1
            j = pboff(iw) + 1
            call bootfep (na,ia,pa0(k),pa1(k),pva(k),wa(k),hwa(k),
     &                    pxa(iw),rt,pres,cfore,hfore,u0,v0)
            call bootfep (nb,ib,pb1(j),pb0(j),pvb(j),wb(j),hwb(j),
     &                    pxb(iw),rt,pres,cback,hback,u1,v1)
            cback = -cback
            ptry(iw) = ptry(iw) + 1
            psum(1,iw) = psum(1,iw) + cfore
            psum(2,iw) = psum(2,iw) + cfore*cfore
            psum(3,iw) = psum(3,iw) + cback
            psum(4,iw) = psum(4,iw) + cback*cback
            val(1) = val(1) + cfore
            val(2) = val(2) + cback
            call bootiter (na,ia,nb,ib,ea(k),eb(j),pcref(iw),rt,rt,
     &                     pbar(1,iw),maxiter,ctrial,tconv)
            if (tconv) then
               call boothbar (na,ia,nb,ib,pa0(k),pa1(k),pb0(j),pb1(j),
     &                        pva(k),pvb(j),ea(k),eb(j),pcref(iw),rt,
     &                        rt,pres,ctrial,htrial)
               pconv(iw) = pconv(iw) + 1
               psum(5,iw) = psum(5,iw) + ctrial
               psum(6,iw) = psum(6,iw) + ctrial*ctrial
               psum(7,iw) = psum(7,iw) + htrial
               psum(8,iw) = psum(8,iw) + htrial*htrial
               val(3) = val(3) + ctrial
               val(4) = val(4) + htrial
            else
               allc = .false.
            end if
  230       continue
         end do
c
c     the direct enthalpy uses the same samples of the end windows
c
         if (wsel(1).gt.0 .and. wsel(nwin).gt.0) then
            u0 = 0.0d0
            do m = 1, wsel(1)
               i = wsamp(woff(1)+draw(woff(1)+m))
               u0 = u0 + rene(sfirst(i)) + pres*rvol(sfirst(i))
            end do
            u1 = 0.0d0
            do m = 1, wsel(nwin)
               i = wsamp(woff(nwin)+draw(woff(nwin)+m))
               u1 = u1 + rene(sfirst(i)) + pres*rvol(sfirst(i))
            end do
            val(5) = u1/dble(wsel(nwin)) - u0/dble(wsel(1))
            tsum(5) = tsum(5) + val(5)
            tsum2(5) = tsum2(5) + val(5)*val(5)
         end if
         if (allf) then
            ntf = ntf + 1
            do k = 1, 2
               tsum(k) = tsum(k) + val(k)
               tsum2(k) = tsum2(k) + val(k)*val(k)
            end do
         end if
         if (allc) then
            ntc = ntc + 1
            do k = 3, 4
               tsum(k) = tsum(k) + val(k)
               tsum2(k) = tsum2(k) + val(k)*val(k)
            end do
         end if
      end do
c
c     find the bootstrap errors of each pair
c
      do iw = 1, npair
         pfail(iw) = 0
         if (pstat(iw) .eq. 0) then
            pfail(iw) = nbst - pconv(iw)
            pfep(2,iw) = barsd (ptry(iw),psum(1,iw),psum(2,iw))
            pfep(4,iw) = barsd (ptry(iw),psum(3,iw),psum(4,iw))
            pbar(3,iw) = 0.0d0
            if (pconv(iw) .ge. 1) then
               pbar(3,iw) = psum(5,iw) / dble(pconv(iw))
            end if
            pbar(4,iw) = barsd (pconv(iw),psum(5,iw),psum(6,iw))
            phs(2,iw) = barsd (pconv(iw),psum(7,iw),psum(8,iw))
         end if
      end do
c
c     open the output file for the free energy values
c
      call lmdaopen ('freeenergy',ifree,freefile)
c
c     print information about the file and the samples used
c
      write (ifree,240)
  240 format (/,' Free Energy from Tinker FEP File :',/)
      write (ifree,250)  fepfile(1:trimtext(fepfile))
  250 format (' FEP File Name :',12x,a)
      write (ifree,260)  nwin
  260 format (' Lambda Windows Found :',7x,i8)
      if (nwhead.gt.0 .and. nwhead.ne.nwin) then
         write (ifree,270)  nwhead
  270    format (' Lambda Windows in Header :',3x,i8,
     &              '  (Mismatch)')
      end if
      write (ifree,290)  temp
  290 format (' Simulation Temperature :',4x,f9.2,' K')
      if (usevol) then
         write (ifree,300)  patm
  300    format (' Simulation Pressure :',7x,f9.2,' Atm')
      end if
      write (ifree,340)  nbst
  340 format (' Bootstrap Trials :',11x,i8)
      if (norphan .gt. 0) then
         write (ifree,350)  norphan
  350    format (' Records without Native Energy :',i6,
     &              '  (Ignored)')
      end if
c
c     print the free energy for each pair of adjacent windows
c
      write (ifree,360)
  360 format (/,' Free Energy by Window Pair :',
     &        //,'  Pair  Lambda A  Lambda B    BAR Iter     +/-',
     &           '    BAR Boot     +/-    Cumul dG',/)
      cumul = 0.0d0
      found = .true.
      do iw = 1, npair
         line = ' '
         write (line(1:26),370)  iw,wlmda(iw),wlmda(iw+1)
  370    format (i6,2f10.4)
         if (pstat(iw) .lt. 0) then
            k = iw
            if (wsel(iw) .ne. 0)  k = iw + 1
            if (wsel(k) .eq. 0) then
               write (line(27:80),380)  k
  380          format (3x,'Not Computed, Window',i4,' has no Samples')
            else
               write (line(27:80),390)
  390          format (3x,'Not Computed, No Neighbor Energies')
            end if
            found = .false.
         else if (pstat(iw) .eq. 1) then
            write (line(27:80),400)
  400       format (3x,'BAR Iteration not Converged')
            found = .false.
         else if (pstat(iw) .eq. 2) then
            write (line(27:80),410)
  410       format (3x,'BAR Failed, No Overlap')
            found = .false.
         else
            write (line(27:46),420)  pbar(1,iw),pbar(2,iw)
  420       format (f12.4,f8.4)
            if (pconv(iw) .ge. 2) then
               write (line(47:66),430)  pbar(3,iw),pbar(4,iw)
  430          format (f12.4,f8.4)
            else
               write (line(47:66),440)
  440          format (9x,'n/a')
            end if
            if (found) then
               cumul = cumul + pbar(1,iw)
               write (line(67:78),450)  cumul
  450          format (f12.4)
            end if
         end if
         write (ifree,460)  line(1:trimtext(line))
  460    format (a)
      end do
c
c     note any samples skipped or bootstrap trials that failed
c
      n = 0
      do iw = 1, npair
         if (pskip(iw).gt.0 .or. pfail(iw).gt.0)  n = n + 1
      end do
      if (n .gt. 0)  write (ifree,470)
  470 format ()
      do iw = 1, npair
         if (pskip(iw) .gt. 0) then
            write (ifree,480)  iw,pskip(iw)
  480       format (' Pair',i4,' :',i8,' Samples Skipped, Missing',
     &                 ' a Neighbor Energy')
         end if
         if (pfail(iw) .gt. 0) then
            write (ifree,490)  iw,pfail(iw),nbst
  490       format (' Pair',i4,' :',i8,' of',i8,' Bootstrap Trials',
     &                 ' did not Converge')
         end if
      end do
c
c     print sampling diagnostics for each pair of windows
c
      write (ifree,500)
  500 format (/,' Sampling Diagnostics by Window Pair :',
     &        //,'  Pair    nA    nB    FEP Fwd    FEP Bwd    Fwd-Bwd',
     &           '  SD dU Fwd  SD dU Bwd',/)
      do iw = 1, npair
         line = ' '
         write (line(1:18),510)  iw,pna(iw),pnb(iw)
  510    format (3i6)
         if (pstat(iw) .lt. 0) then
            write (line(19:80),520)
  520       format (3x,'Not Computed')
         else if (pstat(iw) .eq. 1) then
            write (line(19:80),530)
  530       format (3x,'BAR Iteration not Converged')
         else if (pstat(iw) .eq. 2) then
            write (line(19:80),540)
  540       format (3x,'BAR Failed, No Overlap')
         else
            write (line(19:73),550)  pfep(1,iw),pfep(3,iw),
     &                               pfep(1,iw)-pfep(3,iw),
     &                               psdu(1,iw),psdu(2,iw)
  550       format (5f11.4)
         end if
         write (ifree,560)  line(1:trimtext(line))
  560    format (a)
      end do
c
c     totals need every pair and the expected number of windows
c
      totals = .true.
      if (nwhead.gt.0 .and. nwhead.ne.nwin)  totals = .false.
      nlist = 0
      do iw = 1, npair
         if (pstat(iw) .ne. 0) then
            totals = .false.
            nlist = nlist + 1
            plist(nlist) = iw
         end if
      end do
      if (.not. totals) then
         if (nlist .gt. 0) then
            write (ifree,570)
  570       format (/,' Totals not Computed; Missing or Failed Pairs :')
            write (ifree,580)  (plist(k),k=1,nlist)
  580       format (5x,15i5)
         else
            write (ifree,590)
  590       format (/,' Totals not Computed; Lambda Windows Differ',
     &                 ' from the File Header')
         end if
      end if
c
c     sum the free energies over pairs, with errors from the joint
c     bootstrap trials over all of the windows
c
      if (totals) then
         do k = 1, 4
            sum(k) = 0.0d0
         end do
         do iw = 1, npair
            sum(1) = sum(1) + pfep(1,iw)
            sum(2) = sum(2) + pfep(3,iw)
            sum(3) = sum(3) + pbar(1,iw)
            sum(4) = sum(4) + phs(1,iw)
         end do
         write (ifree,600)
  600    format (/,' Free Energy Totals :',/)
         write (ifree,610)  sum(1),barsd(ntf,tsum(1),tsum2(1))
  610    format (' Free Energy via Forward FEP',9x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         write (ifree,620)  sum(2),barsd(ntf,tsum(2),tsum2(2))
  620    format (' Free Energy via Backward FEP',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (ntc .ge. 2) then
            write (ifree,630)  sum(3),barsd(ntc,tsum(3),tsum2(3))
  630       format (' Free Energy via BAR Iteration',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
            write (ifree,640)  tsum(3)/dble(ntc),
     &                        barsd(ntc,tsum(3),tsum2(3))
  640       format (' Free Energy via BAR Bootstrap',7x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (ifree,650)  sum(3)
  650       format (' Free Energy via BAR Iteration',7x,f12.4,
     &                 ' Kcal/mol, Error not Available')
            write (ifree,660)
  660       format (' Free Energy via BAR Bootstrap',13x,
     &                 'not Available')
         end if
c
c     find the enthalpy directly from the end window energies,
c     and sum the BAR enthalpy and entropy over the pairs
c
         u0 = 0.0d0
         do m = 1, wsel(1)
            i = wsamp(woff(1)+m)
            u0 = u0 + rene(sfirst(i)) + pres*rvol(sfirst(i))
         end do
         u1 = 0.0d0
         do m = 1, wsel(nwin)
            i = wsamp(woff(nwin)+m)
            u1 = u1 + rene(sfirst(i)) + pres*rvol(sfirst(i))
         end do
         hdir = u1/dble(wsel(nwin)) - u0/dble(wsel(1))
         write (ifree,670)
  670    format (/,' Enthalpy and Entropy Totals :',/)
         write (ifree,680)  hdir,barsd(nbst,tsum(5),tsum2(5))
  680    format (' Enthalpy via Direct Estimate',8x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (ntc .ge. 2) then
            write (ifree,690)  sum(4),barsd(ntc,tsum(4),tsum2(4))
  690       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         else
            write (ifree,700)  sum(4)
  700       format (' Enthalpy via BAR Estimate',11x,f12.4,
     &                 ' Kcal/mol, Error not Available')
         end if
         tds = 0.0d0
         do iw = 1, npair
            tds = tds + phs(3,iw)
         end do
         write (ifree,710)  tds
  710    format (' BAR Estimate -T*dS Value',12x,f12.4,' Kcal/mol')
      end if
c
c     close the output file and report its name
c
      close (unit=ifree)
      write (iout,720)  freefile(1:trimtext(freefile))
  720 format (/,' Free Energy Written To :  ',a)
c
c     perform deallocation of some local arrays
c
      deallocate (rstep)
      deallocate (rsim)
      deallocate (rtrl)
      deallocate (rene)
      deallocate (rvol)
      deallocate (sfirst)
      deallocate (slast)
      deallocate (swin)
      deallocate (wlmda)
      deallocate (wsel)
      deallocate (woff)
      deallocate (wsamp)
      deallocate (draw)
      deallocate (mapa)
      deallocate (mapb)
      deallocate (pa0)
      deallocate (pa1)
      deallocate (pva)
      deallocate (pb0)
      deallocate (pb1)
      deallocate (pvb)
      deallocate (wa)
      deallocate (hwa)
      deallocate (ea)
      deallocate (wb)
      deallocate (hwb)
      deallocate (eb)
      deallocate (ia)
      deallocate (ib)
      deallocate (pstat)
      deallocate (pna)
      deallocate (pnb)
      deallocate (paoff)
      deallocate (pboff)
      deallocate (pskip)
      deallocate (pfail)
      deallocate (pconv)
      deallocate (ptry)
      deallocate (plist)
      deallocate (pfep)
      deallocate (pbar)
      deallocate (psdu)
      deallocate (phs)
      deallocate (pxa)
      deallocate (pxb)
      deallocate (pcref)
      deallocate (psum)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ticalc  --  free energy from Tinker TI file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ticalc" reads the block averaged lambda derivatives of a Tinker
c     TI file written by dynamics over a schedule of lambda windows,
c     integrates the mean derivative of each window over lambda by
c     the trapezoid rule, and writes the free energy of each window
c     and the totals to a file
c
c     note the windows are found from the lambda values of the block
c     averages and put in order of increasing lambda for either
c     direction of the simulation, so the values are for lambda
c     0 --> 1; the error of each window mean is the standard error of
c     its block averages, which are taken to be independent, and the
c     error of the full sample total follows from the trapezoid
c     weights; the bootstrap trials resample the block averages of
c     every window at once
c
c
      subroutine ticalc (tifile)
      use iounit
      implicit none
      integer i,j,k,m,n
      integer iti,ifree
      integer nrow,nwin,nwhead
      integer nstepavg,nsingle
      integer nbst,nbstin,ibst
      integer barnbst
      integer freeunit
      integer trimtext
      integer, allocatable :: rwin(:)
      integer, allocatable :: wcnt(:)
      integer, allocatable :: woff(:)
      integer, allocatable :: wrow(:)
      real*8 tol,lmda,dedl
      real*8 cumul,random
      real*8 dg,dgvar,dgse
      real*8 bsum,bsum2
      real*8 barsd,stdev
      real*8 val,wsum,wsum2
      real*8, allocatable :: rlmda(:)
      real*8, allocatable :: rdedl(:)
      real*8, allocatable :: wlmda(:)
      real*8, allocatable :: wmean(:)
      real*8, allocatable :: wse(:)
      real*8, allocatable :: wgt(:)
      logical exist,header
      logical span,totals
      character*240 record
      character*240 string
      character*240 freefile
      character*(*) tifile
      external random
c
c
c     get any requested number of bootstrap trials
c
      nbstin = 0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=10,end=10)  nbstin
   10 continue
      if (nbstin .lt. 0)  nbstin = 0
      if (nbstin .eq. 1)  nbstin = 2
c
c     read the header values and count the block averages
c
      iti = freeunit ()
      open (unit=iti,file=tifile,status='old')
      rewind (unit=iti)
      nwhead = 0
      nstepavg = 0
      nrow = 0
      header = .true.
   20 continue
      read (iti,30,err=40,end=40)  record
   30 format (a240)
      if (header) then
         k = index(record,'Number of Lambda Windows')
         if (k .ne. 0)  read (record(k+24:240),*,err=20,end=20)  nwhead
         k = index(record,'Steps per Block Average')
         if (k .ne. 0)
     &      read (record(k+23:240),*,err=20,end=20)  nstepavg
         if (index(record,'dE/dL') .ne. 0)  header = .false.
      else if (record .ne. ' ') then
         nrow = nrow + 1
      end if
      goto 20
   40 continue
      if (header .or. nrow.eq.0) then
         close (unit=iti)
         write (iout,50)
   50    format (/,' TICALC  --  No Block Averages Found in the',
     &              ' TI File')
         call fatal
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (rlmda(nrow))
      allocate (rdedl(nrow))
      allocate (rwin(nrow))
      allocate (wrow(nrow))
      allocate (wlmda(nrow))
c
c     read the block averages, ignoring a final record that was cut
c     short; a whole record holds 63 columns as written by "saveti"
c
      rewind (unit=iti)
   60 continue
      read (iti,70,err=100,end=100)  record
   70 format (a240)
      if (index(record,'dE/dL') .eq. 0)  goto 60
      m = 0
      n = 0
   80 continue
      read (iti,90,err=100,end=100)  record
   90 format (a240)
      if (record .eq. ' ')  goto 80
      m = m + 1
      if (trimtext(record) .lt. 63) then
         if (m .eq. nrow)  goto 100
         goto 110
      end if
      read (record,*,err=110,end=110)  k,lmda,dedl
      n = n + 1
      rlmda(n) = lmda
      rdedl(n) = dedl
      if (m .lt. nrow)  goto 80
      goto 100
  110 continue
      close (unit=iti)
      write (iout,120)  record(1:trimtext(record))
  120 format (/,' TICALC  --  Unreadable Record in the TI File :',
     &        //,' ',a)
      call fatal
  100 continue
      close (unit=iti)
      nrow = n
c
c     find the windows from the lambda values of the block averages
c
      tol = 1.0d-6
      nwin = 0
      do i = 1, nrow
         call addlmdawin (rlmda(i),tol,nwin,wlmda)
      end do
      if (nwin .lt. 2) then
         write (iout,130)
  130    format (/,' TICALC  --  Fewer than Two Lambda Windows',
     &              ' have Block Averages')
         deallocate (rlmda)
         deallocate (rdedl)
         deallocate (rwin)
         deallocate (wrow)
         deallocate (wlmda)
         return
      end if
c
c     sort the windows into order of increasing lambda
c
      call sortlmdawin (nwin,wlmda)
c
c     perform dynamic allocation of some local arrays
c
      allocate (wcnt(nwin))
      allocate (woff(nwin+1))
      allocate (wmean(nwin))
      allocate (wse(nwin))
      allocate (wgt(nwin))
c
c     assign each block average to its window, count those of each
c     and list the block averages of each window in file order
c
      call grouplmdawin (nrow,rlmda,tol,nwin,wlmda,rwin,wcnt,
     &                   woff,wrow)
c
c     find the mean derivative of each window and its standard error
c
      nsingle = 0
      do j = 1, nwin
         n = wcnt(j)
         wsum = 0.0d0
         wsum2 = 0.0d0
         do k = 1, n
            val = rdedl(wrow(woff(j)+k))
            wsum = wsum + val
            wsum2 = wsum2 + val*val
         end do
         wmean(j) = wsum / dble(n)
         stdev = barsd (n,wsum,wsum2)
         wse(j) = stdev / sqrt(dble(n))
         if (n .eq. 1)  nsingle = nsingle + 1
      end do
c
c     set the trapezoid weight of each window over lambda
c
      do j = 1, nwin
         wgt(j) = 0.0d0
         if (j .gt. 1)  wgt(j) = wgt(j) + 0.5d0*(wlmda(j)-wlmda(j-1))
         if (j .lt. nwin)  wgt(j) = wgt(j) + 0.5d0*(wlmda(j+1)-wlmda(j))
      end do
c
c     totals need windows from lambda 0 to 1 and the expected number
c
      span = (abs(wlmda(1)).lt.tol .and. abs(wlmda(nwin)-1.0d0).lt.tol)
      totals = span
      if (nwhead.gt.0 .and. nwhead.ne.nwin)  totals = .false.
c
c     find the full sample free energy and its error
c
      dg = 0.0d0
      dgvar = 0.0d0
      do j = 1, nwin
         dg = dg + wgt(j)*wmean(j)
         dgvar = dgvar + (wgt(j)*wse(j))**2
      end do
      dgse = sqrt(dgvar)
c
c     apply the trapezoid rule to bootstrap trials that resample the
c     block averages of each window
c
      nbst = 0
      bsum = 0.0d0
      bsum2 = 0.0d0
      if (totals) then
         nbst = nbstin
         if (nbst .eq. 0)  nbst = barnbst (nrow)
         do ibst = 1, nbst
            cumul = 0.0d0
            do j = 1, nwin
               n = wcnt(j)
               wsum = 0.0d0
               do k = 1, n
                  m = int(dble(n)*random()) + 1
                  if (m .gt. n)  m = n
                  wsum = wsum + rdedl(wrow(woff(j)+m))
               end do
               cumul = cumul + wgt(j)*wsum/dble(n)
            end do
            bsum = bsum + cumul
            bsum2 = bsum2 + cumul*cumul
         end do
      end if
c
c     open the output file for the free energy values
c
      call lmdaopen ('freeenergy',ifree,freefile)
c
c     write information about the file and the block averages used
c
      write (ifree,150)
  150 format (/,' Free Energy from Tinker TI File :',/)
      write (ifree,160)  tifile(1:trimtext(tifile))
  160 format (' TI File Name :',13x,a)
      write (ifree,170)  nwin
  170 format (' Lambda Windows Found :',7x,i8)
      if (nwhead.gt.0 .and. nwhead.ne.nwin) then
         write (ifree,180)  nwhead
  180    format (' Lambda Windows in Header :',3x,i8,
     &              '  (Mismatch)')
      end if
      if (nstepavg .gt. 0) then
         write (ifree,190)  nstepavg
  190    format (' Steps per Block Average :',4x,i8)
      end if
      write (ifree,200)  nrow
  200 format (' Block Averages Used :',8x,i8)
      if (nbst .gt. 0) then
         write (ifree,210)  nbst
  210    format (' Bootstrap Trials :',11x,i8)
      end if
c
c     write the mean derivative and free energy for each window
c
      write (ifree,220)
  220 format (/,' Free Energy by Lambda Window :',
     &        //,2x,'Window',4x,'Lambda',4x,'Blocks',4x,'dU/dL Mean',
     &           7x,'+/-',4x,'Cumul dG',/)
      cumul = 0.0d0
      do j = 1, nwin
         if (j .gt. 1)  cumul = cumul + 0.5d0*(wlmda(j)-wlmda(j-1))
     &                                  *(wmean(j)+wmean(j-1))
         if (wcnt(j) .ge. 2) then
            write (ifree,230)  j,wlmda(j),wcnt(j),wmean(j),wse(j),cumul
  230       format (i8,f10.4,i10,f14.4,f10.4,f12.4)
         else
            write (ifree,240)  j,wlmda(j),wcnt(j),wmean(j),cumul
  240       format (i8,f10.4,i10,f14.4,7x,'n/a',f12.4)
         end if
      end do
      if (nsingle .gt. 0) then
         write (ifree,250)  nsingle
  250    format (/,' Windows with a Single Block Average and no',
     &              ' Error :',i5)
      end if
c
c     write the totals, or the reason they were not computed
c
      if (.not. totals) then
         if (.not. span) then
            write (ifree,260)
  260       format (/,' Totals not Computed; Lambda Windows do not',
     &                 ' Span 0 to 1')
         else
            write (ifree,270)
  270       format (/,' Totals not Computed; Lambda Windows Differ',
     &                 ' from the File Header')
         end if
      else
         write (ifree,280)
  280    format (/,' Free Energy Totals :',/)
         write (ifree,290)  dg,dgse
  290    format (' Free Energy via TI Full Sample',6x,f12.4,
     &              ' +/-',f9.4,' Kcal/mol')
         if (nbst .ge. 2) then
            write (ifree,300)  bsum/dble(nbst),barsd(nbst,bsum,bsum2)
  300       format (' Free Energy via TI Bootstrap',8x,f12.4,
     &                 ' +/-',f9.4,' Kcal/mol')
         end if
      end if
c
c     close the output file and report its name
c
      close (unit=ifree)
      write (iout,310)  freefile(1:trimtext(freefile))
  310 format (/,' Free Energy Written To :  ',a)
c
c     perform deallocation of some local arrays
c
      deallocate (rlmda)
      deallocate (rdedl)
      deallocate (rwin)
      deallocate (wrow)
      deallocate (wlmda)
      deallocate (wcnt)
      deallocate (woff)
      deallocate (wmean)
      deallocate (wse)
      deallocate (wgt)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine addlmdawin  --  add a new lambda window  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "addlmdawin" appends the lambda value "lmda" to the list of
c     windows unless it lies within "tol" of a window already listed
c
c
      subroutine addlmdawin (lmda,tol,nwin,wl)
      implicit none
      integer k,nwin
      real*8 lmda,tol
      real*8 wl(*)
c
c
      do k = 1, nwin
         if (abs(wl(k)-lmda) .lt. tol)  return
      end do
      nwin = nwin + 1
      wl(nwin) = lmda
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine sortlmdawin  --  sort the lambda windows  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "sortlmdawin" puts the lambda windows into order of increasing
c     lambda by an insertion sort
c
c
      subroutine sortlmdawin (nwin,wl)
      implicit none
      integer i,j,nwin
      real*8 lmda
      real*8 wl(*)
c
c
      do i = 2, nwin
         lmda = wl(i)
         j = i - 1
         do while (j .ge. 1)
            if (wl(j) .le. lmda)  goto 10
            wl(j+1) = wl(j)
            j = j - 1
         end do
   10    continue
         wl(j+1) = lmda
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine grouplmdawin  --  group values by window  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "grouplmdawin" assigns each of "n" lambda values to the window
c     within "tol" of it, counts the values of each window, and lists
c     the values of window k in file order as list(off(k)+1) through
c     list(off(k+1))
c
c
      subroutine grouplmdawin (n,vals,tol,nwin,wl,win,cnt,off,list)
      implicit none
      integer i,k,n,nwin
      integer win(*)
      integer cnt(*)
      integer off(*)
      integer list(*)
      real*8 tol
      real*8 vals(*)
      real*8 wl(*)
c
c
c     assign each value to its window and count those of each
c
      do k = 1, nwin
         cnt(k) = 0
      end do
      do i = 1, n
         do k = 1, nwin
            if (abs(wl(k)-vals(i)) .lt. tol)  win(i) = k
         end do
         cnt(win(i)) = cnt(win(i)) + 1
      end do
c
c     list the values of each window in file order
c
      off(1) = 0
      do k = 1, nwin
         off(k+1) = off(k) + cnt(k)
         cnt(k) = 0
      end do
      do i = 1, n
         k = win(i)
         cnt(k) = cnt(k) + 1
         list(off(k)+cnt(k)) = i
      end do
      return
      end
