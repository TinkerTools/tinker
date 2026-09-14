c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eabfbias -- apply adaptive biasing force  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eabfbias" evaluates the adaptive biasing force free energy at
c     the current lambda and subtracts it from the energy; the bias
c     depends on lambda alone, so the Cartesian gradient and virial
c     are unchanged; it has no side effects on the history or the
c     lambda particle, so it may be called from a Monte Carlo barostat
c     trial, and it saves the bias derivative for a later "eabfdyn"
c     call in the same step
c
c
      subroutine eabfbias
      use dlmda
      use energi
      implicit none
      real*8 eabflmda,dfdl
c
c
c     free energy and its lambda derivative from the mean force
c
      call efreelmda (eabflmda,dfdl)
      esum = esum - eabflmda
      lmdavbias = -eabflmda
      lmdadfdl = dfdl
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine eabfdyn -- adaptive biasing force sampling  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "eabfdyn" propagates the lambda particle under the adaptive
c     biasing force and records the lambda and dU/dlambda average of
c     each sample interval, reusing the ost propagate, equilibrate and
c     average phases
c
c
      subroutine eabfdyn
      use dlmda
      use mutant
      implicit none
      integer isamp,istep
      integer nskip
      real*8 efreetot
c
c
c     increment lmdastep step counter
c
      lmdastep = lmdastep + 1
c
c     remove the running mean force saved by eabfbias from the
c     unbiased lambda derivative summed by lmdachain
c
      lmdaddgdl = lmdadfdl
      deffdl = dedl - lmdaddgdl
c
c     save all values in the interval, but average only after the
c     propagation and equilibration phases
c
      istep = mod(lmdastep,lmdaintv)
      if (istep .eq. 0) then
         isamp = lmdaintv
      else
         isamp = istep
      end if
      lmdallist(isamp) = lambda
      lmdaflist(isamp) = dedl
c
c     record the averaging phase mean every lmdaintv steps
c
      if (istep .eq. 0) then
         nskip = lmdanpa + lmdanpb
         call avgstd (lmdallist,nskip+1,lmdanpc,lmdaavg,
     &                lmdastd)
         call avgstd (lmdaflist,nskip+1,lmdanpc,dedlavg,dedlstd)
c
c     every interval is kept, since rejecting the noisy intervals
c     would bias the conditional mean of dU/dlambda
c
         nlmdahist = nlmdahist + 1
         if (nlmdahist .gt. sizelmdahist)  call resizeabfhist
         lmdaihist(nlmdahist) = lmdastep
         lmdalhist(nlmdahist) = lmdaavg
         lmdafhist(nlmdahist) = dedlavg
         call addabfhist (nlmdahist)
         lmdadeltag = efreetot()
      end if
c
c     propagate the lambda particle
c
      if (isamp .le. lmdanpa)  call lmdalangevin
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine addabfhist -- add one abf interval sample  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "addabfhist" adds one saved interval average of dU/dlambda to
c     the running sum and count of its lambda bin, and sets the bin
c     bias force to the mean force of the bin
c
c
      subroutine addabfhist (ihist)
      use dlmda
      implicit none
      integer ihist
      integer ilmda
      integer lmdabin
c
c
c     accumulate the sample and update the bin mean force
c
      ilmda = lmdabin(lmdalhist(ihist))
      lmdafsum(ilmda) = lmdafsum(ilmda) + lmdafhist(ihist)
      lmdafwt(ilmda) = lmdafwt(ilmda) + 1.0d0
      lmdafmean(ilmda) = lmdafsum(ilmda) / lmdafwt(ilmda)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine buildabfkernel -- build the abf mean force  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "buildabfkernel" rebuilds the abf mean force at every lambda bin
c     from the full saved history of interval samples
c
c
      subroutine buildabfkernel
      use dlmda
      implicit none
      integer i
      integer ihist
c
c
c     zero out the mean force and its accumulators
c
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
      end do
c
c     loop over saved interval samples
c
      do ihist = 1, nlmdahist
         call addabfhist (ihist)
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine resizeabfhist -- resize abf sample history  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "resizeabfhist" doubles the storage for saved abf interval
c     samples while preserving the previously saved samples
c
c
      subroutine resizeabfhist
      use dlmda
      implicit none
      integer i
      integer oldsize
      integer newsize
      integer, allocatable :: lmdaihist0(:)
      real*8, allocatable :: lmdalhist0(:)
      real*8, allocatable :: lmdafhist0(:)
c
c
c     save the old samples and set the new size
c
      oldsize = sizelmdahist
      newsize = 2 * oldsize
      allocate (lmdaihist0(oldsize))
      allocate (lmdalhist0(oldsize))
      allocate (lmdafhist0(oldsize))
      do i = 1, oldsize
         lmdaihist0(i) = lmdaihist(i)
         lmdalhist0(i) = lmdalhist(i)
         lmdafhist0(i) = lmdafhist(i)
      end do
c
c     allocate the resized arrays and restore the old samples
c
      deallocate (lmdaihist)
      deallocate (lmdalhist)
      deallocate (lmdafhist)
      sizelmdahist = newsize
      allocate (lmdaihist(sizelmdahist))
      allocate (lmdalhist(sizelmdahist))
      allocate (lmdafhist(sizelmdahist))
      do i = 1, oldsize
         lmdaihist(i) = lmdaihist0(i)
         lmdalhist(i) = lmdalhist0(i)
         lmdafhist(i) = lmdafhist0(i)
      end do
      do i = oldsize+1, newsize
         lmdaihist(i) = 0
         lmdalhist(i) = 0.0d0
         lmdafhist(i) = 0.0d0
      end do
      deallocate (lmdaihist0)
      deallocate (lmdalhist0)
      deallocate (lmdafhist0)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine initabffile -- open the abf history file  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "initabffile" opens the external file that will hold the abf
c     bias history; an existing history is read back and continued
c     at an interval boundary, with the settings of the current run
c     taking precedence, while otherwise a new file is started
c
c
      subroutine initabffile
      use bath
      use dlmda
      use files
      use iounit
      implicit none
      integer ihis
      integer nlmda0
      integer lmdaintv0
      integer freeunit
      integer trimtext
      real*8 kelvin0
      real*8 lmdamass0
      real*8 lmdafric0
      real*8 lmdadt0
      logical exist
c
c
c     return if abf has not been initialized
c
      if (.not. use_abf)  return
c
c     continue from the latest existing abf history file
c
      lmdasavefile = filename(1:leng)//'.abf'
      call version (lmdasavefile,'old')
      inquire (file=lmdasavefile,exist=exist)
      if (exist) then
         nlmda0 = nlmda
         lmdaintv0 = lmdaintv
         kelvin0 = kelvin
         lmdamass0 = lmdamass
         lmdafric0 = lmdafric
         lmdadt0 = lmdadt
         call rdabf
         if (nlmda.ne.nlmda0 .or. lmdaintv.ne.lmdaintv0) then
            write (iout,10)  lmdasavefile(1:trimtext(lmdasavefile))
   10       format (/,' INITABFFILE  --  LAMBDA-NBIN and',
     &                 ' OSTHIST-INTERVAL must match the ABF',
     &                 ' history being continued',
     &              /,'                  File Name :  ',a)
            call fatal
         end if
c
c     the settings of the current run win over the stored ones
c
         kelvin = kelvin0
         lmdamass = lmdamass0
         lmdafric = lmdafric0
         lmdadt = lmdadt0
c
c     resume at an interval boundary, since the samples of a partial
c     interval are not saved in the history file
c
         lmdastep = (lmdastep/lmdaintv) * lmdaintv
c
c     rewrite the history from the samples just read, so that rows
c     appended after the last header update are dropped before new
c     samples are appended
c
         ihis = freeunit ()
         open (unit=ihis,file=lmdasavefile,status='replace',
     &         access='stream',form='formatted')
         call prtabf (ihis)
         close (unit=ihis)
         nlmdasave = nlmdahist
         write (iout,20)  lmdasavefile(1:trimtext(lmdasavefile))
   20    format (/,' ABF  --  Bias History Continued From  ',a)
         return
      end if
c
c     otherwise open a new file holding the current history
c
      ihis = freeunit ()
      call version (lmdasavefile,'new')
      open (unit=ihis,file=lmdasavefile,status='new',
     &      access='stream',form='formatted')
      call prtabf (ihis)
      close (unit=ihis)
      nlmdasave = nlmdahist
c
c     report the name of the file holding the bias history
c
      write (iout,30)  lmdasavefile(1:trimtext(lmdasavefile))
   30 format (/,' ABF  --  Bias History Written To  ',a)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine saveabf -- output abf history information  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "saveabf" appends to the external .abf file the interval samples
c     recorded since the previous call, so a running simulation can
c     be monitored; the fixed-size header is updated in place
c
c
      subroutine saveabf
      use dlmda
      implicit none
      integer ihis
      integer freeunit
c
c
c     return if abf has not been initialized
c
      if (.not. use_abf)  return
      if (.not. allocated(lmdaihist))  return
c
c     append the samples recorded since the previous call before
c     updating the header, so an interrupted write leaves the
c     previous sample count as the last complete checkpoint
c
      ihis = freeunit ()
      if (nlmdahist .gt. nlmdasave) then
         open (unit=ihis,file=lmdasavefile,status='old',
     &         access='stream',form='formatted',position='append')
         call prtabfhist (ihis,nlmdasave+1,nlmdahist)
         close (unit=ihis)
      end if
      open (unit=ihis,file=lmdasavefile,status='old',
     &      access='stream',form='unformatted',position='rewind')
      call updbiashead (ihis,abftitle,abflabel)
      close (unit=ihis)
      nlmdasave = nlmdahist
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine rdabf -- input abf history information  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "rdabf" reads a saved abf history from the external .abf file
c     into the abf sample arrays, then rebuilds the mean force of every
c     lambda bin; used both to continue an abf simulation and for
c     analysis
c
c
      subroutine rdabf
      use dlmda
      use iounit
      implicit none
      integer ihis
      integer freeunit
      integer trimtext
      logical exist
      character*240 abffile
      character*240 title
c
c
c     return unless an abf history file is present
c
      abffile = lmdasavefile
      inquire (file=abffile,exist=exist)
      if (.not. exist)  return
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      rewind (unit=ihis)
c
c     read the header, which must be that of an abf history
c
      call rdbiashead (ihis,abffile,title)
      if (index(title,abftitle(1:trimtext(abftitle))) .eq. 0) then
         close (unit=ihis)
         write (iout,10)  abffile(1:trimtext(abffile))
   10    format (/,' RDABF  --  File is not an ABF History File',
     &           /,'           File Name :  ',a)
         call fatal
      end if
c
c     read the interval samples and rebuild the mean force
c
      call rdabfhist (ihis,abffile)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine rdabfhist -- input abf interval samples  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "rdabfhist" reads the interval samples that follow the header of
c     an open .abf file, then closes the file and rebuilds the mean
c     force of every lambda bin
c
c
      subroutine rdabfhist (ihis,abffile)
      use dlmda
      use inform
      use iounit
      implicit none
      integer i,ihis
      integer ihist
      integer trimtext
      real*8 efreetot
      character*240 record
      character*(*) abffile
c
c
c     reallocate the abf arrays to match the history file
c
      if (allocated(lmdaihist))  deallocate (lmdaihist)
      if (allocated(lmdalhist))  deallocate (lmdalhist)
      if (allocated(lmdafhist))  deallocate (lmdafhist)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(lmdaflist))  deallocate (lmdaflist)
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      allocate (lmdaihist(sizelmdahist))
      allocate (lmdalhist(sizelmdahist))
      allocate (lmdafhist(sizelmdahist))
      allocate (lmdallist(lmdaintv))
      allocate (lmdaflist(lmdaintv))
      allocate (lmdafmean(nlmda))
      allocate (lmdafsum(nlmda))
      allocate (lmdafwt(nlmda))
      do i = 1, sizelmdahist
         lmdaihist(i) = 0
         lmdalhist(i) = 0.0d0
         lmdafhist(i) = 0.0d0
      end do
      do i = 1, lmdaintv
         lmdallist(i) = 0.0d0
         lmdaflist(i) = 0.0d0
      end do
c
c     read the saved interval samples
c
      do ihist = 1, nlmdahist
         read (ihis,20,err=90,end=90)  record
         read (record,*,err=90,end=90)  lmdaihist(ihist),
     &      lmdalhist(ihist),lmdafhist(ihist)
      end do
   20 format (a240)
      close (unit=ihis)
c
c     rebuild the mean force from the saved interval samples
c
      call buildabfkernel
      lmdadeltag = efreetot()
      nlmdasave = nlmdahist
      if (debug) then
         write (iout,30)  abffile(1:trimtext(abffile))
   30    format (/,' Reading ABF Bias from :  ',a)
      end if
      return
c
c     malformed history file
c
   90 continue
      close (unit=ihis)
      write (iout,40)  abffile(1:trimtext(abffile))
   40 format (/,' RDABF  --  Error while Reading ABF History File',
     &        /,'           File Name :  ',a)
      call fatal
      end
c
c
c     ##################################################
c     ##                                              ##
c     ##  subroutine prtabf -- output of abf history  ##
c     ##                                              ##
c     ##################################################
c
c
c     "prtabf" writes the abf history header followed by all of the
c     interval samples recorded so far
c
c
      subroutine prtabf (ihis)
      use dlmda
      implicit none
      integer ihis
c
c
c     write the header followed by all current history entries
c
      call prtbiashead (ihis,abftitle,abflabel)
      call prtabfhist (ihis,1,nlmdahist)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine prtabfhist  --  output abf sample range  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "prtabfhist" writes a requested range of saved abf interval
c     samples as the step, lambda and dU/dlambda average of each
c
c
      subroutine prtabfhist (ihis,ifirst,ilast)
      use dlmda
      implicit none
      integer ihis
      integer ifirst
      integer ilast
      integer ihist
c
c
c     write saved interval samples
c
      do ihist = ifirst, ilast
         write (ihis,10)  lmdaihist(ihist),lmdalhist(ihist),
     &                    lmdafhist(ihist)
      end do
   10 format (i12,2d26.16)
      return
      end
