c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eost  --  OST kernel and history tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eost" checks OST helper routines against
c     deterministic inputs and expected outputs
c
c
      subroutine test_eost
      use bath
      use units
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_eost')
c
c
      if (skiptest(tname,'ost'))  return
      call initial
      kelvin = 300.0d0
      call test_eost_index
      call test_eost_resize
      call test_eost_buildindex
      call test_eost_ensure
      call test_eost_gkernels
      call test_eost_fkernel
      call test_eost_kernelbuilds
      call test_eost_eginterpolate
      call test_eost_avgstd
      call test_eost_histstat
      call test_eost_drift
      call test_eost_depcriteria
      call test_eost_efkernel
      call test_eost_vkernelmax
      call test_eost_tempering
      call test_eost_ostdyn
      call test_eost_ostlocal
      call test_eost_ostphase
      call test_eost_ostgate
      call test_eost_save
      call test_eost_meta
      call test_eost_metaimage
      call test_eost_metadyn
      call test_eost_metatemper
      call test_eost_temperkeys
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_save  --  restart append test  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_save" checks that repeated saves append only new
c     histories while updating the fixed-size restart header in place
c
c
      subroutine test_eost_save
      use bath
      use dlmda
      use files
      use ost
      implicit none
      integer i,ihis
      integer ios
      integer leng0
      integer nline
      integer nsave0
      integer size1,size2,size3
      integer freeunit
      real*8 brutevkmax
      real*8 kelvin0
      logical exist
      character*240 filename0
      character*240 ostfile
      character*240 record
c
c     create a small deterministic history and a temporary restart
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 321.0d0
      filename = 'tinkertest-saveost'
      leng = len_trim(filename)
      ostfile = filename(1:leng)//'.ost'
      inquire (file=ostfile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=ostfile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetost (3,5,3)
      use_ost = .true.
      nlmdasave = 0
      lmdastep = 10
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      lmdamass = 1.0d0
      lmdafric = 1.0d0
      lmdadt = 0.001d0
      call sethist (1,0.1d0,1.0d0,2.0d0,0.01d0,1.0d0)
      call sethist (2,0.2d0,2.0d0,3.0d0,0.01d0,1.0d0)
      call sethist (3,0.3d0,3.0d0,4.0d0,0.01d0,1.0d0)
c
c     open the history file the way "dynamic" does before sampling
c
      call initostfile
c
c     append two histories, then save again without adding history
c
      nlmdahist = 1
      call saveost
      inquire (file=ostfile,size=size1)
      nlmdahist = 3
      lmdastep = 30
      call saveost
      inquire (file=ostfile,size=size2)
      lmdastep = 31
      call saveost
      inquire (file=ostfile,size=size3)
      call assert_logical (size2.gt.size1,.true.,
     &                     'saveost appends new history')
      call assert_int (size3,size2,'saveost does not duplicate history')
c
c     the per-gaussian heights survive a round trip, and the derived
c     running maxima are rebuilt from the restored history
c
      kelvin = 111.0d0
      call rdost
      call assert_real (kelvin,321.0d0,1.0d-12,
     &                  'rdost restores the temperature')
      call assert_int (nlmdahist,3,'rdost restores the history count')
      call assert_real (osthhist(1),2.0d0,1.0d-12,
     &                  'rdost restores height 1')
      call assert_real (osthhist(2),3.0d0,1.0d-12,
     &                  'rdost restores height 2')
      call assert_real (osthhist(3),4.0d0,1.0d-12,
     &                  'rdost restores height 3')
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'rdost rebuilds the running maximum')
      end do
c
c     check the refreshed header and the sequential history records
c
      ihis = freeunit ()
      open (unit=ihis,file=ostfile,status='old')
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      read (ihis,'(a)')  record
      read (record,*)  i,i,i,i,i,nsave0,i
      do i = 1, 5
         read (ihis,'(a)')  record
      end do
      nline = 0
      do
         read (ihis,'(a)',iostat=ios)  record
         if (ios .ne. 0)  exit
         nline = nline + 1
      end do
      close (unit=ihis,status='delete')
      call assert_int (nsave0,3,'saveost updates history count')
      call assert_int (nline,3,'saveost history record count')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_index  --  packed index tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_index" checks conversion between lambda
c     grid coordinates and packed histogram indices
c
c
      subroutine test_eost_index
      implicit none
      integer i,j,k
c
c
c     test ij_to_k and k_to_ij
c
      call ij_to_k (3,4,7,k)
      call assert_int (k,24,
     &                 'ij_to_k input i=3 j=4 nrow=7')
      call k_to_ij (24,7,i,j)
      call assert_int (i,3,
     &                 'k_to_ij input k=24 nrow=7 i')
      call assert_int (j,4,
     &                 'k_to_ij input k=24 nrow=7 j')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_resize  --  history resize tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_resize" checks OST history resizing preserves
c     old entries and initializes new storage
c
c
      subroutine test_eost_resize
      use dlmda
      use ost
      implicit none
      integer i
c
c
c     set small history arrays and resize them
c
      call resetost (5,5,2)
      nlmdahist = 2
      do i = 1, 2
         osthist(i) = 10 + i
         ostnext(i) = i - 1
         lmdalhist(i) = 0.25d0 * dble(i)
         lmdafhist(i) = -3.0d0 + 2.0d0*dble(i)
         osthhist(i) = 1.0d0 + dble(i)
         ostwlhist(i) = 0.25d0
         ostwfhist(i) = 1.0d0
      end do
      call resizeosthist
c
c     check old data are preserved and new slots are initialized
c
      call assert_int (sizelmdahist,4,
     &                 'resizeosthist size input size=2')
      do i = 1, 4
         if (i .le. 2) then
      call assert_int (osthist(i),10+i,
     &                 'resizeosthist preserve osthist')
      call assert_int (ostnext(i),i-1,
     &                 'resizeosthist preserve ostnext')
      call assert_real (lmdalhist(i),0.25d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve lmdalhist')
      call assert_real (lmdafhist(i),-3.0d0+2.0d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve lmdafhist')
      call assert_real (osthhist(i),1.0d0+dble(i),1.0d-12,
     &                  'resizeosthist preserve osthhist')
      call assert_real (ostwlhist(i),0.25d0,1.0d-12,
     &                  'resizeosthist preserve ostwlhist')
      call assert_real (ostwfhist(i),1.0d0,1.0d-12,
     &                  'resizeosthist preserve ostwfhist')
         else
      call assert_int (osthist(i),0,
     &                 'resizeosthist init osthist')
      call assert_int (ostnext(i),0,
     &                 'resizeosthist init ostnext')
      call assert_real (lmdalhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init lmdalhist')
      call assert_real (lmdafhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init lmdafhist')
      call assert_real (osthhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init osthhist')
      call assert_real (ostwlhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostwlhist')
      call assert_real (ostwfhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostwfhist')
         end if
      end do
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_eost_buildindex  --  index tests  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_eost_buildindex" checks linked-list construction
c     for gaussian history bins
c
c
      subroutine test_eost_buildindex
      use dlmda
      use ost
      implicit none
      integer k,ilmda,iflmda
c
c
c     put two gaussians in the same bin and one in another bin
c
      call resetost (5,5,3)
      nlmdahist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0,wlmda,wflmda)
      call buildostindex
c
c     expected linked list head is newest entry in the shared bin
c
      ilmda = 3
      iflmda = 3
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call assert_int (osthist(1),k,
     &                 'buildostindex packed shared bin')
      call assert_int (osthead(ilmda,iflmda),2,
     &                 'buildostindex head shared bin')
      call assert_int (ostnext(2),1,
     &                 'buildostindex next shared bin')
      call assert_int (ostnext(1),0,
     &                 'buildostindex tail shared bin')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_ensure  --  flambda resize tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_ensure" checks flambda grid expansion,
c     kernel preservation and index rebuilding
c
c
      subroutine test_eost_ensure
      use dlmda
      use ost
      implicit none
      integer k
      integer nold
      integer oldfli0
      integer ilmda,iflmda
c
c
c     high-side expansion preserves fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (2000.0d0)
      call assert_int (nflmda,2105,
     &                 'ensureflambda high nflmda')
      call assert_int (fli0,3,
     &                 'ensureflambda high fli0')
      call assert_real (gkernel(2,3),7.0d0,1.0d-12,
     &                  'ensureflambda high preserve gkernel')
c
c     low-side expansion shifts fli0 and old gkernel values
c
      call resetost (3,5,1)
      gkernel(2,3) = 7.0d0
      call ensureflambda (-2000.0d0)
      call assert_int (nflmda,2105,
     &                 'ensureflambda low nflmda')
      call assert_int (fli0,2103,
     &                 'ensureflambda low fli0')
      call assert_real (gkernel(2,2103),7.0d0,1.0d-12,
     &                  'ensureflambda low shifted gkernel')
c
c     low-side expansion also rebuilds osthead, ostnext and osthist
c
      call resetost (3,5,2)
      nlmdahist = 2
      call sethist (1,0.50d0,0.0d0,1.0d0,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0,wlmda,wflmda)
      call buildostindex
      call ensureflambda (-2000.0d0)
c
c     after low-side expansion, fli0 shifts from 3 to 2103, so
c     flambda=0.0 now belongs to iflmda=2103
c
      ilmda = 2
      iflmda = 2103
      call ij_to_k (ilmda,iflmda,nlmda,k)
      call assert_int (osthist(1),k,
     &                 'ensureflambda rebuild osthist(1)')
      call assert_int (osthist(2),k,
     &                 'ensureflambda rebuild osthist(2)')
      call assert_int (osthead(ilmda,iflmda),2,
     &                 'ensureflambda rebuild osthead')
      call assert_int (ostnext(2),1,
     &                 'ensureflambda rebuild ostnext head')
      call assert_int (ostnext(1),0,
     &                 'ensureflambda rebuild ostnext tail')
c
c     the resize copies four separate kernels, so seed all of them
c     with values that identify both the array and the bin, then
c     require every entry of every array to land where it belongs
c
      call resetost (3,5,1)
      nold = nflmda
      oldfli0 = fli0
      call seedkernels
      call ensureflambda (2000.0d0)
      call checkkernels ('ensureflambda high',nold,fli0-oldfli0)
      call resetost (3,5,1)
      nold = nflmda
      oldfli0 = fli0
      call seedkernels
      call ensureflambda (-2000.0d0)
      call checkkernels ('ensureflambda low',nold,fli0-oldfli0)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_eost_gkernels  --  g kernel tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_eost_gkernels" checks gaussian kernel updates,
c     mirroring, interpolation and lookup behavior
c
c
      subroutine test_eost_gkernels
      use dlmda
      use math
      use mutant
      use ost
      implicit none
      real*8 egbias,dgdl,dgdfl
      real*8 expected
      real*8 height
      real*8 targetl,targetf
c
c
c     choose height so the normalized gaussian prefactor is one
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     addgkernelhist from a single interior gaussian
c
      call addgkernelhist (1)
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'addgkernelhist center value')
      expected = exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'addgkernelhist off-center value')
c
c     addgkernelhist includes left-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addgkernelhist (1)
      call assert_real (gkernel(1,3),2.0d0,1.0d-12,
     &                  'addgkernelhist left mirror center')
      expected = 2.0d0 * exp(-0.5d0)
      call assert_real (gkernel(2,3),expected,1.0d-12,
     &                  'addgkernelhist left mirror neighbor')
c
c     addgkernelhist includes right-boundary mirror image
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      call addgkernelhist (1)
      call assert_real (gkernel(5,3),2.0d0,1.0d-12,
     &                  'addgkernelhist right mirror center')
      expected = 2.0d0 * exp(-0.5d0)
      call assert_real (gkernel(4,3),expected,1.0d-12,
     &                  'addgkernelhist right mirror neighbor')
c
c     rebuild the original interior gaussian for later checks
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
c
c     buildgkernel rebuilds the same full grid
c
      call buildgkernel
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'buildgkernel center value')
      expected = exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'buildgkernel off-center value')
c
c     updategkernel adds only the newest gaussian to current grid
c
      gkernel = 0.0d0
      call updategkernel
      call assert_real (gkernel(3,3),1.0d0,1.0d-12,
     &                  'updategkernel center value')
c
c     adding a taller gaussian updates only from the new history entry
c
      nlmdahist = 2
      call sethist (2,0.75d0,1.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
      expected = 2.0d0 + exp(-1.0d0)
      call assert_real (gkernel(4,4),expected,1.0d-12,
     &                  'updategkernel added new center')
      expected = 1.0d0 + 2.0d0*exp(-1.0d0)
      call assert_real (gkernel(3,3),expected,1.0d-12,
     &                  'updategkernel keeps old plus new overlap')
c
c     egkernel evaluates the continuous gaussian and derivatives
c
      lambda = 0.75d0
      dedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 2.0d0 + exp(-1.0d0)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel bias lambda=.75 flambda=1')
      call assert_real (dgdl,-4.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel dgdl lambda=.75 flambda=1')
      call assert_real (dgdfl,-exp(-1.0d0),1.0d-12,
     &                  'egkernel dgdfl lambda=.75 flambda=1')
c
c     two gaussians in the same bin are both followed by ostnext
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 2
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.50d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel same-bin linked bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel same-bin linked dgdl')
      call assert_real (dgdfl,0.0d0,1.0d-12,
     &                  'egkernel same-bin linked dgdfl')
c
c     multiple bins are found through osthead lookup, including
c     two gaussians in one bin and one gaussian in another bin
c
      call resetost (5,5,4)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0*height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0*height,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      lambda = 0.75d0
      dedl = 1.0d0
      call egkernel (egbias,dgdl,dgdfl)
      expected = 3.0d0 + 3.0d0*exp(-1.0d0)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel multi-bin mixed bias')
      call assert_real (dgdl,-12.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel multi-bin mixed dgdl')
      call assert_real (dgdfl,-3.0d0*exp(-1.0d0),1.0d-12,
     &                  'egkernel multi-bin mixed dgdfl')
c
c     left endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.0d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel left mirror bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel left mirror dgdl')
c
c     right endpoint includes both real and mirror gaussian images
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 1.0d0
      dedl = 0.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,2.0d0,1.0d-12,
     &                  'egkernel right mirror bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel right mirror dgdl')
c
c     wide histogram widths are used for both grid and continuous bias
c
      call resetost (41,81,3)
      wlhist = 0.05d0
      wfhist = 10.0d0
      maxwlhist = wlhist
      maxwfhist = wfhist
      height = 2.0d0 * pi * wlhist * wfhist
      targetl = 0.525d0
      targetf = 5.0d0
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlhist,wfhist)
      call buildostindex
      call buildgkernel
      expected = exp(-0.25d0)
      call assert_real (gkernel(22,46),expected,1.0d-12,
     &                  'buildgkernel wide hist width')
      lambda = targetl
      dedl = targetf
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,expected,1.0d-12,
     &                  'egkernel wide hist width bias')
      call assert_real (dgdl,-10.0d0*expected,1.0d-12,
     &                  'egkernel wide hist width dgdl')
      call assert_real (dgdfl,-0.05d0*expected,1.0d-12,
     &                  'egkernel wide hist width dgdfl')
c
c     if dU/dlambda is outside the current flambda grid,
c     egkernel returns zero
c
      call resetost (5,5,3)
      height = 2.0d0 * pi * wlmda * wflmda
      nlmdahist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      lambda = 0.50d0
      dedl = 100.0d0
      call egkernel (egbias,dgdl,dgdfl)
      call assert_real (egbias,0.0d0,1.0d-12,
     &                  'egkernel outside flambda bias')
      call assert_real (dgdl,0.0d0,1.0d-12,
     &                  'egkernel outside flambda dgdl')
      call assert_real (dgdfl,0.0d0,1.0d-12,
     &                  'egkernel outside flambda dgdfl')
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine test_eost_fkernel  --  f kernel tests  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "test_eost_fkernel" checks mean-force kernel
c     construction from gaussian history data
c
c
      subroutine test_eost_fkernel
      use bath
      use dlmda
      use math
      use ost
      use units
      implicit none
      real*8 rt
      real*8 expected
      real*8 height
      real*8 wratio
      real*8 w1,w2,w3,w4,w5
c
c     two nonzero gkernel weights at flambda=-1 and +1
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = log(2.0d0) * rt
      gkernel(3,4) = log(4.0d0) * rt
      vkernelmax(3) = log(4.0d0) * rt
      call buildfkernel
      expected = 1.0d0 / 3.0d0
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel weighted mean force')
      call assert_real (lmdafmean(1),0.0d0,1.0d-12,
     &                  'buildfkernel empty row')
c
c     build gkernel incrementally from multiple gaussians using
c     updategkernel, then construct the f kernel
c
      call resetost (5,5,4)
      rt = gasconst * kelvin
      height = 2.0d0 * pi * wlmda * wflmda
c
      nlmdahist = 1
      call sethist (1,0.25d0,-1.0d0,1.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nlmdahist = 2
      call sethist (2,0.50d0, 0.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nlmdahist = 3
      call sethist (3,0.75d0, 1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      call buildfkernel
c
c     expected weighted mean at lambda=0.5
c
      w1 = exp(gkernel(3,1)/rt)
      w2 = exp(gkernel(3,2)/rt)
      w3 = exp(gkernel(3,3)/rt)
      w4 = exp(gkernel(3,4)/rt)
      w5 = exp(gkernel(3,5)/rt)
      
      expected = (-2.d0*w1 - w2 + w4 + 2.d0*w5)
     &     /(w1+w2+w3+w4+w5)
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel from updated gkernel')
c
c     a deeply filled bin overflows exp(g/kT) unless the largest bias
c     in the row is factored out; unshifted this returns a NaN
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = 500.0d0
      gkernel(3,4) = 499.0d0
      vkernelmax(3) = 500.0d0
      call buildfkernel
      wratio = exp((499.0d0-500.0d0)/rt)
      expected = (-1.0d0 + wratio) / (1.0d0 + wratio)
      call assert_real (lmdafmean(3),expected,1.0d-12,
     &                  'buildfkernel large bias no overflow')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_kernelbuilds  --  build tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_kernelbuilds" checks full kernel rebuilds
c     and incremental updates against reference arrays
c
c
      subroutine test_eost_kernelbuilds
      use bath
      use dlmda
      use math
      use ost
      use units
      implicit none
      integer i,j
      integer ihist
      integer nhist
      real*8 rt
      real*8 flmda
      real*8 partfunc
      real*8 fsum
      real*8 weight
      real*8 vmaxref
      real*8 height
      real*8 fmanual(9)
      real*8 fref(9)
      real*8 fsumref(9)
      real*8 gref(9,9)
      real*8 pfref(9)
c
c
c     build a mixed history with overlapping gaussians, endpoint
c     mirror images, and multiple gaussian widths
c
      call resetost (9,9,8)
      rt = gasconst * kelvin
      nhist = 6
      nlmdahist = nhist
      height = 2.0d0 * pi * wlmda * wflmda
      call sethist (1,0.00d0, 0.0d0,0.7d0*height,wlmda,wflmda)
      call sethist (2,0.25d0,-1.0d0,1.1d0*height,wlmda,wflmda)
      call sethist (3,0.50d0, 0.0d0,1.6d0*height,wlmda,wflmda)
      call sethist (4,0.50d0, 0.0d0,0.4d0*height,wlmda,wflmda)
      call sethist (5,0.75d0, 1.0d0,2.3d0*height,wlmda,wflmda)
      call sethist (6,1.00d0, 0.0d0,0.9d0*height,wlmda,wflmda)
      ostwlhist(5) = 2.0d0 * wlmda
      ostwfhist(5) = 2.0d0 * wflmda
      call buildostindex
c
c     buildgkernel must clear stale values before rebuilding
c
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = -123.0d0
         end do
      end do
      call buildgkernel
      call buildfkernel
      do i = 1, nlmda
         fref(i) = lmdafmean(i)
         do j = 1, nflmda
            gref(i,j) = gkernel(i,j)
         end do
      end do
      call assert_real (gkernel(1,1),0.0d0,1.0d-12,
     &                  'buildgkernel clears outside support')
c
c     compute the buildfkernel result independently from gkernel
c
      do i = 1, nlmda
         partfunc = 0.0d0
         fsum = 0.0d0
c
c     the accumulators factor out the largest bias in the row, so the
c     reference has to be built with the same vkernelmax shift
c
         vmaxref = 0.0d0
         do j = 1, nflmda
            if (gref(i,j) .ne. 0.0d0)
     &         vmaxref = max(vmaxref,gref(i,j))
         end do
         do j = 1, nflmda
            if (gref(i,j) .ne. 0.0d0) then
               flmda = dble(j-fli0) * wflmda
               weight = exp((gref(i,j)-vmaxref)/rt)
               fsum = fsum + flmda*weight
               partfunc = partfunc + weight
            end if
         end do
         if (partfunc .eq. 0.0d0) then
            fmanual(i) = 0.0d0
         else
            fmanual(i) = fsum / partfunc
         end if
         fsumref(i) = fsum
         pfref(i) = partfunc
      end do
      call assert_array1 (lmdafmean,fmanual,nlmda,1.0d-12,
     &                    'buildfkernel all-row weighted mean')
c
c     buildkernels must produce the same gkernel, lmdafmean and
c     free energy accumulators as the old full rebuild path
c
      do i = 1, nlmda
         lmdafmean(i) = -123.0d0
         lmdafsum(i) = -123.0d0
         lmdafwt(i) = -123.0d0
         do j = 1, nflmda
            gkernel(i,j) = -123.0d0
         end do
      end do
      call buildkernels
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    'buildkernels gkernel reference')
      call assert_array1 (lmdafmean,fref,nlmda,1.0d-12,
     &                    'buildkernels lmdafmean reference')
      call assert_array1 (lmdafsum,fsumref,nlmda,1.0d-12,
     &                    'buildkernels fsum reference')
      call assert_array1 (lmdafwt,pfref,nlmda,1.0d-12,
     &                    'buildkernels partfunc reference')
c
c     incrementally updating kernels one history at a time must match
c     a full buildkernels rebuild
c
      call resetost (9,9,8)
      height = 2.0d0 * pi * wlmda * wflmda
      do ihist = 1, nhist
         nlmdahist = ihist
         if (ihist .eq. 1) then
            call sethist (ihist,0.00d0,0.0d0,0.7d0*height,wlmda,wflmda)
         else if (ihist .eq. 2) then
            call sethist (ihist,0.25d0,-1.0d0,1.1d0*height,wlmda,wflmda)
         else if (ihist .eq. 3) then
            call sethist (ihist,0.50d0,0.0d0,1.6d0*height,wlmda,wflmda)
         else if (ihist .eq. 4) then
            call sethist (ihist,0.50d0,0.0d0,0.4d0*height,wlmda,wflmda)
         else if (ihist .eq. 5) then
            call sethist (ihist,0.75d0,1.0d0,2.3d0*height,wlmda,wflmda)
            ostwlhist(ihist) = 2.0d0 * wlmda
            ostwfhist(ihist) = 2.0d0 * wflmda
         else
            call sethist (ihist,1.00d0,0.0d0,0.9d0*height,wlmda,wflmda)
         end if
         call buildostindex
         call updatekernels
      end do
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    'updatekernels gkernel reference')
      call assert_array1 (lmdafmean,fref,nlmda,1.0d-12,
     &                    'updatekernels lmdafmean reference')
      call assert_array1 (lmdafsum,fsumref,nlmda,1.0d-12,
     &                    'updatekernels fsum reference')
      call assert_array1 (lmdafwt,pfref,nlmda,1.0d-12,
     &                    'updatekernels partfunc reference')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_avgstd  --  OST average tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_avgstd" checks interval averages and
c     standard deviations for sampled lambda values
c
c
      subroutine test_eost_avgstd
      use dlmda
      implicit none
      integer i
      real*8 stdref
c
c
c     average all saved interval samples after equilibration prefix
c
      call resetost (5,5,1)
      lmdaintv = 6
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 4
      do i = 1, lmdaintv
         lmdallist(i) = dble(i)
         lmdaflist(i) = 2.0d0*dble(i)
      end do
      call avgstd (lmdallist,lmdanpa+lmdanpb+1,lmdanpc,
     &             lmdaavg,lmdastd)
      call avgstd (lmdaflist,lmdanpa+lmdanpb+1,lmdanpc,
     &             dedlavg,dedlstd)
      stdref = sqrt(1.25d0)
      call assert_real (lmdaavg,4.5d0,1.0d-12,
     &                  'avgstd configurable lambda average')
      call assert_real (dedlavg,9.0d0,1.0d-12,
     &                  'avgstd configurable dE/dl average')
      call assert_real (lmdastd,stdref,1.0d-12,
     &                  'avgstd configurable lambda std')
      call assert_real (dedlstd,2.0d0*stdref,1.0d-12,
     &                  'avgstd configurable dE/dl std')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eost_eginterpolate  --  g interp tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eost_eginterpolate" checks interpolated gaussian
c     bias values and derivatives
c
c
      subroutine test_eost_eginterpolate
      use dlmda
      use math
      use mutant
      use ost
      implicit none
      real*8 egbias0,dgdl0,dgdfl0
      real*8 egbias1,dgdl1,dgdfl1
      real*8 height
      real*8 sigl,sigf
c
c
c     grid-point interpolation should reproduce the analytic gaussian
c     sum and derivatives from egkernel exactly to roundoff
c
      call resetost (9,9,4)
      sigl = 2.0d0 * wlmda
      sigf = 2.0d0 * wflmda
      height = 2.0d0 * pi * sigl * sigf
      oststdev = 4.0d0
      nlmdahist = 3
      call sethist (1,0.25d0,-1.0d0,1.1d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.6d0*height,sigl,sigf)
      call sethist (3,0.75d0, 1.0d0,2.3d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      lambda = 0.50d0
      dedl = 0.0d0
      ostinterpol = .false.
      call egkernel (egbias0,dgdl0,dgdfl0)
      call egkernelinterpolate (egbias1,dgdl1,dgdfl1)
      call assert_real (egbias1,egbias0,1.0d-12,
     &                  'egkernel interpolate grid bias')
      call assert_real (dgdl1,dgdl0,1.0d-12,
     &                  'egkernel interpolate grid dgdl')
      call assert_real (dgdfl1,dgdfl0,1.0d-12,
     &                  'egkernel interpolate grid dgdfl')
c
c     off-grid interpolation is approximate; wide gaussians should
c     keep the value and derivative errors small
c
      call resetost (17,17,4)
      sigl = 4.0d0 * wlmda
      sigf = 4.0d0 * wflmda
      height = 2.0d0 * pi * sigl * sigf
      oststdev = 4.0d0
      nlmdahist = 3
      call sethist (1,0.25d0,-2.0d0,0.8d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.2d0*height,sigl,sigf)
      call sethist (3,0.75d0, 2.0d0,1.6d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      lambda = 0.53125d0
      dedl = 0.25d0
      ostinterpol = .false.
      call egkernel (egbias0,dgdl0,dgdfl0)
      call egkernelinterpolate (egbias1,dgdl1,dgdfl1)
      call assert_real (egbias1,egbias0,1.0d-3,
     &                  'egkernel interpolate offgrid bias')
      call assert_real (dgdl1,dgdl0,2.0d-2,
     &                  'egkernel interpolate offgrid dgdl')
      call assert_real (dgdfl1,dgdfl0,2.0d-2,
     &                  'egkernel interpolate offgrid dgdfl')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_efkernel  --  free-energy tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_efkernel" checks DeltaG interpolation and
c     total free-energy integration from f kernels
c
c
      subroutine test_eost_efkernel
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 eostlmda,dfdl
      real*8 expected
      real*8 efreetot
c
c
c     use lmdafmean(lambda)=lambda and integrate to lambda=0.375
c
      call resetost (5,5,1)
      do i = 1, nlmda
         lmdafmean(i) = dble(i-1) * wlmda
      end do
      lambda = 0.375d0
      call efreelmda (eostlmda,dfdl)
      expected = 0.5d0 * lambda * lambda
      call assert_real (eostlmda,expected,1.0d-12,
     &                  'efreelmda DeltaG lambda=.375')
      call assert_real (dfdl,lambda,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda')
c
c     use lmdafmean(lambda)=1+lambda so the endpoint mean forces are
c     nonzero and distinct from each other
c
      call resetost (5,5,1)
      do i = 1, nlmda
         lmdafmean(i) = 1.0d0 + dble(i-1)*wlmda
      end do
c
c     at and below lambda = 0 the free energy is zero and the
c     derivative comes from the first lambda bin
c
      lambda = 0.0d0
      call efreelmda (eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efreelmda DeltaG lambda=0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda lambda=0')
      lambda = -0.25d0
      call efreelmda (eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efreelmda DeltaG lambda below 0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda below 0')
c
c     at lambda = 1 the last interval is integrated in full, and
c     beyond lambda = 1 the loop falls through to the same result
c
      lambda = 1.0d0
      call efreelmda (eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efreelmda DeltaG lambda=1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda lambda=1')
      lambda = 1.25d0
      call efreelmda (eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efreelmda DeltaG lambda above 1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efreelmda dDeltaG/dlambda above 1')
c
c     the trapezoid rule is exact for a linear mean force, so the
c     total must match the efreelmda value at lambda = 1
c
      call assert_real (efreetot(),1.5d0,1.0d-12,
     &                  'efreetot linear lmdafmean')
c
c     a constant mean force integrates to itself over unit lambda
c
      do i = 1, nlmda
         lmdafmean(i) = 2.5d0
      end do
      call assert_real (efreetot(),2.5d0,1.0d-12,
     &                  'efreetot constant lmdafmean')
c
c     an all-zero mean force gives no free energy change
c
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
      end do
      call assert_real (efreetot(),0.0d0,1.0d-12,
     &                  'efreetot zero lmdafmean')
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_eost_meta  --  metadynamics tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_eost_meta" checks metadynamics bias evaluation,
c     free-energy difference and history resizing
c
c
      subroutine test_eost_meta
      use math
      use ost
      implicit none
      real*8 vbias,dvdl
      real*8 pref
      real*8 expected
      real*8 metadeltag
c
c
c     one normalized 1D gaussian centered at lambda=0.5
c
      call resetost (5,5,1)
      call resetmeta (2)
      nmetahist = 1
      metalhist(1) = 0.5d0
      metahhist(1) = 2.0d0
      metawhist(1) = 0.25d0
      pref = metahhist(1) / (metawhist(1)*sqrt(2.0d0*pi))
c
c     a gaussian centered midway is symmetric in its two images
c
      call emetabias (0.5d0,vbias,dvdl)
      expected = pref * (1.0d0 + 2.0d0*exp(-8.0d0))
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias center value')
      call assert_real (dvdl,0.0d0,1.0d-12,
     &                  'emetabias center derivative')
      call emetabias (0.75d0,vbias,dvdl)
      expected = pref * (exp(-0.5d0)+exp(-12.5d0)+exp(-4.5d0))
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias off-center value')
      expected = -16.0d0 * pref
     &              * (0.25d0*exp(-0.5d0) + 1.25d0*exp(-12.5d0)
     &                 - 0.75d0*exp(-4.5d0))
      call assert_real (dvdl,expected,1.0d-12,
     &                  'emetabias off-center derivative')
c
c     symmetric gaussian has zero endpoint free energy difference
c
      call assert_real (metadeltag(),0.0d0,1.0d-12,
     &                  'metadeltag symmetric gaussian')
c
c     resizing preserves old metadynamics history and zeros new slots
c
      nmetahist = 2
      metalhist(2) = 0.25d0
      metahhist(2) = 3.0d0
      metawhist(2) = 0.125d0
      call resizemeta
      call assert_int (sizemetahist,4,
     &                 'resizemeta size input size=2')
      call assert_real (metalhist(2),0.25d0,1.0d-12,
     &                  'resizemeta preserve metalhist(2)')
      call assert_real (metalhist(3),0.0d0,1.0d-12,
     &                  'resizemeta init metalhist(3)')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_histstat  --  interval stat test  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_histstat" checks the average, deviation and fitted
c     drift that histstat computes over the averaging phase
c
c
      subroutine test_eost_histstat
      use dlmda
      use ost
      implicit none
      integer i
      real*8 stdref
c
c
c     average the samples following the equilibration prefix
c
      call resetost (5,5,1)
      lmdaintv = 6
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 4
      do i = 1, lmdaintv
         lmdallist(i) = dble(i)
         lmdaflist(i) = 2.0d0*dble(i)
      end do
      call histstat (lmdallist,lmdaavg,lmdastd,ostlambdaslp)
      call histstat (lmdaflist,dedlavg,dedlstd,ostdedlslp)
      stdref = sqrt(1.25d0)
      call assert_real (lmdaavg,4.5d0,1.0d-12,
     &                  'histstat lambda average')
      call assert_real (dedlavg,9.0d0,1.0d-12,
     &                  'histstat dE/dl average')
      call assert_real (lmdastd,stdref,1.0d-12,
     &                  'histstat lambda deviation')
      call assert_real (dedlstd,2.0d0*stdref,1.0d-12,
     &                  'histstat dE/dl deviation')
c
c     the fitted drift keeps the scale of each ramp
c
      call assert_real (ostlambdaslp,1.0d0,1.0d-12,
     &                  'histstat lambda slope')
      call assert_real (ostdedlslp,2.0d0,1.0d-12,
     &                  'histstat dE/dl slope')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_drift  --  interval drift tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_drift" checks the fitted drift of flat, ramped,
c     folded and offset series over the averaging phase
c
c
      subroutine test_eost_drift
      use dlmda
      use ost
      implicit none
      integer i
      real*8 v(8)
      data v / 4.0d0,3.0d0,2.0d0,1.0d0,1.0d0,2.0d0,3.0d0,4.0d0 /
c
c
c     a flat series has no drift
c
      call resetost (5,5,1)
      lmdaintv = 8
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 8
      do i = 1, lmdaintv
         lmdallist(i) = 7.0d0
      end do
      call histstat (lmdallist,lmdaavg,lmdastd,ostlambdaslp)
      call assert_real (lmdaavg,7.0d0,1.0d-12,
     &                  'histstat flat average')
      call assert_real (ostlambdaslp,0.0d0,1.0d-12,
     &                  'histstat flat slope')
c
c     a decreasing ramp keeps its change per sample
c
      do i = 1, lmdaintv
         lmdallist(i) = -0.5d0*dble(i-1)
      end do
      call histstat (lmdallist,lmdaavg,lmdastd,ostlambdaslp)
      call assert_real (ostlambdaslp,-0.5d0,1.0d-12,
     &                  'histstat ramp slope')
c
c     a folded series has no net drift
c
      do i = 1, lmdaintv
         lmdallist(i) = v(i)
      end do
      call histstat (lmdallist,lmdaavg,lmdastd,ostlambdaslp)
      call assert_real (ostlambdaslp,0.0d0,1.0d-12,
     &                  'histstat folded slope')
c
c     a large offset must not swamp a small drift
c
      do i = 1, lmdaintv
         lmdallist(i) = 5000.0d0 + 1.0d-6*dble(i-1)
      end do
      call histstat (lmdallist,lmdaavg,lmdastd,ostlambdaslp)
      call assert_real (ostlambdaslp,1.0d-6,1.0d-9,
     &                  'histstat offset slope')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_eost_depcriteria  --  deposit gate tests  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_eost_depcriteria" checks the deposition gate against a
c     tolerance built from an absolute and a relative part
c
c
      subroutine test_eost_depcriteria
      use ost
      implicit none
      logical depcriteria
c
c
c     the tolerance is ostcvstd plus ostcvrat times the average
c
      ostcvstd = 10.0d0
      ostcvrat = 0.2d0
      call assert_logical (depcriteria(0.0d0,9.9d0),.true.,
     &                     'depcriteria inside absolute tolerance')
      call assert_logical (depcriteria(0.0d0,10.0d0),.false.,
     &                     'depcriteria at absolute tolerance')
      call assert_logical (depcriteria(50.0d0,19.9d0),.true.,
     &                     'depcriteria inside relative tolerance')
      call assert_logical (depcriteria(50.0d0,20.0d0),.false.,
     &                     'depcriteria at relative tolerance')
      call assert_logical (depcriteria(-50.0d0,19.9d0),.true.,
     &                     'depcriteria relative tolerance is signless')
      call assert_logical (depcriteria(-50.0d0,20.0d0),.false.,
     &                     'depcriteria negative average at tolerance')
c
c     a vanishing tolerance rejects every interval
c
      ostcvstd = 0.0d0
      ostcvrat = 0.0d0
      call assert_logical (depcriteria(1.0d0,0.0d0),.false.,
     &                     'depcriteria zero tolerance')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_vkernelmax  --  bias level test  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_vkernelmax" checks the running maximum of the g
c     kernel along flambda and the global path bias level derived
c     from it, over every path that fills the kernel
c
c
      subroutine test_eost_vkernelmax
      use bath
      use dlmda
      use ost
      implicit none
      integer i
      real*8 vmin
      real*8 brutevkmax
      real*8 ostvminimax
      real*8, allocatable :: gsave(:,:)
c
c
c     two saved sources spread over the whole kernel
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      nlmdahist = 2
      call sethist (1,0.25d0,0.0d0,1.0d0,0.25d0,1.0d0)
      call sethist (2,0.75d0,1.0d0,2.0d0,0.25d0,1.0d0)
      call buildostindex
      call buildkernels
      vmin = brutevkmax (1)
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'buildkernels running maximum')
         vmin = min(vmin,brutevkmax(i))
      end do
      call assert_real (ostvminimax(),vmin,1.0d-12,
     &                  'ostvminimax over lambda')
      call assert_logical (vmin.gt.0.0d0,.true.,
     &                     'ostvminimax reaches every lambda bin')
c
c     the incremental update path stays exact
c
      nlmdahist = 3
      call sethist (3,0.5d0,-1.0d0,1.5d0,0.25d0,1.0d0)
      call buildostindex
      call updatekernels
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'updatekernels running maximum')
      end do
c
c     the g kernel only path reproduces both kernel and maximum
c
      allocate (gsave(nlmda,nflmda))
      gsave = gkernel
      call buildgkernel
      call assert_array2 (gkernel,gsave,nlmda,nflmda,1.0d-12,
     &                    'buildgkernel reproduces gkernel')
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'buildgkernel running maximum')
      end do
      deallocate (gsave)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function brutevkmax  --  direct kernel column maximum  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "brutevkmax" scans the g kernel directly for the largest
c     value at one lambda bin, ignoring the running maximum
c
c
      function brutevkmax (ilmda)
      use ost
      implicit none
      integer ilmda,j
      real*8 brutevkmax
c
c
      brutevkmax = 0.0d0
      do j = 1, nflmda
         brutevkmax = max(brutevkmax,gkernel(ilmda,j))
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_tempering  --  height decay test  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_tempering" checks the deposited gaussian height with
c     the global and local tempering factors disabled, below their
c     thresholds, active alone and active together, and that the
c     height stays positive and never exceeds the full height
c
c
      subroutine test_eost_tempering
      use bath
      use ost
      use units
      implicit none
      integer i,j,k
      real*8 rt
      real*8 h,prev
      real*8 h1,h2
      real*8 hg,hl
      real*8 vg,dl
      real*8 temperedheight
      real*8 vstar(4)
      real*8 delta(5)
      real*8 vgrid(6)
      real*8 dgrid(6)
      real*8 gpair(2,3)
      real*8 hgrid(6,6)
      data vstar / 1.5d0,2.0d0,3.0d0,5.0d0 /
      data delta / 0.0d0,1.0d0,1.5d0,3.0d0,5.0d0 /
      data vgrid / 0.0d0,0.5d0,1.0d0,2.0d0,5.0d0,20.0d0 /
      data dgrid / 0.0d0,0.5d0,1.0d0,2.0d0,5.0d0,20.0d0 /
      data gpair / 1.0d0,1.0d0,2.0d0,0.5d0,0.5d0,2.0d0 /
c
c
c     an untempered run deposits at the full height
c
      kelvin = 300.0d0
      call resetost (5,5,1)
      hbias = 1.0d-5
      rt = gasconst * kelvin
      use_ostgtemp = .false.
      use_ostltemp = .false.
      ostgthresh = 1.0d0
      ostgtempgamma = 1.0d0
      ostlthresh = 1.0d0
      ostltempgamma = 1.0d0
      call assert_real (temperedheight(0.0d0,0.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled empty path')
      call assert_real (temperedheight(50.0d0,50.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled filled path')
      call assert_real (temperedheight(50.0d0,500.0d0),hbias,1.0d-18,
     &                  'temperedheight disabled uneven path')
c
c     at or below both thresholds the height is still untempered
c
      use_ostgtemp = .true.
      use_ostltemp = .true.
      call assert_real (temperedheight(0.0d0,0.0d0),hbias,1.0d-18,
     &                  'temperedheight below both thresholds')
      call assert_real (temperedheight(1.0d0,1.0d0),hbias,1.0d-18,
     &                  'temperedheight at global threshold')
      call assert_real (temperedheight(1.0d0,2.0d0),hbias,1.0d-18,
     &                  'temperedheight at both thresholds')
      call assert_real (temperedheight(0.5d0,1.5d0),hbias,1.0d-18,
     &                  'temperedheight at local threshold')
c
c     the global factor alone decays with the path bias level and
c     ignores how far the deposit bin is ahead of that level
c
      use_ostltemp = .false.
      prev = hbias
      do i = 1, 4
         h = temperedheight(vstar(i),vstar(i))
         call assert_real (h,hbias*exp(-(vstar(i)-1.0d0)/rt),1.0d-18,
     &                     'temperedheight global above threshold')
         call assert_real (temperedheight(vstar(i),vstar(i)+50.0d0),
     &                     h,1.0d-18,'temperedheight global ignores '//
     &                     'the local excess')
         call assert_logical (h.lt.prev,.true.,
     &                        'temperedheight global decays')
         prev = h
      end do
c
c     a larger global tempering factor decays more slowly
c
      ostgtempgamma = 1.0d0
      h1 = temperedheight(3.0d0,3.0d0)
      ostgtempgamma = 2.0d0
      h2 = temperedheight(3.0d0,3.0d0)
      call assert_logical (h2.gt.h1,.true.,
     &                     'temperedheight larger gamma decays slower')
      call assert_real (h2,hbias*exp(-2.0d0/(2.0d0*rt)),1.0d-18,
     &                  'temperedheight gamma scaling')
c
c     the local factor alone depends only on the excess of the
c     deposit bin over the path bias level
c
      use_ostgtemp = .false.
      use_ostltemp = .true.
      ostgtempgamma = 1.0d0
      prev = hbias
      do i = 1, 5
         h = temperedheight(50.0d0,50.0d0+delta(i))
         hl = hbias * exp(-max(0.0d0,delta(i)-1.0d0)/rt)
         call assert_real (h,hl,1.0d-18,
     &                     'temperedheight local excess')
         call assert_real (temperedheight(0.0d0,delta(i)),h,1.0d-18,
     &                     'temperedheight local shift invariance')
         call assert_logical (h.le.prev,.true.,
     &                        'temperedheight local decays')
         prev = h
      end do
      call assert_real (temperedheight(50.0d0,50.0d0),hbias,1.0d-18,
     &                  'temperedheight least filled bin untempered')
c
c     both factors multiply, so with equal settings the path bias
c     level cancels and only the deposit bin bias level remains
c
      use_ostgtemp = .true.
      use_ostltemp = .true.
      ostgthresh = 1.0d0
      ostlthresh = 1.0d0
      ostgtempgamma = 2.0d0
      ostltempgamma = 2.0d0
      h = temperedheight(3.0d0,6.0d0)
      hg = exp(-2.0d0/(2.0d0*rt))
      hl = exp(-2.0d0/(2.0d0*rt))
      call assert_real (h,hbias*hg*hl,1.0d-18,
     &                  'temperedheight factors multiply')
      call assert_real (h,hbias*exp((2.0d0-6.0d0)/(2.0d0*rt)),1.0d-18,
     &                  'temperedheight equal settings')
      ostltempgamma = 0.5d0
      hg = exp(-2.0d0/(2.0d0*rt))
      hl = exp(-2.0d0/(0.5d0*rt))
      call assert_real (temperedheight(3.0d0,6.0d0),hbias*hg*hl,
     &                  1.0d-18,'temperedheight unequal gammas')
c
c     over a sweep of levels the height stays positive, never exceeds
c     the full height, and never grows with either level
c
      ostgthresh = 1.0d0
      ostlthresh = 1.0d0
      do k = 1, 3
         ostgtempgamma = gpair(1,k)
         ostltempgamma = gpair(2,k)
         do i = 1, 6
            do j = 1, 6
               vg = vgrid(i)
               dl = dgrid(j)
               h = temperedheight(vg,vg+dl)
               hgrid(i,j) = h
               call assert_logical (h.gt.0.0d0 .and. h.le.hbias,.true.,
     &                              'temperedheight bounded height')
               if (i .gt. 1) then
                  call assert_logical (h.le.hgrid(i-1,j),.true.,
     &                        'temperedheight nonincreasing in path')
               end if
               if (j .gt. 1) then
                  call assert_logical (h.le.hgrid(i,j-1),.true.,
     &                        'temperedheight nonincreasing in bin')
               end if
            end do
         end do
      end do
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_metaimage  --  meta image tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_metaimage" checks that the metadynamics bias sums
c     the same three reflected lambda images as the g kernel
c
c
      subroutine test_eost_metaimage
      use math
      use ost
      implicit none
      integer m,t
      real*8 vbias,dvdl
      real*8 vref,dref
      real*8 delta,bias
      real*8 pref,sig2
      real*8 single,far
      real*8 src(3)
      real*8 lam(4)
      data lam / 0.0d0,0.1d0,0.5d0,1.0d0 /
c
c
c     one gaussian near the lambda = 0 wall
c
      call resetost (5,5,1)
      call resetmeta (2)
      nmetahist = 1
      metalhist(1) = 0.1d0
      metahhist(1) = 2.0d0
      metawhist(1) = 0.25d0
      pref = metahhist(1) / (metawhist(1)*sqrt(2.0d0*pi))
      sig2 = metawhist(1) * metawhist(1)
      src(1) = 0.1d0
      src(2) = -0.1d0
      src(3) = 1.9d0
c
c     the bias and derivative match a direct three image sum
c
      do t = 1, 4
         vref = 0.0d0
         dref = 0.0d0
         do m = 1, 3
            delta = lam(t) - src(m)
            bias = pref * exp(-0.5d0*delta*delta/sig2)
            vref = vref + bias
            dref = dref - delta*bias/sig2
         end do
         call emetabias (lam(t),vbias,dvdl)
         call assert_real (vbias,vref,1.0d-12,
     &                     'emetabias image sum value')
         call assert_real (dvdl,dref,1.0d-12,
     &                     'emetabias image sum derivative')
      end do
c
c     at the wall the nearby image doubles the bias
c
      single = pref * exp(-0.5d0*0.01d0/sig2)
      far = pref * exp(-0.5d0*1.9d0*1.9d0/sig2)
      call emetabias (0.0d0,vbias,dvdl)
      call assert_real (vbias,2.0d0*single+far,1.0d-12,
     &                  'emetabias doubles at the lambda wall')
      call assert_logical (vbias.gt.single,.true.,
     &                     'emetabias wall raises the bias')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eost_metadyn  --  meta deposit tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eost_metadyn" drives one full deposit interval and
c     checks the gaussian that emetadyn stores
c
c
      subroutine test_eost_metadyn
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 avgref
      real*8 lam(4)
      data lam / 0.1d0,0.2d0,0.4d0,0.6d0 /
c
c
c     a frozen lambda particle keeps the sampled values controlled
c
      call resetost (5,5,1)
      call resetmeta (2)
      lmdaintv = 4
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 2
      hbias = 2.0d0
      wlmda = 0.25d0
      dedl = 0.0d0
      lmdadt = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = lam(istep)
         call emetadyn
         if (istep .lt. lmdaintv) then
            call assert_int (nmetahist,0,
     &                       'emetadyn waits for the interval end')
         end if
      end do
c
c     only the samples after the equilibration prefix are averaged
c
      avgref = (lam(3)+lam(4)) / dble(lmdanpc)
      call assert_int (nmetahist,1,'emetadyn deposits one gaussian')
      call assert_real (metalhist(1),avgref,1.0d-12,
     &                  'emetadyn gaussian center')
      call assert_real (metahhist(1),hbias,1.0d-12,
     &                  'emetadyn gaussian height')
      call assert_real (metawhist(1),wlmda,1.0d-12,
     &                  'emetadyn gaussian width')
      call assert_int (metaihist(1),lmdaintv,
     &                 'emetadyn gaussian step stamp')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_eost_metatemper  --  meta decay tests  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_eost_metatemper" drives several deposit intervals with
c     tempering enabled and checks the deposited heights and the
c     accumulated metadynamics grid
c
c
      subroutine test_eost_metatemper
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer il,istep,k
      integer ndep
      parameter (ndep=5)
      real*8 lmda
      real*8 vd,dd,vi,di
      real*8 refvstar
      real*8 temperedheight
c
c
c     deposit repeatedly at a fixed lambda with tempering on
c
      kelvin = 300.0d0
      call resetost (5,5,1)
      call resetmeta (8)
      lmdaintv = 4
      lmdanpa = 1
      lmdanpb = 1
      lmdanpc = 2
      hbias = 2.0d0
      dedl = 0.0d0
      lmdadt = 0.0d0
      use_ostgtemp = .true.
      ostgthresh = 0.5d0
      ostgtempgamma = 1.0d0
      lmdastep = 0
      do istep = 1, ndep*lmdaintv
         lambda = 0.5d0
         call emetadyn
      end do
      call assert_int (nmetahist,ndep,'emetadyn deposit count')
c
c     the first deposit sees an empty bias, so it is untempered
c
      call assert_real (metahhist(1),hbias,1.0d-12,
     &                  'emetadyn first height untempered')
      call assert_logical (refvstar(1).gt.ostgthresh,.true.,
     &                     'emetadyn crosses the threshold')
c
c     every later height follows the pre-deposit bias level
c
      do k = 2, ndep
         call assert_real (metahhist(k),
     &                     temperedheight(refvstar(k-1),refvstar(k-1)),
     &                     1.0d-12,'emetadyn tempered height')
         call assert_logical (metahhist(k).lt.metahhist(k-1),.true.,
     &                        'emetadyn heights decay')
      end do
c
c     the grid matches the direct sum at each bin center, and the
c     interpolation reproduces it exactly at those nodes
c
      do il = 1, nlmda
         lmda = dble(il-1) * wlmda
         ostinterpol = .false.
         call emetabias (lmda,vd,dd)
         call assert_real (vmetagrid(il),vd,1.0d-12,
     &                     'addmetagrid matches the direct sum')
         call assert_real (dvmetagrid(il),dd,1.0d-12,
     &                     'addmetagrid matches the direct slope')
         ostinterpol = .true.
         call emetabias (lmda,vi,di)
         call assert_real (vi,vd,1.0d-12,
     &                     'emetabiasinterpolate at a node')
         call assert_real (di,dd,1.0d-12,
     &                     'emetabiasinterpolate slope at a node')
      end do
      ostinterpol = .false.
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  function refvstar  --  direct metadynamics bias level  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "refvstar" sums the first "upto" saved metadynamics gaussians
c     directly and returns the smallest bias over the lambda bins
c
c
      function refvstar (upto)
      use dlmda
      use math
      use ost
      implicit none
      integer upto
      integer il,j,m
      real*8 refvstar
      real*8 lambda,delta
      real*8 sig,sig2,pref
      real*8 v
      real*8 src(3)
c
c
      refvstar = 0.0d0
      do il = 1, nlmda
         lambda = dble(il-1) * wlmda
         v = 0.0d0
         do j = 1, upto
            sig = metawhist(j)
            sig2 = sig * sig
            pref = metahhist(j) / (sig*sqrt(2.0d0*pi))
            src(1) = metalhist(j)
            src(2) = -metalhist(j)
            src(3) = 2.0d0 - metalhist(j)
            do m = 1, 3
               delta = lambda - src(m)
               v = v + pref*exp(-0.5d0*delta*delta/sig2)
            end do
         end do
         if (il .eq. 1) then
            refvstar = v
         else
            refvstar = min(refvstar,v)
         end if
      end do
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_eost_ostdyn  --  ost deposit tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_eost_ostdyn" drives eostdyn over two deposit intervals
c     and checks that a settled interval deposits while an unsettled
c     one is rejected
c
c
      subroutine test_eost_ostdyn
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 eostsave
c
c
c     a settled interval deposits one gaussian at the interval end
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      ostcvstd = 1.0d0
      ostcvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.0d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         call eostdyn
         if (istep .lt. lmdaintv) then
            call assert_int (nlmdahist,0,
     &                       'eostdyn waits for the interval end')
         end if
      end do
      call assert_int (nlmdahist,1,
     &                 'eostdyn deposits a settled interval')
      call assert_int (lmdaihist(1),lmdaintv,'eostdyn stamps the step')
      call assert_real (lmdalhist(1),0.5d0,1.0d-12,
     &                  'eostdyn gaussian lambda center')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'eostdyn gaussian flambda center')
      call assert_real (osthhist(1),hbias,1.0d-12,
     &                  'eostdyn untempered gaussian height')
      call assert_real (ostwlhist(1),wlhist,1.0d-12,
     &                  'eostdyn gaussian lambda width')
c
c     an unsettled interval is rejected and changes nothing
c
      eostsave = lmdadeltag
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call eostdyn
      end do
      call assert_real (dedlavg,6.0d0,1.0d-12,
     &                  'eostdyn unsettled interval average')
      call assert_real (dedlstd,5.0d0,1.0d-12,
     &                  'eostdyn unsettled interval deviation')
      call assert_int (nlmdahist,1,'eostdyn rejects an unsettled '//
     &                 'interval')
      call assert_real (lmdadeltag,eostsave,1.0d-12,
     &                  'eostdyn rejection leaves the free energy')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostlocal  --  local tempering test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostlocal" deposits into an unevenly filled kernel
c     with both tempering factors on, and checks that eostdyn takes
c     the height from the pre-deposit bias levels, that the least
c     filled lambda bin deposits at the global height alone, and
c     that growing the flambda grid keeps the bin bias levels
c
c
      subroutine test_eost_ostlocal
      use bath
      use dlmda
      use mutant
      use ost
      use units
      implicit none
      integer i,istep
      integer imax,imin
      real*8 rt,gmin,gl
      real*8 hglobal
      real*8 brutevkmax
      real*8 ostvminimax
      real*8 temperedheight
c
c
c     fill the kernel much higher near lambda of zero
c
      kelvin = 300.0d0
      rt = gasconst * kelvin
      call resetost (5,5,8)
      oststdev = 4.0d0
      nlmdahist = 2
      call sethist (1,0.0d0,0.0d0,5.0d0,0.25d0,1.0d0)
      call sethist (2,0.75d0,0.0d0,1.0d0,0.25d0,1.0d0)
      call buildostindex
      call buildkernels
c
c     settle each deposit interval with both tempering factors on
c
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      ostcvstd = 1.0d0
      ostcvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.0d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      use_ostgtemp = .true.
      use_ostltemp = .true.
      ostgthresh = 0.1d0
      ostlthresh = 0.1d0
      ostgtempgamma = 1.0d0
      ostltempgamma = 1.0d0
c
c     a deposit in the most filled bin is tempered by both factors
c
      imax = 1
      do i = 2, nlmda
         if (brutevkmax(i) .gt. brutevkmax(imax))  imax = i
      end do
      gmin = ostvminimax ()
      gl = vkernelmax(imax)
      call assert_logical (gmin.gt.ostgthresh,.true.,
     &                     'eostdyn global factor active')
      call assert_logical (gl-gmin.gt.ostlthresh,.true.,
     &                     'eostdyn local factor active')
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = dble(imax-1) * wlmda
         dedl = 0.0d0
         call eostdyn
      end do
      hglobal = hbias * exp(-(gmin-ostgthresh)/rt)
      call assert_int (nlmdahist,3,'eostdyn deposits in the full bin')
      call assert_real (osthhist(3),temperedheight(gmin,gl),1.0d-12,
     &                  'eostdyn height from pre-deposit levels')
      call assert_logical (osthhist(3).lt.hglobal,.true.,
     &                     'eostdyn full bin below global height')
c
c     a deposit in the least filled bin sees no local excess
c
      imin = 1
      do i = 2, nlmda
         if (vkernelmax(i) .lt. vkernelmax(imin))  imin = i
      end do
      gmin = ostvminimax ()
      do istep = 1, lmdaintv
         lambda = dble(imin-1) * wlmda
         dedl = 0.0d0
         call eostdyn
      end do
      hglobal = hbias * exp(-max(0.0d0,gmin-ostgthresh)/rt)
      call assert_int (nlmdahist,4,'eostdyn deposits in the least bin')
      call assert_real (osthhist(4),hglobal,1.0d-12,
     &                  'eostdyn least filled bin global height')
c
c     growing the flambda grid keeps the running bin maxima
c
      call ensureflambda (500.0d0)
      do i = 1, nlmda
         call assert_real (vkernelmax(i),brutevkmax(i),1.0d-12,
     &                     'ensureflambda keeps bin bias levels')
      end do
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_temperkeys  --  tempering keywords  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_temperkeys" checks that the tempering keywords set
c     their flags, thresholds and tempering factors, fill a missing
c     value with its default, replace a zero by the default and a
c     negative value by its magnitude
c
c
      subroutine test_eost_temperkeys
      use dlmda
      use keys
      use mutant
      use ost
      implicit none
c
c
c     both keywords set their flags and both values
c
      use_ost = .false.
      use_meta = .false.
      lambda = 0.5d0
      if (allocated(keyline))  deallocate (keyline)
      allocate (keyline(2))
      nkey = 2
      keyline(1) = 'OST-TEMPER-GLOBAL 2.0 3.0'
      keyline(2) = 'ost-temper-local 0.5 0.25'
      call mutate_ost
      call assert_logical (use_ostgtemp,.true.,
     &                     'OST-TEMPER-GLOBAL sets its flag')
      call assert_logical (use_ostltemp,.true.,
     &                     'OST-TEMPER-LOCAL sets its flag')
      call assert_real (ostgthresh,2.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL threshold')
      call assert_real (ostgtempgamma,3.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL tempering factor')
      call assert_real (ostlthresh,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL threshold')
      call assert_real (ostltempgamma,0.25d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL tempering factor')
c
c     a missing tempering factor keeps its default
c
      nkey = 1
      keyline(1) = 'OST-TEMPER-LOCAL 0.5'
      call mutate_ost
      call assert_logical (use_ostgtemp,.false.,
     &                     'OST-TEMPER-LOCAL leaves global off')
      call assert_logical (use_ostltemp,.true.,
     &                     'OST-TEMPER-LOCAL one value sets its flag')
      call assert_real (ostlthresh,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL one value threshold')
      call assert_real (ostltempgamma,1.0d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL default tempering factor')
c
c     zero values take the defaults
c
      keyline(1) = 'OST-TEMPER-GLOBAL 0 0'
      call mutate_ost
      call assert_real (ostgthresh,1.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL zero threshold')
      call assert_real (ostgtempgamma,1.0d0,1.0d-12,
     &                  'OST-TEMPER-GLOBAL zero tempering factor')
c
c     negative values take their magnitude
c
      keyline(1) = 'OST-TEMPER-LOCAL -2.0 -0.5'
      call mutate_ost
      call assert_real (ostlthresh,2.0d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL negative threshold')
      call assert_real (ostltempgamma,0.5d0,1.0d-12,
     &                  'OST-TEMPER-LOCAL negative tempering factor')
c
c     no keywords and the removed keyword leave tempering off
c
      nkey = 0
      call mutate_ost
      call assert_logical (use_ostgtemp .or. use_ostltemp,.false.,
     &                     'tempering off without keywords')
      nkey = 1
      keyline(1) = 'OST-TEMPER'
      call mutate_ost
      call assert_logical (use_ostgtemp,.false.,
     &                     'OST-TEMPER no longer enables tempering')
c
c     restore an empty keyword list and a clean OST state
c
      nkey = 0
      deallocate (keyline)
      call resetost (5,5,1)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostphase  --  phase split tests     ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostphase" checks that setlmdaphase divides the deposit
c     interval into propagation, equilibration and averaging phases,
c     and that the clamps keep a propagation step and enough samples
c     to average without ever losing a sample from the interval
c
c
      subroutine test_eost_ostphase
      use dlmda
      implicit none
c
c
c     the requested ratios divide the interval by truncation
c
      lmdaintv = 10
      lmdaparatio = 0.3d0
      lmdapbratio = 0.3d0
      call setlmdaphase
      call assert_int (lmdanpa,3,'setlmdaphase propagation phase')
      call assert_int (lmdanpb,3,'setlmdaphase equilibration phase')
      call assert_int (lmdanpc,4,'setlmdaphase averaging phase')
      call assert_real (lmdapcratio,0.4d0,1.0d-12,
     &                  'setlmdaphase leftover ratio')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans the whole interval')
c
c     a zero propagation ratio still keeps one propagation step
c
      lmdaintv = 10
      lmdaparatio = 0.0d0
      lmdapbratio = 0.3d0
      call setlmdaphase
      call assert_int (lmdanpa,1,'setlmdaphase keeps one propagation')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans a clamped interval')
c
c     a crowded interval gives back samples to the averaging phase,
c     taking them from the equilibration phase first
c
      lmdaintv = 10
      lmdaparatio = 0.4d0
      lmdapbratio = 0.5d0
      call setlmdaphase
      call assert_int (lmdanpc,2,'setlmdaphase restores the average')
      call assert_int (lmdanpa,4,'setlmdaphase spares the propagation')
      call assert_int (lmdanpb,4,'setlmdaphase trims the equilibration')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans a crowded interval')
c
c     the shortest usable interval still holds all three phases
c
      lmdaintv = 3
      lmdaparatio = 0.4d0
      lmdapbratio = 0.4d0
      call setlmdaphase
      call assert_int (lmdanpa,1,'setlmdaphase minimum propagation')
      call assert_int (lmdanpc,2,'setlmdaphase minimum average')
      call assert_int (lmdanpa+lmdanpb+lmdanpc,lmdaintv,
     &                 'setlmdaphase spans the minimum interval')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eost_ostgate  --  lambda freezing tests  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eost_ostgate" drives eostdyn over one deposit interval and
c     checks that the lambda particle moves only during the leading
c     propagation phase, that lambda is then held exactly fixed, and
c     that the deposited gaussian sits on that frozen lambda
c
c
      subroutine test_eost_ostgate
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      integer istep
      real*8 lam(6)
      real*8 frozen
c
c
c     drive an interval with a deterministic frictionless lambda
c     particle, so that any lambda motion comes from the gate alone
c
      kelvin = 300.0d0
      call resetost (5,5,4)
      lmdaintv = 6
      lmdanpa = 2
      lmdanpb = 2
      lmdanpc = 2
      ostcvstd = 1.0d0
      ostcvrat = 0.0d0
      hbias = 1.0d0
      lmdadt = 0.1d0
      lmdamass = 1.0d0
      lmdafric = 0.0d0
      lmdatheta = 0.25d0 * 3.14159265358979323846d0
      lmdavtheta = 0.0d0
      lambda = 0.5d0
      fastkernel = .true.
      d2edl2 = 0.0d0
      ostbdgdl = 0.0d0
      ostbdgdfl = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         dedl = 1.0d0
         call eostdyn
         lam(istep) = lambda
      end do
c
c     the particle moves while the interval is in its first phase
c
      call assert_logical (lam(1).ne.lam(2),.true.,
     &                     'eostdyn propagates during phase a')
c
c     lambda is then bit identical for the rest of the interval
c
      frozen = lam(2)
      do istep = 3, lmdaintv
         call assert_real (lam(istep),frozen,0.0d0,
     &                     'eostdyn holds lambda after phase a')
      end do
c
c     the averaged lambda is the frozen value, not a smear, so the
c     gaussian is deposited exactly on it
c
      call assert_int (nlmdahist,1,'eostdyn deposits a frozen interval')
      call assert_real (lmdaavg,frozen,0.0d0,
     &                  'eostdyn averages the frozen lambda')
      call assert_real (lmdalhist(1),frozen,0.0d0,
     &                  'eostdyn centers the gaussian on frozen lambda')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'eostdyn centers the gaussian on flambda')
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine resetost  --  reset OST test state  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "resetost" allocates OST arrays and initializes scalar
c     state to deterministic unit-test defaults
c
c
      subroutine resetost (nl,nf,nhist)
      use dlmda
      use mutant
      use ost
      implicit none
      integer nl,nf,nhist
      integer i,j
c
c
c     clear any previous allocation
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(lmdaihist))  deallocate (lmdaihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(lmdaflist))  deallocate (lmdaflist)
      if (allocated(lmdalhist))  deallocate (lmdalhist)
      if (allocated(lmdafhist))  deallocate (lmdafhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      if (allocated(vkernelmax))  deallocate (vkernelmax)
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
c
c     set scalar state
c
      nlmda = nl
      nflmda = nf
      fli0 = (nflmda + 1) / 2
      wlmda = 1.0d0 / dble(nlmda-1)
      wflmda = 1.0d0
      wlmda2 = 0.5d0 * wlmda
      wflmda2 = 0.5d0 * wflmda
      wlhist = 0.005d0
      wfhist = 1.0d0
      maxwlhist = wlhist
      maxwfhist = wfhist
      nlmdahist = 0
      nlmdasave = 0
      sizelmdahist = nhist
      lmdathmap = 'SIN'
      lmdathalpha = 0.999999999d0
      lmdaintv = 10
      lmdanpa = 3
      lmdanpb = 3
      lmdanpc = 4
      lambda = 0.0d0
      lmdaavg = 0.0d0
      lmdastd = 0.0d0
      ostlambdaslp = 0.0d0
      dedl = 0.0d0
      dedlavg = 0.0d0
      dedlstd = 0.0d0
      ostdedlslp = 0.0d0
      deffdl = 0.0d0
      ostcvdif = 0.0d0
      ostcvrat = 0.0d0
      ostcvslp = 0.0d0
      ostcvstd = 0.0d0
      use_ostgtemp = .false.
      use_ostltemp = .false.
      ostgthresh = 0.0d0
      ostgtempgamma = 1.0d0
      ostlthresh = 0.0d0
      ostltempgamma = 1.0d0
      plmdamap = 'QNT'
      elmdamap = 'QNT'
      vlmdamap = 'QNT'
      plmdaexp = 1
      elmdaexp = 1
      vlmdaexp = 1
      plmdainvn = 1
      elmdainvn = 1
      vlmdainvn = 1
      plmdainveps = 0.0d0
      elmdainveps = 0.0d0
      vlmdainveps = 0.0d0
      lmdaparatio = 0.3d0
      lmdapbratio = 0.3d0
      lmdapcratio = 0.4d0
      hbias = 0.0d0
      lmdadeltag = 0.0d0
      oststdev = 1.0d0
      ostinterpol = .false.
      fastkernel = .false.
c
c     allocate arrays
c
      allocate (osthist(sizelmdahist))
      allocate (lmdaihist(sizelmdahist))
      allocate (osthead(nlmda,nflmda))
      allocate (ostnext(sizelmdahist))
      allocate (lmdallist(lmdaintv))
      allocate (lmdaflist(lmdaintv))
      allocate (lmdalhist(sizelmdahist))
      allocate (lmdafhist(sizelmdahist))
      allocate (osthhist(sizelmdahist))
      allocate (ostwlhist(sizelmdahist))
      allocate (ostwfhist(sizelmdahist))
      allocate (lmdafmean(nlmda))
      allocate (lmdafsum(nlmda))
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      allocate (lmdafwt(nlmda))
      allocate (vkernelmax(nlmda))
c
c     initialize arrays
c
      do i = 1, sizelmdahist
         osthist(i) = 0
         lmdaihist(i) = 0
         ostnext(i) = 0
         lmdalhist(i) = 0.0d0
         lmdafhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      do i = 1, lmdaintv
         lmdallist(i) = 0.0d0
         lmdaflist(i) = 0.0d0
      end do
      do i = 1, nlmda
         lmdafmean(i) = 0.0d0
         lmdafsum(i) = 0.0d0
         lmdafwt(i) = 0.0d0
         vkernelmax(i) = 0.0d0
         do j = 1, nflmda
            osthead(i,j) = 0
            gfkernel(i,j) = 0.0d0
            gkernel(i,j) = 0.0d0
            glfkernel(i,j) = 0.0d0
            glkernel(i,j) = 0.0d0
         end do
      end do
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine resetmeta  --  reset meta history  ##
c     ##                                                ##
c     ####################################################
c
c
c     "resetmeta" allocates metadynamics history arrays
c     and clears their contents for a test case
c
c
      subroutine resetmeta (nhist)
      use dlmda
      use ost
      implicit none
      integer nhist
      integer i
c
c
c     clear and allocate metadynamics history arrays
c
      if (allocated(metaihist))  deallocate (metaihist)
      if (allocated(metalhist))  deallocate (metalhist)
      if (allocated(metahhist))  deallocate (metahhist)
      if (allocated(metawhist))  deallocate (metawhist)
      if (allocated(vmetagrid))  deallocate (vmetagrid)
      if (allocated(dvmetagrid))  deallocate (dvmetagrid)
      sizemetahist = nhist
      nmetahist = 0
      allocate (metaihist(sizemetahist))
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      allocate (vmetagrid(nlmda))
      allocate (dvmetagrid(nlmda))
      do i = 1, sizemetahist
         metaihist(i) = 0
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
      end do
      do i = 1, nlmda
         vmetagrid(i) = 0.0d0
         dvmetagrid(i) = 0.0d0
      end do
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine seedkernels  --  seed kernel tags  ##
c     ##                                                ##
c     ####################################################
c
c
c     "seedkernels" fills each flambda-dependent kernel
c     with values that encode its array and bin indices
c
c
      subroutine seedkernels
      use dlmda
      use ost
      implicit none
      integer i,j
      real*8 kerneltag
c
c
      do i = 1, nlmda
         do j = 1, nflmda
            gkernel(i,j) = kerneltag (1,i,j)
            gfkernel(i,j) = kerneltag (2,i,j)
            glkernel(i,j) = kerneltag (3,i,j)
            glfkernel(i,j) = kerneltag (4,i,j)
         end do
      end do
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  function kerneltag  --  encode a kernel bin  ##
c     ##                                               ##
c     ###################################################
c
c
c     "kerneltag" returns a unique numeric tag for a
c     kernel array and lambda/flambda bin pair
c
c
      function kerneltag (ikern,i,j)
      implicit none
      integer ikern,i,j
      real*8 kerneltag
c
c
      kerneltag = 1000.0d0*dble(ikern) + 10.0d0*dble(i) + dble(j)
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine checkkernels  --  check resized kernels  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "checkkernels" compares resized kernels with reference
c     arrays whose seeded block has shifted by offset bins
c
c
      subroutine checkkernels (label,nold,offset)
      use dlmda
      use ost
      implicit none
      integer i,j
      integer nold,offset
      real*8 kerneltag
      real*8, allocatable :: gref(:,:)
      real*8, allocatable :: gfref(:,:)
      real*8, allocatable :: glref(:,:)
      real*8, allocatable :: glfref(:,:)
      character*(*) label
c
c
c     bins outside the copied block must be zero after the resize
c
      allocate (gref(nlmda,nflmda))
      allocate (gfref(nlmda,nflmda))
      allocate (glref(nlmda,nflmda))
      allocate (glfref(nlmda,nflmda))
      do i = 1, nlmda
         do j = 1, nflmda
            gref(i,j) = 0.0d0
            gfref(i,j) = 0.0d0
            glref(i,j) = 0.0d0
            glfref(i,j) = 0.0d0
         end do
      end do
c
c     the seeded block keeps its values but shifts by offset bins
c
      do i = 1, nlmda
         do j = 1, nold
            gref(i,j+offset) = kerneltag (1,i,j)
            gfref(i,j+offset) = kerneltag (2,i,j)
            glref(i,j+offset) = kerneltag (3,i,j)
            glfref(i,j+offset) = kerneltag (4,i,j)
         end do
      end do
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    label//' gkernel')
      call assert_array2 (gfkernel,gfref,nlmda,nflmda,1.0d-12,
     &                    label//' gfkernel')
      call assert_array2 (glkernel,glref,nlmda,nflmda,1.0d-12,
     &                    label//' glkernel')
      call assert_array2 (glfkernel,glfref,nlmda,nflmda,1.0d-12,
     &                    label//' glfkernel')
      deallocate (gref)
      deallocate (gfref)
      deallocate (glref)
      deallocate (glfref)
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine sethist  --  set one history entry  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "sethist" stores one gaussian history entry and its
c     packed lookup bin in the OST arrays
c
c
      subroutine sethist (ihist,lambda,flmda,height,sigl,sigf)
      use dlmda
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda,k
      integer lmdabin,flambdabin
      real*8 lambda,flmda
      real*8 height,sigl,sigf
c
c
c     save real gaussian center, parameters and packed lookup bin
c
      ilmda = lmdabin(lambda)
      iflmda = flambdabin(flmda)
      call ij_to_k (ilmda,iflmda,nlmda,k)
      osthist(ihist) = k
      lmdalhist(ihist) = lambda
      lmdafhist(ihist) = flmda
      osthhist(ihist) = height
      ostwlhist(ihist) = sigl
      ostwfhist(ihist) = sigf
      maxwlhist = max(maxwlhist,sigl)
      maxwfhist = max(maxwfhist,sigf)
      ostnext(ihist) = 0
      return
      end
