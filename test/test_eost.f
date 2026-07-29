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
      call test_eost_efkernel
      call test_eost_save
      call test_eost_meta
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
      logical exist
      character*240 filename0
      character*240 ostfile
      character*240 record
c
c     create a small deterministic history and a temporary restart
c
      filename0 = filename
      leng0 = leng
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
      nosthistsave = 0
      iost = 10
      osttheta = 0.0d0
      ostvtheta = 0.0d0
      ostmass = 1.0d0
      ostfriction = 1.0d0
      ostdt = 0.001d0
      call sethist (1,0.1d0,1.0d0,2.0d0,0.01d0,1.0d0)
      call sethist (2,0.2d0,2.0d0,3.0d0,0.01d0,1.0d0)
      call sethist (3,0.3d0,3.0d0,4.0d0,0.01d0,1.0d0)
c
c     append two histories, then save again without adding history
c
      nosthist = 1
      call saveost
      inquire (file=ostfile,size=size1)
      nosthist = 3
      iost = 30
      call saveost
      inquire (file=ostfile,size=size2)
      iost = 31
      call saveost
      inquire (file=ostfile,size=size3)
      call assert_logical (size2.gt.size1,.true.,
     &                     'saveost appends new history')
      call assert_int (size3,size2,'saveost does not duplicate history')
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
      use ost
      implicit none
      integer i
c
c
c     set small history arrays and resize them
c
      call resetost (5,5,2)
      nosthist = 2
      do i = 1, 2
         osthist(i) = 10 + i
         ostnext(i) = i - 1
         ostlhist(i) = 0.25d0 * dble(i)
         ostfhist(i) = -3.0d0 + 2.0d0*dble(i)
         osthhist(i) = 1.0d0 + dble(i)
         ostwlhist(i) = 0.25d0
         ostwfhist(i) = 1.0d0
      end do
      call resizeosthist
c
c     check old data are preserved and new slots are initialized
c
      call assert_int (sizeosthist,4,
     &                 'resizeosthist size input size=2')
      do i = 1, 4
         if (i .le. 2) then
      call assert_int (osthist(i),10+i,
     &                 'resizeosthist preserve osthist')
      call assert_int (ostnext(i),i-1,
     &                 'resizeosthist preserve ostnext')
      call assert_real (ostlhist(i),0.25d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve ostlhist')
      call assert_real (ostfhist(i),-3.0d0+2.0d0*dble(i),1.0d-12,
     &                  'resizeosthist preserve ostfhist')
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
      call assert_real (ostlhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostlhist')
      call assert_real (ostfhist(i),0.0d0,1.0d-12,
     &                  'resizeosthist init ostfhist')
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
      use ost
      implicit none
      integer k,ilmda,iflmda
c
c
c     put two gaussians in the same bin and one in another bin
c
      call resetost (5,5,3)
      nosthist = 3
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
      nosthist = 2
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
      nosthist = 1
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
      nosthist = 1
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
      nosthist = 1
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
      nosthist = 1
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
      nosthist = 2
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
      ostlambda = 0.75d0
      ostdedl = 1.0d0
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
      nosthist = 2
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.50d0
      ostdedl = 0.0d0
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
      nosthist = 3
      call sethist (1,0.50d0,0.0d0,1.0d0*height,wlmda,wflmda)
      call sethist (2,0.50d0,0.0d0,2.0d0*height,wlmda,wflmda)
      call sethist (3,0.75d0,1.0d0,3.0d0*height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.75d0
      ostdedl = 1.0d0
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
      nosthist = 1
      call sethist (1,0.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.0d0
      ostdedl = 0.0d0
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
      nosthist = 1
      call sethist (1,1.0d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 1.0d0
      ostdedl = 0.0d0
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
      nosthist = 1
      call sethist (1,0.50d0,0.0d0,height,wlhist,wfhist)
      call buildostindex
      call buildgkernel
      expected = exp(-0.25d0)
      call assert_real (gkernel(22,46),expected,1.0d-12,
     &                  'buildgkernel wide hist width')
      ostlambda = targetl
      ostdedl = targetf
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
      nosthist = 1
      call sethist (1,0.50d0,0.0d0,height,wlmda,wflmda)
      call buildostindex
      ostlambda = 0.50d0
      ostdedl = 100.0d0
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
      use math
      use ost
      use units
      implicit none
      real*8 rt
      real*8 expected
      real*8 height
      real*8 w1,w2,w3,w4,w5
c
c     two nonzero gkernel weights at flambda=-1 and +1
c
      call resetost (5,5,1)
      rt = gasconst * kelvin
      gkernel(3,2) = log(2.0d0) * rt
      gkernel(3,4) = log(4.0d0) * rt
      call buildfkernel
      expected = 1.0d0 / 3.0d0
      call assert_real (fkernel(3),expected,1.0d-12,
     &                  'buildfkernel weighted mean force')
      call assert_real (fkernel(1),0.0d0,1.0d-12,
     &                  'buildfkernel empty row')
c
c     build gkernel incrementally from multiple gaussians using
c     updategkernel, then construct the f kernel
c
      call resetost (5,5,4)
      rt = gasconst * kelvin
      height = 2.0d0 * pi * wlmda * wflmda
c
      nosthist = 1
      call sethist (1,0.25d0,-1.0d0,1.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nosthist = 2
      call sethist (2,0.50d0, 0.0d0,2.0d0*height,wlmda,wflmda)
      call buildostindex
      call updategkernel
c
      nosthist = 3
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
      call assert_real (fkernel(3),expected,1.0d-12,
     &                  'buildfkernel from updated gkernel')
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
      nosthist = nhist
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
         fref(i) = fkernel(i)
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
         do j = 1, nflmda
            if (gref(i,j) .ne. 0.0d0) then
               flmda = dble(j-fli0) * wflmda
               weight = exp(gref(i,j)/rt)
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
      call assert_array1 (fkernel,fmanual,nlmda,1.0d-12,
     &                    'buildfkernel all-row weighted mean')
c
c     buildkernels must produce the same gkernel, fkernel and
c     free energy accumulators as the old full rebuild path
c
      do i = 1, nlmda
         fkernel(i) = -123.0d0
         fsumkernel(i) = -123.0d0
         pfkernel(i) = -123.0d0
         do j = 1, nflmda
            gkernel(i,j) = -123.0d0
         end do
      end do
      call buildkernels
      call assert_array2 (gkernel,gref,nlmda,nflmda,1.0d-12,
     &                    'buildkernels gkernel reference')
      call assert_array1 (fkernel,fref,nlmda,1.0d-12,
     &                    'buildkernels fkernel reference')
      call assert_array1 (fsumkernel,fsumref,nlmda,1.0d-12,
     &                    'buildkernels fsum reference')
      call assert_array1 (pfkernel,pfref,nlmda,1.0d-12,
     &                    'buildkernels partfunc reference')
c
c     incrementally updating kernels one history at a time must match
c     a full buildkernels rebuild
c
      call resetost (9,9,8)
      height = 2.0d0 * pi * wlmda * wflmda
      do ihist = 1, nhist
         nosthist = ihist
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
      call assert_array1 (fkernel,fref,nlmda,1.0d-12,
     &                    'updatekernels fkernel reference')
      call assert_array1 (fsumkernel,fsumref,nlmda,1.0d-12,
     &                    'updatekernels fsum reference')
      call assert_array1 (pfkernel,pfref,nlmda,1.0d-12,
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
      use ost
      implicit none
      integer i
      real*8 stdref
c
c
c     average all saved interval samples after equilibration prefix
c
      call resetost (5,5,1)
      iosthist = 6
      ostnequil = 2
      ostnavg = 4
      do i = 1, iosthist
         ostllist(i) = dble(i)
         ostflist(i) = 2.0d0*dble(i)
      end do
      call avgstd (ostllist,ostnequil+1,ostnavg,
     &             ostlambdaavg,ostlambdastd)
      call avgstd (ostflist,ostnequil+1,ostnavg,
     &             ostdedlavg,ostdedlstd)
      stdref = sqrt(1.25d0)
      call assert_real (ostlambdaavg,4.5d0,1.0d-12,
     &                  'avgstd configurable lambda average')
      call assert_real (ostdedlavg,9.0d0,1.0d-12,
     &                  'avgstd configurable dE/dl average')
      call assert_real (ostlambdastd,stdref,1.0d-12,
     &                  'avgstd configurable lambda std')
      call assert_real (ostdedlstd,2.0d0*stdref,1.0d-12,
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
      use math
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
      nosthist = 3
      call sethist (1,0.25d0,-1.0d0,1.1d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.6d0*height,sigl,sigf)
      call sethist (3,0.75d0, 1.0d0,2.3d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      ostlambda = 0.50d0
      ostdedl = 0.0d0
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
      nosthist = 3
      call sethist (1,0.25d0,-2.0d0,0.8d0*height,sigl,sigf)
      call sethist (2,0.50d0, 0.0d0,1.2d0*height,sigl,sigf)
      call sethist (3,0.75d0, 2.0d0,1.6d0*height,sigl,sigf)
      call buildostindex
      call buildkernels
      ostlambda = 0.53125d0
      ostdedl = 0.25d0
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
      use ost
      implicit none
      integer i
      real*8 eostlmda,dfdl
      real*8 expected
      real*8 etotfkernel
c
c
c     use fkernel(lambda)=lambda and integrate to lambda=0.375
c
      call resetost (5,5,1)
      do i = 1, nlmda
         fkernel(i) = dble(i-1) * wlmda
      end do
      ostlambda = 0.375d0
      call efkernel (eostlmda,dfdl)
      expected = 0.5d0 * ostlambda * ostlambda
      call assert_real (eostlmda,expected,1.0d-12,
     &                  'efkernel DeltaG lambda=.375')
      call assert_real (dfdl,ostlambda,1.0d-12,
     &                  'efkernel dDeltaG/dlambda')
c
c     use fkernel(lambda)=1+lambda so the endpoint mean forces are
c     nonzero and distinct from each other
c
      call resetost (5,5,1)
      do i = 1, nlmda
         fkernel(i) = 1.0d0 + dble(i-1)*wlmda
      end do
c
c     at and below lambda = 0 the free energy is zero and the
c     derivative comes from the first lambda bin
c
      ostlambda = 0.0d0
      call efkernel (eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efkernel DeltaG lambda=0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efkernel dDeltaG/dlambda lambda=0')
      ostlambda = -0.25d0
      call efkernel (eostlmda,dfdl)
      call assert_real (eostlmda,0.0d0,1.0d-12,
     &                  'efkernel DeltaG lambda below 0')
      call assert_real (dfdl,1.0d0,1.0d-12,
     &                  'efkernel dDeltaG/dlambda below 0')
c
c     at lambda = 1 the last interval is integrated in full, and
c     beyond lambda = 1 the loop falls through to the same result
c
      ostlambda = 1.0d0
      call efkernel (eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efkernel DeltaG lambda=1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efkernel dDeltaG/dlambda lambda=1')
      ostlambda = 1.25d0
      call efkernel (eostlmda,dfdl)
      call assert_real (eostlmda,1.5d0,1.0d-12,
     &                  'efkernel DeltaG lambda above 1')
      call assert_real (dfdl,2.0d0,1.0d-12,
     &                  'efkernel dDeltaG/dlambda above 1')
c
c     the trapezoid rule is exact for a linear mean force, so the
c     total must match the efkernel value at lambda = 1
c
      call assert_real (etotfkernel(),1.5d0,1.0d-12,
     &                  'etotfkernel linear fkernel')
c
c     a constant mean force integrates to itself over unit lambda
c
      do i = 1, nlmda
         fkernel(i) = 2.5d0
      end do
      call assert_real (etotfkernel(),2.5d0,1.0d-12,
     &                  'etotfkernel constant fkernel')
c
c     an all-zero mean force gives no free energy change
c
      do i = 1, nlmda
         fkernel(i) = 0.0d0
      end do
      call assert_real (etotfkernel(),0.0d0,1.0d-12,
     &                  'etotfkernel zero fkernel')
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
      call emetabias (0.5d0,vbias,dvdl)
      call assert_real (vbias,pref,1.0d-12,
     &                  'emetabias center value')
      call assert_real (dvdl,0.0d0,1.0d-12,
     &                  'emetabias center derivative')
      call emetabias (0.75d0,vbias,dvdl)
      expected = pref * exp(-0.5d0)
      call assert_real (vbias,expected,1.0d-12,
     &                  'emetabias off-center value')
      call assert_real (dvdl,-4.0d0*expected,1.0d-12,
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
      if (allocated(ostihist))  deallocate (ostihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(ostllist))  deallocate (ostllist)
      if (allocated(ostflist))  deallocate (ostflist)
      if (allocated(ostlhist))  deallocate (ostlhist)
      if (allocated(ostfhist))  deallocate (ostfhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(fkernel))  deallocate (fkernel)
      if (allocated(fsumkernel))  deallocate (fsumkernel)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(pfkernel))  deallocate (pfkernel)
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
      nosthist = 0
      nosthistsave = 0
      sizeosthist = nhist
      iosthist = 10
      ostnequil = 5
      ostnavg = 5
      ostlambda = 0.0d0
      ostlambdaavg = 0.0d0
      ostlambdastd = 0.0d0
      ostdedl = 0.0d0
      ostdedlavg = 0.0d0
      ostdedlstd = 0.0d0
      deffdl = 0.0d0
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
      osteqratio = 0.5d0
      hbias = 0.0d0
      eosttot = 0.0d0
      oststdev = 1.0d0
      ostinterpol = .false.
c
c     allocate arrays
c
      allocate (osthist(sizeosthist))
      allocate (ostihist(sizeosthist))
      allocate (osthead(nlmda,nflmda))
      allocate (ostnext(sizeosthist))
      allocate (ostllist(iosthist))
      allocate (ostflist(iosthist))
      allocate (ostlhist(sizeosthist))
      allocate (ostfhist(sizeosthist))
      allocate (osthhist(sizeosthist))
      allocate (ostwlhist(sizeosthist))
      allocate (ostwfhist(sizeosthist))
      allocate (fkernel(nlmda))
      allocate (fsumkernel(nlmda))
      allocate (gfkernel(nlmda,nflmda))
      allocate (gkernel(nlmda,nflmda))
      allocate (glfkernel(nlmda,nflmda))
      allocate (glkernel(nlmda,nflmda))
      allocate (pfkernel(nlmda))
c
c     initialize arrays
c
      do i = 1, sizeosthist
         osthist(i) = 0
         ostihist(i) = 0
         ostnext(i) = 0
         ostlhist(i) = 0.0d0
         ostfhist(i) = 0.0d0
         osthhist(i) = 0.0d0
         ostwlhist(i) = 0.0d0
         ostwfhist(i) = 0.0d0
      end do
      do i = 1, iosthist
         ostllist(i) = 0.0d0
         ostflist(i) = 0.0d0
      end do
      do i = 1, nlmda
         fkernel(i) = 0.0d0
         fsumkernel(i) = 0.0d0
         pfkernel(i) = 0.0d0
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
      sizemetahist = nhist
      nmetahist = 0
      allocate (metaihist(sizemetahist))
      allocate (metalhist(sizemetahist))
      allocate (metahhist(sizemetahist))
      allocate (metawhist(sizemetahist))
      do i = 1, sizemetahist
         metaihist(i) = 0
         metalhist(i) = 0.0d0
         metahhist(i) = 0.0d0
         metawhist(i) = 0.0d0
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
      use ost
      implicit none
      integer ihist
      integer ilmda,iflmda,k
      integer lambdabin,flambdabin
      real*8 lambda,flmda
      real*8 height,sigl,sigf
c
c
c     save real gaussian center, parameters and packed lookup bin
c
      ilmda = lambdabin(lambda)
      iflmda = flambdabin(flmda)
      call ij_to_k (ilmda,iflmda,nlmda,k)
      osthist(ihist) = k
      ostlhist(ihist) = lambda
      ostfhist(ihist) = flmda
      osthhist(ihist) = height
      ostwlhist(ihist) = sigl
      ostwfhist(ihist) = sigf
      maxwlhist = max(maxwlhist,sigl)
      maxwfhist = max(maxwfhist,sigf)
      ostnext(ihist) = 0
      return
      end
