c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf  --  adaptive biasing force tests  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf" checks the adaptive biasing force bias, the interval
c     sampling built on the ost history, and the abf history restart
c
c
      subroutine test_eabf
      use bath
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_eabf')
c
c
      if (skiptest(tname,'abf'))  return
      call initial
      kelvin = 300.0d0
      call test_eabf_bias
      call test_eabf_dyn
      call test_eabf_mean
      call test_eabf_gate
      call test_eabf_resize
      call test_eabf_restart
      call test_eabf_orphan
      call final
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eabf_bias  --  lambda-only bias test  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eabf_bias" checks that eabfbias removes the integrated
c     mean force from the energy, saves its slope for the lambda
c     particle, and leaves the virial untouched
c
c
      subroutine test_eabf_bias
      use dlmda
      use energi
      use mutant
      use virial
      implicit none
      integer i,j
c
c
c     seed a linear mean force f(lambda) = 2 + 4*lambda on the grid
c
      call resetabf (5,4)
      use_abf = .true.
      do i = 1, nlmda
         lmdafmean(i) = 2.0d0 + 4.0d0*dble(i-1)*wlmda
      end do
      do i = 1, 3
         do j = 1, 3
            vir(j,i) = 1.0d0
         end do
      end do
      lambda = 0.3d0
      esum = 1.0d0
      call eabfbias
c
c     the free energy 2*lambda + 2*lambda**2 leaves the energy and
c     its slope 2 + 4*lambda is saved for the lambda particle
c
      call assert_real (esum,1.0d0-0.78d0,1.0d-12,
     &                  'eabfbias removes the free energy')
      call assert_real (lmdavbias,-0.78d0,1.0d-12,
     &                  'eabfbias saves the bias energy')
      call assert_real (lmdadfdl,3.2d0,1.0d-12,
     &                  'eabfbias saves the mean force')
      call assert_real (vir(2,1),1.0d0,0.0d0,
     &                  'eabfbias leaves the virial')
      use_abf = .false.
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eabf_dyn  --  interval sample tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eabf_dyn" drives eabfdyn over two sample intervals and
c     checks the sample, bin accumulators and free energy of
c     a settled interval, and that a noisy interval is still kept
c
c
      subroutine test_eabf_dyn
      use dlmda
      use mutant
      implicit none
      integer istep
c
c
c     a settled interval records its lambda and dU/dlambda average
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.5d0
      lmdastep = 0
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 3.0d0
         call eabfdyn
         if (istep .lt. lmdaintv) then
            call assert_int (nlmdahist,0,
     &                       'eabfdyn waits for the interval end')
         end if
      end do
      call assert_real (deffdl,2.5d0,1.0d-12,
     &                  'eabfdyn removes the mean force')
      call assert_int (nlmdahist,1,'eabfdyn records one interval')
      call assert_int (lmdaihist(1),lmdaintv,'eabfdyn stamps the step')
      call assert_real (lmdalhist(1),0.5d0,0.0d0,
     &                  'eabfdyn sample lambda')
      call assert_real (lmdafhist(1),3.0d0,1.0d-12,
     &                  'eabfdyn sample dU/dlambda')
      call assert_real (lmdafwt(3),1.0d0,0.0d0,'eabfdyn bin count')
      call assert_real (lmdafmean(3),3.0d0,1.0d-12,'eabfdyn mean force')
      call assert_real (lmdadeltag,0.75d0,1.0d-12,
     &                  'eabfdyn free energy estimate')
c
c     a noisy interval is still recorded, since abf has no gate
c
      do istep = 1, lmdaintv
         lambda = 0.5d0
         dedl = 1.0d0
         if (mod(istep,2) .eq. 0)  dedl = 11.0d0
         call eabfdyn
      end do
      call assert_int (nlmdahist,2,'eabfdyn keeps a noisy interval')
      call assert_real (lmdafwt(3),2.0d0,0.0d0,
     &                  'eabfdyn noisy bin count')
      call assert_real (lmdafmean(3),4.5d0,1.0d-12,
     &                  'eabfdyn noisy mean force')
      use_abf = .false.
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_eabf_mean  --  bin mean force tests  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_eabf_mean" checks that each sample sets the bin bias force
c     to the running mean of its bin, and that a rebuild from the
c     history reproduces the accumulated bins
c
c
      subroutine test_eabf_mean
      use dlmda
      implicit none
      integer i
      real*8 ref
c
c
c     the bias force follows the running mean of the bin samples
c
      call resetabf (5,8)
      do i = 1, 6
         lmdalhist(i) = 0.25d0
         lmdafhist(i) = dble(i)
         nlmdahist = i
         call addabfhist (i)
         ref = 0.5d0 * dble(i+1)
         call assert_real (lmdafmean(2),ref,1.0d-12,
     &                     'addabfhist sets the running mean force')
      end do
c
c     a rebuild from the history reproduces the accumulated bins
c
      call buildabfkernel
      call assert_real (lmdafwt(2),6.0d0,0.0d0,
     &                  'buildabfkernel bin count')
      call assert_real (lmdafsum(2),21.0d0,1.0d-12,
     &                  'buildabfkernel bin sum')
      call assert_real (lmdafmean(2),3.5d0,1.0d-12,
     &                  'buildabfkernel bin mean force')
      call assert_real (lmdafmean(1),0.0d0,0.0d0,
     &                  'buildabfkernel empty bin')
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_eabf_gate  --  lambda freezing tests  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_eabf_gate" drives eabfdyn over one sample interval and
c     checks that the lambda particle moves only during the leading
c     propagation phase and that the sample sits on the frozen lambda
c
c
      subroutine test_eabf_gate
      use dlmda
      use mutant
      implicit none
      integer istep
      real*8 lam(6)
      real*8 frozen
c
c
c     drive an interval with a deterministic frictionless lambda
c     particle, so that any lambda motion comes from the gate alone
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 6
      lmdanpa = 2
      lmdanpb = 2
      lmdanpc = 2
      lmdadt = 0.1d0
      lmdamass = 1.0d0
      lmdafric = 0.0d0
      lmdatheta = 0.25d0 * 3.14159265358979323846d0
      lmdavtheta = 0.0d0
      lambda = 0.5d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, lmdaintv
         dedl = 1.0d0
         call eabfdyn
         lam(istep) = lambda
      end do
c
c     the particle moves only while the interval is in phase a
c
      call assert_logical (lam(1).ne.lam(2),.true.,
     &                     'eabfdyn propagates during phase a')
      frozen = lam(2)
      do istep = 3, lmdaintv
         call assert_real (lam(istep),frozen,0.0d0,
     &                     'eabfdyn holds lambda after phase a')
      end do
c
c     the sample is recorded exactly on the frozen lambda
c
      call assert_int (nlmdahist,1,'eabfdyn samples a frozen interval')
      call assert_real (lmdalhist(1),frozen,0.0d0,
     &                  'eabfdyn sample on the frozen lambda')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'eabfdyn sample dU/dlambda on the gate')
      use_abf = .false.
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf_resize  --  sample history growth  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf_resize" records more interval samples than the abf
c     history holds and checks that the history doubles while every
c     sample and the bin mean force are kept
c
c
      subroutine test_eabf_resize
      use dlmda
      use mutant
      implicit none
      integer istep
c
c
c     three samples overflow a history sized for two
c
      call resetabf (5,2)
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      do istep = 1, 3*lmdaintv
         lambda = 0.25d0
         dedl = dble((istep-1)/lmdaintv + 1)
         call eabfdyn
      end do
      call assert_int (nlmdahist,3,'resizeabfhist sample count')
      call assert_int (sizelmdahist,4,'resizeabfhist doubles the size')
      call assert_int (lmdaihist(1),4,'resizeabfhist keeps step 1')
      call assert_int (lmdaihist(3),12,'resizeabfhist keeps step 3')
      call assert_real (lmdafhist(1),1.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 1')
      call assert_real (lmdafhist(2),2.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 2')
      call assert_real (lmdafhist(3),3.0d0,1.0d-12,
     &                  'resizeabfhist keeps sample 3')
      call assert_real (lmdalhist(3),0.25d0,0.0d0,
     &                  'resizeabfhist keeps the sample lambda')
      call assert_real (lmdafmean(2),2.0d0,1.0d-12,
     &                  'resizeabfhist bin mean force')
      use_abf = .false.
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_eabf_restart  --  history restart test  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_eabf_restart" writes an abf history, checks that repeated
c     saves append only new samples, then continues the history in a
c     fresh run and checks the restored bins, lambda particle and
c     interval boundary, with the run settings winning over the file
c
c
      subroutine test_eabf_restart
      use bath
      use dlmda
      use files
      use mutant
      use ost
      implicit none
      integer i,ihis
      integer iline
      integer istep
      integer leng0
      integer size1,size2,size3
      integer freeunit
      integer trimtext
      real*8 kelvin0
      real*8 eostref
      real*8 fref(5)
      real*8 sref(5)
      real*8 pref(5)
      logical exist
      logical hastitle
      logical haslabel
      character*240 filename0
      character*240 abffile
      character*240 record
c
c
c     start an abf history in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-saveabf'
      leng = len_trim(filename)
      abffile = filename(1:leng)//'.abf'
      inquire (file=abffile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=abffile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetabf (5,4)
      use_ost = .false.
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdamass = 2.0d0
      lmdafric = 0.5d0
      lmdatheta = 0.3d0
      lmdavtheta = 0.1d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      call initabffile
      inquire (file=abffile,exist=exist)
      call assert_logical (exist,.true.,'initabffile creates the file')
c
c     the history title marks the file as an abf history
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      read (ihis,10)  record
   10 format (a240)
      hastitle = (index(record,abftitle(1:trimtext(abftitle))) .gt. 0)
      do iline = 2, 8
         read (ihis,10)  record
      end do
      close (unit=ihis)
      haslabel = (index(record,abflabel(1:trimtext(abflabel))) .gt. 0)
      call assert_logical (hastitle,.true.,
     &                     'initabffile writes the abf title')
      call assert_logical (haslabel,.true.,
     &                     'initabffile writes the sample label')
c
c     record two intervals, then save twice with no new samples
c
      do istep = 1, 2*lmdaintv
         lambda = 0.25d0
         dedl = 1.0d0
         if (istep .gt. lmdaintv) then
            lambda = 0.75d0
            dedl = 5.0d0
         end if
         call eabfdyn
      end do
      call saveabf
      inquire (file=abffile,size=size1)
      call saveabf
      inquire (file=abffile,size=size2)
      call assert_int (size2,size1,'saveabf does not duplicate')
c
c     a third interval and part of a fourth, then save again
c
      do istep = 1, lmdaintv+2
         lambda = 0.75d0
         dedl = 7.0d0
         call eabfdyn
      end do
      call saveabf
      inquire (file=abffile,size=size3)
      call assert_logical (size3.gt.size2,.true.,
     &                     'saveabf appends new samples')
      call assert_int (nlmdahist,3,'eabfdyn restart sample count')
      call assert_int (lmdastep,14,'eabfdyn restart partial step')
      do i = 1, 5
         fref(i) = lmdafmean(i)
         sref(i) = lmdafsum(i)
         pref(i) = lmdafwt(i)
      end do
      eostref = lmdadeltag
c
c     a fresh run with matching grid and interval continues the file
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      kelvin = 310.0d0
      lmdamass = 9.0d0
      lmdafric = 0.2d0
      lmdadt = 0.002d0
      lmdatheta = 0.0d0
      lmdavtheta = 0.0d0
      call initabffile
      call assert_logical (lmdasavefile.eq.abffile,.true.,
     &                     'initabffile continues the same file')
      call assert_int (nlmdahist,3,'initabffile restores the samples')
      call assert_int (nlmdasave,3,'initabffile restores the save')
      call assert_int (lmdastep,12,'initabffile resumes on an interval')
      call assert_array1 (lmdafmean,fref,5,1.0d-12,
     &                    'initabffile restores the mean force')
      call assert_array1 (lmdafsum,sref,5,1.0d-12,
     &                    'initabffile restores the bin sums')
      call assert_array1 (lmdafwt,pref,5,0.0d0,
     &                    'initabffile restores the bin counts')
      call assert_real (lmdadeltag,eostref,1.0d-12,
     &                  'initabffile restores the free energy')
      call assert_real (lmdatheta,0.3d0,1.0d-14,
     &                  'initabffile restores the lambda theta')
      call assert_real (lmdavtheta,0.1d0,1.0d-14,
     &                  'initabffile restores the theta velocity')
      call assert_real (kelvin,310.0d0,0.0d0,
     &                  'initabffile keeps the run temperature')
      call assert_real (lmdamass,9.0d0,0.0d0,
     &                  'initabffile keeps the run lambda mass')
      call assert_real (lmdadt,0.002d0,0.0d0,
     &                  'initabffile keeps the run lambda step')
      call assert_logical (allocated(osthist),.false.,
     &                     'rdabf leaves the ost history unallocated')
      call assert_logical (allocated(gkernel),.false.,
     &                     'rdabf leaves the g kernel unallocated')
c
c     remove the temporary file and restore the global state
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      close (unit=ihis,status='delete')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_abf = .false.
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_eabf_orphan  --  stray row restart test  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_eabf_orphan" appends a sample row past the header count,
c     as left by a run stopped between appending samples and updating
c     the header, and checks that continuing the history drops it
c
c
      subroutine test_eabf_orphan
      use bath
      use dlmda
      use files
      use mutant
      implicit none
      integer ihis
      integer istep
      integer leng0
      integer freeunit
      real*8 kelvin0
      logical exist
      character*240 filename0
      character*240 abffile
c
c
c     start an abf history with two intervals in a temporary file
c
      filename0 = filename
      leng0 = leng
      kelvin0 = kelvin
      kelvin = 300.0d0
      filename = 'tinkertest-orphanabf'
      leng = len_trim(filename)
      abffile = filename(1:leng)//'.abf'
      inquire (file=abffile,exist=exist)
      if (exist) then
         ihis = freeunit ()
         open (unit=ihis,file=abffile,status='old')
         close (unit=ihis,status='delete')
      end if
      call resetabf (5,4)
      use_ost = .false.
      use_abf = .true.
      lmdaintv = 4
      lmdanpa = 0
      lmdanpb = 0
      lmdanpc = 4
      lmdadt = 0.0d0
      lmdadfdl = 0.0d0
      lmdastep = 0
      call initabffile
      do istep = 1, 2*lmdaintv
         lambda = 0.25d0
         dedl = 1.0d0
         call eabfdyn
      end do
      call saveabf
c
c     append a stray row past the header count
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old',position='append')
      write (ihis,10)  999,0.9d0,1.0d6
   10 format (i12,2d26.16)
      close (unit=ihis)
c
c     continuing the history drops the stray row before new samples
c
      call resetabf (5,4)
      use_abf = .true.
      lmdaintv = 4
      kelvin = 300.0d0
      call initabffile
      call assert_int (nlmdahist,2,'initabffile skips the stray row')
      do istep = 1, lmdaintv
         lambda = 0.75d0
         dedl = 5.0d0
         call eabfdyn
      end do
      call saveabf
c
c     a fresh read holds the kept and new samples and no stray bin
c
      call resetabf (5,4)
      use_abf = .true.
      lmdasavefile = abffile
      call rdabf
      call assert_int (nlmdahist,3,'rdabf reads the kept samples')
      call assert_int (lmdaihist(3),12,'rdabf reads the new step')
      call assert_real (lmdalhist(3),0.75d0,0.0d0,
     &                  'rdabf reads the new sample lambda')
      call assert_real (lmdafhist(3),5.0d0,1.0d-12,
     &                  'rdabf reads the new sample average')
      call assert_real (lmdafwt(5),0.0d0,0.0d0,
     &                  'rdabf builds no bin from the stray row')
c
c     remove the temporary file and restore the global state
c
      ihis = freeunit ()
      open (unit=ihis,file=abffile,status='old')
      close (unit=ihis,status='delete')
      filename = filename0
      leng = leng0
      kelvin = kelvin0
      use_abf = .false.
      return
      end
c
c
c     #####################################################
c     ##                                                 ##
c     ##  subroutine resetabf  --  reset ABF test state  ##
c     ##                                                 ##
c     #####################################################
c
c
c     "resetabf" sets up the ost test state and then frees every array
c     that the adaptive biasing force method must never touch, so any
c     stray access from an abf routine fails the test
c
c
      subroutine resetabf (nl,nhist)
      use ost
      implicit none
      integer nl,nhist
c
c
c     start from the full deterministic ost test state
c
      call resetost (nl,5,nhist)
c
c     free the arrays used only by ost and metadynamics
c
      deallocate (osthist)
      deallocate (ostnext)
      deallocate (osthead)
      deallocate (osthhist)
      deallocate (ostwlhist)
      deallocate (ostwfhist)
      deallocate (gfkernel)
      deallocate (gkernel)
      deallocate (glfkernel)
      deallocate (glkernel)
      deallocate (vkernelmax)
      deallocate (ostlmdaavgbin)
      deallocate (ostlmdaslpbin)
      deallocate (ostlmdastdbin)
      deallocate (ostdedlavgbin)
      deallocate (ostdedlslpbin)
      deallocate (ostdedlstdbin)
      return
      end
