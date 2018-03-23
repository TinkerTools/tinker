c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2018  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine kopenmm  --  assign OpenMM parameter values  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "kopenmm" assigns OpenMM parameters and options to be passed
c     from the Tinker interface into the OpenMM library
c
c
      subroutine kopenmm
      use keys
      use openmm
      implicit none
      integer i,next
      logical found
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     initialize the OpenMM platform, CUDA device and precision
c
      ommPlatform = 'CUDA     '
      cudaDevice = '                '
      cudaPrecision = 'MIXED '
c
c     get the CUDA device setting from an environmental variable
c
      found = .false.
      call getenv ('cuda_device',record)
      if (record(1:1) .ne. ' ')  found = .true.
      if (.not. found) then
         call getenv ('CUDA_DEVICE',record)
         if (record(1:1) .ne. ' ')  found = .true.
      end if
      if (.not. found) then
         cudaDevice = '0               '
      else
         cudaDevice = record(1:16)
      end if
c
c     check for keywords containing any altered parameters
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:16) .eq. 'OPENMM-PLATFORM ') then
            call gettext (record,ommPlatform,next)
            call upcase (ommPlatform)
         else if (keyword(1:12) .eq. 'CUDA-DEVICE ') then
            call gettext (record,cudaDevice,next)
            call upcase (cudaDevice)
         else if (keyword(1:15) .eq. 'CUDA-PRECISION ') then
            call gettext (record,cudaPrecision,next)
            call upcase (cudaPrecision)
         end if
      end do
      return
      end
