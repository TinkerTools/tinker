c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdcd  --  output of DCD-format coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdcd" writes out a set of Cartesian coordinates to a
c     coordinates file in CHARMM DCD binary format compatible
c     with the VMD visualization software and other packages
c
c     note the format used is based on the "dcdplugin.c" code from
c     the NAMD and VMD programs, and tutorial 4.1 from the software
c     package GENESIS: Generalized-Ensemble Simulation System
c
c     variables and parameters:
c
c     header     type of data (CORD=coordinates, VELD=velocities)
c     nframe     number of frames stored in the DCD file
c     nprev      number of previous integration steps
c     ncrdsav    frequency in steps for saving coordinate frames
c     nstep      number of integration steps in the total run
c     nvelsav    frequency of coordinate saves with velocity data
c     ndfree     number of degrees of freedom for the system
c     nfixat     number of fixed atoms for the system
c     usebox     flag for periodic boundaries (1=true, 0=false)
c     use4d      flag for 4D trajectory (1=true, 0=false)
c     usefq      flag for fluctuating charges (1=true, 0=false)
c     merged     result of merge without checks (1=true, 0=false)
c     vcharmm    version of CHARMM software for compatibility
c
c     in general a value of zero for any of the above indicates that
c     the particular feature is unused
c
c
      subroutine prtdcd (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use titles
      implicit none
      integer i,idcd
      integer zero,one
      integer nframe,nprev
      integer ncrdsav,nstep
      integer nvelsav,ndfree
      integer nfixat,usebox
      integer use4d,usefq
      integer merged,vcharmm
      integer ntitle
      real*4 tdelta
      logical opened,first
      character*4 header
      character*240 dcdfile
c
c
c     open the output unit if not already done
c
      inquire (unit=idcd,opened=opened)
      if (.not. opened) then
         dcdfile = filename(1:leng)//'.dcd'
         call version (dcdfile,'new')
         open (unit=idcd,file=dcdfile,form='unformatted',status='new')
      end if
c
c     write header info along with title and number of atoms
c
      if (first) then
         first = .false.
         zero = 0
         one = 1
         header = 'CORD'
         nframe = zero
         nprev = zero
         ncrdsav = zero
         nstep = zero
         nvelsav = zero
         ndfree = zero
         nfixat = zero
         tdelta = 0.0
         usebox = zero
         if (use_bounds)  usebox = one
         use4d = zero
         usefq = zero
         merged = zero
         vcharmm = 24
         ntitle = one
         write (idcd)  header,nframe,nprev,ncrdsav,nstep,
     &                 nvelsav,zero,zero,ndfree,nfixat,
     &                 tdelta,usebox,use4d,usefq,merged,
     &                 zero,zero,zero,zero,zero,vcharmm
         write (idcd)  ntitle,title(1:80)
         write (idcd)  n
      end if
c
c     append the lattice values based on header flag value
c
      if (use_bounds) then
         write (idcd)  xbox,gamma_cos,ybox,beta_cos,alpha_cos,zbox
      end if
c
c     append the atomic coordinates along each axis in turn
c
      write (idcd)  (real(x(i)),i=1,n)
      write (idcd)  (real(y(i)),i=1,n)
      write (idcd)  (real(z(i)),i=1,n)
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
