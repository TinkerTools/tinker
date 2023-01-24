c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtarc  --  output of a coordinates archive  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtarc" writes a Cartesian coordinates archive as either
c     a formatted or binary disk file
c
c
      subroutine prtarc (iarc,first)
      use output
      implicit none
      integer iarc
      logical first
c
c
c     write archive file as either formatted or binary file
c
      if (archive) then
         call prtarcf (iarc)
      else if (binary) then
         call prtarcb (iarc,first)
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prtarcf  --  output of Tinker archive file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtarcf" writes out a set of Cartesian coordinates for
c     all active atoms in the Tinker XYZ archive format
c
c
      subroutine prtarcf (iarc)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use titles
      use usage
      implicit none
      integer i,j,k,iarc
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 arcfile
c
c
c     open output unit if not already done
c
      inquire (unit=iarc,opened=opened)
      if (.not. opened) then
         arcfile = filename(1:leng)//'.arc'
         call version (arcfile,'new')
         open (unit=iarc,file=arcfile,status='new')
      end if
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,x(i),y(i),z(i))
         crdmax = max(crdmax,x(i),y(i),z(i))
      end do
      crdsiz = 6
      if (crdmin .le. -1000.0d0)  crdsiz = 7
      if (crdmax .ge. 10000.0d0)  crdsiz = 7
      if (crdmin .le. -10000.0d0)  crdsiz = 8
      if (crdmax .ge. 100000.0d0)  crdsiz = 8
      crdsiz = crdsiz + max(6,digits)
      size = 0
      call numeral (crdsiz,crdc,size)
      if (digits .le. 6) then
         digc = '6 '
      else if (digits .le. 8) then
         digc = '8'
      else
         digc = '10'
      end if
c
c     write out the number of atoms and the title
c
      if (ltitle .eq. 0) then
         fstr = '('//atmc//')'
         write (iarc,fstr(1:4))  nuse
      else
         fstr = '('//atmc//',2x,a)'
         write (iarc,fstr(1:9))  nuse,title(1:ltitle)
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (iarc,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the coordinate line for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      do i = 1, n
         if (use(i)) then
            k = n12(i)
            if (k .eq. 0) then
               write (iarc,fstr)  iuse(i),name(i),x(i),y(i),z(i),
     &                            type(i)
            else
               write (iarc,fstr)  iuse(i),name(i),x(i),y(i),z(i),
     &                            type(i),(iuse(i12(j,i)),j=1,k)
            end if
         end if
      end do
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iarc)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine prtarcb  --  output binary-format archive file  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtarcb" writes out a set of Cartesian coordinates for all
c     active atoms in the CHARMM DCD binary format
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
      subroutine prtarcb (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use titles
      use usage
      implicit none
      integer i,k,idcd
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
         write (idcd)  nuse
      end if
c
c     append the lattice values based on header flag value
c
      if (use_bounds) then
         write (idcd)  xbox,gamma_cos,ybox,beta_cos,alpha_cos,zbox
      end if
c
c     remove unused atoms from the coordinates to be output
c
      if (nuse .ne. n) then
         k = 0
         do i = 1, n
            if (use(i)) then
               k = k + 1
               x(k) = x(i)
               y(k) = y(i)
               z(k) = z(i)
            end if
         end do
      end if
c
c     append the atomic coordinates along each axis in turn
c
      write (idcd)  (real(x(i)),i=1,nuse)
      write (idcd)  (real(y(i)),i=1,nuse)
      write (idcd)  (real(z(i)),i=1,nuse)
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
