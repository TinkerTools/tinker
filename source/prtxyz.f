c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtxyz  --  output of XYZ atomic coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtxyz" writes out a set of Cartesian atomic coordinates
c     to an external disk file in Tinker XYZ format
c
c
      subroutine prtxyz (ixyz)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use output
      use titles
      implicit none
      integer i,j,k,ixyz
      integer ii
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 xyzfile
c
c
c     open the output unit if not already done
c
      inquire (unit=ixyz,opened=opened)
      if (.not. opened) then
         xyzfile = filename(1:leng)//'.xyz'
         call version (xyzfile,'new')
         open (unit=ixyz,file=xyzfile,status='new')
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
      if (.not. onlysave) then
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (ixyz,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (ixyz,fstr(1:9))  n,title(1:ltitle)
         end if
      else
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (ixyz,fstr(1:4))  nonly
         else
            fstr = '('//atmc//',2x,a)'
            write (ixyz,fstr(1:9))  nonly,title(1:ltitle)
         end if
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (ixyz,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the atomic coordinates for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      if (.not. onlysave) then
         do i = 1, n
            k = n12(i)
            if (k .eq. 0) then
               write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i)
            else
               write (ixyz,fstr)  i,name(i),x(i),y(i),z(i),type(i),
     &                            (i12(j,i),j=1,k)
            end if
         end do
      else
         do ii = 1, nonly
            i = ionly(ii)
            k = n12(i)
            if (k .eq. 0) then
               write (ixyz,fstr)  ii,name(i),x(i),y(i),z(i),type(i)
            else
               write (ixyz,fstr)  ii,name(i),x(i),y(i),z(i),type(i),
     &                            (ionlyinv(i12(j,i)),j=1,k)
            end if
         end do
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=ixyz)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine prtvec3  --  output of atomic 3D vectors  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "prtvec3" writes out a set of atomic vectors
c     to an external disk file in Tinker XYZ format
c
c
      subroutine prtvec3 (iunit,mode)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use output
      use titles
      implicit none
      integer i,j,k
      integer ii,iunit
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*3 mode
      character*25 fstr
      character*240 outputfile
c
c
c     open the output unit if not already done
c
      inquire (unit=iunit,opened=opened)
      if (.not. opened) then
         if (mode .eq. 'XYZ') then
            outputfile = filename(1:leng)//'.xyz'
         else if (mode .eq. 'VEL') then
            outputfile = filename(1:leng)//'.vel'
         else if (mode .eq. 'FRC') then
            outputfile = filename(1:leng)//'.frc'
         else if (mode .eq. 'UIN') then
            outputfile = filename(1:leng)//'.uind'
         else if (mode .eq. 'UST') then
            outputfile = filename(1:leng)//'.ustc'
         else if (mode .eq. 'UCH') then
            outputfile = filename(1:leng)//'.uchg'
         else if (mode .eq. 'UDR') then
            outputfile = filename(1:leng)//'.udir'
         else if (mode .eq. 'DEF') then
            outputfile = filename(1:leng)//'.def'
         else if (mode .eq. 'TEF') then
            outputfile = filename(1:leng)//'.tef'
         end if
         call version (outputfile,'new')
         open (unit=iunit,file=outputfile,status='new')
      end if
c
c     copy the vector values to a common array
c
      call copyvec3 (print3n,mode)
c
c     check for large systems needing extended formatting
c
      atmc = 'i6'
      if (n .ge. 100000)  atmc = 'i7'
      if (n .ge. 1000000)  atmc = 'i8'
      crdmin = 0.0d0
      crdmax = 0.0d0
      do i = 1, n
         crdmin = min(crdmin,print3n(1,i),print3n(2,i),print3n(3,i))
         crdmax = max(crdmax,print3n(1,i),print3n(2,i),print3n(3,i))
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
      if (.not. onlysave) then
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (iunit,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (iunit,fstr(1:9))  n,title(1:ltitle)
         end if
      else
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (iunit,fstr(1:4))  nonly
         else
            fstr = '('//atmc//',2x,a)'
            write (iunit,fstr(1:9))  nonly,title(1:ltitle)
         end if
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (iunit,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the atomic coordinates for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      if (.not. onlysave) then
         do i = 1, n
            k = n12(i)
            if (k .eq. 0) then
               write (iunit,fstr)  i,name(i),(print3n(j,i),j=1,3),
     &                             type(i)
            else
               write (iunit,fstr)  i,name(i),(print3n(j,i),j=1,3),
     &                             type(i),(i12(j,i),j=1,k)
            end if
         end do
      else
         do ii = 1, nonly
            i = ionly(ii)
            k = n12(i)
            if (k .eq. 0) then
               write (iunit,fstr)  ii,name(i),(print3n(j,i),j=1,3),
     &                             type(i)
            else
               write (iunit,fstr)  ii,name(i),(print3n(j,i),j=1,3),
     &                             type(i),(ionlyinv(i12(j,i)),j=1,k)
            end if
         end do
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iunit)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine copyvec3  --  copy vector values to print3n  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "copyvec3" copies the values of the requested vector to the
c     common array "print3n"
c
c
      subroutine copyvec3 (print3n,mode)
      use atoms
      use boxes
      use deriv
      use expol
      use moldyn
      use mpole
      use polar
      use polpot
      use units
      implicit none
      integer i,j,k
      real*8 c
      real*8 xi,yi,zi
      real*8 print3n(3,*)
      character*3 mode
c
c
c     copy the requested vector values to the common array "print3n"
c
      if (mode .eq. 'XYZ') then
         do i = 1, n
            print3n(1,i) = x(i)
            print3n(2,i) = y(i)
            print3n(3,i) = z(i)
         end do
      else if (mode .eq. 'VEL') then
         do i = 1, n
            print3n(1,i) = v(1,i)
            print3n(2,i) = v(2,i)
            print3n(3,i) = v(3,i)
         end do
      else if (mode .eq. 'FRC') then
         do i = 1, n
            print3n(1,i) = -desum(1,i)
            print3n(2,i) = -desum(2,i)
            print3n(3,i) = -desum(3,i)
         end do
      else if (mode .eq. 'UIN') then
         do i = 1, n
            print3n(1,i) = debye*uind(1,i)
            print3n(2,i) = debye*uind(2,i)
            print3n(3,i) = debye*uind(3,i)
         end do
      else if (mode .eq. 'UST') then
         do i = 1, n
            print3n(1,i) = debye*rpole(2,i)
            print3n(2,i) = debye*rpole(3,i)
            print3n(3,i) = debye*rpole(4,i)
         end do
      else if (mode .eq. 'UCH') then
         do i = 1, n
            c = rpole(1,i)
            print3n(1,i) = c*debye*(x(i)-xcenter)
            print3n(2,i) = c*debye*(y(i)-ycenter)
            print3n(3,i) = c*debye*(z(i)-zcenter)
         end do
      else if (mode .eq. 'UDR') then
         if (.not. use_expol) then
            do i = 1, n
               do j = 1, 3
                  print3n(j,i) = debye*udir(j,i)
               end do
            end do
         else
            do i = 1, n
               do j = 1, 3
                  print3n(j,i) = 0.0d0
                  do k = 1, 3
                     print3n(j,i) = print3n(j,i)
     &                             + udir(k,i)*polinv(j,k,i)
                  end do
                  print3n(j,i) = debye*print3n(j,i)
               end do
            end do
         end if
      else if (mode .eq. 'DEF') then
         do i = 1, n
            c = elefield/polarity(i)
            do j = 1, 3
               print3n(j,i) = c * udir(j,i)
            end do
         end do
      else if (mode .eq. 'TEF') then
         if (.not. use_expol) then
            do i = 1, n
               c = elefield/polarity(i)
               do j = 1, 3
                  print3n(j,i) = c * uind(j,i)
               end do
            end do
         else
            do i = 1, n
               c = elefield/polarity(i)
               do j = 1, 3
                  print3n(j,i) = 0.0d0
                  do k = 1, 3
                     print3n(j,i) = print3n(j,i)
     &                             + uind(k,i)*polscale(j,k,i)
                  end do
                  print3n(j,i) = c * print3n(j,i)
               end do
            end do
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtdcd  --  output of DCD atomic coordinates  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtdcd" writes out a set of Cartesian atomic coordinates to
c     a file in CHARMM DCD binary format compatible with the VMD
c     visualization software and other packages
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
      use output
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
         ncrdsav = one
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
         if (.not. onlysave) then
            write (idcd)  n
         else
            write (idcd)  nonly
         end if
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
      if (.not. onlysave) then
         write (idcd)  (real(x(i)),i=1,n)
         write (idcd)  (real(y(i)),i=1,n)
         write (idcd)  (real(z(i)),i=1,n)
      else
         write (idcd)  (real(x(ionly(i))),i=1,nonly)
         write (idcd)  (real(y(ionly(i))),i=1,nonly)
         write (idcd)  (real(z(ionly(i))),i=1,nonly)
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine prtdcdv3  --  output of DCD atomic 3D vectors  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "prtdcdV3" writes out a set of atomic 3D vectors to
c     a file in CHARMM DCD binary format compatible with the VMD
c     visualization software and other packages
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
      subroutine prtdcdv3 (idcd,first,mode)
      use atoms
      use bound
      use boxes
      use files
      use output
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
      character*3 mode
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
      if (.not. opened) then
         if (mode .eq. 'XYZ') then
            dcdfile = filename(1:leng)//'.dcd'
         else if (mode .eq. 'VEL') then
            dcdfile = filename(1:leng)//'.dcdv'
         else if (mode .eq. 'FRC') then
            dcdfile = filename(1:leng)//'.dcdf'
         else if (mode .eq. 'UIN') then
            dcdfile = filename(1:leng)//'.dcdui'
         else if (mode .eq. 'UST') then
            dcdfile = filename(1:leng)//'.dcdus'
         else if (mode .eq. 'UCH') then
            dcdfile = filename(1:leng)//'.dcduc'
         else if (mode .eq. 'UDR') then
            dcdfile = filename(1:leng)//'.dcdud'
         else if (mode .eq. 'DEF') then
            dcdfile = filename(1:leng)//'.dcdde'
         else if (mode .eq. 'TEF') then
            dcdfile = filename(1:leng)//'.dcdte'
         end if
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
         ncrdsav = one
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
         if (.not. onlysave) then
            write (idcd)  n
         else
            write (idcd)  nonly
         end if
      end if
c
c     append the lattice values based on header flag value
c
      if (use_bounds) then
         write (idcd)  xbox,gamma_cos,ybox,beta_cos,alpha_cos,zbox
      end if
c
c     copy the vector values to a common array
c
      call copyvec3 (print3n,mode)
c
c     append the atomic coordinates along each axis in turn
c
      if (.not. onlysave) then
         write (idcd)  (real(print3n(1,i)),i=1,n)
         write (idcd)  (real(print3n(2,i)),i=1,n)
         write (idcd)  (real(print3n(3,i)),i=1,n)
      else
         write (idcd)  (real(print3n(1,ionly(i))),i=1,nonly)
         write (idcd)  (real(print3n(2,ionly(i))),i=1,nonly)
         write (idcd)  (real(print3n(3,ionly(i))),i=1,nonly)
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end