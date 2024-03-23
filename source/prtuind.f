c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2023  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine prtuind  --  output of atomic induced dipoles  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "prtuind" writes out a set of induced dipole components
c     to an external disk file in Tinker XYZ format
c
c
      subroutine prtuind (iind)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use polar
      use output
      use titles
      use units
      implicit none
      integer i,j,k,iind
      integer ii
      integer size,crdsiz
      real*8 crdmin,crdmax
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 indfile
c
c
c     open the output unit if not already done
c
      inquire (unit=iind,opened=opened)
      if (.not. opened) then
         indfile = filename(1:leng)//'.uind'
         call version (indfile,'new')
         open (unit=iind,file=indfile,status='new')
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
            write (iind,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (iind,fstr(1:9))  n,title(1:ltitle)
         end if
      else
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (iind,fstr(1:4))  nonly
         else
            fstr = '('//atmc//',2x,a)'
            write (iind,fstr(1:9))  nonly,title(1:ltitle)
         end if
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (iind,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the induced dipole components for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      if (.not. onlysave) then
         do i = 1, n
            k = n12(i)
            if (k .eq. 0) then
               write (iind,fstr)  i,name(i),(debye*uind(j,i),j=1,3),
     &                            type(i)
            else
               write (iind,fstr)  i,name(i),(debye*uind(j,i),j=1,3),
     &                            type(i),(i12(j,i),j=1,k)
            end if
         end do
      else
         do ii = 1, nonly
            i = ionly(ii)
            k = n12(i)
            if (k .eq. 0) then
               write (iind,fstr) ii,name(i),(debye*uind(j,i),j=1,3),
     &                           type(i)
            else
               write (iind,fstr) ii,name(i),(debye*uind(j,i),j=1,3),
     &                           type(i), (ionlyinv(i12(j,i)),j=1,k)
            end if
         end do
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=iind)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine prtdcdu  --  output of DCD induced dipoles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "prtdcdu" writes out a set of induced dipole components to
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
      subroutine prtdcdu (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use polar
      use output
      use titles
      use units
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
         dcdfile = filename(1:leng)//'.dcdu'
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
c     append the induced dipoles along each axis in turn
c
      if (.not. onlysave) then
         write (idcd)  (real(debye*uind(1,i)),i=1,n)
         write (idcd)  (real(debye*uind(2,i)),i=1,n)
         write (idcd)  (real(debye*uind(3,i)),i=1,n)
      else
         write (idcd)  (real(debye*uind(1,ionly(i))),i=1,nonly)
         write (idcd)  (real(debye*uind(2,ionly(i))),i=1,nonly)
         write (idcd)  (real(debye*uind(3,ionly(i))),i=1,nonly)
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine prtustc  --  output of atomic static dipoles  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "prtustc" writes out a set of static dipole components
c     to an external disk file in Tinker XYZ format
c
c
      subroutine prtustc (istc)
      use atomid
      use atoms
      use bound
      use boxes
      use couple
      use files
      use inform
      use mpole
      use output
      use titles
      use units
      implicit none
      integer i,j,k,istc
      integer ii
      integer size,crdsiz
      real*8 crdmin,crdmax
      real*8 c,xd,yd,zd
      logical opened
      character*2 atmc
      character*2 crdc
      character*2 digc
      character*25 fstr
      character*240 stcfile
c
c
c     open the output unit if not already done
c
      inquire (unit=istc,opened=opened)
      if (.not. opened) then
         stcfile = filename(1:leng)//'.ustc'
         call version (stcfile,'new')
         open (unit=istc,file=stcfile,status='new')
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
            write (istc,fstr(1:4))  n
         else
            fstr = '('//atmc//',2x,a)'
            write (istc,fstr(1:9))  n,title(1:ltitle)
         end if
      else
         if (ltitle .eq. 0) then
            fstr = '('//atmc//')'
            write (istc,fstr(1:4))  nonly
         else
            fstr = '('//atmc//',2x,a)'
            write (istc,fstr(1:9))  nonly,title(1:ltitle)
         end if
      end if
c
c     write out the periodic cell lengths and angles
c
      if (use_bounds) then
         fstr = '(1x,6f'//crdc//'.'//digc//')'
         write (istc,fstr)  xbox,ybox,zbox,alpha,beta,gamma
      end if
c
c     write out the static dipole components for each atom
c
      fstr = '('//atmc//',2x,a3,3f'//crdc//
     &          '.'//digc//',i6,8'//atmc//')'
      if (.not. onlysave) then
         do i = 1, n
            c = rpole(1,i)
            xd = (x(i)*c + rpole(2,i)) * debye
            yd = (y(i)*c + rpole(3,i)) * debye
            zd = (z(i)*c + rpole(4,i)) * debye
            k = n12(i)
            if (k .eq. 0) then
               write (istc,fstr)  i,name(i),xd,yd,zd,type(i)
            else
               write (istc,fstr)  i,name(i),xd,yd,zd,type(i),
     &                            (i12(j,i),j=1,k)
            end if
         end do
      else
         do ii = 1, nonly
            i = ionly(ii)
            c = rpole(1,i)
            xd = (x(i)*c + rpole(2,i)) * debye
            yd = (y(i)*c + rpole(3,i)) * debye
            zd = (z(i)*c + rpole(4,i)) * debye
            k = n12(i)
            if (k .eq. 0) then
               write (istc,fstr) ii,name(i),xd,yd,zd,type(i)
            else
               write (istc,fstr) ii,name(i),xd,yd,zd,type(i),
     &                           (ionlyinv(i12(j,i)),j=1,k)
            end if
         end do
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=istc)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine prtdcdd  --  output of DCD static dipoles  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "prtdcdd" writes out a set of static dipole components to
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
      subroutine prtdcdd (idcd,first)
      use atoms
      use bound
      use boxes
      use files
      use mpole
      use output
      use titles
      use units
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
         dcdfile = filename(1:leng)//'.dcdd'
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
c     append the static dipoles along each axis in turn
c
      if (.not. onlysave) then
         write (idcd) (real(debye*(x(i)*rpole(1,i) + rpole(2,i))),i=1,n)
         write (idcd) (real(debye*(y(i)*rpole(1,i) + rpole(3,i))),i=1,n)
         write (idcd) (real(debye*(z(i)*rpole(1,i) + rpole(4,i))),i=1,n)
      else
         write (idcd) (real(debye*(x(ionly(i))*rpole(1,ionly(i))
     &                      + rpole(2,ionly(i)))),i=1,nonly)
         write (idcd) (real(debye*(y(ionly(i))*rpole(1,ionly(i))
     &                      + rpole(3,ionly(i)))),i=1,nonly)
         write (idcd) (real(debye*(z(ionly(i))*rpole(1,ionly(i))
     &                      + rpole(4,ionly(i)))),i=1,nonly)
      end if
c
c     close the output unit if opened by this routine
c
      if (.not. opened)  close (unit=idcd)
      return
      end
