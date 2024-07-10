c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine optsave  --  save optimization info and results  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "optsave" is used by the optimizers to write imtermediate
c     coordinates and other relevant information; also checks for
c     user requested termination of an optimization
c
c
      subroutine optsave (ncycle,f,xx)
      use atomid
      use atoms
      use bound
      use deriv
      use files
      use group
      use iounit
      use math
      use mpole
      use omega
      use output
      use polar
      use potent
      use scales
      use socket
      use titles
      use units
      use usage
      use zcoord
      implicit none
      integer i,j,ii,lext
      integer ixyz,ifrc
      integer iind,iend
      integer ncycle,nvar
      integer freeunit
      integer trimtext
      real*8 f,xx(*)
      logical exist,first
      character*7 ext
      character*240 xyzfile
      character*240 frcfile
      character*240 indfile
      character*240 endfile
c
c
c     nothing to do if coordinate type is undefined
c
      if (coordtype .eq. 'NONE')  return
c
c     check scaling factors for optimization parameters
c
      if (.not. set_scale) then
         set_scale = .true.
         if (coordtype .eq. 'CARTESIAN') then
            if (.not. allocated(scale))  allocate (scale(3*n))
            do i = 1, 3*n
               scale(i) = 1.0d0
            end do
         else if (coordtype .eq. 'INTERNAL') then
            if (.not. allocated(scale))  allocate (scale(nomega))
            do i = 1, nomega
               scale(i) = 1.0d0
            end do
         else if (coordtype .eq. 'RIGIDBODY') then
            if (.not. allocated(scale))  allocate (scale(6*ngrp))
            do i = 1, 6*ngrp
               scale(i) = 1.0d0
            end do
         end if
      end if
c
c     convert optimization parameters to atomic coordinates
c
      if (coordtype .eq. 'CARTESIAN') then
         nvar = 0
         do ii = 1, nuse
            i = iuse(ii)
            nvar = nvar + 1
            x(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            y(i) = xx(nvar) / scale(nvar)
            nvar = nvar + 1
            z(i) = xx(nvar) / scale(nvar)
         end do
      else if (coordtype .eq. 'INTERNAL') then
         do i = 1, nomega
            dihed(i) = xx(i) / scale(i)
            ztors(zline(i)) = dihed(i) * radian
         end do
      end if
c
c     move stray molecules into periodic box if desired
c
      if (coordtype .eq. 'CARTESIAN') then
         if (use_bounds)  call bounds
      end if
c
c     save coordinates to archive or numbered structure file
c
      ixyz = freeunit ()
      if (cyclesave) then
         if (dcdsave) then
            xyzfile = filename(1:leng)
            call suffix (xyzfile,'dcd','old')
            inquire (file=xyzfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=ixyz,file=xyzfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=ixyz,file=xyzfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcd (ixyz,first)
         else if (arcsave) then
            xyzfile = filename(1:leng)
            call suffix (xyzfile,'arc','old')
            inquire (file=xyzfile,exist=exist)
            if (exist) then
               call openend (ixyz,xyzfile)
            else
               open (unit=ixyz,file=xyzfile,status='new')
            end if
         else
            lext = 3
            call numeral (ncycle,ext,lext)
            xyzfile = filename(1:leng)//'.'//ext(1:lext)
            call version (xyzfile,'new')
            open (unit=ixyz,file=xyzfile,status='new')
         end if
      else
         xyzfile = outfile
         call version (xyzfile,'old')
         open (unit=ixyz,file=xyzfile,status='old')
         rewind (unit=ixyz)
      end if
c
c     update intermediate file with desired coordinate type
c
      if (coordtype .eq. 'CARTESIAN') then
         if (.not. dcdsave)  call prtxyz (ixyz)
      else if (coordtype .eq. 'INTERNAL') then
         call prtint (ixyz)
      else if (coordtype .eq. 'RIGIDBODY') then
         call prtxyz (ixyz)
      end if
      close (unit=ixyz)
c
c     save the force vector components for the current step
c
      if (frcsave .and. coordtype.eq.'CARTESIAN') then
         ifrc = freeunit ()
         if (cyclesave) then
            if (dcdsave) then
               frcfile = filename(1:leng)
               call suffix (frcfile,'dcdf','old')
               inquire (file=frcfile,exist=exist)
               if (exist) then
                  first = .false.
                  open (unit=ifrc,file=frcfile,form='unformatted',
     &                     status='old',position='append')
               else
                  first = .true.
                  open (unit=ifrc,file=frcfile,form='unformatted',
     &                     status='new')
               end if
               call prtdcdf (ifrc,first)
            else if (arcsave) then
               frcfile = filename(1:leng)
               call suffix (frcfile,'frc','old')
               inquire (file=frcfile,exist=exist)
               if (exist) then
                  call openend (ifrc,frcfile)
               else
                  open (unit=ifrc,file=frcfile,status='new')
               end if
               call prtfrc (ifrc)
            else
               frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
               call version (frcfile,'new')
               open (unit=ifrc,file=frcfile,status='new')
               call prtfrc (ifrc)
            end if
         else
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            call version (frcfile,'old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               open (unit=ifrc,file=frcfile,status='old')
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
            rewind (unit=ifrc)
            call prtfrc (ifrc)
         end if
         close (unit=ifrc)
         write (iout,10)  frcfile(1:trimtext(frcfile))
   10    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar .and. coordtype.eq.'CARTESIAN') then
         iind = freeunit ()
         if (cyclesave) then
            if (dcdsave) then
               indfile = filename(1:leng)
               call suffix (indfile,'dcdu','old')
               inquire (file=indfile,exist=exist)
               if (exist) then
                  first = .false.
                  open (unit=iind,file=indfile,form='unformatted',
     &                     status='old',position='append')
               else
                  first = .true.
                  open (unit=iind,file=indfile,form='unformatted',
     &                     status='new')
               end if
               call prtdcdu (iind,first)
            else if (arcsave) then
               indfile = filename(1:leng)
               call suffix (indfile,'uind','old')
               inquire (file=indfile,exist=exist)
               if (exist) then
                  call openend (iind,indfile)
               else
                  open (unit=iind,file=indfile,status='new')
               end if
               call prtuind (iind)
            else
               indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
               call version (indfile,'new')
               open (unit=iind,file=indfile,status='new')
               call prtuind (iind)
            end if
         else
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            call version (indfile,'old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               open (unit=iind,file=indfile,status='old')
            else
               open (unit=iind,file=indfile,status='new')
            end if
            rewind (unit=iind)
            call prtuind (iind)
         end if
         close (unit=iind)
         write (iout,20)  indfile(1:trimtext(indfile))
   20    format (' Induced Dipole File',10x,a)
      end if
c
c     send data via external socket communication if desired
c
      if (.not.sktstart .or. use_socket) then
         if (coordtype .eq. 'INTERNAL')  call makexyz
         call sktopt (ncycle,f)
      end if
c
c     test for requested termination of the optimization
c
      endfile = 'tinker.end'
      inquire (file=endfile,exist=exist)
      if (.not. exist) then
         endfile = filename(1:leng)//'.end'
         inquire (file=endfile,exist=exist)
         if (exist) then
            iend = freeunit ()
            open (unit=iend,file=endfile,status='old')
            close (unit=iend,status='delete')
         end if
      end if
      if (exist) then
         write (iout,30)
   30    format (/,' OPTSAVE  --  Optimization Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
      return
      end
