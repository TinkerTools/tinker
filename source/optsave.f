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
      integer iopt,ifrc
      integer iind,iend
      integer ncycle,nvar
      integer freeunit
      integer trimtext
      real*8 f,xx(*)
      logical exist,first
      character*7 ext
      character*240 optfile
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
      iopt = freeunit ()
      if (cyclesave) then
         if (dcdsave) then
            optfile = filename(1:leng)
            call suffix (optfile,'dcd','old')
            inquire (file=optfile,exist=exist)
            if (exist) then
               first = .false.
               open (unit=iopt,file=optfile,form='unformatted',
     &                  status='old',position='append')
            else
               first = .true.
               open (unit=iopt,file=optfile,form='unformatted',
     &                  status='new')
            end if
            call prtdcd (iopt,first)
         else if (arcsave) then
            optfile = filename(1:leng)
            call suffix (optfile,'arc','old')
            inquire (file=optfile,exist=exist)
            if (exist) then
               call openend (iopt,optfile)
            else
               open (unit=iopt,file=optfile,status='new')
            end if
         else
            lext = 3
            call numeral (ncycle,ext,lext)
            optfile = filename(1:leng)//'.'//ext(1:lext)
            call version (optfile,'new')
            open (unit=iopt,file=optfile,status='new')
         end if
      else
         optfile = outfile
         call version (optfile,'old')
         open (unit=iopt,file=optfile,status='old')
         rewind (unit=iopt)
      end if
c
c     update intermediate file with desired coordinate type
c
      if (coordtype .eq. 'CARTESIAN') then
         if (.not. dcdsave)  call prtxyz (iopt)
      else if (coordtype .eq. 'INTERNAL') then
         call prtint (iopt)
      else if (coordtype .eq. 'RIGIDBODY') then
         call prtxyz (iopt)
      end if
      close (unit=iopt)
c
c     save the force vector components for the current step
c
      if (frcsave .and. coordtype.eq.'CARTESIAN') then
         ifrc = freeunit ()
         if (cyclesave) then
            if (arcsave) then
               frcfile = filename(1:leng)
               call suffix (frcfile,'frc','old')
               inquire (file=frcfile,exist=exist)
               if (exist) then
                  call openend (ifrc,frcfile)
               else
                  open (unit=ifrc,file=frcfile,status='new')
               end if
            else
               frcfile = filename(1:leng)//'.'//ext(1:lext)//'f'
               call version (frcfile,'new')
               open (unit=ifrc,file=frcfile,status='new')
            end if
         else
            frcfile = filename(1:leng)
            call suffix (frcfile,'frc','old')
            inquire (file=frcfile,exist=exist)
            if (exist) then
               open (unit=ifrc,file=frcfile,status='old')
               rewind (unit=ifrc)
            else
               open (unit=ifrc,file=frcfile,status='new')
            end if
         end if
         write (ifrc,10)  n,title(1:ltitle)
   10    format (i6,2x,a)
         do i = 1, n
            write (ifrc,20)  i,name(i),(-desum(j,i),j=1,3)
   20       format (i6,2x,a3,3x,d13.6,3x,d13.6,3x,d13.6)
         end do
         close (unit=ifrc)
         write (iout,30)  frcfile(1:trimtext(frcfile))
   30    format (' Force Vector File',11x,a)
      end if
c
c     save the current induced dipole moment at each site
c
      if (uindsave .and. use_polar .and. coordtype.eq.'CARTESIAN') then
         iind = freeunit ()
         if (cyclesave) then
            if (arcsave) then
               indfile = filename(1:leng)
               call suffix (indfile,'uind','old')
               inquire (file=indfile,exist=exist)
               if (exist) then
                  call openend (iind,indfile)
               else
                  open (unit=iind,file=indfile,status='new')
               end if
            else
               indfile = filename(1:leng)//'.'//ext(1:lext)//'u'
               call version (indfile,'new')
               open (unit=iind,file=indfile,status='new')
            end if
         else
            indfile = filename(1:leng)
            call suffix (indfile,'uind','old')
            inquire (file=indfile,exist=exist)
            if (exist) then
               open (unit=iind,file=indfile,status='old')
               rewind (unit=iind)
            else
               open (unit=iind,file=indfile,status='new')
            end if
         end if
         write (iind,40)  n,title(1:ltitle)
   40    format (i6,2x,a)
         do ii = 1, npole
            i = ipole(ii)
            if (polarity(i) .ne. 0.0d0) then
               write (iind,50)  i,name(i),(debye*uind(j,i),j=1,3)
   50          format (i6,2x,a3,3f12.6)
            end if
         end do
         close (unit=iind)
         write (iout,60)  indfile(1:trimtext(indfile))
   60    format (' Induced Dipole File',10x,a)
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
         write (iout,70)
   70    format (/,' OPTSAVE  --  Optimization Calculation Ending',
     &              ' due to User Request')
         call fatal
      end if
      return
      end
