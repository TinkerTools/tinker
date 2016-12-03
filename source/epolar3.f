c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2015  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine epolar3  --  induced dipole energy & analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "epolar3" calculates the dipole polarizability interaction
c     from the induced dipoles times the electric field
c
c
      subroutine epolar3
      use sizes
      use action
      use analyz
      use atoms
      use atomid
      use boxes
      use chgpot
      use energi
      use ewald
      use inform
      use iounit
      use limits
      use math
      use mpole
      use polar
      use potent
      use units
      implicit none
      integer i,j,ii
      integer nepo
      real*8 e,epo,f,fi
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 dix,diy,diz
      real*8 uix,uiy,uiz
      real*8, allocatable :: aepo(:)
      logical header,huge
c
c
c     zero out the total polarization energy and partitioning
c
      nep = 0
      ep = 0.0d0
      do i = 1, n
         aep(i) = 0.0d0
      end do
      if (npole .eq. 0)  return
c
c     check the sign of multipole components at chiral sites
c
      if (.not. use_mpole)  call chkpole
c
c     rotate the multipole components into the global frame
c
      if (.not. use_mpole)  call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce
c
c     perform dynamic allocation of some local arrays
c
      allocate (aepo(n))
c
c     set the energy conversion factor
c
      f = -0.5d0 * electric / dielec
c
c     print header information if debug output was requested
c
      header = .true.
      if (debug) then
         header = .false.
         write (iout,10)
   10    format (/,' Individual Dipole Polarization Interactions :',
     &           //,' Type',9x,'Atom Name',24x,'Alpha',8x,'Energy',/)
      end if
c
c     initialize local variables for OpenMP calculation
c
      epo = 0.0d0
      nepo = nep
      do i = 1, n
         aepo(i) = aep(i)
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private)
!$OMP& shared(npole,ipole,polarity,f,uind,udirp,name,verbose,
!$OMP& debug,header,iout,epo,nepo,aepo)
!$OMP DO reduction(+:epo,nepo,aepo) schedule(guided)
c
c     get polarization energy via induced dipoles times field
c
      do i = 1, npole
         if (polarity(i) .ne. 0.0d0) then
            ii = ipole(i)
            fi = f / polarity(i)
            e = 0.0d0
            do j = 1, 3
               e = e + fi*uind(j,i)*udirp(j,i)
            end do
            nepo = nepo + 1
            epo = epo + e
            aepo(ii) = aepo(ii) + e
c
c     print a message if the energy for this site is large
c
            huge = (abs(e) .gt. 100.0d0)
            if (debug .or. (verbose.and.huge)) then
               if (header) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Individual Dipole Polarization',
     &                       ' Interactions :',
     &                    //,' Type',9x,'Atom Name',24x,'Alpha',
     &                       8x,'Energy',/)
               end if
               write (iout,30)  ii,name(ii),polarity(i),e
   30          format (' Polar',5x,i7,'-',a3,16x,2f14.4)
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
c
c     add local copies to global variables for OpenMP calculation
c
      ep = ep + epo
      nep = nepo
      do i = 1, n
         aep(i) = aepo(i)
      end do
c
c     compute the cell dipole boundary correction term
c
      if (use_ewald) then
         if (boundary .eq. 'VACUUM') then
            f = electric / dielec
            xd = 0.0d0
            yd = 0.0d0
            zd = 0.0d0
            xu = 0.0d0
            yu = 0.0d0
            zu = 0.0d0
            do i = 1, npole
               ii = ipole(i)
               dix = rpole(2,i)
               diy = rpole(3,i)
               diz = rpole(4,i)
               uix = uind(1,i)
               uiy = uind(2,i)
               uiz = uind(3,i)
               xd = xd + dix + rpole(1,i)*x(ii)
               yd = yd + diy + rpole(1,i)*y(ii)
               zd = zd + diz + rpole(1,i)*z(ii)
               xu = xu + uix
               yu = yu + uiy
               zu = zu + uiz
            end do
            e = (2.0d0/3.0d0) * f * (pi/volbox) * (xd*xu+yd*yu+zd*zu)
            nep = nep + 1
            ep = ep + e
            do i = 1, npole
               ii = ipole(i)
               aep(ii) = aep(ii) + e/dble(npole)
            end do
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (aepo)
      return
      end
