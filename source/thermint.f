c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine settisched  --  build the lambda schedule  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "settisched" fills the list of lambda values visited by the
c     thermodynamic integration windows, either from the explicit
c     "TI-WINDOW" keyword lines already parsed into "tilmdalist"
c     or, by default, as the "tinbin" evenly spaced values that
c     descend from one to zero; the fraction of the run spent in
c     each window is resolved into "tifraclist" summing to one
c
c
      subroutine settisched (ntiwin,tinbinset)
      use iounit
      use thrmint
      implicit none
      integer i
      integer ntiwin
      integer nunspec
      logical tinbinset
      logical up,down
      real*8 fsum,frem
      real*8, allocatable :: tlist(:)
c
c
c     with no explicit windows, space them evenly over [0,1]
c
      if (ntiwin .eq. 0) then
         if (tinbin .lt. 2) then
            write (iout,10)
   10       format (/,' SETTISCHED  --  TI-NBIN must be at least 2')
            call fatal
         end if
         if (allocated(tilmdalist))  deallocate (tilmdalist)
         if (allocated(tifraclist))  deallocate (tifraclist)
         allocate (tilmdalist(tinbin))
         allocate (tifraclist(tinbin))
         do i = 1, tinbin
            tilmdalist(i) = 1.0d0 - dble(i-1)/dble(tinbin-1)
            tifraclist(i) = 1.0d0 / dble(tinbin)
         end do
         tilmdalist(tinbin) = 0.0d0
      else
c
c     an explicit schedule sets the number of windows itself
c
         if (tinbinset) then
            write (iout,20)
   20       format (/,' SETTISCHED  --  TI-NBIN and TI-WINDOW cannot',
     &                 ' both be used')
            call fatal
         end if
         tinbin = ntiwin
c
c     compact the parsed values down to the number of windows
c
         allocate (tlist(tinbin))
         do i = 1, tinbin
            tlist(i) = tilmdalist(i)
         end do
         deallocate (tilmdalist)
         allocate (tilmdalist(tinbin))
         do i = 1, tinbin
            tilmdalist(i) = tlist(i)
         end do
         do i = 1, tinbin
            tlist(i) = tifraclist(i)
         end do
         deallocate (tifraclist)
         allocate (tifraclist(tinbin))
         do i = 1, tinbin
            tifraclist(i) = tlist(i)
         end do
         deallocate (tlist)
c
c     a time share that was given must be a positive fraction
c
         fsum = 0.0d0
         nunspec = 0
         do i = 1, tinbin
            if (tifraclist(i) .lt. 0.0d0) then
               nunspec = nunspec + 1
            else if (tifraclist(i) .eq. 0.0d0) then
               write (iout,50)  i
   50          format (/,' SETTISCHED  --  TI-WINDOW',i5,' was given',
     &                    ' a time fraction of zero')
               call fatal
            else
               fsum = fsum + tifraclist(i)
            end if
         end do
c
c     share whatever time is left among the windows that did not
c     ask for a fraction of their own
c
         if (nunspec .gt. 0) then
            frem = 1.0d0 - fsum
            if (frem .le. 0.0d0) then
               write (iout,60)  fsum
   60          format (/,' SETTISCHED  --  TI-WINDOW fractions total',
     &                    f12.6,' leaving no time for the windows',
     &                 /,'                 without an explicit value')
               call fatal
            end if
            do i = 1, tinbin
               if (tifraclist(i) .lt. 0.0d0) then
                  tifraclist(i) = frem / dble(nunspec)
               end if
            end do
         end if
c
c     rescale the time shares so that they span the whole run
c
         fsum = 0.0d0
         do i = 1, tinbin
            fsum = fsum + tifraclist(i)
         end do
         do i = 1, tinbin
            tifraclist(i) = tifraclist(i) / fsum
         end do
c
c     each window must sit within the physical lambda range
c
         do i = 1, tinbin
            if (tilmdalist(i).lt.0.0d0 .or.
     &          tilmdalist(i).gt.1.0d0) then
               write (iout,30)  i,tilmdalist(i)
   30          format (/,' SETTISCHED  --  TI-WINDOW',i5,' value',
     &                    f12.6,' is outside [0,1]')
               call fatal
            end if
         end do
c
c     the schedule must walk in one direction without repeating
c
         if (tinbin .ge. 2) then
            up = .true.
            down = .true.
            do i = 2, tinbin
               if (tilmdalist(i) .le. tilmdalist(i-1))  up = .false.
               if (tilmdalist(i) .ge. tilmdalist(i-1))  down = .false.
            end do
            if (.not.up .and. .not.down) then
               write (iout,40)
   40          format (/,' SETTISCHED  --  TI-WINDOW values must',
     &                    ' increase or decrease monotonically')
               call fatal
            end if
         end if
      end if
c
c     start the schedule at its first window
c
      tibin = 1
      tilmda = tilmdalist(1)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine inittidyn  --  set up lambda window layout  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "inittidyn" divides a dynamics run of "nstep" steps among the
c     lambda windows according to their requested time fractions,
c     sizes the block average accumulators, and puts the main lambda
c     at the start of the schedule
c
c
      subroutine inittidyn (nstep)
      use dlmda
      use thrmint
      implicit none
      integer i
      integer nstep
      real*8 cum
c
c
c     give each window the share of the run that it asked for, and
c     pin the last boundary so the whole trajectory is covered
c
      if (allocated(tiwinend))  deallocate (tiwinend)
      allocate (tiwinend(tinbin))
      cum = 0.0d0
      do i = 1, tinbin
         cum = cum + tifraclist(i)
         tiwinend(i) = nint(cum*dble(nstep))
      end do
      tiwinend(tinbin) = nstep
c
c     size the accumulators to the schedule and map the sublambdas
c
      call settiblocks
      call mapsublmda (tilmda)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine settiblocks  --  size the block accumulators  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "settiblocks" counts the block averages the window boundaries
c     in "tiwinend" can hold, allocates the recording arrays to that
c     exact length, and rewinds the schedule to its first window
c
c
      subroutine settiblocks
      use iounit
      use thrmint
      implicit none
      integer i
      integer nw,ne,nb
      integer istart
c
c
c     count the blocks each window can hold; a window too short for
c     a complete block still runs, but records nothing
c
      tinbtot = 0
      istart = 0
      do i = 1, tinbin
         nw = tiwinend(i) - istart
         if (nw .lt. 1) then
            write (iout,10)  i,tilmdalist(i)
   10       format (/,' SETTIBLOCKS  --  TI-WINDOW',i5,' at lambda',
     &                 f12.6,' was given no dynamics steps')
            call fatal
         end if
         ne = int(dble(nw) * tieqratio)
         nb = (nw-ne) / tinstepavg
         if (nb .eq. 0) then
            write (iout,20)  i,tilmdalist(i),tifraclist(i)
   20       format (/,' SETTIBLOCKS  --  TI-WINDOW',i5,' at lambda',
     &                 f12.6,' with fraction',f12.6,
     &              /,'                  is shorter than TI-NSTEPAVG',
     &                 ' and will record no samples')
         end if
         tinbtot = tinbtot + nb
         istart = tiwinend(i)
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(tilmdahist))  deallocate (tilmdahist)
      if (allocated(tilmdadedl))  deallocate (tilmdadedl)
      if (allocated(tilmdadedlstd))  deallocate (tilmdadedlstd)
      allocate (tilmdahist(max(1,tinbtot)))
      allocate (tilmdadedl(max(1,tinbtot)))
      allocate (tilmdadedlstd(max(1,tinbtot)))
c
c     zero out the block averages recorded along the schedule
c
      do i = 1, max(1,tinbtot)
         tilmdahist(i) = 0.0d0
         tilmdadedl(i) = 0.0d0
         tilmdadedlstd(i) = 0.0d0
      end do
      tinbcount = 0
      tinbsave = 0
c
c     start the schedule at its first window
c
      tibin = 1
      tilmda = tilmdalist(1)
      call settiwindow
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine settiwindow  --  size the current TI window  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "settiwindow" sets the step count, equilibration length and
c     block capacity of the lambda window that is currently active
c
c
      subroutine settiwindow
      use thrmint
      implicit none
c
c
c     measure the current window against the preceding boundary
c
      tiwindow = tiwinend(tibin)
      if (tibin .gt. 1)  tiwindow = tiwinend(tibin) - tiwinend(tibin-1)
      tinequil = int(dble(tiwindow) * tieqratio)
      tinblock = (tiwindow-tinequil) / tinstepavg
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prttihead  --  start the block average file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prttihead" creates the external file that collects the block
c     averaged lambda derivatives and writes the fixed header that
c     describes the window layout and the lambda schedule; the block
c     averages themselves are appended later by "saveti"
c
c
      subroutine prttihead
      use files
      use iounit
      use thrmint
      implicit none
      integer w
      integer iti
      integer freeunit
      integer trimtext
c
c
c     open a new file, keeping any output from a previous run
c
      iti = freeunit ()
      tifile = filename(1:leng)//'.ti'
      call version (tifile,'new')
      open (unit=iti,file=tifile,status='new')
c
c     write a header describing the window and block layout
c
      write (iti,10)
   10 format ('# tinker thermodynamic integration')
      write (iti,20)  tinbin,tinstepavg,tieqratio,tiwinend(tinbin),
     &                tinbtot
   20 format ('# tinbin ',i0,' tinstepavg ',i0,' tieqratio ',f12.6,
     &           ' nstep ',i0,' tinbtot ',i0)
c
c     record the full schedule, so the file still describes every
c     window when an interrupted run leaves some of them empty
c
      do w = 1, tinbin
         write (iti,30)  w,tilmdalist(w),tifraclist(w),tiwinend(w)
   30    format ('# schedule ',i0,1x,f12.8,1x,f12.8,1x,i0)
      end do
      write (iti,40)
   40 format ('# index lambda dedl dedlstd')
      close (unit=iti)
c
c     report the name of the file holding the block averages
c
      write (iout,50)  tifile(1:trimtext(tifile))
   50 format (/,' TI  --  dU/dlambda Block Averages Written To  ',a)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine etidyn  --  thermodynamic integration sampling  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "etidyn" collects the lambda derivative of the potential energy
c     at each dynamics step, averages it into blocks over the
c     production portion of the current lambda window, and advances
c     the window once its final step is reached
c
c
      subroutine etidyn (istep)
      use dlmda
      use thrmint
      implicit none
      integer istep
      integer tistep
      integer tiprod
      integer tistart
      real*8 avg,std
c
c
c     nothing is left to sample once the schedule has run out
c
      if (tibin .gt. tinbin)  return
c
c     find the position of this step within the current window
c
      tistart = 0
      if (tibin .gt. 1)  tistart = tiwinend(tibin-1)
      tistep = istep - tistart
c
c     nothing is stored while the window is equilibrating
c
      if (tistep .gt. tinequil) then
         tiprod = tistep - tinequil
         tidedllist(mod(tiprod-1,tinstepavg)+1) = dedl
c
c     reduce a full block into its average and deviation, keeping
c     the lambda that produced it alongside the block itself
c
         if (mod(tiprod,tinstepavg) .eq. 0) then
            call avgstd (tidedllist,1,tinstepavg,avg,std)
            if (tinbcount .lt. tinbtot) then
               tinbcount = tinbcount + 1
               tilmdahist(tinbcount) = tilmda
               tilmdadedl(tinbcount) = avg
               tilmdadedlstd(tinbcount) = std
            end if
         end if
      end if
c
c     move on to the next lambda window at the window boundary
c
      if (istep .eq. tiwinend(tibin))  call tischedule
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine tischedule  --  advance the lambda window  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "tischedule" moves the main lambda to the value of the next
c     window in the schedule; past the final window the lambda is
c     left where it is, and "etidyn" stops collecting samples
c
c
      subroutine tischedule
      use thrmint
      implicit none
c
c
c     take the next lambda and refresh the sublambda values
c
      tibin = tibin + 1
      if (tibin .le. tinbin) then
         tilmda = tilmdalist(tibin)
         call settiwindow
         call mapsublmda (tilmda)
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine saveti  --  append new TI block averages  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "saveti" appends to the external file the block averaged
c     lambda derivatives and their deviations recorded since the
c     previous call, so a running simulation can be monitored as
c     the lambda windows are completed
c
c
      subroutine saveti
      use dlmda
      use thrmint
      implicit none
      integer i
      integer iti
      integer freeunit
c
c
c     return if thermodynamic integration was never initialized
c
      if (.not. use_ti)  return
      if (.not. allocated(tilmdadedl))  return
c
c     skip the file entirely when no block has completed since the
c     last time the averages were written out
c
      if (tinbcount .le. tinbsave)  return
c
c     append the block averages recorded since the previous call
c
      iti = freeunit ()
      open (unit=iti,file=tifile,status='old',position='append')
      do i = tinbsave+1, tinbcount
         write (iti,10)  i,tilmdahist(i),tilmdadedl(i),
     &                   tilmdadedlstd(i)
   10    format (i6,1x,f12.8,1p,1x,e20.10,1x,e20.10)
      end do
      tinbsave = tinbcount
      close (unit=iti)
      return
      end
