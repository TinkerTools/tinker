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
c     descend from one to zero
c
c
      subroutine settisched (ntiwin,tinbinset)
      use iounit
      use thrmint
      implicit none
      integer i
      integer ntiwin
      logical tinbinset
      logical up,down
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
         allocate (tilmdalist(tinbin))
         do i = 1, tinbin
            tilmdalist(i) = 1.0d0 - dble(i-1)/dble(tinbin-1)
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
         deallocate (tlist)
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
c     "inittidyn" divides a dynamics run of "nstep" steps into equal
c     lambda windows, sizes the block average accumulators, and puts
c     the main lambda at the start of the schedule
c
c
      subroutine inittidyn (nstep)
      use dlmda
      use iounit
      use thrmint
      implicit none
      integer i,j
      integer nstep
c
c
c     split the trajectory evenly among the requested lambda windows
c
      tiwindow = nstep / tinbin
      if (tiwindow .lt. 1) then
         write (iout,10)
   10    format (/,' INITTIDYN  --  Fewer Dynamics Steps than',
     &              ' TI-NBIN Windows')
         call fatal
      end if
c
c     discard the leading fraction of each window as equilibration
c
      tinequil = int(dble(tiwindow) * tieqratio)
      if (tiwindow-tinequil .lt. tinstepavg) then
         write (iout,20)
   20    format (/,' INITTIDYN  --  Production Block is Shorter',
     &              ' than TI-NSTEPAVG')
         call fatal
      end if
c
c     each window holds the same number of complete blocks
c
      tinblock = (tiwindow-tinequil) / tinstepavg
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(tinbcount))  deallocate (tinbcount)
      if (allocated(tilmdadedl))  deallocate (tilmdadedl)
      if (allocated(tilmdadedlstd))  deallocate (tilmdadedlstd)
      allocate (tinbcount(tinbin))
      allocate (tilmdadedl(tinbin,tinblock))
      allocate (tilmdadedlstd(tinbin,tinblock))
c
c     zero out the block averages for each lambda window
c
      do i = 1, tinbin
         tinbcount(i) = 0
         do j = 1, tinblock
            tilmdadedl(i,j) = 0.0d0
            tilmdadedlstd(i,j) = 0.0d0
         end do
      end do
c
c     start the schedule at the first window and map the sublambdas
c
      tibin = 1
      tilmda = tilmdalist(1)
      call mapsublmda (tilmda)
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
      real*8 avg,std
c
c
c     a trailing partial window is left unsampled when the step
c     count is not an exact multiple of the number of windows
c
      if (tibin .gt. tinbin)  return
c
c     find the position of this step within the current window
c
      tistep = mod(istep-1,tiwindow) + 1
c
c     nothing is stored while the window is equilibrating
c
      if (tistep .gt. tinequil) then
         tiprod = tistep - tinequil
         tidedllist(mod(tiprod-1,tinstepavg)+1) = dedl
c
c     reduce a full block into its average and deviation
c
         if (mod(tiprod,tinstepavg) .eq. 0) then
            call avgstd (tidedllist,1,tinstepavg,avg,std)
            if (tinbcount(tibin) .lt. tinblock) then
               tinbcount(tibin) = tinbcount(tibin) + 1
               tilmdadedl(tibin,tinbcount(tibin)) = avg
               tilmdadedlstd(tibin,tinbcount(tibin)) = std
            end if
         end if
      end if
c
c     move on to the next lambda window at the window boundary
c
      if (tistep .eq. tiwindow)  call tischedule
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
         call mapsublmda (tilmda)
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine tiprint  --  output of block averaged dU/dl  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "tiprint" writes the block averaged lambda derivatives and
c     their deviations to an external file, one row per block, as a
c     whitespace delimited table
c
c
      subroutine tiprint
      use dlmda
      use files
      use iounit
      use thrmint
      implicit none
      integer b,w
      integer iti
      integer freeunit
      integer trimtext
      real*8 lam
      character*240 tifile
c
c
c     return if thermodynamic integration was never initialized
c
      if (.not. use_ti)  return
      if (.not. allocated(tilmdadedl))  return
c
c     open a new file to hold the block averages
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
      write (iti,20)  tinbin,tinstepavg,tieqratio,tiwindow,tinequil
   20 format ('# tinbin ',i0,' tinstepavg ',i0,' tieqratio ',f12.6,
     &           ' tiwindow ',i0,' tinequil ',i0)
c
c     record the full schedule, so the file still describes every
c     window when an interrupted run leaves some of them empty
c
      do w = 1, tinbin
         write (iti,30)  w,tilmdalist(w)
   30    format ('# schedule ',i0,1x,f12.8)
      end do
      write (iti,40)
   40 format ('# window lambda block dedl dedlstd')
c
c     write the block averages for each lambda window in turn
c
      do w = 1, tinbin
         lam = tilmdalist(w)
         do b = 1, tinbcount(w)
            write (iti,50)  w,lam,b,tilmdadedl(w,b),
     &                      tilmdadedlstd(w,b)
   50       format (i6,1x,f12.8,1x,i6,1p,1x,e20.10,1x,e20.10)
         end do
      end do
      close (unit=iti)
c
c     report the name of the file holding the block averages
c
      write (iout,60)  tifile(1:trimtext(tifile))
   60 format (/,' TI  --  dU/dlambda Block Averages Written To  ',a)
      return
      end
