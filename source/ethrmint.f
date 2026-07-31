c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
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
      tibin = 0
      tilmda = 1.0d0
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
      if (tibin .ge. tinbin)  return
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
            if (tinbcount(tibin+1) .lt. tinblock) then
               tinbcount(tibin+1) = tinbcount(tibin+1) + 1
               tilmdadedl(tibin+1,tinbcount(tibin+1)) = avg
               tilmdadedlstd(tibin+1,tinbcount(tibin+1)) = std
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
c     "tischedule" steps the main lambda to the value of the next
c     lambda window; both endpoints are sampled, so the tinbin
c     windows span [0,1] in tinbin-1 decrements
c
c
      subroutine tischedule
      use thrmint
      implicit none
c
c
c     decrement the main lambda and refresh the sublambda values
c
      tibin = tibin + 1
      tilmda = 1.0d0 - dble(tibin)/dble(tinbin-1)
      if (tilmda .lt. 0.0d0)  tilmda = 0.0d0
      call mapsublmda (tilmda)
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
      write (iti,30)
   30 format ('# window lambda block dedl dedlstd')
c
c     write the block averages for each lambda window in turn
c
      do w = 1, tinbin
         lam = 1.0d0 - dble(w-1)/dble(tinbin-1)
         if (lam .lt. 0.0d0)  lam = 0.0d0
         do b = 1, tinbcount(w)
            write (iti,40)  w-1,lam,b-1,tilmdadedl(w,b),
     &                      tilmdadedlstd(w,b)
   40       format (i6,1x,f12.8,1x,i6,1p,1x,e20.10,1x,e20.10)
         end do
      end do
      close (unit=iti)
c
c     report the name of the file holding the block averages
c
      write (iout,50)  tifile(1:trimtext(tifile))
   50 format (/,' TI  --  dU/dlambda Block Averages Written To  ',a)
      return
      end
