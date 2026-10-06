c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine inittidyn  --  set up TI window averaging  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "inittidyn" divides a dynamics run of "nstep" steps among the
c     lambda windows according to their requested time fractions,
c     sizes the block average accumulators, and puts the main lambda
c     at the start of the schedule
c
c
      subroutine inittidyn (nstep)
      implicit none
      integer nstep
c
c
c     lay out the lambda windows over the length of the run
c
      call initlmdawin (nstep)
c
c     size the accumulators to the schedule
c
      call settiblocks
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
c     in "lmdawinend" can hold, and allocates the recording arrays
c     to that exact length
c
c
      subroutine settiblocks
      use dlmda
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
      do i = 1, nlmdawin
         nw = lmdawinend(i) - istart
         ne = int(dble(nw) * lmdawinratio)
         nb = (nw-ne) / tinstepavg
         if (nb .eq. 0) then
            write (iout,10)  i,lmdawinlist(i),lmdawinfrac(i)
   10       format (/,' SETTIBLOCKS  --  LAMBDA-WINDOW',i5,
     &                 ' at lambda',f12.6,' with fraction',f12.6,
     &              /,'                  is shorter than TI-NSTEPAVG',
     &                 ' and will record no samples')
         end if
         tinbtot = tinbtot + nb
         istart = lmdawinend(i)
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
c     describes the window and block layout; the block averages
c     themselves are appended later by "saveti"
c
c
      subroutine prttihead
      use dlmda
      use files
      use iounit
      use thrmint
      implicit none
      integer iti
      integer freeunit
      integer trimtext
c
c
c     open a new file, keeping any output from a previous run
c
      iti = freeunit ()
      lmdasavefile = filename(1:leng)//'.ti'
      call version (lmdasavefile,'new')
      open (unit=iti,file=lmdasavefile,status='new')
c
c     start the file with the standard Tinker banner message
c
      call promo (iti)
c
c     write a header describing the window and block layout
c
      write (iti,10)
   10 format (/,' Thermodynamic Integration Parameters :')
      write (iti,20)  nlmdawin,tinstepavg,lmdawinratio,
     &                lmdawinend(nlmdawin),tinbtot
   20 format (/,' Number of Lambda Windows',9x,i10,
     &        /,' Steps per Block Average',10x,i10,
     &        /,' Equilibration Ratio Value',10x,f8.3,
     &        /,' Number of Dynamics Steps',9x,i10,
     &        /,' Total Number of Block Averages',3x,i10)
c
c     label the columns of the block averages appended by "saveti"
c
      write (iti,30)
   30 format (/,' Block Averaged Lambda Derivatives :')
      write (iti,40)
   40 format (/,4x,'Index',8x,'Lambda',15x,'dE/dL',15x,'StDev',/)
      close (unit=iti)
c
c     report the name of the file holding the block averages
c
      write (iout,50)  lmdasavefile(1:trimtext(lmdasavefile))
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
c     at a dynamics step and averages it into blocks over the
c     production portion of the current lambda window
c
c
      subroutine etidyn (istep)
      use dlmda
      use mutant
      use thrmint
      implicit none
      integer istep
      integer iprod
      real*8 avg,std
c
c
c     nothing is stored while the window is equilibrating, or once
c     the schedule has run out
c
      call lmdawinphase (istep,iprod)
      if (iprod .eq. 0)  return
      tidedllist(mod(iprod-1,tinstepavg)+1) = dedl
c
c     reduce a full block into its average and deviation, keeping
c     the lambda that produced it alongside the block itself
c
      if (mod(iprod,tinstepavg) .eq. 0) then
         call avgstd (tidedllist,1,tinstepavg,avg,std)
         if (tinbcount .lt. tinbtot) then
            tinbcount = tinbcount + 1
            tilmdahist(tinbcount) = lambda
            tilmdadedl(tinbcount) = avg
            tilmdadedlstd(tinbcount) = std
         end if
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
      open (unit=iti,file=lmdasavefile,status='old',position='append')
      do i = tinbsave+1, tinbcount
         write (iti,10)  i,tilmdahist(i),tilmdadedl(i),
     &                   tilmdadedlstd(i)
   10    format (i9,f14.8,1p,e20.10,e20.10)
      end do
      tinbsave = tinbcount
      close (unit=iti)
      return
      end
