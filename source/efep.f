c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine efeptrial  --  energy at trial lambda values  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "efeptrial" computes the potential energy of the current
c     structure at each of a list of trial values of the main lambda
c     and at the native main lambda, as needed for free energy
c     perturbation between neighboring lambda windows
c
c     the native lambda is evaluated last, which reinstalls all of
c     the native lambda dependent parameters, and the induced dipoles
c     are returned to the values they had on entry; prediction of the
c     induced dipoles is assumed to be off, as is enforced for this
c     method by "predict", so no predictor history is changed
c
c
      subroutine efeptrial (ntrial,ltrial,etrial,enative)
      use atoms
      use mutant
      use polar
      use potent
      implicit none
      integer i
      integer ntrial
      real*8 energy
      real*8 enative
      real*8 lnative
      real*8 ltrial(*)
      real*8 etrial(*)
      real*8, allocatable :: udirold(:,:)
      real*8, allocatable :: udirpold(:,:)
      real*8, allocatable :: uindold(:,:)
      real*8, allocatable :: uinpold(:,:)
      logical savepol
c
c
c     perform dynamic allocation of some local arrays
c
      savepol = use_polar
      if (savepol) then
         allocate (udirold(3,n))
         allocate (udirpold(3,n))
         allocate (uindold(3,n))
         allocate (uinpold(3,n))
      end if
c
c     save the native main lambda and the native induced dipoles
c
      lnative = lambda
      if (savepol) then
         udirold = udir
         udirpold = udirp
         uindold = uind
         uinpold = uinp
      end if
c
c     get the potential energy at each of the trial lambda values
c
      do i = 1, ntrial
         lambda = ltrial(i)
         etrial(i) = energy ()
      end do
c
c     get the native energy last to reinstall the native parameters
c
      lambda = lnative
      enative = energy ()
c
c     restore the native induced dipoles
c
      if (savepol) then
         udir = udirold
         udirp = udirpold
         uind = uindold
         uinp = uinpold
      end if
c
c     perform deallocation of some local arrays
c
      if (savepol) then
         deallocate (udirold)
         deallocate (udirpold)
         deallocate (uindold)
         deallocate (uinpold)
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine initfepdyn  --  set up FEP window sampling  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "initfepdyn" divides a dynamics run of "nstep" steps among the
c     lambda windows, counts the trial energies that will be found
c     over the production portion of each window, and sizes the
c     arrays that record them
c
c
      subroutine initfepdyn (nstep)
      use dlmda
      use fep
      use iounit
      implicit none
      integer i
      integer nstep
      integer nw,ne,ns
c
c
c     lay out the lambda windows over the length of the run
c
      call initlmdawin (nstep)
c
c     count the samples of each window; a sample holds the native
c     energy and that of each neighbor, of which the first and
c     the last window have only one
c
      nfeptot = 0
      do i = 1, nlmdawin
         call lmdawinsize (i,nw,ne)
         ns = (nw-ne) / fepintv
         if (ns .eq. 0) then
            write (iout,10)  i,lmdawinlist(i),lmdawinfrac(i)
   10       format (/,' INITFEPDYN  --  LAMBDA-WINDOW',i5,
     &                 ' at lambda',f12.6,' with fraction',f12.6,
     &              /,'                 is shorter than',
     &                 ' LAMBDA-INTERVAL and will record no samples')
         end if
         if (i.eq.1 .or. i.eq.nlmdawin) then
            nfeptot = nfeptot + 2*ns
         else
            nfeptot = nfeptot + 3*ns
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(fepstep))  deallocate (fepstep)
      if (allocated(fepene))  deallocate (fepene)
      if (allocated(feplmda))  deallocate (feplmda)
      if (allocated(feptrial))  deallocate (feptrial)
      allocate (fepstep(max(1,nfeptot)))
      allocate (fepene(max(1,nfeptot)))
      allocate (feplmda(max(1,nfeptot)))
      allocate (feptrial(max(1,nfeptot)))
c
c     zero out the energies recorded along the schedule
c
      do i = 1, max(1,nfeptot)
         fepstep(i) = 0
         fepene(i) = 0.0d0
         feplmda(i) = 0.0d0
         feptrial(i) = 0.0d0
      end do
      nfep = 0
      nfepsave = 0
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine prtfephead  --  start the trial energy file  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "prtfephead" creates the external file that collects the trial
c     lambda energies and writes the fixed header that describes the
c     lambda windows and the sampling; the energies themselves are
c     appended later by "savefep"
c
c
      subroutine prtfephead
      use bath
      use dlmda
      use fep
      use files
      use iounit
      implicit none
      integer i
      integer ifep
      integer freeunit
      integer trimtext
c
c
c     open a new file, keeping any output from a previous run
c
      ifep = freeunit ()
      lmdasavefile = filename(1:leng)//'.fep'
      call version (lmdasavefile,'new')
      open (unit=ifep,file=lmdasavefile,status='new')
c
c     start the file with the standard Tinker banner message
c
      call promo (ifep)
c
c     write a header describing the windows and the sampling
c
      write (ifep,10)
   10 format (/,' Free Energy Perturbation Parameters :')
      write (ifep,20)  nlmdawin,fepintv,lmdawinratio,
     &                 lmdawinend(nlmdawin),nfeptot,kelvin
   20 format (/,' Number of Lambda Windows',9x,i10,
     &        /,' Steps between Samples',12x,i10,
     &        /,' Equilibration Ratio Value',10x,f8.3,
     &        /,' Number of Dynamics Steps',9x,i10,
     &        /,' Total Number of Energy Values',4x,i10,
     &        /,' Simulation Temperature',13x,f8.2)
c
c     list the main lambda and the final step of each window
c
      write (ifep,30)
   30 format (/,' Lambda Windows :',
     &        //,4x,'Window',8x,'Lambda',5x,'Last Step',/)
      do i = 1, nlmdawin
         write (ifep,40)  i,lmdawinlist(i),lmdawinend(i)
   40    format (i10,f14.8,i14)
      end do
c
c     label the columns of the energies appended by "savefep"
c
      write (ifep,50)
   50 format (/,' Potential Energies at Trial Lambda Values :')
      write (ifep,60)
   60 format (/,7x,'Step',4x,'Sim Lambda',2x,'Trial Lambda',
     &           14x,'Energy',/)
      close (unit=ifep)
c
c     report the name of the file holding the trial energies
c
      write (iout,70)  lmdasavefile(1:trimtext(lmdasavefile))
   70 format (/,' FEP  --  Trial Lambda Energies Written To  ',a)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine efepdyn  --  free energy perturbation sampling  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "efepdyn" finds the potential energy at the main lambda of the
c     current window and at those of its neighboring windows, at
c     regular intervals over the production portion of the window,
c     and records the values for output
c
c
      subroutine efepdyn (istep)
      use dlmda
      use fep
      use mutant
      implicit none
      integer i
      integer istep
      integer iprod
      integer ntrial
      real*8 enative
      real*8 ltrial(2)
      real*8 etrial(2)
c
c
c     nothing is found while the window is equilibrating, once the
c     schedule has run out, or between the sampling intervals
c
      call lmdawinphase (istep,iprod)
      if (iprod .eq. 0)  return
      if (mod(iprod,fepintv) .ne. 0)  return
c
c     the trial values are those of the windows before and after
c
      ntrial = 0
      if (lmdawin .gt. 1) then
         ntrial = ntrial + 1
         ltrial(ntrial) = lmdawinlist(lmdawin-1)
      end if
      if (lmdawin .lt. nlmdawin) then
         ntrial = ntrial + 1
         ltrial(ntrial) = lmdawinlist(lmdawin+1)
      end if
      if (nfep+ntrial+1 .gt. nfeptot)  return
c
c     get the energy at the native and at the trial lambda values
c
      call efeptrial (ntrial,ltrial,etrial,enative)
c
c     record the native energy, then that of each trial lambda
c
      nfep = nfep + 1
      fepstep(nfep) = istep
      fepene(nfep) = enative
      feplmda(nfep) = lambda
      feptrial(nfep) = lambda
      do i = 1, ntrial
         nfep = nfep + 1
         fepstep(nfep) = istep
         fepene(nfep) = etrial(i)
         feplmda(nfep) = lambda
         feptrial(nfep) = ltrial(i)
      end do
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine savefep  --  append new trial lambda energies  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "savefep" appends to the external file the trial lambda
c     energies recorded since the previous call, so that a running
c     simulation can be monitored as the lambda windows are completed
c
c
      subroutine savefep
      use dlmda
      use fep
      implicit none
      integer i
      integer ifep
      integer freeunit
c
c
c     return if free energy perturbation was never initialized
c
      if (.not. use_fep)  return
      if (.not. allocated(fepene))  return
c
c     skip the file entirely when no energy has been recorded since
c     the last time the values were written out
c
      if (nfep .le. nfepsave)  return
c
c     append the energies recorded since the previous call
c
      ifep = freeunit ()
      open (unit=ifep,file=lmdasavefile,status='old',position='append')
      do i = nfepsave+1, nfep
         write (ifep,10)  fepstep(i),feplmda(i),feptrial(i),fepene(i)
   10    format (i11,2f14.8,f20.8)
      end do
      nfepsave = nfep
      close (unit=ifep)
      return
      end
