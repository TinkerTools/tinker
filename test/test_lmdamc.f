c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine test_lmdamc  --  Monte Carlo lambda moves  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "test_lmdamc" checks the Monte Carlo propagation of the main
c     lambda for the adaptive bias methods; the keywords, the phase
c     and frequency of the trials, the size and reflection of a move,
c     the energy and the induced dipoles kept by an accepted move and
c     restored by a rejected one, and short dynamics runs
c
c
      subroutine test_lmdamc
      implicit none
      logical skiptest
      character*(*) tname
      parameter (tname='test_lmdamc')
c
c
      if (skiptest(tname,'lmdamc'))  return
      call test_lmdamc_keys
      call test_lmdamc_particle
      call test_lmdamc_move ('abf')
      call test_lmdamc_move ('ost')
      call test_lmdamc_move ('meta')
      call test_lmdamc_trial
      call test_lmdamc_deposit
      call test_lmdamc_stale ('abf')
      call test_lmdamc_stale ('ost')
      call test_lmdamc_stalemeta
      call test_lmdamc_nopol
      call test_lmdamc_run
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_lmdamc_keys  --  keyword settings  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_lmdamc_keys" checks the defaults of the sample interval
c     and of the Monte Carlo lambda move, that the lambda particle is
c     the default propagator, and the values read from the keywords
c
c
      subroutine test_lmdamc_keys
      use dlmda
      use uprior
      implicit none
c
c
c     an adaptive bias run moves the lambda particle by default
c
      call pushdir ('file/mutate')
      call loadfix_keyadd ('ionwat','179_ionwat_ast_l05.key',
     &                     'lambda-mode abf')
      call assert_logical (use_lmdamc,.false.,'lmdamc default off')
      call assert_int (lmdaintv,100,'lmdamc default interval')
      call assert_int (lmdanpa,30,'lmdamc default propagate steps')
      call assert_int (lmdanpb,30,'lmdamc default hold steps')
      call assert_int (lmdanpc,40,'lmdamc default average steps')
      call assert_int (lmdamctrial,10,'lmdamc default trial')
      call assert_real (lmdamcmove,0.02d0,0.0d0,'lmdamc default move')
      call final
c
c     the keywords turn on the moves and set their rate and size
c
      call loadfix_keyadd ('ionwat','lmdamc_ionwat.key',
     &                     'lambda-mode abf')
      call assert_logical (use_lmdamc,.true.,'lmdamc keyword on')
      call assert_int (lmdaintv,4,'lmdamc keyword interval')
      call assert_int (lmdanpa,2,'lmdamc keyword propagate steps')
      call assert_int (lmdamctrial,1,'lmdamc keyword trial')
      call assert_real (lmdamcmove,0.05d0,0.0d0,'lmdamc keyword move')
      call assert_int (nlmdamctry,0,'lmdamc starts with no trials')
      call assert_int (nlmdamcacc,0,'lmdamc starts with no accepts')
      call final
c
c     a move larger than the lambda range is cut back to the range
c
      call loadfix_keyadd ('ionwat','lmdamc_ionwat.key',
     &                     'lambda-mc-move -3.0')
      call assert_real (lmdamcmove,1.0d0,0.0d0,'lmdamc move capped')
      call final
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_lmdamc_particle  --  particle left idle  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_lmdamc_particle" checks that the interval sampling no
c     longer propagates the lambda particle once the lambda is moved
c     by Monte Carlo, while the samples are still deposited
c
c
      subroutine test_lmdamc_particle
      use bath
      use dlmda
      use mutant
      implicit none
      integer istep
      real*8 theta0
c
c
      call pushdir ('file/mutate')
      call loadfix_keyadd ('ionwat','lmdamc_ionwat.key',
     &                     'lambda-mode abf')
      kelvin = 298.0d0
      lmdastep = 0
      lambda = 0.5d0
      call lmdathetainv (lambda,lmdatheta)
      theta0 = lmdatheta
      lmdavtheta = 0.0d0
      lmdadfdl = 0.0d0
      do istep = 1, 2*lmdaintv
         dedl = 40.0d0
         call elmdadyn
      end do
      call assert_real (lambda,0.5d0,0.0d0,'lmdamc particle lambda')
      call assert_real (lmdatheta,theta0,0.0d0,'lmdamc particle theta')
      call assert_real (lmdavtheta,0.0d0,0.0d0,
     &                  'lmdamc particle velocity')
      call assert_int (nlmdahist,2,'lmdamc particle still deposits')
c
c     the same steps move the lambda once the particle is in charge
c
      use_lmdamc = .false.
      do istep = 1, lmdaintv
         dedl = 40.0d0
         call elmdadyn
      end do
      call assert_logical (lambda.ne.0.5d0,.true.,
     &                     'lmdamc particle moves when propagated')
      kelvin = 0.0d0
      call final
      call popdir
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine mcload  --  set up a Monte Carlo system  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "mcload" sets up the ion and water fixture for the adaptive
c     bias method named by "mode" with Monte Carlo lambda moves, as
c     the "dynamic" program would, and adds a bias so that the bias
c     change of a move is not zero; "epot" returns the energy of the
c     current structure including that bias
c
c
      subroutine mcload (mode,epot)
      use atoms
      use dlmda
      use mutant
      use ost
      implicit none
      integer i
      real*8 epot
      real*8, allocatable :: derivs(:,:)
      character*(*) mode
c
c
      call loadfix_keyadd ('ionwat','lmdamc_ionwat.key',
     &                     'lambda-mode '//mode)
      use_ostdyn = use_ost
      use_metadyn = use_meta
      use_abfdyn = use_abf
      call sd_mdinit (298.0d0)
      allocate (derivs(3,n))
c
c     give the abf mean force a slope, or deposit one gaussian
c
      if (mode .eq. 'abf') then
         do i = 1, nlmda
            lmdafmean(i) = 30.0d0 * dble(i-1) * wlmda
         end do
      else if (mode .eq. 'ost') then
         call gradient (epot,derivs)
         hbias = 2.0d0
         lmdaavg = lambda
         dedlavg = dedl
         call ostdeposit
      else
         hbias = 2.0d0
         lmdaavg = lambda
         call metadeposit
      end if
      call gradient (epot,derivs)
      deallocate (derivs)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine mcenergy  --  energy with the lambda bias  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "mcenergy" finds the potential energy plus the adaptive bias
c     at the current lambda, independently of the dynamics step
c
c
      subroutine mcenergy (mode,epot)
      use atoms
      use dlmda
      use mutant
      use virial
      implicit none
      real*8 epot
      real*8 energy
      real*8 vbias,dvdl
      real*8 virsave(3,3)
      real*8, allocatable :: derivs(:,:)
      character*(*) mode
c
c
      if (mode .eq. 'abf') then
         call efreelmda (lambda,vbias,dvdl)
         epot = energy () - vbias
      else if (mode .eq. 'ost') then
         allocate (derivs(3,n))
         virsave = vir
         call gradient (epot,derivs)
         vir = virsave
         deallocate (derivs)
      else
         call emetabias (lambda,vbias,dvdl)
         epot = energy () + vbias
      end if
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_lmdamc_move  --  accepted and rejected  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_lmdamc_move" forces a rejected and then an accepted
c     Monte Carlo lambda move for one adaptive bias method; the
c     rejected move must return the lambda, the energy, the bias, the
c     induced dipoles and the installed multipole parameters to what
c     they were, and the accepted move must leave the energy, the
c     induced dipoles and the multipoles of the new lambda; neither
c     may touch the induced dipole predictor, which is turned off
c
c
      subroutine test_lmdamc_move (mode)
      use atoms
      use deriv
      use dlmda
      use ielscf
      use mpole
      use mutant
      use polar
      use uprior
      use virial
      implicit none
      real*8 epot,e0,eref
      real*8 lold,lmap,dldth
      real*8 dedl0,vbias0
      real*8 vir0(3,3)
      real*8, allocatable :: uind0(:,:)
      real*8, allocatable :: uinp0(:,:)
      real*8, allocatable :: udir0(:,:)
      real*8, allocatable :: udirp0(:,:)
      real*8, allocatable :: pole0(:,:)
      real*8, allocatable :: rpole0(:,:)
      real*8, allocatable :: polar0(:)
      real*8, allocatable :: desum0(:,:)
      character*(*) mode
      character*11 tag
c
c
      tag = 'lmdamc '//mode
      call pushdir ('file/mutate')
      call mcload (mode,e0)
c
c     the predictor that the keyfile asks for is turned off
c
      call assert_logical (use_pred,.false.,trim(tag)//' no predictor')
      call assert_logical (use_ielscf,.false.,trim(tag)//' no iel scf')
      call assert_int (maxualt,0,trim(tag)//' no predictor storage')
      call assert_int (nualt,0,trim(tag)//' no predictor history')
c
c     save the native state of the step that was just evaluated
c
      allocate (uind0(3,n))
      allocate (uinp0(3,n))
      allocate (udir0(3,n))
      allocate (udirp0(3,n))
      allocate (pole0(maxpole,n))
      allocate (rpole0(maxpole,n))
      allocate (polar0(n))
      allocate (desum0(3,n))
      uind0 = uind
      uinp0 = uinp
      udir0 = udir
      udirp0 = udirp
      pole0 = pole
      rpole0 = rpole
      polar0 = polarity
      desum0 = desum
      vir0 = vir
      dedl0 = dedl
      vbias0 = lmdavbias
      lold = lambda
c
c     an old energy far below the trial energy rejects the move
c
      lmdastep = 1
      epot = -1.0d30
      call lmdamonte (epot)
      call assert_int (nlmdamctry,1,trim(tag)//' reject trial count')
      call assert_int (nlmdamcacc,0,trim(tag)//' reject accept count')
      call assert_real (epot,-1.0d30,0.0d0,trim(tag)//' reject epot')
      call assert_real (lambda,lold,0.0d0,trim(tag)//' reject lambda')
      call assert_real (dedl,dedl0,0.0d0,trim(tag)//' reject dedl')
      call assert_real (lmdavbias,vbias0,0.0d0,
     &                    trim(tag)//' reject bias')
      call assert_array2 (uind,uind0,3,n,0.0d0,
     &                    trim(tag)//' reject uind')
      call assert_array2 (uinp,uinp0,3,n,0.0d0,
     &                    trim(tag)//' reject uinp')
      call assert_array2 (udir,udir0,3,n,0.0d0,
     &                    trim(tag)//' reject udir')
      call assert_array2 (udirp,udirp0,3,n,0.0d0,
     &                    trim(tag)//' reject udirp')
      call assert_array2 (pole,pole0,maxpole,n,0.0d0,
     &                    trim(tag)//' reject pole')
      call assert_array2 (rpole,rpole0,maxpole,n,0.0d0,
     &                    trim(tag)//' reject rpole')
      call assert_array1 (polarity,polar0,n,0.0d0,
     &                    trim(tag)//' reject polarity')
      call assert_array2 (vir,vir0,3,3,0.0d0,
     &                    trim(tag)//' reject virial')
      call assert_array2 (desum,desum0,3,n,0.0d0,
     &                    trim(tag)//' reject desum')
      call assert_int (nualt,0,trim(tag)//' reject predictor history')
c
c     the native energy is still that of the step after the rejection
c
      call mcenergy (mode,eref)
      call assert_real (eref,e0,1.0d-6,
     &                    trim(tag)//' reject native energy')
c
c     an old energy far above the trial energy accepts the move
c
      lmdastep = 1
      epot = e0 + 1.0d6
      call lmdamonte (epot)
      call assert_int (nlmdamctry,2,trim(tag)//' accept trial count')
      call assert_int (nlmdamcacc,1,trim(tag)//' accept accept count')
      call assert_logical (lambda.ne.lold,.true.,
     &                     trim(tag)//' accept changes lambda')
      call assert_logical (abs(lambda-lold).le.lmdamcmove,.true.,
     &                     trim(tag)//' accept within move')
      call assert_array2 (vir,vir0,3,3,0.0d0,
     &                    trim(tag)//' accept virial')
      call lmdathetamap (lmdatheta,lmap,dldth)
      call assert_real (lmap,lambda,1.0d-10,trim(tag)//' accept theta')
      call assert_real (lmdavtheta,0.0d0,0.0d0,
     &                  trim(tag)//' accept theta velocity')
      call assert_int (nualt,0,trim(tag)//' accept predictor history')
c
c     the accepted energy, induced dipoles and multipoles are those
c     of an independent evaluation at the new lambda
c
      uind0 = uind
      uinp0 = uinp
      rpole0 = rpole
      pole0 = pole
      call mcenergy (mode,eref)
      call assert_real (epot,eref,1.0d-6,trim(tag)//' accept energy')
      call assert_array2 (uind,uind0,3,n,1.0d-8,
     &                    trim(tag)//' accept uind')
      call assert_array2 (uinp,uinp0,3,n,1.0d-8,
     &                    trim(tag)//' accept uinp')
      call assert_array2 (pole,pole0,maxpole,n,0.0d0,
     &                    trim(tag)//' accept pole')
      call assert_array2 (rpole,rpole0,maxpole,n,0.0d0,
     &                    trim(tag)//' accept rpole')
      deallocate (uind0)
      deallocate (uinp0)
      deallocate (udir0)
      deallocate (udirp0)
      deallocate (pole0)
      deallocate (rpole0)
      deallocate (polar0)
      deallocate (desum0)
      call sd_final
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_lmdamc_trial  --  phase, rate and walls  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_lmdamc_trial" checks that ost and abf try a lambda move
c     only in the propagation phase of the sample interval while
c     metadynamics tries at any step, that the mean trial rate is
c     honored, and that a move from either end of the lambda range is
c     reflected back into the range
c
c
      subroutine test_lmdamc_trial
      use bath
      use dlmda
      use mutant
      implicit none
      integer i
      real*8 epot,e0
      real*8 lmin,lmax
c
c
c     no trial is made once the propagation phase is over
c
      call pushdir ('file/mutate')
      call mcload ('abf',e0)
      do i = 1, 10
         lmdastep = 3
         epot = e0 + 1.0d6
         call lmdamonte (epot)
         lmdastep = 4
         epot = e0 + 1.0d6
         call lmdamonte (epot)
      end do
      call assert_int (nlmdamctry,0,'lmdamc abf holds lambda')
      call assert_real (lambda,0.5d0,0.0d0,'lmdamc abf held lambda')
      call assert_real (epot,e0+1.0d6,0.0d0,'lmdamc abf held epot')
c
c     a long mean time between trials leaves the lambda alone
c
      lmdamctrial = 1000000000
      do i = 1, 10
         lmdastep = 1
         epot = e0 + 1.0d6
         call lmdamonte (epot)
      end do
      call assert_int (nlmdamctry,0,'lmdamc rare trials')
      lmdamctrial = 1
c
c     moves from the lambda end points are reflected into the range
c
      lmin = 1.0d0
      lmax = 0.0d0
      do i = 1, 10
         lambda = 0.0d0
         lmdastep = 1
         epot = 1.0d30
         call lmdamonte (epot)
         lmin = min(lmin,lambda)
         lmax = max(lmax,lambda)
      end do
      call assert_int (nlmdamcacc,10,'lmdamc low wall accepts')
      call assert_logical (lmin.ge.0.0d0,.true.,'lmdamc low wall min')
      call assert_logical (lmax.le.lmdamcmove,.true.,
     &                     'lmdamc low wall max')
      lmin = 1.0d0
      lmax = 0.0d0
      do i = 1, 10
         lambda = 1.0d0
         lmdastep = 1
         epot = 1.0d30
         call lmdamonte (epot)
         lmin = min(lmin,lambda)
         lmax = max(lmax,lambda)
      end do
      call assert_int (nlmdamcacc,20,'lmdamc high wall accepts')
      call assert_logical (lmax.le.1.0d0,.true.,'lmdamc high wall max')
      call assert_logical (lmin.ge.1.0d0-lmdamcmove,.true.,
     &                     'lmdamc high wall min')
      call sd_final
c
c     metadynamics tries a move in every phase of the interval
c
      call mcload ('meta',e0)
      do i = 1, 5
         lmdastep = 3
         epot = -1.0d30
         call lmdamonte (epot)
         lmdastep = 4
         epot = -1.0d30
         call lmdamonte (epot)
      end do
      call assert_int (nlmdamctry,10,'lmdamc meta tries any step')
      call assert_int (nlmdamcacc,0,'lmdamc meta forced rejects')
      call sd_final
      call popdir
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_lmdamc_deposit  --  bias added before move  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_lmdamc_deposit" checks that a lambda move tried on the
c     step that deposits a metadynamics gaussian takes its old bias
c     from the history with that gaussian in place; the energy of the
c     step was found before the deposit, and a move judged against it
c     would see the whole of the new gaussian as an energy barrier
c
c
      subroutine test_lmdamc_deposit
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 epot,e0,eref
      real*8 lold
c
c
c     a tall gaussian lands on the current lambda after the forces
c     of the step were found, as at the end of a sample interval
c
      call pushdir ('file/mutate')
      call mcload ('meta',e0)
      hbias = 1000.0d0
      lmdaavg = lambda
      call metadeposit
c
c     any move within two widths of the peak is then downhill
c
      lmdamcmove = 2.0d0 * wlmda
      lold = lambda
      lmdastep = 1
      epot = e0
      call lmdamonte (epot)
      call assert_int (nlmdamctry,1,'lmdamc deposit trial made')
      call assert_int (nlmdamcacc,1,'lmdamc deposit move accepted')
      call assert_logical (lambda.ne.lold,.true.,
     &                     'lmdamc deposit changes lambda')
      call mcenergy ('meta',eref)
      call assert_real (epot,eref,1.0d-6,'lmdamc deposit energy')
      call sd_final
      call popdir
      return
      end
c
c
c     #######################################################
c     ##                                                   ##
c     ##  subroutine mctrial  --  one forced lambda trial  ##
c     ##                                                   ##
c     #######################################################
c
c
c     "mctrial" tries one Monte Carlo lambda move from the lambda
c     "lmda" on step "kstep" of the sample interval, handing over
c     "efake" as the energy of the step, and reports whether the move
c     was accepted along with the energy that comes back
c
c
      subroutine mctrial (lmda,kstep,efake,epot,took)
      use dlmda
      use mutant
      implicit none
      integer kstep
      integer nacc
      real*8 lmda,efake,epot
      logical took
c
c
      lambda = lmda
      lmdastep = kstep
      epot = efake
      nacc = nlmdamcacc
      call lmdamonte (epot)
      took = (nlmdamcacc .gt. nacc)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_lmdamc_stale  --  outdated step energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_lmdamc_stale" checks for abf or ost that the old energy
c     of a lambda move is found again when the energy of the step no
c     longer describes the system, after a scaling barostat, a rigid
c     body step or, for ost, a deposit on the same step, and that the
c     energy of the step is trusted and nothing is found otherwise
c
c     a large constant mean force makes the energy with its bias
c     rise steeply with lambda, so that at a very low temperature
c     every move from lambda of one is downhill and accepted, and
c     every move from lambda of zero is uphill and rejected, as long
c     as the old energy is the true one; a false energy handed over
c     for the step shows whether it was used
c
c
      subroutine test_lmdamc_stale (mode)
      use atoms
      use bath
      use deriv
      use dlmda
      use mdstuf
      use mutant
      use polar
      use virial
      implicit none
      integer i
      real*8 epot,e0,eref
      real*8 dedl0,vbias0
      real*8 vir0(3,3)
      real*8, allocatable :: uind0(:,:)
      real*8, allocatable :: udir0(:,:)
      real*8, allocatable :: desum0(:,:)
      logical took
      character*(*) mode
      character*11 tag
c
c
      tag = 'lmdamc '//mode
      call pushdir ('file/mutate')
      call mcload (mode,e0)
      do i = 1, nlmda
         lmdafmean(i) = -1.0d5
      end do
      kelvin = 1.0d-6
      allocate (uind0(3,n))
      allocate (udir0(3,n))
      allocate (desum0(3,n))
c
c     after a scaling barostat a downhill move is accepted whatever
c     energy the step hands over, and it returns the true new energy
c
      isobaric = .true.
      barostat = 'BUSSI'
      call mctrial (1.0d0,1,-1.0d30,epot,took)
      call assert_logical (took,.true.,
     &                     trim(tag)//' scaled downhill accepted')
      call assert_logical (lambda.lt.1.0d0,.true.,
     &                     trim(tag)//' scaled downhill moves')
      call mcenergy (mode,eref)
      call assert_real (epot,eref,1.0d-6,
     &                  trim(tag)//' scaled accepted energy')
c
c     an uphill move is rejected, and the rejection leaves the energy
c     of the step and the state that the step left untouched
c
      lambda = 0.0d0
      uind0 = uind
      udir0 = udir
      desum0 = desum
      vir0 = vir
      dedl0 = dedl
      vbias0 = lmdavbias
      call mctrial (0.0d0,1,1.0d30,epot,took)
      call assert_logical (took,.false.,
     &                     trim(tag)//' scaled uphill rejected')
      call assert_real (lambda,0.0d0,0.0d0,
     &                  trim(tag)//' scaled reject lambda')
      call assert_real (epot,1.0d30,0.0d0,
     &                  trim(tag)//' scaled reject epot')
      call assert_real (dedl,dedl0,0.0d0,
     &                  trim(tag)//' scaled reject dedl')
      call assert_real (lmdavbias,vbias0,0.0d0,
     &                  trim(tag)//' scaled reject bias')
      call assert_array2 (uind,uind0,3,n,0.0d0,
     &                    trim(tag)//' scaled reject uind')
      call assert_array2 (udir,udir0,3,n,0.0d0,
     &                    trim(tag)//' scaled reject udir')
      call assert_array2 (vir,vir0,3,3,0.0d0,
     &                    trim(tag)//' scaled reject virial')
      call assert_array2 (desum,desum0,3,n,0.0d0,
     &                    trim(tag)//' scaled reject desum')
c
c     with nothing moved the energy of the step is trusted, so the
c     same false energies now decide both moves the other way
c
      isobaric = .false.
      call mctrial (1.0d0,1,-1.0d30,epot,took)
      call assert_logical (took,.false.,
     &                     trim(tag)//' unmoved low energy trusted')
      call mctrial (0.0d0,1,1.0d30,epot,took)
      call assert_logical (took,.true.,
     &                     trim(tag)//' unmoved high energy trusted')
c
c     a rigid body step finds its forces before it moves the bodies
c
      integrate = 'RIGIDBODY'
      call mctrial (1.0d0,1,-1.0d30,epot,took)
      call assert_logical (took,.true.,
     &                     trim(tag)//' rigid body energy found again')
      integrate = 'BEEMAN'
c
c     on a step that deposits, ost finds its energy again while abf
c     only has to take its bias from the history as it now stands
c
      lmdanpa = lmdaintv
      call mctrial (1.0d0,lmdaintv,-1.0d30,epot,took)
      call assert_logical (took,mode.eq.'ost',
     &                     trim(tag)//' deposit step old energy')
      deallocate (uind0)
      deallocate (udir0)
      deallocate (desum0)
      call sd_final
      call popdir
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_lmdamc_stalemeta  --  stale energy, meta  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_lmdamc_stalemeta" makes the same check for metadynamics,
c     where a tall gaussian centered on the current lambda makes every
c     nearby move downhill once the old energy is the true one
c
c
      subroutine test_lmdamc_stalemeta
      use bath
      use dlmda
      use mutant
      use ost
      implicit none
      real*8 epot,e0,eref
      real*8 lpeak
      logical took
c
c
      call pushdir ('file/mutate')
      call mcload ('meta',e0)
      lpeak = lambda
      hbias = 1000.0d0
      lmdaavg = lpeak
      call metadeposit
      lmdamcmove = 2.0d0 * wlmda
c
c     after a scaling barostat the old energy is found again
c
      isobaric = .true.
      barostat = 'BUSSI'
      call mctrial (lpeak,3,-1.0d30,epot,took)
      call assert_logical (took,.true.,
     &                     'lmdamc meta scaled downhill accepted')
      call mcenergy ('meta',eref)
      call assert_real (epot,eref,1.0d-6,
     &                  'lmdamc meta scaled accepted energy')
c
c     with nothing moved the energy of the step is trusted
c
      isobaric = .false.
      call mctrial (lpeak,3,-1.0d30,epot,took)
      call assert_logical (took,.false.,
     &                     'lmdamc meta unmoved low energy trusted')
      call sd_final
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_lmdamc_nopol  --  fixed multipole model  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_lmdamc_nopol" checks that a rejected lambda move leaves
c     the global frame multipoles as the step left them for a model
c     with no polarization, where the multipole routine with a lambda
c     derivative leaves the unscaled multipoles in place and nothing
c     later in the step installs the scaled ones
c
c
      subroutine test_lmdamc_nopol
      use atoms
      use dlmda
      use mpole
      use mutant
      use output
      implicit none
      real*8 epot,e0
      real*8, allocatable :: derivs(:,:)
      real*8, allocatable :: rpole0(:,:)
      logical took
c
c
      call pushdir ('file/mutate')
      call loadfix_keyadd ('ionwat','lmdamc_nopol.key',
     &                     'lambda-mode abf')
      use_abfdyn = use_abf
      call sd_mdinit (298.0d0)
      allocate (derivs(3,n))
      allocate (rpole0(maxpole,n))
      call gradient (e0,derivs)
      rpole0 = rpole
c
c     an old energy far below the trial energy rejects the move;
c     the multipoles are kept only when the output that reads them
c     has been asked for
c
      ustcsave = .true.
      call mctrial (lambda,1,-1.0d30,epot,took)
      ustcsave = .false.
      call assert_logical (took,.false.,'lmdamc nopol move rejected')
      call assert_array2 (rpole,rpole0,maxpole,n,0.0d0,
     &                    'lmdamc nopol reject rpole')
      deallocate (derivs)
      deallocate (rpole0)
      call sd_final
      call popdir
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_lmdamc_run  --  short dynamics runs  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_lmdamc_run" runs short dynamics trajectories with Monte
c     Carlo lambda moves for each adaptive bias method, at constant
c     volume and under the Monte Carlo and the scaling barostats, and
c     checks that a run without a thermostat is refused
c
c
      subroutine test_lmdamc_run
      implicit none
c
c
      call mcrun ('abf','ionwat 60 1.0 0.02 2 298','montecarlo',
     &            'abf','nvt')
      call mcrun ('meta','ionwat 60 1.0 0.02 2 298','montecarlo',
     &            'meta','nvt')
      call mcrun ('ost','ionwat 60 1.0 0.02 2 298','montecarlo',
     &            'ost','nvt')
      call mcrun ('ost','ionwat 60 1.0 0.02 4 298 1.0','montecarlo',
     &            'ost','npt montecarlo')
      call mcrun ('abf','ionwat 60 1.0 0.02 4 298 1.0','bussi',
     &            'abf','npt bussi')
      call mcrun ('meta','ionwat 60 1.0 0.02 4 298 1.0','bussi',
     &            'meta','npt bussi')
      call mcrun ('ost','ionwat 60 1.0 0.02 4 298 1.0','berendsen',
     &            'ost','npt berendsen')
      call mcrun ('abf','ionwat 10 1.0 0.02 1','montecarlo',
     &            'none','nve')
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine mcrun  --  check one dynamics run  ##
c     ##                                                ##
c     ####################################################
c
c
c     "mcrun" makes one dynamics run in a scratch directory with
c     Monte Carlo lambda moves, and checks that it ends normally,
c     reports the acceptance of the moves and writes the history file
c     with the extension "ext"; an extension of "none" marks a run
c     without a thermostat, which has to stop with a message; "baro"
c     names the barostat put in the keyfile and "ens" the ensemble
c     in the label of each check
c
c
      subroutine mcrun (mode,args,baro,ext,ens)
      implicit none
      integer ist
      logical exist
      character*(*) mode
      character*(*) args,baro,ext,ens
      character*40 tag
      character*512 cmd
c
c
c     copy the fixture, pointing its key at the parameter file
c
      tag = 'lmdamc run '//trim(mode)//' '//trim(ens)
      call pushdir ('file/mutate')
      cmd = 'rm -rf ../lmdamc_run ; mkdir -p ../lmdamc_run ; '//
     &      'cp ionwat.xyz ../lmdamc_run/ ; '//
     &      'cp lmdamc_ionwat.key ../lmdamc_run/ionwat.key'
      call execute_command_line (cmd)
      call popdir
      call pushdir ('file/lmdamc_run')
      call tnist_append ('ionwat.key','lambda-mode '//mode)
      call tnist_append ('ionwat.key','lambda-mc-trial 2')
      call tnist_append ('ionwat.key','barostat '//baro)
      call tnist_append ('ionwat.key','volume-trial 3')
      call tnist_append ('ionwat.key','randomseed 123456789')
      call run_prog ('dynamic',args,'out.txt',ist)
      if (ist .ne. -1) then
         if (ext .eq. 'none') then
            cmd = 'grep -q "Require a Constant Temperature" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' refused')
         else
            call assert_int (ist,0,trim(tag)//' status')
            cmd = 'grep -q "Lambda MC Accept" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,0,trim(tag)//' acceptance line')
            cmd = 'grep -q "Illegal Value" out.txt'
            call execute_command_line (cmd,exitstat=ist)
            call assert_int (ist,1,trim(tag)//' finite energy')
            inquire (file='ionwat.'//trim(ext),exist=exist)
            call assert_logical (exist,.true.,trim(tag)//' history')
         end if
      end if
      call popdir
      call tnist_clean ('lmdamc_run')
      return
      end
