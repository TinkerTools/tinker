c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2009 by Chuanjie Wu & Jay William Ponder  ##
c     ##  COPYRIGHT (C) 2026 by Moses Chung & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate  --  set parameters for hybrid system  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate" constructs the hybrid hamiltonian for a specified
c     initial state, final state and mutation parameter "lambda"
c
c     note torsional and most electrostatics terms apply "lambda"
c     by directly scaling parameters, while vdw and repulsion energy
c     terms use soft core functions from the references cited below
c
c     literature references:
c
c     T. Steinbrecher, D. L. Mobley and D. A. Case, "Nonlinear Scaling
c     Schemes for Lennard-Jones Interactions in Free Energy
c     Calculations", Journal of Chemical Physics, 127, 214108 (2007)
c
c     D. Jiao, P. A. Golubkov, T. A. Darden and P. Ren, "Calculation
c     of Protein-Ligand Binding Free Energy by Using a Polarizable
c     Potential", PNAS, 105, 6290-6295 (2008)
c
c
      subroutine mutate
      use atoms
      use bndstr
      use dlmda
      use inform
      use iounit
      use katoms
      use keys
      use mutant
      use potent
      implicit none
      integer i,j,k,ihyb
      integer it0,it1
      integer igrp
      integer next,size
      integer ntbnd
      integer, allocatable :: list(:)
      integer, allocatable :: itbnd(:,:)
      character*4 leg(3)
      character*20 keyword
      character*240 record
      character*240 string
      data leg  / 'LIG2', 'VDWM', 'LIG1' /
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(imut))  deallocate (imut)
      if (allocated(type0))  deallocate (type0)
      if (allocated(class0))  deallocate (class0)
      if (allocated(type1))  deallocate (type1)
      if (allocated(class1))  deallocate (class1)
      if (allocated(mutg))  deallocate (mutg)
      if (allocated(fluxsc))  deallocate (fluxsc)
      if (allocated(pdltflx))  deallocate (pdltflx)
      allocate (imut(n))
      allocate (type0(n))
      allocate (class0(n))
      allocate (type1(n))
      allocate (class1(n))
      allocate (mutg(n))
      allocate (fluxsc(n))
      allocate (pdltflx(n))
c
c     perform dynamic allocation of some local arrays
c
      size = 40
      allocate (list(size))
      allocate (itbnd(2,nbond))
c
c     set defaults for lambda perturbation scaling values
c
      lambda = 1.0d0
      elambda = 1.0d0
      vlambda = 1.0d0
      tlambda = 1.0d0
      use_mainlmda = .false.
c
c     set defaults for lambda scaling for lambda derivatives
c
      setelambda = .false.
      setplambda = .false.
      setvlambda = .false.
      plambda = 1.0d0
c
c     set defaults for vdw coupling type and soft core vdw
c
      vcouple = 0
      scexp = 5.0d0
      scalphav = 0.7d0
      scexpm = 0.0d0
      scexpp = 2.0d0
c
c     set defaults for soft core Pauli repulsion
c
      scalphar = 1.0d0
      scexpr = 5.0d0
      scexprm = 3.0d0
      scexprp = 1.0d0
c
c     set defaults for soft core damped dispersion
c
      scrdsp = 2.0d0
      scexpd = 4.0d0
      scexpdq = 1.0d0
c
c     zero out number of hybrid atoms and mutated torsions
c
      nmut = 0
      nmutb = 0
      use_rel = .false.
      do i = 1, n
         mutg(i) = 0
         fluxsc(i) = 1.0d0
         pdltflx(i) = 0.0d0
      end do
      ntbnd = 0
      do i = 1, nbond
         itbnd(1,i) = 0
         itbnd(2,i) = 0
      end do
c
c     search keywords for free energy perturbation options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  lambda
            use_mainlmda = .true.
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  elambda
            setelambda = .true.
         else if (keyword(1:11) .eq. 'POL-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  plambda
            setplambda = .true.
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  vlambda
            setvlambda = .true.
         else if (keyword(1:13) .eq. 'VDW-SOFTCORE ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  scalphav,scexp,
     &                                     scexpm,scexpp
         else if (keyword(1:15) .eq. 'REPEL-SOFTCORE ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  scalphar,scexpr,
     &                                     scexprm,scexprp
         else if (keyword(1:14) .eq. 'DISP-SOFTCORE ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  scrdsp,scexpd,scexpdq
         else if (keyword(1:12) .eq. 'TORS-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  tlambda
         else if (keyword(1:15) .eq. 'VDW-ANNIHILATE ') then
            vcouple = 1
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:240)
            read (string,*,err=30,end=30)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mutg(ihyb) = 1
            type0(nmut) = it0
            type1(nmut) = it1
            class0(nmut) = atmcls(it0)
            class1(nmut) = atmcls(it1)
         else if (keyword(1:7).eq.'LIGAND ' .or.
     &            keyword(1:8).eq.'LIGAND1 ' .or.
     &            keyword(1:8).eq.'LIGAND2 ') then
            if (keyword(1:8) .eq. 'LIGAND2 ') then
               igrp = 2
            else
               igrp = 1
            end if
            do k = 1, size
               list(k) = 0
            end do
            string = record(next:240)
            read (string,*,err=10,end=10)  (list(k),k=1,size)
   10       continue
            k = 1
            do while (list(k) .ne. 0)
               if (list(k).gt.0 .and. list(k).le.n) then
                  j = list(k)
                  call setligand (j,igrp)
                  k = k + 1
               else
                  do j = max(1,abs(list(k))), min(n,abs(list(k+1)))
                     call setligand (j,igrp)
                  end do
                  k = k + 2
               end if
            end do
         else if (keyword(1:15) .eq. 'ROTATABLE-BOND ') then
            do k = 1, size
               list(k) = 0
            end do
            string = record(next:240)
            read (string,*,err=20,end=20)  (list(k),k=1,size)
   20       continue
            k = 1
            do while (list(k) .ne. 0)
               ntbnd = ntbnd + 1
               itbnd(1,ntbnd) = list(k)
               itbnd(2,ntbnd) = list(k+1)
               k = k + 2
            end do
         end if
   30    continue
      end do
c
c     keep the main lambda and the sublambdas within the physical
c     lambda range
c
      lambda = min(1.0d0,max(0.0d0,lambda))
      elambda = min(1.0d0,max(0.0d0,elambda))
      plambda = min(1.0d0,max(0.0d0,plambda))
      vlambda = min(1.0d0,max(0.0d0,vlambda))
c
c     set plambda to elambda if no values given
c
      if (.not. setplambda)  plambda = elambda
c
c     turn off hybrid potentials if no sites are mutated
c
      use_mutate = .true.
      if (nmut .eq. 0)  use_mutate = .false.
c
c     set the options for each flavor of the lambda calculation
c
      call mutate_dlmda
      call mutate_ost
      call mutate_abf
      call mutate_meta
      call mutate_window
      call mutate_ti
      call mutate_fep
      call mutate_check
c
c     map the active main lambda and install its dependent parameters
c
      call refreshsublmda
c
c     scale electrostatic parameter values based on lambda
c
      if (.not.use_rel .and.
     &    elambda.ge.0.0d0 .and. elambda.lt.1.0d0) then
         call altelec
      end if
c
c     scale torsional parameter values based on lambda
c
      if (.not.use_rel .and.
     &    tlambda.ge.0.0d0 .and. tlambda.lt.1.0d0) then
         if (ntbnd .ne. 0)  call alttors (ntbnd,itbnd)
      end if
c
c     scale implicit solvation parameter values based on lambda
c
      if (.not.use_rel .and.
     &    elambda.ge.0.0d0 .and. elambda.lt.1.0d0) then
         call altsolv
      end if
c
c     write status of current hybrid potential lambda values
c
      if (use_mutate .and. .not.silent) then
         write (iout,40)
   40    format (/,' Free Energy Perturbation Parameters :')
         write (iout,50)  nmut,vlambda,elambda,plambda,tlambda
   50    format (/,' Number of FEP Hybrid Atoms',9x,i8,
     &           /,' van der Waals Lambda Value',9x,f8.3,
     &           /,' Electrostatics Lambda Value',8x,f8.3,
     &           /,' Polarization Lambda Value',10x,f8.3,
     &           /,' Torsion Angle Lambda Value',9x,f8.3)
c
c     report the mode chosen along each axis of the calculation
c
         if (use_rel) then
            write (iout,60)
   60       format (/,' Free Energy Mode',12x,'Staged Relative')
            do i = 1, 3
               if (use_relstg(i)) then
                  write (iout,70)  leg(i),relstglo(i),relstghi(i)
   70             format (' Staged Leg ',a4,' Range',6x,2f8.3)
               end if
            end do
         else
            write (iout,80)
   80       format (/,' Free Energy Mode',19x,'Absolute')
         end if
         if (use_ost) then
            write (iout,90)
   90       format (' Sampling Mode',27x,'OST')
         else if (use_meta) then
            write (iout,100)
  100       format (' Sampling Mode',18x,'Metadynamics')
         else if (use_ti) then
            write (iout,110)
  110       format (' Sampling Mode',28x,'TI')
         else if (use_abf) then
            write (iout,115)
  115       format (' Sampling Mode',27x,'ABF')
         else if (use_fep) then
            write (iout,117)
  117       format (' Sampling Mode',27x,'FEP')
         else
            write (iout,120)
  120       format (' Sampling Mode',18x,'Fixed Lambda')
         end if
         if (use_epdt) then
            write (iout,150)
  150       format (' Polarization Topology',18x,'Dual')
         else
            write (iout,160)
  160       format (' Polarization Topology',16x,'Single')
         end if
         if (use_mainlmda) then
            write (iout,190)  lambda
  190       format (' Main Lambda Value',18x,f8.3)
         end if
         if (use_rel) then
            write (iout,200)  nmut-nmutb,nmutb
  200       format (/,' Relative Dual Topology Active :',
     &              /,' Number of Ligand1 Atoms',12x,i8,
     &              /,' Number of Ligand2 Atoms',12x,i8)
         end if
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (list)
      deallocate (itbnd)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate_dlmda  --  lambda derivative settings  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate_dlmda" sets the options governing the lambda derivative,
c     which method drives the main lambda, the dual topology end state
c     treatment and the mapping from the main lambda to the individual
c     electrostatic, polarization and van der Waals sublambdas
c
c
      subroutine mutate_dlmda
      use atoms
      use bndstr
      use charge
      use chgpen
      use dipole
      use dlmda
      use expol
      use iounit
      use keys
      use mplpot
      use mpole
      use mutant
      use polar
      use polpot
      use potent
      implicit none
      integer i,j,k
      integer next
      real*8 temp
      real*8 lo,hi
      character*4 legword
      character*8 propword
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     flag for use of lambda derivative
c
      use_abf = .false.
      use_abfdyn = .false.
      use_d2lmda = .false.
      use_dlmda = .false.
      use_elmdamap = .false.
      use_epdt = .false.
      use_fep = .false.
      use_lmdamc = .false.
      use_meta = .false.
      use_metadyn = .false.
      use_ost = .false.
      use_ostdyn = .false.
      use_prst = .false.
      use_plmdamap = .false.
      use_ti = .false.
      use_vlmdamap = .false.
c
c     set defaults describing the flavor of the lambda calculation
c
      lmdasampmode = 'NONE'
c
c     set defaults for dual topology
c
      epdtexp = 1
c
c     set defaults for the staged relative free energy schedule
c
      relstage = 'VDWM'
      do i = 1, 3
         use_relstg(i) = .false.
         relstglo(i) = 0.0d0
         relstghi(i) = 1.0d0
      end do
c
c     set default mapping from main lambda to sublambda
c
      qntelmda1 = 1.0d0
      qntelmda0 = 0.0d0
      qntplmda1 = 1.0d0
      qntplmda0 = 0.0d0
      qntvlmda1 = 1.0d0
      qntvlmda0 = 0.0d0
      elmdamap = 'QNT'
      plmdamap = 'QNT'
      vlmdamap = 'QNT'
      elmdaexp = 1
      plmdaexp = 1
      vlmdaexp = 1
      elmdainvn = 4
      plmdainvn = 4
      vlmdainvn = 4
      elmdainveps = 0.3d0
      plmdainveps = 0.3d0
      vlmdainveps = 0.3d0
      elmdaapmn = 12
      plmdaapmn = 12
      vlmdaapmn = 12
      elmdaapmrho = 4.0d0
      plmdaapmrho = 4.0d0
      vlmdaapmrho = 4.0d0
c
c     a sublambda not driven by the main lambda has an identity chain
c
      deldlmda = 1.0d0
      dpldlmda = 1.0d0
      dvldlmda = 1.0d0
      d2eldlmda2 = 0.0d0
      d2pldlmda2 = 0.0d0
      d2vldlmda2 = 0.0d0
c
c     the lambda derivative is undefined until a gradient is taken
c
      dedl = 0.0d0
      d2edl2 = 0.0d0
c
c     set default lambda sample interval and its phase ratios
c
      lmdastep = 0
      lmdaintv = 100
      lmdaparatio = 0.3d0
      lmdapbratio = 0.3d0
      nlmdasave = 0
c
c     set defaults for the lambda particle propagation
c
      lmdaavg = 0.0d0
      lmdastd = 0.0d0
      dedlavg = 0.0d0
      dedlstd = 0.0d0
      deffdl = 0.0d0
      lmdavtheta = 0.0d0
      lmdamass = 25.0d0
      lmdafric = 0.01d0
      lmdadt = 0.001d0
      lmdathmap = 'TRI'
      lmdathalpha = 0.999999999d0
c
c     set defaults for the Monte Carlo lambda trial moves
c
      lmdamctrial = 10
      lmdamcmove = 0.02d0
      nlmdamctry = 0
      nlmdamcacc = 0
c
c     set default lambda bins and lambda bias free energy
c
      nlmda = 201
      lmdaddgdl = 0.0d0
      lmdadeltag = 0.0d0
c
c     set default convergence gate for accepting interval samples
c
      use_lmdacv = .false.
      lmdacvstd = 50.0d0
      lmdacvrat = 0.5d0
c
c     search keywords for lambda derivative options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:13) .eq. 'LAMBDA-DERIV ') then
            use_dlmda = .true.
         else if (keyword(1:14) .eq. 'LAMBDA-DERIV2 ') then
            use_dlmda = .true.
            use_d2lmda = .true.
         else if (keyword(1:12) .eq. 'LAMBDA-MODE ') then
            call getword (record,lmdasampmode,next)
            call upcase (lmdasampmode)
         else if (keyword(1:16) .eq. 'LAMBDA-INTERVAL ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdaintv
         else if (keyword(1:15) .eq. 'LAMBDA-EQRATIO ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdaparatio,lmdapbratio
         else if (keyword(1:10) .eq. 'LAMBDA-DT ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdadt
         else if (keyword(1:12) .eq. 'LAMBDA-MASS ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdamass
         else if (keyword(1:16) .eq. 'LAMBDA-FRICTION ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdafric
         else if (keyword(1:17) .eq. 'LAMBDA-THETA-MAP ') then
            call getword (record,lmdathmap,next)
            call upcase (lmdathmap)
         else if (keyword(1:19) .eq. 'LAMBDA-THETA-ALPHA ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdathalpha
         else if (keyword(1:17) .eq. 'LAMBDA-PROPAGATE ') then
            call getword (record,propword,next)
            call upcase (propword)
            if (propword .eq. 'MC') then
               use_lmdamc = .true.
            else if (propword .eq. 'LANGEVIN') then
               use_lmdamc = .false.
            else
               write (iout,3)
    3          format (/,' MUTATE_DLMDA  --  LAMBDA-PROPAGATE takes',
     &                    ' either LANGEVIN or MC')
               call fatal
            end if
         else if (keyword(1:16) .eq. 'LAMBDA-MC-TRIAL ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdamctrial
         else if (keyword(1:15) .eq. 'LAMBDA-MC-MOVE ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdamcmove
         else if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  nlmda
         else if (keyword(1:15) .eq. 'LAMBDA-CONVCRI ') then
            use_lmdacv = .true.
            string = record(next:240)
            read (string,*,err=10,end=10)  lmdacvstd,lmdacvrat
         else if (keyword(1:13) .eq. 'POL-DUALTOPO ') then
            use_epdt = .true.
         else if (keyword(1:17) .eq. 'POL-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  epdtexp
         else if (keyword(1:15) .eq. 'ELE-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  qntelmda0, qntelmda1
         else if (keyword(1:15) .eq. 'POL-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  qntplmda0, qntplmda1
         else if (keyword(1:15) .eq. 'VDW-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  qntvlmda0, qntvlmda1
         else if (keyword(1:13) .eq. 'ELE-LMDA-MAP ') then
            use_elmdamap = .true.
            call getword (record,elmdamap,next)
            call upcase (elmdamap)
         else if (keyword(1:13) .eq. 'POL-LMDA-MAP ') then
            use_plmdamap = .true.
            call getword (record,plmdamap,next)
            call upcase (plmdamap)
         else if (keyword(1:13) .eq. 'VDW-LMDA-MAP ') then
            use_vlmdamap = .true.
            call getword (record,vlmdamap,next)
            call upcase (vlmdamap)
         else if (keyword(1:13) .eq. 'ELE-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  elmdaexp
         else if (keyword(1:13) .eq. 'POL-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  plmdaexp
         else if (keyword(1:13) .eq. 'VDW-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  vlmdaexp
         else if (keyword(1:15) .eq. 'ELE-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  elmdainvn
         else if (keyword(1:15) .eq. 'POL-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  plmdainvn
         else if (keyword(1:15) .eq. 'VDW-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  vlmdainvn
         else if (keyword(1:17) .eq. 'ELE-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  elmdainveps
         else if (keyword(1:17) .eq. 'POL-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  plmdainveps
         else if (keyword(1:17) .eq. 'VDW-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  vlmdainveps
         else if (keyword(1:10) .eq. 'REL-STAGE ') then
            use_rel = .true.
            use_mainlmda = .true.
            call getword (record,legword,next)
            call upcase (legword)
            k = 0
            if (legword .eq. 'LIG2')  k = 1
            if (legword .eq. 'VDWM')  k = 2
            if (legword .eq. 'LIG1')  k = 3
            if (k .eq. 0) then
               write (iout,2)
    2          format (/,' MUTATE_DLMDA  --  REL-STAGE requires the',
     &                    ' leg to be named; use LIG2 to discharge',
     &                    ' ligand 2, VDWM to morph van der Waals, or',
     &                    ' LIG1 to charge ligand 1')
               call fatal
            end if
c
c     a leg spans the whole main lambda unless it names both bounds
c
            use_relstg(k) = .true.
            relstglo(k) = 0.0d0
            relstghi(k) = 1.0d0
            string = record(next:240)
            read (string,*,err=5,end=5)  lo,hi
            relstglo(k) = lo
            relstghi(k) = hi
            goto 10
    5       continue
            read (string,*,err=10,end=10)  lo
            write (iout,7)
    7       format (/,' MUTATE_DLMDA  --  REL-STAGE takes both bounds',
     &                 ' of the main lambda range of a leg, or neither')
            call fatal
         end if
   10    continue
      end do
c
c     keep an odd lambda grid and define the lambda bin width
c
      if (nlmda .lt. 3)  nlmda = 3
      if (mod(nlmda,2) .eq. 0)  nlmda = nlmda + 1
      wlmda = 1.0d0 / dble(nlmda-1)
      wlmda2 = 0.5d0 * wlmda
c
c     keep the convergence tolerances nonnegative
c
      lmdacvstd = abs(lmdacvstd)
      lmdacvrat = abs(lmdacvrat)
c
c     split the sample interval into its propagation, equilibration
c     and averaging phases
c
      if (lmdaintv .lt. 1)  lmdaintv = 1
      if (lmdaparatio .lt. 0.0d0)  lmdaparatio = 0.0d0
      if (lmdapbratio .lt. 0.0d0)  lmdapbratio = 0.0d0
      call setlmdaphase
c
c     keep the Monte Carlo lambda move within a single reflection
c
      if (lmdamctrial .lt. 1)  lmdamctrial = 1
      lmdamcmove = min(1.0d0,abs(lmdamcmove))
c
c     validate the theta map of the lambda particle
c
      if (lmdathmap.ne.'SIN' .and. lmdathmap.ne.'TRI') then
         lmdathmap = 'TRI'
      end if
      if (lmdathalpha.le.0.0d0 .or. lmdathalpha.ge.1.0d0) then
         lmdathalpha = 0.999999999d0
      end if
c
c     start the lambda particle from the current main lambda
c
      call lmdathetainv (lambda,lmdatheta)
c
c     the lambda sampling mode sets its method flag, and every
c     sampling method moves a main lambda by its derivative
c
      if (lmdasampmode .eq. 'OST')  use_ost = .true.
      if (lmdasampmode .eq. 'META')  use_meta = .true.
      if (lmdasampmode .eq. 'TI')  use_ti = .true.
      if (lmdasampmode .eq. 'ABF')  use_abf = .true.
      if (lmdasampmode .eq. 'FEP')  use_fep = .true.
      if (use_ost .or. use_meta .or. use_ti .or. use_abf) then
         use_dlmda = .true.
         use_mainlmda = .true.
      end if
c
c     free energy perturbation moves a main lambda from window to
c     window, but makes no use of the lambda derivative
c
      if (use_fep)  use_mainlmda = .true.
c
c     only orthogonal space tempering consumes the second, force and
c     virial lambda derivatives, which LAMBDA-DERIV2 also asks for
c
      if (use_ost)  use_d2lmda = .true.
c
c     a main lambda drives sublambdas that name a map
c
      if (use_mainlmda .and. use_rel) then
c
c     the staged relative schedule maps every sublambda itself
c
         use_elmdamap = .true.
         use_plmdamap = .true.
         use_vlmdamap = .true.
      else if (use_mainlmda) then
         if (.not. (use_elmdamap .or. use_plmdamap .or.
     &              use_vlmdamap)) then
            if (use_dlmda) then
               write (iout,20)
   20          format (/,' MUTATE_DLMDA  --  A Lambda Derivative',
     &                    ' requires an explicit map for each driven',
     &                    ' sublambda; add the ELE-LMDA-MAP,',
     &                    ' POL-LMDA-MAP or VDW-LMDA-MAP keywords')
               call fatal
            end if
            if (use_fep) then
               write (iout,25)
   25          format (/,' MUTATE_DLMDA  --  Free Energy Perturbation',
     &                    ' requires an explicit map for each driven',
     &                    ' sublambda; add the ELE-LMDA-MAP,',
     &                    ' POL-LMDA-MAP or VDW-LMDA-MAP keywords')
               call fatal
            end if
            use_elmdamap = .true.
            use_plmdamap = .true.
            use_vlmdamap = .true.
            elmdamap = 'EXP'
            plmdamap = 'EXP'
            vlmdamap = 'EXP'
         end if
c
c     a sublambda the main lambda does not drive is held at its value
c     and leaves the chain rule
c
         if (.not. use_elmdamap)  deldlmda = 0.0d0
         if (.not. use_plmdamap)  dpldlmda = 0.0d0
         if (.not. use_vlmdamap)  dvldlmda = 0.0d0
      end if
c
c     set the terms that carry a lambda derivative
c
      call setdlmdaterms
c
c     second, force and virial lambda derivatives of polarization
c     need the dual topology path; the first alone uses one state
c
      if (use_d2lmda .and. use_pdlmda)  use_epdt = .true.
c
c     validate mapping schemes from main lambda to sublambdas
c
      if (elmdamap.ne.'QNT' .and. elmdamap.ne.'EXP'
     &       .and. elmdamap.ne.'INV'
     &       .and. elmdamap.ne.'APM') then
         elmdamap = 'QNT'
      end if
      if (plmdamap.ne.'QNT' .and. plmdamap.ne.'EXP'
     &       .and. plmdamap.ne.'INV'
     &       .and. plmdamap.ne.'APM') then
         plmdamap = 'QNT'
      end if
      if (vlmdamap.ne.'QNT' .and. vlmdamap.ne.'EXP'
     &       .and. vlmdamap.ne.'INV'
     &       .and. vlmdamap.ne.'APM') then
         vlmdamap = 'QNT'
      end if
      if (epdtexp .lt. 1)  epdtexp = 1
      if (elmdaexp .lt. 1)  elmdaexp = 1
      if (plmdaexp .lt. 1)  plmdaexp = 1
      if (vlmdaexp .lt. 1)  vlmdaexp = 1
      if (elmdainvn .lt. 1)  elmdainvn = 1
      if (plmdainvn .lt. 1)  plmdainvn = 1
      if (vlmdainvn .lt. 1)  vlmdainvn = 1
      if (elmdainveps .lt. 0.0d0)  elmdainveps = -elmdainveps
      if (plmdainveps .lt. 0.0d0)  plmdainveps = -plmdainveps
      if (vlmdainveps .lt. 0.0d0)  vlmdainveps = -vlmdainveps
      if (elmdaapmn .lt. 2)  elmdaapmn = 2
      if (plmdaapmn .lt. 2)  plmdaapmn = 2
      if (vlmdaapmn .lt. 2)  vlmdaapmn = 2
      if (elmdaapmrho .lt. 1.0d0)  elmdaapmrho = 1.0d0
      if (plmdaapmrho .lt. 1.0d0)  plmdaapmrho = 1.0d0
      if (vlmdaapmrho .lt. 1.0d0)  vlmdaapmrho = 1.0d0
      elmdaapmrho = min(elmdaapmrho,dble(elmdaapmn)-0.001d0)
      plmdaapmrho = min(plmdaapmrho,dble(plmdaapmn)-0.001d0)
      vlmdaapmrho = min(vlmdaapmrho,dble(vlmdaapmn)-0.001d0)
c
c     check sublambda intervals are in [0,1] and ordered
c
      if (qntelmda0 .lt. 0.0d0)  qntelmda0 = 0.0d0
      if (qntelmda1 .lt. 0.0d0)  qntelmda1 = 0.0d0
      if (qntplmda0 .lt. 0.0d0)  qntplmda0 = 0.0d0
      if (qntplmda1 .lt. 0.0d0)  qntplmda1 = 0.0d0
      if (qntvlmda0 .lt. 0.0d0)  qntvlmda0 = 0.0d0
      if (qntvlmda1 .lt. 0.0d0)  qntvlmda1 = 0.0d0
      if (qntelmda0 .gt. 1.0d0)  qntelmda0 = 1.0d0
      if (qntelmda1 .gt. 1.0d0)  qntelmda1 = 1.0d0
      if (qntplmda0 .gt. 1.0d0)  qntplmda0 = 1.0d0
      if (qntplmda1 .gt. 1.0d0)  qntplmda1 = 1.0d0
      if (qntvlmda0 .gt. 1.0d0)  qntvlmda0 = 1.0d0
      if (qntvlmda1 .gt. 1.0d0)  qntvlmda1 = 1.0d0
      if (qntelmda1 .lt. qntelmda0) then
         temp = qntelmda0
         qntelmda0 = qntelmda1
         qntelmda1 = temp
      end if
      if (qntplmda1 .lt. qntplmda0) then
         temp = qntplmda0
         qntplmda0 = qntplmda1
         qntplmda1 = temp
      end if
      if (qntvlmda1 .lt. qntvlmda0) then
         temp = qntvlmda0
         qntvlmda0 = qntvlmda1
         qntvlmda1 = temp
      end if
c
c     single topology evaluates polarization from one parameter state
c     at plambda, independently of the electrostatic lambda state
c
      use_prst = (use_mutate .and. .not.use_epdt .and. use_polar)
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(pchgorig))  deallocate (pchgorig)
      if (allocated(bdplorig))  deallocate (bdplorig)
      if (allocated(poleorig))  deallocate (poleorig)
      if (allocated(pcoreorig))  deallocate (pcoreorig)
      if (allocated(pvalorig))  deallocate (pvalorig)
      if (allocated(polarityorig))  deallocate (polarityorig)
      if (allocated(douindorig))  deallocate (douindorig)
      if (allocated(prepeporig))  deallocate (prepeporig)
      if (allocated(dpscdl))  deallocate (dpscdl)
      allocate (pchgorig(n))
      allocate (bdplorig(nbond))
      allocate (poleorig(maxpole,n))
      allocate (pcoreorig(n))
      allocate (pvalorig(n))
      allocate (polarityorig(n))
      allocate (douindorig(n))
      allocate (prepeporig(n))
      allocate (dpscdl(3,3,n))
c
c     copy original parameters for lambda derivative calculations
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            pchgorig(k) = pchg(k)
         end do
      end if
      if (use_dipole) then
         do i = 1, ndipole
            bdplorig(i) = bdpl(i)
         end do
      end if
      if (use_mpole .or. use_polar) then
         do i = 1, npole
            k = ipole(i)
            do j = 1, 13
               poleorig(j,k) = pole(j,k)
            end do
            if (use_chgpen) then
               pcoreorig(k) = pcore(k)
               pvalorig(k) = pval(k)
            end if
         end do
      end if
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            polarityorig(k) = polarity(k)
            douindorig(k) = douind(k)
            if (use_expol)  prepeporig(k) = prepep(k)
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setdlmdaterms  --  terms with a lambda deriv  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setdlmdaterms" decides which of the three terms carries a lambda
c     derivative, and so which of them has to be routed to its "empole4"
c     flavored energy routine rather than the plain gradient one
c
c     a term qualifies when the main lambda drives its sublambda through
c     a map, except on a staged relative schedule, where a term that no
c     declared leg drives has its sublambda pinned to a constant; a
c     pinned sublambda has a flat chain rule, so its term has nothing
c     for the lambda derivative to sample and the plain routine gives
c     the same answer for less work
c
c     "gradient" zeroes every lambda derivative accumulator before it
c     dispatches, so a term left out here reports exact zeros
c
c
      subroutine setdlmdaterms
      use dlmda
      use mutant
      use potent
      implicit none
c
c
      use_edlmda = use_dlmda .and. use_elmdamap .and. use_mutate
      use_pdlmda = use_dlmda .and. use_plmdamap .and. use_mutate
      use_vdlmda = use_dlmda .and. use_vlmdamap .and. use_mutate
      if (use_rel) then
         if (.not. (use_relstg(1) .or. use_relstg(3))) then
            use_edlmda = .false.
            use_pdlmda = .false.
         end if
         if (.not. use_relstg(2))  use_vdlmda = .false.
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mutate_ost  --  orthogonal space tempering set  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mutate_ost" sets the flambda grid and the gaussian widths,
c     heights and tempering used by orthogonal space tempering, then
c     allocates the histogram and the bias kernels when the method is
c     active; the lambda grid, sample interval, convergence gate and
c     lambda particle parameters are set by "mutate_dlmda"
c
c
      subroutine mutate_ost
      use dlmda
      use keys
      use ost
      implicit none
      integer i,k
      integer next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default ost bias derivative
c
      ostdgdl = 0.0d0
c
c     set defaults for tempering of the deposited gaussian heights
c
      use_ostgtemp = .false.
      use_ostltemp = .false.
      ostgthresh = 1.0d0
      ostgtempgamma = 1.0d0
      ostlthresh = 1.0d0
      ostltempgamma = 1.0d0
c
c     set default ost flambda bin values
c
      nflmda = 1001
      wflmda = 1.0d0
      wlhist = 0.005d0
      wfhist = 1.0d0
      fli0 = (nflmda + 1) / 2 + (nflmda - 1) / 4
      hbias = 0.00001d0
      oststdev = 4.0d0
      fastkernel = .true.
      ostinterpol = .false.
c
c     search keywords for orthogonal space tempering options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:14) .eq. 'FLAMBDA-WIDTH ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  wflmda
         else if (keyword(1:7) .eq. 'WLHIST ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  wlhist
         else if (keyword(1:7) .eq. 'WFHIST ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  wfhist
         else if (keyword(1:11) .eq. 'OST-STDDEV ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  oststdev
         else if (keyword(1:16) .eq. 'OST-INTERPOLATE ') then
            ostinterpol = .true.
            fastkernel = .true.
         else if (keyword(1:6) .eq. 'HBIAS ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  hbias
         else if (keyword(1:18) .eq. 'OST-TEMPER-GLOBAL ') then
            use_ostgtemp = .true.
            string = record(next:240)
            read (string,*,err=10,end=10)  ostgthresh,ostgtempgamma
         else if (keyword(1:17) .eq. 'OST-TEMPER-LOCAL ') then
            use_ostltemp = .true.
            string = record(next:240)
            read (string,*,err=10,end=10)  ostlthresh,ostltempgamma
         end if
   10    continue
      end do
c
c     define flambda width and range
c
      wflmda2 = 0.5d0 * wflmda
      fli0 = (nflmda + 1) / 2 + (nflmda - 1) / 4
      if (wlhist .lt. 0.0d0) then
         wlhist = -wlhist
      else if (wlhist .eq. 0.0d0) then
         wlhist = 0.005d0
      end if
      if (wfhist .lt. 0.0d0) then
         wfhist = -wfhist
      else if (wfhist .eq. 0.0d0) then
         wfhist = 1.0d0
      end if
      if (hbias .lt. 0.0d0)  hbias = -hbias
      if (ostgthresh .lt. 0.0d0) then
         ostgthresh = -ostgthresh
      else if (ostgthresh .eq. 0.0d0) then
         ostgthresh = 1.0d0
      end if
      if (ostgtempgamma .lt. 0.0d0) then
         ostgtempgamma = -ostgtempgamma
      else if (ostgtempgamma .eq. 0.0d0) then
         ostgtempgamma = 1.0d0
      end if
      if (ostlthresh .lt. 0.0d0) then
         ostlthresh = -ostlthresh
      else if (ostlthresh .eq. 0.0d0) then
         ostlthresh = 1.0d0
      end if
      if (ostltempgamma .lt. 0.0d0) then
         ostltempgamma = -ostltempgamma
      else if (ostltempgamma .eq. 0.0d0) then
         ostltempgamma = 1.0d0
      end if
      maxwlhist = wlhist
      maxwfhist = wfhist
c
c     free the histogram and kernel arrays of any earlier setup,
c     so that only the arrays of the chosen method remain
c
      if (allocated(osthhist))  deallocate (osthhist)
      if (allocated(osthist))  deallocate (osthist)
      if (allocated(lmdaihist))  deallocate (lmdaihist)
      if (allocated(osthead))  deallocate (osthead)
      if (allocated(ostnext))  deallocate (ostnext)
      if (allocated(lmdallist))  deallocate (lmdallist)
      if (allocated(lmdaflist))  deallocate (lmdaflist)
      if (allocated(lmdalhist))  deallocate (lmdalhist)
      if (allocated(lmdafhist))  deallocate (lmdafhist)
      if (allocated(ostwlhist))  deallocate (ostwlhist)
      if (allocated(ostwfhist))  deallocate (ostwfhist)
      if (allocated(lmdafmean))  deallocate (lmdafmean)
      if (allocated(lmdafsum))  deallocate (lmdafsum)
      if (allocated(gfkernel))  deallocate (gfkernel)
      if (allocated(gkernel))  deallocate (gkernel)
      if (allocated(glfkernel))  deallocate (glfkernel)
      if (allocated(glkernel))  deallocate (glkernel)
      if (allocated(lmdafwt))  deallocate (lmdafwt)
      if (allocated(vkernelmax))  deallocate (vkernelmax)
c
c     allocate ost histogram and kernels
c
      if (use_ost) then
         sizelmdahist = 10000
         nlmdahist = 0
         allocate (osthhist(sizelmdahist))
         allocate (osthist(sizelmdahist))
         allocate (lmdaihist(sizelmdahist))
         allocate (osthead(nlmda,nflmda))
         allocate (ostnext(sizelmdahist))
         allocate (lmdallist(lmdaintv))
         allocate (lmdaflist(lmdaintv))
         allocate (lmdalhist(sizelmdahist))
         allocate (lmdafhist(sizelmdahist))
         allocate (ostwlhist(sizelmdahist))
         allocate (ostwfhist(sizelmdahist))
         allocate (lmdafmean(nlmda))
         allocate (lmdafsum(nlmda))
         allocate (gfkernel(nlmda,nflmda))
         allocate (gkernel(nlmda,nflmda))
         allocate (glfkernel(nlmda,nflmda))
         allocate (glkernel(nlmda,nflmda))
         allocate (lmdafwt(nlmda))
         allocate (vkernelmax(nlmda))
c
c     initialize ost histogram and kernels
c
         do i = 1, nlmda
            lmdafmean(i) = 0.0d0
            lmdafsum(i) = 0.0d0
            lmdafwt(i) = 0.0d0
            vkernelmax(i) = 0.0d0
            do k = 1, nflmda
               gfkernel(i,k) = 0.0d0
               gkernel(i,k) = 0.0d0
               glfkernel(i,k) = 0.0d0
               glkernel(i,k) = 0.0d0
               osthead(i,k) = 0
            end do
         end do
         do i = 1, lmdaintv
            lmdallist(i) = 0.0d0
            lmdaflist(i) = 0.0d0
         end do
         do i = 1, sizelmdahist
            osthist(i) = 0
            lmdaihist(i) = 0
            ostnext(i) = 0
            lmdalhist(i) = 0.0d0
            lmdafhist(i) = 0.0d0
            osthhist(i) = 0.0d0
            ostwlhist(i) = 0.0d0
            ostwfhist(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine mutate_abf  --  adaptive biasing force arrays  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "mutate_abf" allocates the interval samples, the interval lists
c     and the lambda bins, the only histogram arrays that the adaptive
c     biasing force method uses
c
c
      subroutine mutate_abf
      use dlmda
      implicit none
      integer i
c
c
c     allocate and zero the abf samples, interval lists and bins
c
      if (use_abf) then
         if (allocated(lmdaihist))  deallocate (lmdaihist)
         if (allocated(lmdalhist))  deallocate (lmdalhist)
         if (allocated(lmdafhist))  deallocate (lmdafhist)
         if (allocated(lmdallist))  deallocate (lmdallist)
         if (allocated(lmdaflist))  deallocate (lmdaflist)
         if (allocated(lmdafmean))  deallocate (lmdafmean)
         if (allocated(lmdafsum))  deallocate (lmdafsum)
         if (allocated(lmdafwt))  deallocate (lmdafwt)
         sizelmdahist = 10000
         nlmdahist = 0
         allocate (lmdaihist(sizelmdahist))
         allocate (lmdalhist(sizelmdahist))
         allocate (lmdafhist(sizelmdahist))
         allocate (lmdallist(lmdaintv))
         allocate (lmdaflist(lmdaintv))
         allocate (lmdafmean(nlmda))
         allocate (lmdafsum(nlmda))
         allocate (lmdafwt(nlmda))
         do i = 1, sizelmdahist
            lmdaihist(i) = 0
            lmdalhist(i) = 0.0d0
            lmdafhist(i) = 0.0d0
         end do
         do i = 1, lmdaintv
            lmdallist(i) = 0.0d0
            lmdaflist(i) = 0.0d0
         end do
         do i = 1, nlmda
            lmdafmean(i) = 0.0d0
            lmdafsum(i) = 0.0d0
            lmdafwt(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine mutate_meta  --  metadynamics parameters  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "mutate_meta" allocates the history of gaussians deposited along
c     the main lambda coordinate and the grid holding the accumulated
c     metadynamics bias
c
c
      subroutine mutate_meta
      use dlmda
      use ost
      implicit none
      integer i
c
c
c     zero out the count of deposited metadynamics gaussians
c
      nmetahist = 0
      nmethistsave = 0
c
c     allocate metadynamics gaussian history
c
      if (use_meta) then
         if (allocated(metalhist))  deallocate (metalhist)
         if (allocated(metahhist))  deallocate (metahhist)
         if (allocated(metawhist))  deallocate (metawhist)
         if (allocated(metaihist))  deallocate (metaihist)
         if (allocated(lmdallist))  deallocate (lmdallist)
         if (allocated(vmetagrid))  deallocate (vmetagrid)
         if (allocated(dvmetagrid))  deallocate (dvmetagrid)
         sizemetahist = 10000
         allocate (metalhist(sizemetahist))
         allocate (metahhist(sizemetahist))
         allocate (metawhist(sizemetahist))
         allocate (metaihist(sizemetahist))
         allocate (lmdallist(lmdaintv))
         allocate (vmetagrid(nlmda))
         allocate (dvmetagrid(nlmda))
         do i = 1, sizemetahist
            metalhist(i) = 0.0d0
            metahhist(i) = 0.0d0
            metawhist(i) = 0.0d0
            metaihist(i) = 0
         end do
         do i = 1, lmdaintv
            lmdallist(i) = 0.0d0
         end do
         do i = 1, nlmda
            vmetagrid(i) = 0.0d0
            dvmetagrid(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine mutate_window  --  lambda window schedule  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "mutate_window" sets the schedule of lambda windows visited in
c     turn by a sequential lambda method, either from explicit
c     "LAMBDA-WINDOW" values or from the "LAMBDA-NWINDOW" count of
c     evenly spaced windows; the first "LAMBDA-EQRATIO" value is the
c     fraction of each window discarded as equilibration
c
c
      subroutine mutate_window
      use dlmda
      use iounit
      use keys
      implicit none
      integer i
      integer next
      integer nwin
      real*8 frac
      logical nwinset
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for the sequential lambda windows
c
      lmdawin = 0
      lmdawineq = 0
      lmdawinlen = 0
      nlmdawin = 21
      lmdawinratio = 0.5d0
c
c     only a sequential lambda method visits the lambda windows
c
      if (allocated(lmdawinlist))  deallocate (lmdawinlist)
      if (allocated(lmdawinfrac))  deallocate (lmdawinfrac)
      if (.not.use_ti .and. .not.use_fep)  return
c
c     size the lambda window schedule to the worst case
c
      allocate (lmdawinlist(max(1,nkey)))
      allocate (lmdawinfrac(max(1,nkey)))
      do i = 1, max(1,nkey)
         lmdawinlist(i) = 0.0d0
         lmdawinfrac(i) = -1.0d0
      end do
      nwin = 0
      nwinset = .false.
c
c     search keywords for lambda window options; a window value
c     that cannot be read is left invalid and is rejected below, while
c     an equilibration ratio that cannot be read keeps its default
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:15) .eq. 'LAMBDA-NWINDOW ') then
            nlmdawin = 0
            nwinset = .true.
            read (string,*,err=10,end=10)  nlmdawin
         else if (keyword(1:14) .eq. 'LAMBDA-WINDOW ') then
            nwin = nwin + 1
            lmdawinlist(nwin) = -1.0d0
            read (string,*,err=10,end=10)  lmdawinlist(nwin)
            read (string,*,err=10,end=10)  lmdawinlist(nwin),frac
            if (.not. (frac .gt. 0.0d0))  frac = 0.0d0
            lmdawinfrac(nwin) = frac
         else if (keyword(1:15) .eq. 'LAMBDA-EQRATIO ') then
            read (string,*,err=10,end=10)  lmdawinratio
         end if
   10    continue
      end do
c
c     some of each window must be left after the equilibration, so
c     an equilibration ratio outside [0,1) is reset to its default
c
      if (.not. (lmdawinratio.ge.0.0d0 .and.
     &           lmdawinratio.lt.1.0d0))  lmdawinratio = 0.5d0
c
c     build the schedule of lambda windows from the values found
c
      call setlmdasched (nwin,nwinset)
c
c     free energy perturbation needs a neighbor for each window
c
      if (use_fep .and. nlmdawin.lt.2) then
         write (iout,30)
   30    format (/,' MUTATE_WINDOW  --  Free Energy Perturbation',
     &              ' needs at least two lambda windows')
         call fatal
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine mutate_ti  --  dU/dlambda block averaging  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "mutate_ti" sets the block size used to average dU/dlambda
c     over the lambda windows of a thermodynamic integration
c
c
      subroutine mutate_ti
      use dlmda
      use keys
      use thrmint
      implicit none
      integer i
      integer next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for thermodynamic integration averaging
c
      tinbcount = 0
      tinbsave = 0
      tinbtot = 0
      tinstepavg = 100
c
c     search keywords for thermodynamic integration options; a
c     value that cannot be read keeps its default
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:12) .eq. 'TI-NSTEPAVG ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  tinstepavg
         end if
   10    continue
      end do
c
c     each block average must hold at least one dynamics step, so
c     a smaller block is reset to its default
c
      if (use_ti) then
         if (tinstepavg .lt. 1)  tinstepavg = 100
c
c     allocate the dU/dlambda block buffer
c
         if (allocated(tidedllist))  deallocate (tidedllist)
         allocate (tidedllist(tinstepavg))
         do i = 1, tinstepavg
            tidedllist(i) = 0.0d0
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine mutate_fep  --  trial energy sample interval  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "mutate_fep" sets the number of dynamics steps between the
c     samples of the trial lambda energies found during free energy
c     perturbation over lambda windows, which is taken from the first
c     "LAMBDA-INTERVAL" value
c
c
      subroutine mutate_fep
      use dlmda
      use fep
      use keys
      implicit none
      integer i
      integer next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for free energy perturbation sampling
c
      nfep = 0
      nfepsave = 0
      nfeptot = 0
      fepintv = 1000
c
c     search keywords for free energy perturbation options; a
c     value that cannot be read keeps its default
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:16) .eq. 'LAMBDA-INTERVAL ') then
            string = record(next:240)
            read (string,*,err=10,end=10)  fepintv
         end if
   10    continue
      end do
c
c     the samples must be at least one dynamics step apart, so a
c     smaller interval is reset to its default
c
      if (fepintv .lt. 1)  fepintv = 1000
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine mutate_check  --  hybrid system consistency  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "mutate_check" tests the free energy options that span more than
c     one of the setup routines, with each routine having already
c     validated the keywords that it owns
c
c
      subroutine mutate_check
      use dlmda
      use dsppot
      use iounit
      use keys
      use limits
      use mplpot
      use mutant
      use polpot
      use potent
      use vdwpot
      implicit none
      integer i,k,next
      integer nstg
      real*8 eps,temp
      real*8 dspcut
      logical polsave,qntstep
      character*20 keyword
      character*240 record
c
c
c     reset any soft core buffered 14-7 parameter outside its allowed
c     range to the default value, where the second lambda derivative
c     diverges at an endpoint unless the lambda exponent is at least
c     two and the one minus lambda shift exponent is one or at least
c     two; the buffering constants keep the pair energy bounded as the
c     shift vanishes, so the shift exponents need no further limit
c
      if (scalphav .lt. 0.0d0)  scalphav = 0.7d0
      if (scexp .lt. 1.0d0)  scexp = 5.0d0
      if (scexpm .lt. 0.0d0)  scexpm = 0.0d0
      if (scexpp .lt. 1.0d0)  scexpp = 2.0d0
      if (use_d2lmda .and. scexp.lt.2.0d0)  scexp = 5.0d0
      if (use_d2lmda .and. scexpp.ne.1.0d0 .and. scexpp.lt.2.0d0)
     &   scexpp = 2.0d0
c
c     reset any soft core Pauli repulsion parameter outside its allowed
c     range to the default value, where the lambda exponent must outrun
c     the vanishing shift so the pair energy stays bounded as lambda
c     goes to zero, and the second lambda derivative diverges at an
c     endpoint unless the exponent is at least two and the one minus
c     lambda shift exponent is one or at least two
c
      if (scalphar .lt. 0.0d0)  scalphar = 1.0d0
      if (scexpr .lt. 1.0d0)  scexpr = 5.0d0
      if (scexprm .lt. 0.0d0)  scexprm = 3.0d0
      if (scexprp .lt. 1.0d0)  scexprp = 1.0d0
      if (use_d2lmda .and. scexpr.lt.2.0d0)  scexpr = 5.0d0
      if (use_d2lmda .and. scexprp.ne.1.0d0 .and. scexprp.lt.2.0d0)
     &   scexprp = 1.0d0
      if (scalphar.gt.0.0d0 .and. scexpr.le.0.5d0*scexprm) then
         scexpr = 5.0d0
         scexprm = 3.0d0
      end if
c
c     reset any soft core damped dispersion parameter outside its
c     allowed range to the default value, where the switch radius must
c     lie inside the dispersion cutoff, or inside the start of its
c     taper when the long range correction fills the taper window, so
c     that reciprocal space and the long range correction see only the
c     plain lambda power, and the second lambda derivative diverges at
c     zero lambda unless the exponent is one or at least two and the
c     sum of the exponents is at least two
c
      dspcut = dispcut
      if (use_dcorr)  dspcut = disptaper
      if (use_dewald)  dspcut = dewaldcut
      if (scrdsp.lt.0.0d0 .or. scrdsp.ge.dspcut)  scrdsp = 2.0d0
      if (scexpd .lt. 1.0d0)  scexpd = 4.0d0
      if (scexpdq .le. 0.0d0)  scexpdq = 1.0d0
      if (use_d2lmda .and. scexpd.ne.1.0d0 .and. scexpd.lt.2.0d0)
     &   scexpd = 4.0d0
      if (use_d2lmda .and. scexpd+scexpdq.lt.2.0d0)  scexpdq = 1.0d0
c
c     the implicit solvation parameters are scaled once at setup and
c     carry no lambda derivative, so they cannot follow a lambda that
c     moves by its derivative or by Monte Carlo
c
      if (use_solv .and. (use_dlmda.or.use_lmdamc)) then
         write (iout,20)
   20    format (/,' MUTATE_CHECK  --  Lambda Dynamics, Lambda',
     &              ' Derivatives and Lambda Monte Carlo are not',
     &              ' Available with Implicit Solvation')
         call fatal
      end if
c
c     single topology currently supports the scalar derivative for
c     mutual polarization with Thole damping, as in AMOEBA and AMOEBA+,
c     or charge penetration damping, as in HIPPO, with or without
c     exchange polarization, in direct and Ewald modes
c
      if (use_pdlmda .and. use_prst) then
         if (poltyp.ne.'MUTUAL' .or.
     &       .not.(use_thole.or.use_chgpen)) then
            write (iout,30)
   30       format (/,' MUTATE_CHECK  --  Single Topology',
     &                 ' Polarization dU/dLambda supports MUTUAL',
     &                 ' polarization with Thole or charge penetration',
     &                 ' damping only, in either direct or Ewald',
     &                 ' mode; use POLARIZATION MUTUAL, or add',
     &                 ' LAMBDA-DERIV2 or POL-DUALTOPO')
            call fatal
         end if
      end if
c
c     every sublambda is mapped from the main lambda, so a lambda
c     derivative has nothing to differentiate without one
c
      if (use_dlmda .and. .not.use_mainlmda) then
         write (iout,50)
   50    format (/,' MUTATE_CHECK  --  A Lambda Derivative requires',
     &              ' a main lambda; add the LAMBDA keyword and a map',
     &              ' for each driven sublambda')
         call fatal
      end if
c
c     a quintic taper over an empty window is a step in the sublambda,
c     whose lambda derivative is undefined at the step
c
      if (use_dlmda .and. use_mainlmda .and. .not.use_rel) then
         qntstep = .false.
         if (use_elmdamap .and. qntelmda1-qntelmda0.lt.1.0d-8 .and.
     &       elmdamap.ne.'EXP' .and. elmdamap.ne.'INV' .and.
     &       elmdamap.ne.'APM')  qntstep = .true.
         if (use_plmdamap .and. qntplmda1-qntplmda0.lt.1.0d-8 .and.
     &       plmdamap.ne.'EXP' .and. plmdamap.ne.'INV' .and.
     &       plmdamap.ne.'APM')  qntstep = .true.
         if (use_vlmdamap .and. qntvlmda1-qntvlmda0.lt.1.0d-8 .and.
     &       vlmdamap.ne.'EXP' .and. vlmdamap.ne.'INV' .and.
     &       vlmdamap.ne.'APM')  qntstep = .true.
         if (qntstep) then
            write (iout,55)
   55       format (/,' MUTATE_CHECK  --  A Lambda Derivative requires',
     &                 ' a sublambda range of nonzero width for each',
     &                 ' QNT map')
            call fatal
         end if
      end if
c
c     the staged relative schedule maps every sublambda on its own,
c     so a sublambda set by its own keyword would be overwritten
c
      if (use_rel .and.
     &    (setelambda .or. setplambda .or. setvlambda)) then
         write (iout,60)
   60    format (/,' MUTATE_CHECK  --  REL-STAGE sets each sublambda',
     &              ' from its own schedule; remove the ELE-LAMBDA,',
     &              ' POL-LAMBDA and VDW-LAMBDA keywords')
         call fatal
      end if
c
c     a main lambda maps each sublambda that names a map, so one set by
c     its own keyword would be overwritten; a sublambda the main lambda
c     does not drive is free to be pinned by its own keyword, which is
c     how a leg holds a term at a fixed coupling state while the main
c     lambda morphs another term
c
      if (use_mainlmda .and. use_elmdamap .and. setelambda) then
         write (iout,70)
   70    format (/,' MUTATE_CHECK  --  LAMBDA drives elambda through',
     &              ' a lambda map; remove the ELE-LAMBDA keyword, or',
     &              ' name a map only for the sublambdas LAMBDA drives')
         call fatal
      end if
      if (use_mainlmda .and. use_plmdamap .and. setplambda) then
         write (iout,80)
   80    format (/,' MUTATE_CHECK  --  LAMBDA drives plambda through',
     &              ' a lambda map; remove the POL-LAMBDA keyword, or',
     &              ' name a map only for the sublambdas LAMBDA drives')
         call fatal
      end if
      if (use_mainlmda .and. use_vlmdamap .and. setvlambda) then
         write (iout,90)
   90    format (/,' MUTATE_CHECK  --  LAMBDA drives vlambda through',
     &              ' a lambda map; remove the VDW-LAMBDA keyword, or',
     &              ' name a map only for the sublambdas LAMBDA drives')
         call fatal
      end if
c
c     the staged relative schedule morphs one ligand into another
c
      if (use_rel .and. nmutb.eq.0) then
         write (iout,100)
  100    format (/,' MUTATE_CHECK  --  REL-STAGE requires a second',
     &              ' ligand group; add the LIGAND2 keyword')
         call fatal
      end if
c
c     the staged legs must follow one another without overlap in the
c     order LIG2, VDWM and LIG1 of rising main lambda; bounds given in
c     reverse are swapped, and two legs that meet share one bound
c
      eps = 1.0d-8
      nstg = 0
      k = 0
      do i = 1, 3
         if (use_relstg(i)) then
            nstg = nstg + 1
            relstglo(i) = min(1.0d0,max(0.0d0,relstglo(i)))
            relstghi(i) = min(1.0d0,max(0.0d0,relstghi(i)))
            if (relstghi(i) .lt. relstglo(i)) then
               temp = relstglo(i)
               relstglo(i) = relstghi(i)
               relstghi(i) = temp
            end if
            if (k .ne. 0) then
               if (relstglo(i) .lt. relstghi(k)-eps) then
                  write (iout,101)
  101             format (/,' MUTATE_CHECK  --  REL-STAGE ranges must',
     &                       ' not overlap and must follow the order',
     &                       ' LIG2, VDWM, LIG1 as the main lambda',
     &                       ' rises')
                  call fatal
               end if
               if (relstglo(i) .lt. relstghi(k)+eps) then
                  relstglo(i) = relstghi(k)
               end if
            end if
            if (relstghi(i) .le. relstglo(i)) then
               write (iout,102)
  102          format (/,' MUTATE_CHECK  --  REL-STAGE range must span',
     &                    ' a nonzero width of the main lambda')
               call fatal
            end if
            k = i
         end if
      end do
c
c     van der Waals has to morph between the two charging legs
c
      if (use_relstg(1) .and. use_relstg(3)
     &       .and. .not.use_relstg(2)) then
         write (iout,103)
  103    format (/,' MUTATE_CHECK  --  REL-STAGE with the LIG2 and',
     &              ' LIG1 legs requires the VDWM leg between them')
         call fatal
      end if
c
c     legs that share a main lambda meet with a flat slope only on
c     the quintic map
c
      if (nstg.gt.1 .and.
     &    (elmdamap.ne.'QNT' .or. vlmdamap.ne.'QNT')) then
         write (iout,104)
  104    format (/,' MUTATE_CHECK  --  REL-STAGE with more than one',
     &              ' leg requires the QNT lambda map; remove the',
     &              ' ELE-LMDA-MAP and VDW-LMDA-MAP keywords')
         call fatal
      end if
c
c     a second ligand group makes the free energy a relative one,
c     which is only computed by the staged relative schedule
c
      if (nmutb.gt.0 .and. .not.use_rel) then
         write (iout,105)
  105    format (/,' MUTATE_CHECK  --  Relative Free Energy Requires',
     &              ' the Staged Schedule; add the REL-STAGE Keyword')
         call fatal
      end if
c
c     only the buffered 14-7 and Lennard-Jones van der Waals forms
c     scale a mutated pair with lambda
c
      if (use_mutate .and. use_vdw .and. vdwtyp.ne.'BUFFERED-14-7'
     &       .and. vdwtyp.ne.'LENNARD-JONES') then
         write (iout,109)
  109    format (/,' MUTATE_CHECK  --  Mutated Atoms are only',
     &              ' Available for the Buffered 14-7 and',
     &              ' Lennard-Jones van der Waals Potentials')
         call fatal
      end if
c
c     only the buffered 14-7 van der Waals potential and the atomic
c     multipole electrostatics carry a lambda derivative
c
      if (use_vdlmda .and. use_vdw .and.
     &    vdwtyp.ne.'BUFFERED-14-7') then
         write (iout,106)
  106    format (/,' MUTATE_CHECK  --  Lambda Derivatives are only',
     &              ' Available for the Buffered 14-7 van der Waals',
     &              ' Potential')
         call fatal
      end if
      if (use_edlmda .and.
     &    (use_charge.or.use_dipole.or.use_chgdpl)) then
         write (iout,107)
  107    format (/,' MUTATE_CHECK  --  Lambda Derivatives are not',
     &              ' Available for Partial Charge or Bond Dipole',
     &              ' Electrostatics')
         call fatal
      end if
c
c     relative free energy handles only the buffered 14-7 form of
c     the van der Waals potential
c
      if (use_rel .and. use_vdw .and. vdwtyp.ne.'BUFFERED-14-7') then
         write (iout,110)
  110    format (/,' MUTATE_CHECK  --  Relative Free Energy is only',
     &              ' Available for the Buffered 14-7 van der Waals',
     &              ' Potential')
         call fatal
      end if
c
c     partial charge and bond dipole electrostatics have no relative
c     free energy treatment
c
      if (use_rel .and. (use_charge .or. use_dipole)) then
         write (iout,111)
  111    format (/,' MUTATE_CHECK  --  Relative Free Energy is not',
     &              ' Available for Partial Charge or Bond Dipole',
     &              ' Electrostatics')
         call fatal
      end if
c
c     the HIPPO repulsion and dispersion terms have no relative free
c     energy treatment
c
      if (use_rel .and. (use_repel .or. use_disp)) then
         write (iout,112)
  112    format (/,' MUTATE_CHECK  --  Relative Free Energy is not',
     &              ' Available for HIPPO Repulsion or Dispersion')
         call fatal
      end if
c
c     implicit solvation parameters are not scaled per ligand group
c
      if (use_rel .and. use_solv) then
         write (iout,113)
  113    format (/,' MUTATE_CHECK  --  Relative Free Energy is not',
     &              ' Available with Implicit Solvation')
         call fatal
      end if
c
c     the ost and abf sample interval must keep a propagation phase,
c     an equilibration phase and samples to average at the fixed
c     lambda, so equilibration ratios that leave too few samples or an
c     interval too short for all three phases are reset to defaults
c
      if (use_ost .or. use_abf) then
         if (lmdaparatio+lmdapbratio .ge. 0.9d0) then
            lmdaparatio = 0.3d0
            lmdapbratio = 0.3d0
         end if
         if (lmdaintv .lt. 3)  lmdaintv = 100
         call setlmdaphase
      end if
c
c     dual topology leaves the induced dipoles of the last end state
c     computed, not of the current lambda, so saving them or anything
c     built from them is refused; the save keywords are read here as
c     "mdinit" and "optinit" parse them only after this check
c
      if (use_epdt .and. use_polar) then
         polsave = .false.
         do i = 1, nkey
            next = 1
            record = keyline(i)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:13).eq.'SAVE-UINDUCE ' .or.
     &          keyword(1:13).eq.'SAVE-UDIRECT ' .or.
     &          keyword(1:13).eq.'SAVE-USYSTEM ' .or.
     &          keyword(1:13).eq.'SAVE-TEFIELD ' .or.
     &          keyword(1:13).eq.'SAVE-DEFIELD ')  polsave = .true.
         end do
         if (polsave) then
            write (iout,140)
  140       format (/,' MUTATE_CHECK  --  Induced Dipoles cannot be',
     &                 ' Saved with Dual Topology Polarization; remove',
     &                 ' the SAVE-UINDUCE, SAVE-UDIRECT, SAVE-USYSTEM,',
     &                 ' SAVE-TEFIELD and SAVE-DEFIELD keywords')
            call fatal
         end if
      end if
c
c     the implicit solvation parameters are scaled once at setup,
c     and so do not follow the main lambda to its trial values
c
      if (use_fep .and. use_solv) then
         write (iout,150)
  150    format (/,' MUTATE_CHECK  --  Free Energy Perturbation is',
     &              ' not Available with Implicit Solvation')
         call fatal
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine setligand  --  register a ligand hybrid atom  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "setligand" appends atom "j" to the list of mutated hybrid atoms
c     as a member of alchemical group "igrp" (1 for the first ligand,
c     2 for the second ligand of a relative dual topology calculation)
c
c
      subroutine setligand (j,igrp)
      use atomid
      use atoms
      use mutant
      implicit none
      integer j,igrp
c
c
      nmut = nmut + 1
      imut(nmut) = j
      mutg(j) = igrp
      if (igrp .eq. 2) then
         type0(nmut) = type(j)
         type1(nmut) = 0
         class0(nmut) = class(j)
         class1(nmut) = 0
         nmutb = nmutb + 1
      else
         type0(nmut) = 0
         type1(nmut) = type(j)
         class0(nmut) = 0
         class1(nmut) = class(j)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altelec  --  mutated electrostatic parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altelec" constructs the electrostatic parameters of every site
c     from their original values and the scale that "grpscale" gives
c     the site's group, so the charging ligand carries "elambda", the
c     other ligand of a staged relative leg is annihilated, and the
c     environment keeps its unscaled values
c
c     the charge flux parameters stay unscaled, and "alterchg" scales
c     the flux change in charge of each site by the same scale, so
c     every charge is linear in the scale of its own group
c
c     the exchange polarization size "prepep" is scaled linearly with
c     the polarizability, so each pairwise overlap carries the product
c     of the two group scales and vanishes for an annihilated site
c
c     note charge transfer electrostatics is not treated by parameter
c     scaling due to the functional form used; instead the potential
c     routines scale each pair by the product of the "grpscale" scales
c     of its two groups at the polarization lambda, as the multipoles
c     are, or under dual topology by the polarization weight between
c     the pair scales of the two endpoint states, as given by
c     "ctlmdasc"; the staged relative legs follow from the group
c     scales, and the lambda derivative is added to that of the
c     polarization term; pairs within the ligand are always
c     annihilated, whatever the van der Waals coupling
c
c
      subroutine altelec
      use charge
      use dipole
      use dlmda
      use mutant
      use potent
      implicit none
      integer i,k
      integer k1,k2
      real*8 emsc(0:2)
      real*8 demsc(0:2)
c
c
c     get the electrostatic scale of the environment and each ligand
c
      call grpscale (elambda,emsc,demsc)
c
c     set scaled parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            pchg(k) = pchgorig(k) * emsc(mutg(k))
            pchg0(k) = pchg(k)
         end do
      end if
c
c     set scaled parameters for bond dipole models
c
      if (use_dipole) then
         do i = 1, ndipole
            k1 = idpl(1,i)
            k2 = idpl(2,i)
            bdpl(i) = bdplorig(i) * emsc(max(mutg(k1),mutg(k2)))
         end do
      end if
c
c     set scaled multipole, polarizability and charge flux parameters
c
      call altgrpsc (emsc,use_mpole.or.use_polar)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine altemdt  --  dual topology end state reset  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "altemdt" switches the electrostatic parameters to the state
c     given by the lambda value "elmda", as needed by the multipole
c     dual topology energy routines
c
c     the charge flux monopoles are recomputed, and the multipoles
c     are checked for chirality inversion and rotated into the global
c     frame, so that the global frame multipoles are left consistent
c     with the requested state for any later energy term
c
c
      subroutine altemdt (elmda)
      use mutant
      use potent
      implicit none
      real*8 elmda
c
c
c     set electrostatic parameters for the requested lambda state
c
      elambda = elmda
      call altelec
      if (use_chgflx)  call alterchg
c
c     get global frame multipoles for the requested lambda state
c
      call chkpole
      call rotpole ('MPOLE')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine altpolr  --  mutated polarization parameters  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "altpolr" constructs the polarization parameters of every site
c     from their original values and the scale that "grpscale" gives
c     the site's group at "plambda", so the charging ligand carries
c     "plambda" and the other ligand of a staged relative leg is
c     annihilated; charge flux and the exchange polarization size are
c     scaled as in "altelec"
c
c
      subroutine altpolr
      use mutant
      use potent
      implicit none
      real*8 sc(0:2)
      real*8 dsc(0:2)
c
c
c     get the polarization scale of the environment and each ligand
c
      call grpscale (plambda,sc,dsc)
c
c     set scaled multipole, polarizability and charge flux parameters
c
      call altgrpsc (sc,use_polar)
c
c     update monopoles for charge flux at the requested lambda state
c
      if (use_chgflx)  call alterchg
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altgrpsc  --  install group scaled parameters  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altgrpsc" constructs the multipole, charge penetration,
c     polarizability, exchange polarization size and charge flux
c     scale of every site from their original values and the scale
c     "sc" of the site's group, where the multipoles are only set if
c     "dopole" is true
c
c
      subroutine altgrpsc (sc,dopole)
      use atoms
      use chgpen
      use dlmda
      use expol
      use mplpot
      use mpole
      use mutant
      use polar
      use polpot
      use potent
      implicit none
      integer i,j,k
      real*8 psc
      real*8 sc(0:2)
      logical dopole
c
c
c     set scaled parameters for atomic multipole models
c
      if (dopole) then
         do i = 1, npole
            k = ipole(i)
            psc = sc(mutg(k))
            do j = 1, 13
               pole(j,k) = poleorig(j,k) * psc
            end do
            mono0(k) = pole(1,k)
            if (use_chgpen) then
               pcore(k) = pcoreorig(k) * psc
               pval(k) = pvalorig(k) * psc
               pval0(k) = pval(k)
            end if
         end do
      end if
c
c     set scaled parameters for atomic polarizability models
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            psc = sc(mutg(k))
            polarity(k) = polarityorig(k) * psc
            douind(k) = douindorig(k)
            if (psc .eq. 0.0d0)  douind(k) = .false.
            if (use_expol)  prepep(k) = prepeporig(k) * psc
         end do
      end if
c
c     set the scale "alterchg" applies to each charge flux change
c
      if (use_chgflx) then
         do i = 1, n
            fluxsc(i) = sc(mutg(i))
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine altepdt  --  polarization lambda state reset  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "altepdt" switches the electrostatic parameters to the state
c     given by the polarization lambda value "plmda", as needed by
c     the dual topology polarization energy routines
c
c     the multipoles are checked for chirality inversion and rotated
c     into the global frame, since the polarization routines do not
c     rotate when the multipole term is in use
c
c
      subroutine altepdt (plmda)
      use mutant
      implicit none
      real*8 plmda
c
c     set polarization parameters for the requested lambda state
c
      plambda = plmda
      call altpolr
c
c     get global frame multipoles for the requested lambda state
c
      call chkpole
      call rotpole ('MPOLE')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altepset  --  polarization lambda state setup  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altepset" readies the electrostatic parameters for polarization
c     at "plambda"; when it matches "elambda", the state installed by
c     "altelec" already holds them, so only the charge flux monopoles
c     and the global frame multipoles are refreshed and "same" tells
c     the caller that no restore is needed afterward
c
c
      subroutine altepset (same)
      use mutant
      use potent
      implicit none
      logical same
c
c
c     reinstall the parameters only when the two lambdas differ
c
      same = (abs(plambda-elambda) .le. 1.0d-12)
      if (same) then
         if (use_chgflx)  call alterchg
         call chkpole
         call rotpole ('MPOLE')
      else
         call altepdt (plambda)
      end if
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine alteprst  --  restore electrostatic lambda  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "alteprst" restores the electrostatic parameters to the state
c     installed by "altelec" for the current value of "elambda", so
c     that a rescale performed by "altepdt" leaves no trace for any
c     later energy term
c
c
      subroutine alteprst
      use mutant
      implicit none
      real*8 elmdaorig
c
c
c     reinstall the electrostatic lambda state, which also restores
c     every parameter that "altpolr" scaled
c
      elmdaorig = elambda
      call altemdt (elmdaorig)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine alttors  --  mutated torsional parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "alttors" constructs mutated torsional parameters based
c     on the lambda mutation parameter "tlambda"
c
c
      subroutine alttors (ntbnd,itbnd)
      use mutant
      use potent
      use tors
      implicit none
      integer i,j
      integer ia,ib,ic,id
      integer kb,kc
      integer ntbnd
      integer itbnd(2,*)
c
c
c     set scaled parameters for specified rotatable bonds
c
      if (use_tors) then
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (mutg(ia).ne.0 .and. mutg(ib).ne.0 .and.
     &          mutg(ic).ne.0 .and. mutg(id).ne.0) then
               do j = 1, ntbnd
                  kb = itbnd(1,j)
                  kc = itbnd(2,j)
                  if ((kb.eq.ib .and. kc.eq.ic) .or.
     &                (kb.eq.ic .and. kc.eq.ib)) then
                     tors1(1,i) = tors1(1,i) * tlambda
                     tors2(1,i) = tors2(1,i) * tlambda
                     tors3(1,i) = tors3(1,i) * tlambda
                     tors4(1,i) = tors4(1,i) * tlambda
                     tors5(1,i) = tors5(1,i) * tlambda
                     tors6(1,i) = tors6(1,i) * tlambda
                  end if
               end do
            end if
         end do
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine altsolv  --  mutated solvation parameters  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "altsolv" constructs mutated implicit solvation parameters
c     based on the lambda mutation parameter "elambda"
c
c
      subroutine altsolv
      use atoms
      use mutant
      use nonpol
      use potent
      use solute
      implicit none
      integer i
c
c
c     set scaled parameters for implicit solvation models
c
      if (use_solv) then
         do i = 1, n
            if (mutg(i) .ne. 0) then
               shct(i) = shct(i) * elambda
               radcav(i) = radcav(i) * elambda
               raddsp(i) = raddsp(i) * elambda
               epsdsp(i) = epsdsp(i) * elambda
               cdsp(i) = cdsp(i) * elambda
            end if
         end do
      end if
      return
      end
