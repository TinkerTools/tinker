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
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(imut))  deallocate (imut)
      if (allocated(type0))  deallocate (type0)
      if (allocated(class0))  deallocate (class0)
      if (allocated(type1))  deallocate (type1)
      if (allocated(class1))  deallocate (class1)
      if (allocated(mut))  deallocate (mut)
      if (allocated(mutg))  deallocate (mutg)
      if (allocated(subon))  deallocate (subon)
      allocate (imut(n))
      allocate (type0(n))
      allocate (class0(n))
      allocate (type1(n))
      allocate (class1(n))
      allocate (mut(n))
      allocate (mutg(n))
      allocate (subon(n))
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
      scalpha = 0.7d0
c
c     zero out number of hybrid atoms and mutated torsions
c
      nmut = 0
      nmutb = 0
      use_rel = .false.
      use_past = .false.
      use_subsys = .false.
      do i = 1, n
         mut(i) = .false.
         mutg(i) = 0
         subon(i) = .true.
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
            read (string,*,err=30)  lambda
            use_mainlmda = .true.
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  elambda
            setelambda = .true.
         else if (keyword(1:11) .eq. 'POL-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  plambda
            setplambda = .true.
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  vlambda
            setvlambda = .true.
         else if (keyword(1:12) .eq. 'TORS-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  tlambda
         else if (keyword(1:15) .eq. 'VDW-ANNIHILATE ') then
            vcouple = 1
         else if (keyword(1:7) .eq. 'MUTATE ') then
            string = record(next:240)
            read (string,*,err=30)  ihyb,it0,it1
            nmut = nmut + 1
            imut(nmut) = ihyb
            mut(ihyb) = .true.
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
c     set plambda to elambda if no values given
c
      if (.not. setplambda)  plambda = elambda
c
c     keep the main lambda within the physical lambda range
c
      if (lambda .lt. 0.0d0)  lambda = 0.0d0
      if (lambda .gt. 1.0d0)  lambda = 1.0d0
c
c     a second ligand group makes the free energy a relative one
c
      use_rel = (nmutb .gt. 0)
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
      call mutate_meta
      call mutate_ti
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
         if (use_relstage) then
            write (iout,60)  relstage
   60       format (/,' Free Energy Mode',12x,'Staged Relative',
     &              /,' Staged Leg',29x,a4)
         else if (use_rel) then
            write (iout,70)
   70       format (/,' Free Energy Mode',19x,'Relative')
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
         else
            write (iout,120)
  120       format (' Sampling Mode',18x,'Fixed Lambda')
         end if
         if (use_emdt) then
            write (iout,130)
  130       format (' Electrostatics Topology',16x,'Dual')
         else
            write (iout,140)
  140       format (' Electrostatics Topology',14x,'Single')
         end if
         if (use_epdt) then
            write (iout,150)
  150       format (' Polarization Topology',18x,'Dual')
         else
            write (iout,160)
  160       format (' Polarization Topology',16x,'Single')
         end if
         if (use_evdt) then
            write (iout,170)
  170       format (' van der Waals Topology',17x,'Dual')
         else
            write (iout,180)
  180       format (' van der Waals Topology',15x,'Single')
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
      use angbnd
      use atoms
      use bndstr
      use cflux
      use charge
      use chgpen
      use dipole
      use dlmda
      use iounit
      use keys
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer next
      real*8 temp
      logical setpolmap
      logical setpolrng
      character*4 legword
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     flag for use of lambda derivative
c
      use_dlmda = .false.
      use_elmdamap = .false.
      use_emdt = .false.
      use_epdt = .false.
      use_evdt = .false.
      use_meta = .false.
      use_metadyn = .false.
      use_ost = .false.
      use_ostdyn = .false.
      use_plmda = .false.
      use_plmdamap = .false.
      use_ti = .false.
      use_vlmdamap = .false.
c
c     set defaults describing the flavor of the lambda calculation
c
      lmdaengymode = 'ABS'
      lmdasampmode = 'NONE'
c
c     set defaults for dual topology
c
      emdtexp = 1
      epdtexp = 1
      evdtexp = 1
c
c     interpolate between the two coupled states unless a staged
c     leg says otherwise
c
      erelst0 = rellig2
      erelst1 = rellig1
      prelst0 = rellig2
      prelst1 = rellig1
      vrelst0 = rellig2
      vrelst1 = rellig1
c
c     set defaults for the staged relative free energy schedule
c
      use_relstage = .false.
      relstage = 'VDWM'
      setpolmap = .false.
      setpolrng = .false.
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
c     search keywords for lambda derivative options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:13) .eq. 'LAMBDA-DERIV ') then
            use_dlmda = .true.
         else if (keyword(1:4) .eq. 'OST ') then
            use_dlmda = .true.
            use_ost = .true.
            use_mainlmda = .true.
            lmdasampmode = 'OST'
         else if (keyword(1:8) .eq. 'METADYN ') then
            use_dlmda = .true.
            use_meta = .true.
            use_mainlmda = .true.
            lmdasampmode = 'META'
         else if (keyword(1:11) .eq. 'THERM-INTG ') then
            use_dlmda = .true.
            use_ti = .true.
            use_mainlmda = .true.
            lmdasampmode = 'TI'
         else if (keyword(1:13) .eq. 'ELE-DUALTOPO ') then
            use_emdt = .true.
         else if (keyword(1:17) .eq. 'ELE-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  emdtexp
         else if (keyword(1:13) .eq. 'POL-DUALTOPO ') then
            use_epdt = .true.
         else if (keyword(1:17) .eq. 'POL-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  epdtexp
         else if (keyword(1:13) .eq. 'VDW-DUALTOPO ') then
            use_evdt = .true.
         else if (keyword(1:17) .eq. 'VDW-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  evdtexp
         else if (keyword(1:15) .eq. 'ELE-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=10)  qntelmda0, qntelmda1
         else if (keyword(1:15) .eq. 'POL-LMDA-RANGE ') then
            setpolrng = .true.
            string = record(next:240)
            read (string,*,err=10)  qntplmda0, qntplmda1
         else if (keyword(1:15) .eq. 'VDW-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=10)  qntvlmda0, qntvlmda1
         else if (keyword(1:13) .eq. 'ELE-LMDA-MAP ') then
            use_elmdamap = .true.
            call getword (record,elmdamap,next)
            call upcase (elmdamap)
         else if (keyword(1:13) .eq. 'POL-LMDA-MAP ') then
            use_plmdamap = .true.
            setpolmap = .true.
            call getword (record,plmdamap,next)
            call upcase (plmdamap)
         else if (keyword(1:13) .eq. 'VDW-LMDA-MAP ') then
            use_vlmdamap = .true.
            call getword (record,vlmdamap,next)
            call upcase (vlmdamap)
         else if (keyword(1:13) .eq. 'ELE-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  elmdaexp
         else if (keyword(1:13) .eq. 'POL-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  plmdaexp
         else if (keyword(1:13) .eq. 'VDW-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=10)  vlmdaexp
         else if (keyword(1:15) .eq. 'ELE-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10)  elmdainvn
         else if (keyword(1:15) .eq. 'POL-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10)  plmdainvn
         else if (keyword(1:15) .eq. 'VDW-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=10)  vlmdainvn
         else if (keyword(1:17) .eq. 'ELE-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10)  elmdainveps
         else if (keyword(1:17) .eq. 'POL-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10)  plmdainveps
         else if (keyword(1:17) .eq. 'VDW-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=10)  vlmdainveps
         else if (keyword(1:10) .eq. 'REL-STAGE ') then
            use_rel = .true.
            use_relstage = .true.
            call getword (record,legword,next)
            call upcase (legword)
            relstage = legword
         end if
   10    continue
      end do
c
c     a main lambda drives sublambdas that name a map
c
      if (use_mainlmda .and. use_relstage) then
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
c     enable dual topology for relative free energy
c
      if (use_rel) then
         lmdaengymode = 'REL'
         use_emdt = .true.
         use_epdt = .true.
         use_evdt = .true.
      end if
c
c     ost requires second and force lambda derivatives
c
      if (use_ost)  use_epdt = .true.
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
      if (emdtexp .lt. 1)  emdtexp = 1
      if (epdtexp .lt. 1)  epdtexp = 1
      if (evdtexp .lt. 1)  evdtexp = 1
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
c     a staged run drives one leg, so the leg must be named; the map it
c     walks, the window of that map and the dual topology exponent of
c     the term it drives are all free, as on any other relative leg
c
      if (use_relstage) then
         if (relstage.ne.'LIG1' .and. relstage.ne.'LIG2'
     &          .and. relstage.ne.'VDWM') then
            write (iout,30)
   30       format (/,' MUTATE_DLMDA  --  REL-STAGE requires the leg',
     &                 ' to be named; use LIG2 to discharge ligand 2,',
     &                 ' VDWM to morph van der Waals, or LIG1 to',
     &                 ' charge ligand 1')
            call fatal
         end if
c
c     polarization stages with the multipoles on its own map, so a map
c     or a window given for it would be silently ignored
c
         if (setpolrng .or. setpolmap) then
            write (iout,40)
   40       format (/,' MUTATE_DLMDA  --  REL-STAGE stages',
     &                 ' polarization with the multipoles; remove the',
     &                 ' POL-LMDA-MAP and POL-LMDA-RANGE keywords')
            call fatal
         end if
      end if
c
c     absolute single topology evaluates polarization from one parameter
c     state at plambda, independently of the electrostatic lambda state
c
      use_past = (use_mutate .and. .not.use_rel .and. .not.use_epdt
     &                .and. use_polar)
      use_plmda = use_past
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(pchgorig))  deallocate (pchgorig)
      if (allocated(pchg0orig))  deallocate (pchg0orig)
      if (allocated(bdplorig))  deallocate (bdplorig)
      if (allocated(poleorig))  deallocate (poleorig)
      if (allocated(pcoreorig))  deallocate (pcoreorig)
      if (allocated(pvalorig))  deallocate (pvalorig)
      if (allocated(pval0orig))  deallocate (pval0orig)
      if (allocated(polarityorig))  deallocate (polarityorig)
      if (allocated(bflxorig))  deallocate (bflxorig)
      if (allocated(aflxorig))  deallocate (aflxorig)
      if (allocated(abflxorig))  deallocate (abflxorig)
      if (allocated(douindorig))  deallocate (douindorig)
      allocate (pchgorig(n))
      allocate (pchg0orig(n))
      allocate (bdplorig(nbond))
      allocate (poleorig(maxpole,n))
      allocate (pcoreorig(n))
      allocate (pvalorig(n))
      allocate (pval0orig(n))
      allocate (polarityorig(n))
      allocate (bflxorig(nbond))
      allocate (aflxorig(2,nangle))
      allocate (abflxorig(2,nangle))
      allocate (douindorig(n))
c
c     copy original parameters for lambda derivative calculations
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            pchgorig(k) = pchg(k)
            pchg0orig(k) = pchg0(k)
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
               pval0orig(k) = pval0(k)
            end if
         end do
      end if
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            polarityorig(k) = polarity(k)
            douindorig(k) = douind(k)
         end do
      end if
      if (use_chgflx) then
         do i = 1, nbond
            bflxorig(i) = bflx(i)
         end do
         do i = 1, nangle
            aflxorig(1,i) = aflx(1,i)
            aflxorig(2,i) = aflx(2,i)
            abflxorig(1,i) = abflx(1,i)
            abflxorig(2,i) = abflx(2,i)
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
c     a map, except on a staged relative leg, which walks one window and
c     pins the sublambdas of the other terms to a constant; a pinned
c     sublambda has a flat chain rule, so its term has nothing for the
c     lambda derivative to sample and the plain routine gives the same
c     answer for less work
c
c     "gradient" zeroes every lambda derivative accumulator before it
c     dispatches, so a term left out here reports exact zeros
c
c
      subroutine setdlmdaterms
      use dlmda
      implicit none
c
c
      use_edlmda = use_dlmda .and. use_elmdamap
      use_pdlmda = use_dlmda .and. use_plmdamap
      use_vdlmda = use_dlmda .and. use_vlmdamap
      if (use_relstage) then
         if (relstage .eq. 'VDWM') then
            use_edlmda = .false.
            use_pdlmda = .false.
         else
            use_vdlmda = .false.
         end if
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
c     "mutate_ost" sets the lambda grid, the gaussian deposit interval,
c     the convergence criteria and the lambda particle parameters used
c     by orthogonal space tempering, then allocates the histogram and
c     the bias kernels when the method is active
c
c
      subroutine mutate_ost
      use dlmda
      use keys
      use math
      use mutant
      use ost
      implicit none
      integer i,k
      integer next
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set default ost update intervals
c
      iost = 0
      iosthist = 10
      ostddgdl = 0.0d0
      ostdgdl = 0.0d0
      ostparatio = 0.3d0
      ostpbratio = 0.3d0
      nosthistsave = 0
c
c     set default criteria for judging convergence of a deposit
c
      ostcvbin = 2
      ostcvdif = 25.0d0
      ostcvrat = 0.1d0
      ostcvslp = 1.0d0
      ostcvstd = 10.0d0
c
c     set defaults for tempering of the deposited gaussian heights
c
      ostemper = .false.
      tempergamma = 1.0d0
      temperthresh = 1.0d0
c
c     set defaults for the lambda particle propagation
c
      ostlambdaavg = 0.0d0
      ostlambdastd = 0.0d0
      ostdedlavg = 0.0d0
      ostdedlstd = 0.0d0
      deffdl = 0.0d0
      osttheta = pi / 2.0d0
      ostvtheta = 0.0d0
      ostmass = 25.0d0
      ostfriction = 0.01d0
      ostdt = 0.001d0
c
c     set default ost lambda bin values
c
      nlmda = 201
      nflmda = 1001
      wflmda = 1.0d0
      wlhist = 0.005d0
      wfhist = 1.0d0
      fli0 = (nflmda + 1) / 2 + (nflmda - 1) / 4
      hbias = 0.00001d0
      oststdev = 4.0d0
      eosttot = 0.0d0
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
         if (keyword(1:17) .eq. 'OSTHIST-INTERVAL ') then
            string = record(next:240)
            read (string,*,err=10)  iosthist
         else if (keyword(1:12) .eq. 'OSTPA-RATIO ') then
            string = record(next:240)
            read (string,*,err=10)  ostparatio
         else if (keyword(1:12) .eq. 'OSTPB-RATIO ') then
            string = record(next:240)
            read (string,*,err=10)  ostpbratio
         else if (keyword(1:8) .eq. 'OST-DT ') then
            string = record(next:240)
            read (string,*,err=10)  ostdt
         else if (keyword(1:9) .eq. 'OST-MASS ') then
            string = record(next:240)
            read (string,*,err=10)  ostmass
         else if (keyword(1:13) .eq. 'OST-FRICTION ') then
            string = record(next:240)
            read (string,*,err=10)  ostfriction
         else if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            string = record(next:240)
            read (string,*,err=10)  nlmda
         else if (keyword(1:14) .eq. 'FLAMBDA-WIDTH ') then
            string = record(next:240)
            read (string,*,err=10)  wflmda
         else if (keyword(1:7) .eq. 'WLHIST ') then
            string = record(next:240)
            read (string,*,err=10)  wlhist
         else if (keyword(1:7) .eq. 'WFHIST ') then
            string = record(next:240)
            read (string,*,err=10)  wfhist
         else if (keyword(1:11) .eq. 'OST-STDDEV ') then
            string = record(next:240)
            read (string,*,err=10)  oststdev
         else if (keyword(1:16) .eq. 'OST-INTERPOLATE ') then
            ostinterpol = .true.
            fastkernel = .true.
         else if (keyword(1:6) .eq. 'HBIAS ') then
            string = record(next:240)
            read (string,*,err=10)  hbias
         else if (keyword(1:13) .eq. 'OST-CONV-BIN ') then
            string = record(next:240)
            read (string,*,err=10)  ostcvbin
         else if (keyword(1:16) .eq. 'OST-CONVCRI-DIF ') then
            string = record(next:240)
            read (string,*,err=10)  ostcvdif
         else if (keyword(1:16) .eq. 'OST-CONVCRI-RAT ') then
            string = record(next:240)
            read (string,*,err=10)  ostcvrat
         else if (keyword(1:16) .eq. 'OST-CONVCRI-SLP ') then
            string = record(next:240)
            read (string,*,err=10)  ostcvslp
         else if (keyword(1:16) .eq. 'OST-CONVCRI-STD ') then
            string = record(next:240)
            read (string,*,err=10)  ostcvstd
         else if (keyword(1:11) .eq. 'OST-TEMPER ') then
            ostemper = .true.
         else if (keyword(1:17) .eq. 'OST-TEMPER-GAMMA ') then
            string = record(next:240)
            read (string,*,err=10)  tempergamma
         else if (keyword(1:18) .eq. 'OST-TEMPER-THRESH ') then
            string = record(next:240)
            read (string,*,err=10)  temperthresh
         end if
   10    continue
      end do
c
c     define lambda width and flambda range
c
      if (nlmda .lt. 3)  nlmda = 3
      if (mod(nlmda,2) .eq. 0)  nlmda = nlmda + 1
      wlmda = 1.0d0 / dble(nlmda-1)
      wlmda2 = 0.5d0 * wlmda
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
      maxwlhist = wlhist
      maxwfhist = wfhist
c
c     split the deposit interval into its propagation, equilibration
c     and averaging phases
c
      if (iosthist .lt. 1)  iosthist = 1
      if (ostparatio .lt. 0.0d0)  ostparatio = 0.0d0
      if (ostpbratio .lt. 0.0d0)  ostpbratio = 0.0d0
      call setostphase
c
c     start the lambda particle from the current main lambda
c
      osttheta = asin(sqrt(lambda))
c
c     allocate the convergence sub-bins
c
      if (ostcvbin .lt. 0)  ostcvbin = 0
      if (use_ost .or. use_meta) then
         if (allocated(ostlmdaavgbin))  deallocate (ostlmdaavgbin)
         if (allocated(ostlmdaslpbin))  deallocate (ostlmdaslpbin)
         if (allocated(ostlmdastdbin))  deallocate (ostlmdastdbin)
         if (allocated(ostdedlavgbin))  deallocate (ostdedlavgbin)
         if (allocated(ostdedlslpbin))  deallocate (ostdedlslpbin)
         if (allocated(ostdedlstdbin))  deallocate (ostdedlstdbin)
         allocate (ostlmdaavgbin(max(ostcvbin,1)))
         allocate (ostlmdaslpbin(max(ostcvbin,1)))
         allocate (ostlmdastdbin(max(ostcvbin,1)))
         allocate (ostdedlavgbin(max(ostcvbin,1)))
         allocate (ostdedlslpbin(max(ostcvbin,1)))
         allocate (ostdedlstdbin(max(ostcvbin,1)))
         do i = 1, max(ostcvbin,1)
            ostlmdaavgbin(i) = 0.0d0
            ostlmdaslpbin(i) = 0.0d0
            ostlmdastdbin(i) = 0.0d0
            ostdedlavgbin(i) = 0.0d0
            ostdedlslpbin(i) = 0.0d0
            ostdedlstdbin(i) = 0.0d0
         end do
      end if
c
c     allocate ost histogram and kernels
c
      if (use_ost) then
         if (allocated(osthhist))  deallocate (osthhist)
         if (allocated(osthist))  deallocate (osthist)
         if (allocated(ostihist))  deallocate (ostihist)
         if (allocated(osthead))  deallocate (osthead)
         if (allocated(ostnext))  deallocate (ostnext)
         if (allocated(ostllist))  deallocate (ostllist)
         if (allocated(ostflist))  deallocate (ostflist)
         if (allocated(ostlhist))  deallocate (ostlhist)
         if (allocated(ostfhist))  deallocate (ostfhist)
         if (allocated(ostwlhist))  deallocate (ostwlhist)
         if (allocated(ostwfhist))  deallocate (ostwfhist)
         if (allocated(fkernel))  deallocate (fkernel)
         if (allocated(fsumkernel))  deallocate (fsumkernel)
         if (allocated(gfkernel))  deallocate (gfkernel)
         if (allocated(gkernel))  deallocate (gkernel)
         if (allocated(glfkernel))  deallocate (glfkernel)
         if (allocated(glkernel))  deallocate (glkernel)
         if (allocated(pfkernel))  deallocate (pfkernel)
         if (allocated(vkernelmax))  deallocate (vkernelmax)
         sizeosthist = 10000
         nosthist = 0
         allocate (osthhist(sizeosthist))
         allocate (osthist(sizeosthist))
         allocate (ostihist(sizeosthist))
         allocate (osthead(nlmda,nflmda))
         allocate (ostnext(sizeosthist))
         allocate (ostllist(iosthist))
         allocate (ostflist(iosthist))
         allocate (ostlhist(sizeosthist))
         allocate (ostfhist(sizeosthist))
         allocate (ostwlhist(sizeosthist))
         allocate (ostwfhist(sizeosthist))
         allocate (fkernel(nlmda))
         allocate (fsumkernel(nlmda))
         allocate (gfkernel(nlmda,nflmda))
         allocate (gkernel(nlmda,nflmda))
         allocate (glfkernel(nlmda,nflmda))
         allocate (glkernel(nlmda,nflmda))
         allocate (pfkernel(nlmda))
         allocate (vkernelmax(nlmda))
c
c     initialize ost histogram and kernels
c
         do i = 1, nlmda
            fkernel(i) = 0.0d0
            fsumkernel(i) = 0.0d0
            pfkernel(i) = 0.0d0
            vkernelmax(i) = 0.0d0
            do k = 1, nflmda
               gfkernel(i,k) = 0.0d0
               gkernel(i,k) = 0.0d0
               glfkernel(i,k) = 0.0d0
               glkernel(i,k) = 0.0d0
               osthead(i,k) = 0
            end do
         end do
         do i = 1, iosthist
            ostllist(i) = 0.0d0
            ostflist(i) = 0.0d0
         end do
         do i = 1, sizeosthist
            osthist(i) = 0
            ostihist(i) = 0
            ostnext(i) = 0
            ostlhist(i) = 0.0d0
            ostfhist(i) = 0.0d0
            osthhist(i) = 0.0d0
            ostwlhist(i) = 0.0d0
            ostwfhist(i) = 0.0d0
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
         if (allocated(ostllist))  deallocate (ostllist)
         if (allocated(vmetagrid))  deallocate (vmetagrid)
         if (allocated(dvmetagrid))  deallocate (dvmetagrid)
         sizemetahist = 10000
         allocate (metalhist(sizemetahist))
         allocate (metahhist(sizemetahist))
         allocate (metawhist(sizemetahist))
         allocate (metaihist(sizemetahist))
         allocate (ostllist(iosthist))
         allocate (vmetagrid(nlmda))
         allocate (dvmetagrid(nlmda))
         do i = 1, sizemetahist
            metalhist(i) = 0.0d0
            metahhist(i) = 0.0d0
            metawhist(i) = 0.0d0
            metaihist(i) = 0
         end do
         do i = 1, iosthist
            ostllist(i) = 0.0d0
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
c     ########################################################
c     ##                                                    ##
c     ##  subroutine mutate_ti  --  lambda window schedule  ##
c     ##                                                    ##
c     ########################################################
c
c
c     "mutate_ti" sets the lambda window schedule and the block size
c     used to average dU/dlambda during a thermodynamic integration,
c     either from explicit "TI-WINDOW" values or from the "TI-NBIN"
c     count of evenly spaced windows
c
c
      subroutine mutate_ti
      use dlmda
      use iounit
      use keys
      use thrmint
      implicit none
      integer i
      integer next
      integer ntiwin
      real*8 frac
      real*8 temp
      logical tinbinset
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     set defaults for thermodynamic integration windows
c
      tibin = 0
      tinbin = 21
      tinblock = 0
      tinbcount = 0
      tinbsave = 0
      tinbtot = 0
      tinequil = 0
      tinstepavg = 100
      tiwindow = 0
      tieqratio = 0.5d0
c
c     size the lambda window schedule to the worst case
c
      if (allocated(tilmdalist))  deallocate (tilmdalist)
      if (allocated(tifraclist))  deallocate (tifraclist)
      allocate (tilmdalist(max(1,nkey)))
      allocate (tifraclist(max(1,nkey)))
      do i = 1, max(1,nkey)
         tilmdalist(i) = 0.0d0
         tifraclist(i) = -1.0d0
      end do
      ntiwin = 0
      tinbinset = .false.
c
c     search keywords for thermodynamic integration options
c
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'TI-NBIN ') then
            string = record(next:240)
            read (string,*,err=20)  tinbin
            tinbinset = .true.
         else if (keyword(1:10) .eq. 'TI-WINDOW ') then
            temp = 0.0d0
            frac = -1.0d0
            string = record(next:240)
            read (string,*,err=10,end=10)  temp,frac
   10       continue
            ntiwin = ntiwin + 1
            tilmdalist(ntiwin) = temp
            tifraclist(ntiwin) = frac
         else if (keyword(1:12) .eq. 'TI-NSTEPAVG ') then
            string = record(next:240)
            read (string,*,err=20)  tinstepavg
         else if (keyword(1:15) .eq. 'TI-EQUIL-RATIO ') then
            string = record(next:240)
            read (string,*,err=20)  tieqratio
         end if
   20    continue
      end do
c
c     the lambda windows must span [0,1] and leave room to average
c
      if (use_ti) then
         if (tinstepavg .lt. 1) then
            write (iout,30)
   30       format (/,' MUTATE_TI  --  TI-NSTEPAVG must be positive')
            call fatal
         end if
         if (tieqratio.lt.0.0d0 .or. tieqratio.ge.1.0d0) then
            write (iout,40)
   40       format (/,' MUTATE_TI  --  TI-EQUIL-RATIO must be',
     &                 ' in [0,1)')
            call fatal
         end if
         call settisched (ntiwin,tinbinset)
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
      use iounit
      use limits
      use mplpot
      use mutant
      use ost
      use polpot
      use potent
      implicit none
c
c
c     only one method can sample the main lambda at a time
c
      if ((use_ost .and. use_meta) .or. (use_ost .and. use_ti)
     &       .or. (use_meta .and. use_ti)) then
         write (iout,10)
   10    format (/,' MUTATE_CHECK  --  Only one of OST, METADYN and',
     &              ' THERM-INTG can be active')
         call fatal
      end if
c
c     ost requires polarization lambda derivatives not yet available
c     for absolute single topology
c
      if (use_ost .and. use_past) then
         write (iout,20)
   20    format (/,' MUTATE_CHECK  --  OST cannot be used with',
     &              ' absolute single topology polarization; add',
     &              ' POL-DUALTOPO or remove OST')
         call fatal
      end if
c
c     absolute single topology currently supports the scalar derivative
c     for mutual Thole polarization in direct and Ewald modes
c
      if (use_pdlmda .and. use_past) then
         if (poltyp.ne.'MUTUAL' .or. .not.use_thole .or. use_chgpen
     &          .or. use_expol .or. use_solv) then
            write (iout,30)
   30       format (/,' MUTATE_CHECK  --  Absolute Single Topology',
     &                 ' Polarization dU/dLambda supports MUTUAL Thole',
     &                 ' polarization only, in either direct or Ewald',
     &                 ' mode; use POLARIZATION MUTUAL without charge',
     &                 ' penetration, exchange polarization or',
     &                 ' implicit solvent')
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
c     the staged relative schedule maps every sublambda on its own,
c     so a sublambda set by its own keyword would be overwritten
c
      if (use_relstage .and.
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
      if (use_relstage .and. nmutb.eq.0) then
         write (iout,100)
  100    format (/,' MUTATE_CHECK  --  REL-STAGE requires a second',
     &              ' ligand group; add the LIGAND2 keyword')
         call fatal
      end if
c
c     relative free energy requires the isolated ligand van der Waals
c     terms that annihilation would remove from the thermodynamic cycle
c
      if (use_rel .and. vcouple.eq.1) then
         write (iout,110)
  110    format (/,' MUTATE_CHECK  --  VDW-ANNIHILATE is not',
     &              ' compatible with relative free energy; the',
     &              ' isolated-ligand van der Waals terms are needed',
     &              ' to preserve the relative thermodynamic cycle;',
     &              ' remove the VDW-ANNIHILATE keyword')
         call fatal
      end if
c
c     the ost deposit interval must keep a propagation phase, an
c     equilibration phase and samples to average at the fixed lambda
c
      if (use_ost .and. ostparatio+ostpbratio.ge.0.9d0) then
         write (iout,120)
  120    format (/,' MUTATE_CHECK  --  OSTPA-RATIO plus OSTPB-RATIO',
     &              ' must be less than 0.9 to leave samples for the',
     &              ' fixed lambda average')
         call fatal
      end if
      if (use_ost .and. iosthist.lt.3) then
         write (iout,130)
  130    format (/,' MUTATE_CHECK  --  OSTHIST-INTERVAL must be at',
     &              ' least 3 to hold a propagation, equilibration',
     &              ' and averaging phase')
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
      mut(j) = .true.
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
c     "altelec" constructs mutated electrostatic parameters based
c     on the lambda mutation parameter "elambda"
c
c     note charge transfer electrostatics is not treated by parameter
c     scaling due to the functional form used, and must be done via
c     modification of pairwise energy terms in the potential routines
c
c
      subroutine altelec
      use angbnd
      use atoms
      use bndstr
      use cflux
      use charge
      use chgpen
      use dipole
      use dlmda
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer k1,k2
      integer ia,ib,ic
c
c
c     set scaled parameters for partial charge models
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            if (mut(k)) then
               pchg(k) = pchgorig(k) * elambda
            end if
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
            if (mut(k1) .or. mut(k2)) then
               bdpl(i) = bdplorig(i) * elambda
            end if
         end do
      end if
c
c     set scaled parameters for atomic multipole models
c
      if (use_mpole .or. use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,k) = poleorig(j,k) * elambda
               end do
               mono0(k) = pole(1,k)
               if (use_chgpen) then
                  pcore(k) = pcoreorig(k) * elambda
                  pval(k) = pvalorig(k) * elambda
                  pval0(k) = pval(k)
               end if
            end if
         end do
      end if
c
c     set scaled parameters for atomic polarizability models
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               polarity(k) = polarityorig(k) * elambda
               douind(k) = douindorig(k)
               if (elambda .eq. 0.0d0)  douind(k) = .false.
            end if
         end do
      end if
c
c     set scaled parameters for bond stretch charge flux
c
      if (use_chgflx) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (mut(ia) .and. mut(ib)) then
               bflx(i) = bflxorig(i) * elambda
            end if
         end do
      end if
c
c     set scaled parameters for angle bend charge flux
c
      if (use_chgflx) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (mut(ia) .and. mut(ib) .and. mut(ic)) then
               aflx(1,i) = aflxorig(1,i) * elambda
               aflx(2,i) = aflxorig(2,i) * elambda
               abflx(1,i) = abflxorig(1,i) * elambda
               abflx(2,i) = abflxorig(2,i) * elambda
            end if
         end do
      end if
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
c     "altpolr" constructs mutated polarization parameters based
c     on the lambda mutation parameter "plambda"
c
c
      subroutine altpolr
      use angbnd
      use bndstr
      use cflux
      use chgpen
      use dlmda
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer k1,k2
      integer ia,ib,ic
c
c
c     set scaled parameters for atomic multipole models
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,k) = poleorig(j,k) * plambda
               end do
               mono0(k) = pole(1,k)
               if (use_chgpen) then
                  pcore(k) = pcoreorig(k) * plambda
                  pval(k) = pvalorig(k) * plambda
                  pval0(k) = pval(k)
               end if
            end if
         end do
      end if
c
c     set scaled parameters for atomic polarizability models
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               polarity(k) = polarityorig(k) * plambda
               douind(k) = douindorig(k)
               if (plambda .eq. 0.0d0)  douind(k) = .false.
            end if
         end do
      end if
c
c     set scaled parameters for bond stretch charge flux
c
      if (use_chgflx) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (mut(ia) .and. mut(ib)) then
               bflx(i) = bflxorig(i) * plambda
            end if
         end do
      end if
c
c     set scaled parameters for angle bend charge flux
c
      if (use_chgflx) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (mut(ia) .and. mut(ib) .and. mut(ic)) then
               aflx(1,i) = aflxorig(1,i) * plambda
               aflx(2,i) = aflxorig(2,i) * plambda
               abflx(1,i) = abflxorig(1,i) * plambda
               abflx(2,i) = abflxorig(2,i) * plambda
            end if
         end do
      end if
c
c     update monopoles for charge flux at the requested lambda state
c
      if (use_chgflx)  call alterchg
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
c     note "altelec" scales the permanent multipoles under "use_mpole"
c     while "altpolr" scales them under "use_polar"; with the multipole
c     term not in use they were never scaled by "altelec" and must be
c     returned to their unscaled values here
c
c
      subroutine alteprst
      use angbnd
      use bndstr
      use cflux
      use chgpen
      use dlmda
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer ia,ib,ic
      real*8 elmdaorig
c
c
c     multipole term in use, so the multipole reset restores every
c     parameter that "altpolr" scaled
c
      if (use_mpole) then
         elmdaorig = elambda
         call altemdt (elmdaorig)
         return
      end if
c
c     return the permanent multipoles to their unscaled values, and
c     the polarizabilities to the electrostatics lambda state
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (mut(k)) then
               do j = 1, 13
                  pole(j,k) = poleorig(j,k)
               end do
               mono0(k) = pole(1,k)
               if (use_chgpen) then
                  pcore(k) = pcoreorig(k)
                  pval(k) = pvalorig(k)
                  pval0(k) = pval(k)
               end if
               polarity(k) = polarityorig(k) * elambda
               douind(k) = douindorig(k)
               if (elambda .eq. 0.0d0)  douind(k) = .false.
            end if
         end do
      end if
c
c     restore scaled parameters for bond stretch charge flux
c
      if (use_chgflx) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (mut(ia) .and. mut(ib)) then
               bflx(i) = bflxorig(i) * elambda
            end if
         end do
      end if
c
c     restore scaled parameters for angle bend charge flux
c
      if (use_chgflx) then
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (mut(ia) .and. mut(ib) .and. mut(ic)) then
               aflx(1,i) = aflxorig(1,i) * elambda
               aflx(2,i) = aflxorig(2,i) * elambda
               abflx(1,i) = abflxorig(1,i) * elambda
               abflx(2,i) = abflxorig(2,i) * elambda
            end if
         end do
      end if
c
c     update monopoles for charge flux, then get global frame
c     multipoles for the restored state
c
      if (use_chgflx)  call alterchg
      call chkpole
      call rotpole ('MPOLE')
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
            if (mut(ia) .and. mut(ib) .and. mut(ic) .and. mut(id)) then
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
            if (mut(i)) then
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
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine submask  --  select relative subsystem atoms  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "submask" flags which atoms are active in the subsystem currently
c     being built for a relative dual topology energy term; group A is
c     included when "la", group B when "lb", and the environment when
c     "le", with all other alchemical atoms zeroed out
c
c
      subroutine submask (la,lb,le)
      use atoms
      use mutant
      implicit none
      integer i
      logical la,lb,le
c
c
      do i = 1, n
         if (mutg(i) .eq. 1) then
            subon(i) = la
         else if (mutg(i) .eq. 2) then
            subon(i) = lb
         else
            subon(i) = le
         end if
      end do
      use_subsys = .not. (la .and. lb .and. le)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine relslot  --  relative subsystem slot lookup  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "relslot" returns the group mask of the subsystem in slot "k"
c     along with whether that subsystem belongs to the coupling states
c     "ist0" and "ist1" holding the two interpolation endpoints
c
c     only five subsystems are reachable through "submask", and the
c     three coupling states of a relative dual topology are sums of
c     them, the decoupled state being the one that the plain relative
c     schedule never needs and the staged schedule is built around,
c
c        slot   la     lb     le     subsystem
c          1    T      F      T      ligand 1 with environment
c          2    F      T      T      ligand 2 with environment
c          3    F      F      T      environment alone
c          4    T      F      F      ligand 1 alone
c          5    F      T      F      ligand 2 alone
c
c        rellig1 = slots 1 and 5 ,   ligand 1 bound, ligand 2 free
c        rellig2 = slots 2 and 4 ,   ligand 2 bound, ligand 1 free
c        relnone = slots 3, 4 and 5 ,  neither ligand bound
c
c
      subroutine relslot (k,ist0,ist1,la,lb,le,in0,in1)
      implicit none
      integer k,ist0,ist1
      logical la,lb,le
      logical in0,in1
      logical subla(5),sublb(5),suble(5)
      logical relmem(5,3)
      save subla,sublb,suble,relmem
      data subla  / .true., .false.,.false.,.true., .false. /
      data sublb  / .false.,.true., .false.,.false.,.true.  /
      data suble  / .true., .true., .true., .false.,.false. /
      data relmem / .true., .false.,.false.,.false.,.true.,
     &              .false.,.true., .false.,.true., .false.,
     &              .false.,.false.,.true., .true., .true.  /
c
c
      la = subla(k)
      lb = sublb(k)
      le = suble(k)
      in0 = relmem(k,ist0)
      in1 = relmem(k,ist1)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine setsubelec  --  subsystem electrostatic state  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "setsubelec" installs the electrostatic parameters for the atom
c     subsystem flagged by "subon", using full original values for the
c     active atoms and zero for the inactive atoms, then refreshes the
c     charge flux monopoles and global frame multipoles so any later
c     energy term is consistent with the requested subsystem
c
c
      subroutine setsubelec
      use angbnd
      use atoms
      use bndstr
      use cflux
      use charge
      use chgpen
      use dipole
      use dlmda
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer k1,k2
      integer ia,ib,ic
c
c
c     partial charge models
c
      if (use_charge) then
         do i = 1, nion
            k = iion(i)
            if (subon(k)) then
               pchg(k) = pchgorig(k)
            else
               pchg(k) = 0.0d0
            end if
            pchg0(k) = pchg(k)
         end do
      end if
c
c     bond dipole models
c
      if (use_dipole) then
         do i = 1, ndipole
            k1 = idpl(1,i)
            k2 = idpl(2,i)
            if (subon(k1) .and. subon(k2)) then
               bdpl(i) = bdplorig(i)
            else
               bdpl(i) = 0.0d0
            end if
         end do
      end if
c
c     atomic multipole models
c
      if (use_mpole) then
         do i = 1, npole
            k = ipole(i)
            if (subon(k)) then
               do j = 1, 13
                  pole(j,k) = poleorig(j,k)
               end do
               if (use_chgpen) then
                  pcore(k) = pcoreorig(k)
                  pval(k) = pvalorig(k)
                  pval0(k) = pval(k)
               end if
            else
               do j = 1, 13
                  pole(j,k) = 0.0d0
               end do
               if (use_chgpen) then
                  pcore(k) = 0.0d0
                  pval(k) = 0.0d0
                  pval0(k) = 0.0d0
               end if
            end if
            mono0(k) = pole(1,k)
         end do
      end if
c
c     atomic polarizability models
c
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (subon(k)) then
               polarity(k) = polarityorig(k)
               douind(k) = douindorig(k)
            else
               polarity(k) = 0.0d0
               douind(k) = .false.
            end if
         end do
      end if
c
c     bond stretch charge flux
c
      if (use_chgflx) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (subon(ia) .and. subon(ib)) then
               bflx(i) = bflxorig(i)
            else
               bflx(i) = 0.0d0
            end if
         end do
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (subon(ia) .and. subon(ib) .and. subon(ic)) then
               aflx(1,i) = aflxorig(1,i)
               aflx(2,i) = aflxorig(2,i)
               abflx(1,i) = abflxorig(1,i)
               abflx(2,i) = abflxorig(2,i)
            else
               aflx(1,i) = 0.0d0
               aflx(2,i) = 0.0d0
               abflx(1,i) = 0.0d0
               abflx(2,i) = 0.0d0
            end if
         end do
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine altemdtsub  --  subsystem multipole end state  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "altemdtsub" switches the electrostatic parameters to the atom
c     subsystem containing group A when "la", group B when "lb", and
c     the environment when "le"; the charge flux monopoles are updated
c     and the global frame multipoles are rebuilt as in "altemdt"
c
c
      subroutine altemdtsub (la,lb,le)
      use potent
      implicit none
      logical la,lb,le
c
c
      call submask (la,lb,le)
      call setsubelec
      if (use_chgflx)  call alterchg
      call chkpole
      call rotpole ('MPOLE')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine altpolrsub  --  subsystem polarization state  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "altpolrsub" installs the polarizability and multipole values for
c     the atom subsystem flagged by "subon", using full original values
c     for active atoms and zero for inactive atoms, for the polarization
c     dual topology subsystem energies
c
c
      subroutine altpolrsub (la,lb,le)
      use angbnd
      use bndstr
      use cflux
      use chgpen
      use dlmda
      use mplpot
      use mpole
      use mutant
      use polar
      use potent
      implicit none
      integer i,j,k
      integer ia,ib,ic
      logical la,lb,le
c
c
      call submask (la,lb,le)
      if (use_polar) then
         do i = 1, npole
            k = ipole(i)
            if (subon(k)) then
               do j = 1, 13
                  pole(j,k) = poleorig(j,k)
               end do
               if (use_chgpen) then
                  pcore(k) = pcoreorig(k)
                  pval(k) = pvalorig(k)
                  pval0(k) = pval(k)
               end if
               polarity(k) = polarityorig(k)
               douind(k) = douindorig(k)
            else
               do j = 1, 13
                  pole(j,k) = 0.0d0
               end do
               if (use_chgpen) then
                  pcore(k) = 0.0d0
                  pval(k) = 0.0d0
                  pval0(k) = 0.0d0
               end if
               polarity(k) = 0.0d0
               douind(k) = .false.
            end if
            mono0(k) = pole(1,k)
         end do
      end if
c
c     set subsystem parameters for charge flux
c
      if (use_chgflx) then
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (subon(ia) .and. subon(ib)) then
               bflx(i) = bflxorig(i)
            else
               bflx(i) = 0.0d0
            end if
         end do
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (subon(ia) .and. subon(ib) .and. subon(ic)) then
               aflx(1,i) = aflxorig(1,i)
               aflx(2,i) = aflxorig(2,i)
               abflx(1,i) = abflxorig(1,i)
               abflx(2,i) = abflxorig(2,i)
            else
               aflx(1,i) = 0.0d0
               aflx(2,i) = 0.0d0
               abflx(1,i) = 0.0d0
               abflx(2,i) = 0.0d0
            end if
         end do
      end if
c
c     update monopoles for charge flux in the requested subsystem
c
      if (use_chgflx)  call alterchg
c
c     get global frame multipoles for the requested subsystem
c
      call chkpole
      call rotpole ('MPOLE')
      return
      end
