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
      use angbnd
      use atomid
      use atoms
      use bndstr
      use cflux
      use charge
      use chgpen
      use dipole
      use dlmda
      use inform
      use iounit
      use katoms
      use keys
      use math
      use mplpot
      use mpole
      use mutant
      use ost
      use polar
      use potent
      implicit none
      integer i,j,k,ihyb
      integer it0,it1
      integer k1,k2
      integer igrp
      integer next,size
      integer ntbnd
      integer, allocatable :: list(:)
      integer, allocatable :: itbnd(:,:)
      real*8 eps
      real*8 temp
      logical setplambda
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
      vlambda = 1.0d0
      elambda = 1.0d0
      tlambda = 1.0d0
c
c     set defaults for lambda scaling for lambda derivatives
c
      setplambda = .false.
      plambda = 1.0d0
c
c     set defaults for vdw coupling type and soft core vdw
c
      vcouple = 0
      scexp = 5.0d0
      scalpha = 0.7d0
c
c     flag for use of lambda derivative
c
      use_dlmda = .false.
      use_emdt = .false.
      use_evdt = .false.
      use_epdt = .false.
      use_plmda = .false.
      use_ost = .false.
      use_ostdyn = .false.
      ostrestart = .false.
      use_meta = .false.
      use_metadyn = .false.
      metarestart = .false.
c
c     set defaults for dual topology
c
      emdtexp = 2
      evdtexp = 2
      epdtexp = 2
c
c     set default ost update intervals
c
      iost = 0
      iosthist = 10
      ostnequil = 5
      ostnavg = 5
      ostddgdl = 0.0d0
      ostdgdl = 0.0d0
      osteqratio = 0.5d0
c
c     set default mapping from main lambda to sublambda
c
      ostplmda1 = 1.0d0
      ostplmda0 = 0.5d0
      ostelmda1 = 0.8d0
      ostelmda0 = 0.3d0
      ostvlmda1 = 0.5d0
      ostvlmda0 = 0.0d0
      ostpmap = 'QNT'
      ostemap = 'QNT'
      ostvmap = 'QNT'
      ostepexp = 1
      ostemexp = 1
      ostevexp = 1
      ostinvepn = 4
      ostinvemn = 4
      ostinvevn = 4
      ostinvepeps = 0.3d0
      ostinvemeps = 0.3d0
      ostinveveps = 0.3d0
      ostlambda = 1.0d0
      ostlambdaavg = 0.0d0
      ostlambdastd = 0.0d0
      ostdedlavg = 0.0d0
      ostdedlstd = 0.0d0
      ostdedl = 0.0d0
      deffdl = 0.0d0
      osttheta = pi / 2.0d0
      ostvtheta = 0.0d0
      ostmass = 25.0d0
      ostfriction = 0.01d0
      ostdt = 0.001d0
c
c     set default flag to compute specific lambda deriv
c
      use_pol4i = .false.
      use_pol4f = .false.
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
      nmetahist = 0
      fastkernel = .true.
      ostinterpol = .false.
c
c     zero out number of hybrid atoms and mutated torsions
c
      nmut = 0
      nmutb = 0
      use_rel = .false.
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
         else if (keyword(1:11) .eq. 'VDW-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  vlambda
         else if (keyword(1:11) .eq. 'ELE-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  elambda
         else if (keyword(1:11) .eq. 'POL-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  plambda
            setplambda = .true.
         else if (keyword(1:12) .eq. 'TORS-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  tlambda
         else if (keyword(1:15) .eq. 'VDW-ANNIHILATE ') then
            vcouple = 1
         else if (keyword(1:13) .eq. 'ELE-DUALTOPO ') then
            use_emdt = .true.
         else if (keyword(1:17) .eq. 'ELE-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  emdtexp
         else if (keyword(1:13) .eq. 'POL-DUALTOPO ') then
            use_epdt = .true.
         else if (keyword(1:17) .eq. 'POL-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  epdtexp
         else if (keyword(1:13) .eq. 'VDW-DUALTOPO ') then
            use_evdt = .true.
         else if (keyword(1:17) .eq. 'VDW-DUALTOPO-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  evdtexp
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
         else if (keyword(1:13) .eq. 'LAMBDA-DERIV ') then
            use_dlmda = .true.
         else if (keyword(1:4) .eq. 'OST ') then
            use_dlmda = .true.
            use_ost = .true.
         else if (keyword(1:8) .eq. 'METADYN ') then
            use_dlmda = .true.
            use_meta = .true.
         else if (keyword(1:17) .eq. 'OSTHIST-INTERVAL ') then
            string = record(next:240)
            read (string,*,err=30)  iosthist
         else if (keyword(1:15) .eq. 'OSTEQUIL-RATIO ') then
            string = record(next:240)
            read (string,*,err=30)  osteqratio
         else if (keyword(1:15) .eq. 'POL-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=30)  ostplmda0, ostplmda1
         else if (keyword(1:15) .eq. 'ELE-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=30)  ostelmda0, ostelmda1
         else if (keyword(1:15) .eq. 'VDW-LMDA-RANGE ') then
            string = record(next:240)
            read (string,*,err=30)  ostvlmda0, ostvlmda1
         else if (keyword(1:13) .eq. 'POL-LMDA-MAP ') then
            call getword (record,ostpmap,next)
            call upcase (ostpmap)
         else if (keyword(1:13) .eq. 'ELE-LMDA-MAP ') then
            call getword (record,ostemap,next)
            call upcase (ostemap)
         else if (keyword(1:13) .eq. 'VDW-LMDA-MAP ') then
            call getword (record,ostvmap,next)
            call upcase (ostvmap)
         else if (keyword(1:13) .eq. 'POL-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  ostepexp
         else if (keyword(1:13) .eq. 'ELE-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  ostemexp
         else if (keyword(1:13) .eq. 'VDW-LMDA-EXP ') then
            string = record(next:240)
            read (string,*,err=30)  ostevexp
         else if (keyword(1:15) .eq. 'POL-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=30)  ostinvepn
         else if (keyword(1:15) .eq. 'ELE-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=30)  ostinvemn
         else if (keyword(1:15) .eq. 'VDW-LMDA-INV-N ') then
            string = record(next:240)
            read (string,*,err=30)  ostinvevn
         else if (keyword(1:17) .eq. 'POL-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=30)  ostinvepeps
         else if (keyword(1:17) .eq. 'ELE-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=30)  ostinvemeps
         else if (keyword(1:17) .eq. 'VDW-LMDA-INV-EPS ') then
            string = record(next:240)
            read (string,*,err=30)  ostinveveps
         else if (keyword(1:11) .eq. 'OST-LAMBDA ') then
            string = record(next:240)
            read (string,*,err=30)  ostlambda
         else if (keyword(1:8) .eq. 'OST-DT ') then
            string = record(next:240)
            read (string,*,err=30)  ostdt
         else if (keyword(1:9) .eq. 'OST-MASS ') then
            string = record(next:240)
            read (string,*,err=30)  ostmass
         else if (keyword(1:13) .eq. 'OST-FRICTION ') then
            string = record(next:240)
            read (string,*,err=30)  ostfriction
         else if (keyword(1:12) .eq. 'LAMBDA-NBIN ') then
            string = record(next:240)
            read (string,*,err=30)  nlmda
         else if (keyword(1:14) .eq. 'FLAMBDA-WIDTH ') then
            string = record(next:240)
            read (string,*,err=30)  wflmda
         else if (keyword(1:7) .eq. 'WLHIST ') then
            string = record(next:240)
            read (string,*,err=30)  wlhist
         else if (keyword(1:7) .eq. 'WFHIST ') then
            string = record(next:240)
            read (string,*,err=30)  wfhist
         else if (keyword(1:11) .eq. 'OST-STDDEV ') then
            string = record(next:240)
            read (string,*,err=30)  oststdev
         else if (keyword(1:16) .eq. 'OST-INTERPOLATE ') then
            ostinterpol = .true.
            fastkernel = .true.
         else if (keyword(1:6) .eq. 'HBIAS ') then
            string = record(next:240)
            read (string,*,err=30)  hbias
         end if
   30    continue
      end do
c
c     only one adaptive lambda bias method can be active at a time
c
      if (use_ost .and. use_meta) then
         write (iout,40)
   40    format (/,' MUTATE  --  OST and METADYN cannot both be active')
         call fatal
      end if
c
c     lambda derivatives of polarization require dual topology
c
      if (use_dlmda)  use_epdt = .true.
      if (use_epdt) then
         use_pol4i = .true.
         use_pol4f = .true.
      end if
c
c     validate mapping schemes from main lambda to sublambdas
c
      if (ostpmap.ne.'QNT' .and. ostpmap.ne.'EXP'
     &       .and. ostpmap.ne.'INV') then
         ostpmap = 'QNT'
      end if
      if (ostemap.ne.'QNT' .and. ostemap.ne.'EXP'
     &       .and. ostemap.ne.'INV') then
         ostemap = 'QNT'
      end if
      if (ostvmap.ne.'QNT' .and. ostvmap.ne.'EXP'
     &       .and. ostvmap.ne.'INV') then
         ostvmap = 'QNT'
      end if
      if (emdtexp .lt. 1) then
         emdtexp = 1
      end if
      if (evdtexp .lt. 1) then
         evdtexp = 1
      end if
      if (ostepexp .lt. 1) then
         ostepexp = 1
      end if
      if (ostemexp .lt. 1) then
         ostemexp = 1
      end if
      if (ostevexp .lt. 1) then
         ostevexp = 1
      end if
      if (ostinvepn .lt. 1) then
         ostinvepn = 1
      end if
      if (ostinvemn .lt. 1) then
         ostinvemn = 1
      end if
      if (ostinvevn .lt. 1) then
         ostinvevn = 1
      end if
      if (ostinvepeps .lt. 0.0d0) then
         ostinvepeps = -ostinvepeps
      end if
      if (ostinvemeps .lt. 0.0d0) then
         ostinvemeps = -ostinvemeps
      end if
      if (ostinveveps .lt. 0.0d0) then
         ostinveveps = -ostinveveps
      end if
c
c     set plambda to elambda if no values given
c
      if (.not. setplambda)  plambda = elambda
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
c     bound ost equilibration ratio to leave at least one sample
c
      if (iosthist .lt. 1)  iosthist = 1
      if (osteqratio .lt. 0.0d0) then
         osteqratio = 0.0d0
      else if (osteqratio .ge. 1.0d0) then
         osteqratio = 1.0d0 - 1.0d0/dble(iosthist)
      end if
      ostnequil = int(osteqratio*dble(iosthist))
      ostnequil = max(0,min(ostnequil,iosthist-1))
      ostnavg = iosthist - ostnequil
c
c     allocate ost histogram and kernels
c
      if (use_ost) then
         if (allocated(osthhist))  deallocate (osthhist)
         if (allocated(osthist))  deallocate (osthist)
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
         sizeosthist = 10000
         nosthist = 0
         allocate (osthhist(sizeosthist))
         allocate (osthist(sizeosthist))
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
c
c     initialize ost histogram and kernels
c
         do i = 1, nlmda
            fkernel(i) = 0.0d0
            fsumkernel(i) = 0.0d0
            pfkernel(i) = 0.0d0
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
            ostnext(i) = 0
            ostlhist(i) = 0.0d0
            ostfhist(i) = 0.0d0
            osthhist(i) = 0.0d0
            ostwlhist(i) = 0.0d0
            ostwfhist(i) = 0.0d0
         end do
      end if
c
c     allocate metadynamics gaussian history
c
      if (use_meta) then
         if (allocated(metalhist))  deallocate (metalhist)
         if (allocated(metahhist))  deallocate (metahhist)
         if (allocated(metawhist))  deallocate (metawhist)
         sizemetahist = 10000
         nmetahist = 0
         allocate (metalhist(sizemetahist))
         allocate (metahhist(sizemetahist))
         allocate (metawhist(sizemetahist))
         do i = 1, sizemetahist
            metalhist(i) = 0.0d0
            metahhist(i) = 0.0d0
            metawhist(i) = 0.0d0
         end do
      end if
c
c     check sublambda intervals are in [0,1] and ordered
c
      if (ostplmda0 .lt. 0.0d0)  ostplmda0 = 0.0d0
      if (ostplmda1 .lt. 0.0d0)  ostplmda1 = 0.0d0
      if (ostelmda0 .lt. 0.0d0)  ostelmda0 = 0.0d0
      if (ostelmda1 .lt. 0.0d0)  ostelmda1 = 0.0d0
      if (ostvlmda0 .lt. 0.0d0)  ostvlmda0 = 0.0d0
      if (ostvlmda1 .lt. 0.0d0)  ostvlmda1 = 0.0d0
      if (ostlambda .lt. 0.0d0)  ostlambda = 0.0d0
      if (ostplmda0 .gt. 1.0d0)  ostplmda0 = 1.0d0
      if (ostplmda1 .gt. 1.0d0)  ostplmda1 = 1.0d0
      if (ostelmda0 .gt. 1.0d0)  ostelmda0 = 1.0d0
      if (ostelmda1 .gt. 1.0d0)  ostelmda1 = 1.0d0
      if (ostvlmda0 .gt. 1.0d0)  ostvlmda0 = 1.0d0
      if (ostvlmda1 .gt. 1.0d0)  ostvlmda1 = 1.0d0
      if (ostlambda .gt. 1.0d0)  ostlambda = 1.0d0
      if (.not. ostrestart .and. .not. metarestart)
     &   osttheta = asin(sqrt(ostlambda))
      if (ostplmda1 .lt. ostplmda0) then
         temp = ostplmda0
         ostplmda0 = ostplmda1
         ostplmda1 = temp
      end if
      if (ostelmda1 .lt. ostelmda0) then
         temp = ostelmda0
         ostelmda0 = ostelmda1
         ostelmda1 = temp
      end if
      if (ostvlmda1 .lt. ostvlmda0) then
         temp = ostvlmda0
         ostvlmda0 = ostvlmda1
         ostvlmda1 = temp
      end if
c
c     get mapping from main lambda to sub-lambdas
c
      if (use_ost .or. use_meta)  call mapsublmda
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
            k1 = idpl(1,i)
            k2 = idpl(2,i)
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
      end if
      if (use_chgflx) then
         do i = 1, nangle
            aflxorig(1,i) = aflx(1,i)
            aflxorig(2,i) = aflx(2,i)
            abflxorig(1,i) = abflx(1,i)
            abflxorig(2,i) = abflx(2,i)
         end do
      end if
c
c     scale electrostatic parameter values based on lambda
c
c     enable two-ligand relative dual topo if a second ligand exists;
c     in that mode the resting parameters are left at their unscaled
c     values and each subsystem state is built on demand by the combiner
c
      use_rel = (nmutb .gt. 0)
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
c     turn off hybrid potentials if no sites are mutated
c
      use_mutate = .true.
      if (nmut .eq. 0)  use_mutate = .false.
c
c     enable rescale of the electrostatic parameters inside the
c     polarization routines when the polarization lambda has been
c     decoupled from the electrostatics lambda; when the two are
c     equal, "altelec" has already installed the correct values
c
      if (use_mutate .and. .not.use_rel .and. .not.use_epdt
     &       .and. use_polar) then
         if (plambda .ne. elambda)  use_plmda = .true.
      end if
c
c     write status of current hybrid potential lambda values
c
      if (use_mutate .and. .not.silent) then
         write (iout,50)
   50    format (/,' Free Energy Perturbation Parameters :')
         write (iout,60)  nmut,vlambda,elambda,plambda,tlambda
   60    format (/,' Number of FEP Hybrid Atoms',9x,i8,
     &           /,' van der Waals Lambda Value',9x,f8.3,
     &           /,' Electrostatics Lambda Value',8x,f8.3,
     &           /,' Polarization Lambda Value',10x,f8.3,
     &           /,' Torsion Angle Lambda Value',9x,f8.3)
         if (use_dlmda) then
            write (iout,70)
   70       format (/,' Lambda Sampling Mode',8x,'Lambda Dynamics')
         else
            write (iout,80)
   80       format (/,' Lambda Sampling Mode',11x,'Fixed Lambda')
         end if
         if (use_plmda) then
            write (iout,90)
   90       format (' Polarization Lambda Decoupled from',
     &              ' Electrostatics Lambda')
         end if
         if (use_rel) then
            write (iout,100)  nmut-nmutb,nmutb
  100       format (/,' Relative Dual Topology Active :',
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
      if (use_mpole) then
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
      use dlmda
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
