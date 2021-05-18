c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine kchgtrn  --  charge transfer term assignment  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "kchgtrn" assigns charge magnitude and damping parameters for
c     charge transfer interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kchgtrn
      use atomid
      use atoms
      use chgpen
      use chgtrn
      use inform
      use iounit
      use kctrn
      use keys
      use mplpot
      use mpole
      use polar
      use polpot
      use potent
      use sizes
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 chtrn,actrn
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing charge transfer parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            k = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  chtrn,actrn
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Charge Transfer',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',13x,'Charge',11x,'Damp',/)
               end if
               if (k .le. maxclass) then
                  ctchg(k) = chtrn
                  ctdmp(k) = actrn
                  if (.not. silent) then
                     write (iout,30)  k,chtrn,actrn
   30                format (6x,i6,7x,f15.4,f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KCHGTRN  --  Too many Charge',
     &                       ' Transfer Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(chgct))  deallocate (chgct)
      if (allocated(dmpct))  deallocate (dmpct)
      allocate (chgct(n))
      allocate (dmpct(n))
c
c     assign the charge transfer charge and alpha parameters 
c
      nct = n
      do i = 1, n
         ic = class(i)
         chgct(i) = ctchg(ic)
         dmpct(i) = ctdmp(ic)
      end do
c
c     process keywords containing atom specific charge transfer
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:7) .eq. 'CHGTRN ') then
            ia = 0
            chtrn = 0.0d0
            actrn = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,chtrn,actrn
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Charge Transfer Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',16x,'Charge',11x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,chtrn,actrn
   60             format (6x,i6,7x,f15.4,f15.4)
               end if
               chgct(ia) = chtrn
               dmpct(ia) = actrn
            end if
   70       continue
         end if
      end do
c
c     remove zero or undefined electrostatic sites from the list
c
      if (use_chgtrn) then
         npole = 0
         ncp = 0
         npolar = 0
         nct = 0
         do i = 1, n
            if (polarity(i) .eq. 0.0d0)  douind(i) = .false.
            if (polsiz(i).ne.0 .or. polarity(i).ne.0.0d0 .or.
     &             chgct(i).ne.0.0d0 .or. dmpct(i).ne.0.0d0) then
               npole = npole + 1
               ipole(npole) = i
               pollist(i) = npole
               zaxis(npole) = zaxis(i)
               xaxis(npole) = xaxis(i)
               yaxis(npole) = yaxis(i)
               polaxe(npole) = polaxe(i)
               do k = 1, maxpole
                  pole(k,npole) = pole(k,i)
               end do
               mono0(npole) = pole(1,i)
               if (palpha(i) .ne. 0.0d0)  ncp = ncp + 1
               pcore(npole) = pcore(i)
               pval(npole) = pval(i)
               pval0(npole) = pval(i)
               palpha(npole) = palpha(i)
               if (polarity(i) .ne. 0.0d0) then
                  npolar = npolar + 1
                  ipolar(npolar) = npole
                  douind(i) = .true.
               end if
               if (dirdamp(i) .ne. 0.0d0)  use_dirdamp = .true.
               polarity(npole) = polarity(i)
               thole(npole) = thole(i)
               dirdamp(npole) = dirdamp(i)
               pdamp(npole) = pdamp(i)
               if (chgct(i).ne.0.0d0 .or. dmpct(i).ne.0.0d0) then
                  nct = nct + 1
               end if
               chgct(npole) = chgct(i)
               dmpct(npole) = dmpct(i)
            end if
         end do
      end if
c
c     test multipoles at chiral sites and invert if necessary
c
      if (use_chgtrn)  call chkpole
c
c     turn off individual electrostatic potentials if not used
c
      if (npole .eq. 0)  use_mpole = .false.
      if (ncp .ne. 0)  use_chgpen = .true.
      if (ncp .ne. 0)  use_thole = .false.
      if (use_dirdamp)  use_thole = .true.
      if (npolar .eq. 0)  use_polar = .false.
      if (nct .eq. 0)  use_chgtrn = .false.
      return
      end
