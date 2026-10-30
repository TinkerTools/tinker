c
c
c     ############################################################
c     ##  COPYRIGHT (C) 2018 by Joshua Rackers & Jay W. Ponder  ##
c     ##                   All Rights Reserved                  ##
c     ############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine kdisp  --  dispersion parameter assignment  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "kdisp" assigns C6 coefficients and damping parameters for
c     dispersion interactions and processes any new or changed
c     values for these parameters
c
c
      subroutine kdisp
      use atomid
      use atoms
      use disp
      use dsppot
      use inform
      use iounit
      use kdsp
      use keys
      use limits
      use potent
      use sizes
      implicit none
      integer i,k,ii
      integer ia,ic,next
      real*8 cs,adsp
      real*8 csum
      logical header
      logical, allocatable :: clsset(:)
      logical, allocatable :: dspset(:)
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing damped dispersion parameters,
c     noting the classes they set
c
      allocate (clsset(maxclass))
      do k = 1, maxclass
         clsset(k) = .false.
      end do
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'DISPERSION ') then
            k = 0
            cs = 0.0d0
            adsp = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  cs,adsp
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Damped Dispersion',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',16x,'C6',12x,'Damp',/)
               end if
               if (k .le. maxclass) then
                  dspsix(k) = cs
                  dspdmp(k) = adsp
                  clsset(k) = .true.
                  if (.not. silent) then
                     write (iout,30)  k,cs,adsp
   30                format (6x,i6,7x,f15.4,f15.4)
                  end if
               else
                  write (iout,40)
   40             format (/,' KDISP  --  Too many Damped',
     &                       ' Dispersion Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(idisp))  deallocate (idisp)
      if (allocated(csix))  deallocate (csix)
      if (allocated(adisp))  deallocate (adisp)
      allocate (idisp(n))
      allocate (csix(n))
      allocate (adisp(n))
c
c     assign the dispersion C6 values and alpha parameters,
c     noting the sites set by a keyword
c
      allocate (dspset(n))
      do i = 1, n
         csix(i) = 0.0d0
         adisp(i) = 0.0d0
         dspset(i) = .false.
         ic = class(i)
         if (ic .ne. 0) then
            csix(i) = dspsix(ic)
            adisp(i) = dspdmp(ic)
            dspset(i) = clsset(ic)
         end if
      end do
c
c     process keywords containing atom specific dispersion parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11) .eq. 'DISPERSION ') then
            ia = 0
            cs = 0.0d0
            adsp = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,cs,adsp
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Dispersion Values for',
     &                       ' Specific Atoms :',
     &                    //,8x,'Atom',19x,'C6',12x,'Damp',/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,cs,adsp
   60             format (6x,i6,7x,f15.4,f15.4)
               end if
               csix(ia) = cs
               adisp(ia) = adsp
               dspset(ia) = .true.
            end if
   70       continue
         end if
      end do
c
c     map the buffered 14-7 vdw tail onto unset dispersion sites
c
      if (use_dhal)  call kdsphal (dspset)
      deallocate (clsset)
      deallocate (dspset)
c 
c     remove zero and undefined dispersion sites from the list
c     
      ndisp = 0
      do i = 1, n
         if (csix(i) .ne. 0.0d0) then 
            ndisp = ndisp + 1
            idisp(ndisp) = i
         end if
      end do
c
c     compute pairwise sum of C6 coefficients needed for PME, which
c     is the square of their sum
c
      csixpr = 0.0d0
      if (use_dewald) then
         csum = 0.0d0
         do ii = 1, ndisp
            i = idisp(ii)
            csum = csum + csix(i)
         end do
         csixpr = csum * csum
      end if
c
c     turn off the dispersion potential if not used
c
      if (ndisp .eq. 0) then
         use_disp = .false.
         use_dhal = .false.
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine kdsphal  --  dispersion C6 from buffered 14-7  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "kdsphal" sets the vdw switching window that hands the buffered
c     14-7 tail to an Ewald dispersion sum, and gives each site without
c     a dispersion value the square root of the C6 whose tail beyond
c     the switch matches that of its own hal pair, using the same
c     numerical integration as the long range correction; sites given
c     a dispersion value, including zero, keep their C6 and damping,
c     but damping must be set for either all or none of the dispersion
c     sites
c
c
      subroutine kdsphal (dspset)
      use atomid
      use atoms
      use disp
      use dsppot
      use inform
      use iounit
      use limits
      use vdw
      use vdwpot
      implicit none
      integer i,j,it,nt
      integer nstep,ndelta
      integer ndamp,nfree
      real*8 r,r2,r4
      real*8 w,dw,e
      real*8 range,rdelta
      real*8 offset
      real*8 ehal,elj
      real*8 rv,rv7,eps
      real*8 ehalpr
      real*8, allocatable :: csixt(:)
      logical header
      logical dspset(*)
      logical, allocatable :: done(:)
c
c
c     set the vdw switching window
c
      dhoff = vdwcut
      dhcut = min(vdwtaper,vdwcut)
c
c     set number of steps and range for numerical integration
c
      nstep = 2
      range = 100.0d0
      ndelta = int(dble(nstep)*(range-dhcut))
      rdelta = (range-dhcut) / dble(ndelta)
      offset = dhcut - 0.5d0*rdelta
c
c     perform dynamic allocation of some local arrays
c
      nt = size(radmin,1)
      allocate (csixt(nt))
      allocate (done(nt))
      do it = 1, nt
         csixt(it) = -1.0d0
         done(it) = .false.
      end do
c
c     find the C6 of each needed vdw type from its own hal pair
c
      header = .true.
      do i = 1, n
         it = jvdw(i)
         if (csix(i).ne.0.0d0 .or. dspset(i) .or. it.le.0)  goto 10
         if (csixt(it) .lt. 0.0d0) then
            rv = radmin(it,it)
            rv7 = rv**7
            eps = epsilon(it,it)
            ehal = 0.0d0
            elj = 0.0d0
            do j = 1, ndelta
               r = offset + dble(j)*rdelta
               r2 = r * r
               r4 = r2 * r2
               e = ehalpr (r,rv,rv7,eps)
               call dhalwt (r,w,dw)
               ehal = ehal + w*e*r2
               elj = elj - w/r4
            end do
c
c     a tail that is not attractive has no C6 equivalent
c
            csixt(it) = 0.0d0
            if (ehal.lt.0.0d0 .and. elj.ne.0.0d0) then
               csixt(it) = sqrt(ehal/elj)
            else if (eps .gt. 0.0d0) then
               write (iout,20)  class(i)
   20          format (/,' KDSPHAL  --  Buffered 14-7 Tail of Atom',
     &                    ' Class',i6,' has no Dispersion Equivalent')
               call fatal
            end if
         end if
         csix(i) = csixt(it)
         if (verbose .and. .not.done(it)) then
            done(it) = .true.
            if (header) then
               header = .false.
               write (iout,30)  dhoff,dhcut
   30          format (/,' Dispersion Values Mapped from Buffered',
     &                    ' 14-7 vdW :',
     &                 //,' vdW Cutoff',f10.4,' Ang and Taper',f10.4,
     &                    ' Ang',
     &                 //,5x,'Atom Class',11x,'Sqrt C6',/)
            end if
            write (iout,40)  class(i),csix(i)
   40       format (6x,i6,7x,f15.4)
         end if
   10    continue
      end do
c
c     damping must be set for all of the dispersion sites or none
c
      ndamp = 0
      nfree = 0
      do i = 1, n
         if (csix(i) .ne. 0.0d0) then
            if (adisp(i) .gt. 0.0d0) then
               ndamp = ndamp + 1
            else
               nfree = nfree + 1
            end if
         end if
      end do
      if (ndamp.ne.0 .and. nfree.ne.0) then
         write (iout,50)
   50    format (/,' KDSPHAL  --  Dispersion Damping must be Set',
     &              ' for All or No Sites with a Buffered 14-7 Tail')
         call fatal
      end if
      dhdamp = (ndamp .ne. 0)
c
c     perform deallocation of some local arrays
c
      deallocate (csixt)
      deallocate (done)
      return
      end
