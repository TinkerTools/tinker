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
      use sizes
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
      implicit none
      integer i,k
      integer ia,ic,next
      real*8 cs,adsp
      real*8 csixi
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing damped dispersion parameters
c
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
c     assign the dispersion C6 values and alpha parameters 
c
      do i = 1, n
         ic = class(i)
         csix(i) = dspsix(ic)
         adisp(i) = dspdmp(ic)
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
            end if
   70       continue
         end if
      end do
c 
c     remove zero and undefined dispersion sites from the list
c     
      ndisp = 0
      do i = 1, n
         if (csix(i) .ne. 0.0d0) then 
            ndisp = ndisp + 1
            idisp(ndisp) = i
            csix(ndisp) = csix(i)
            adisp(ndisp) = adisp(i)
         end if
      end do
c
c     compute pairwise sum of C6 coefficients needed for PME
c
      csixpr = 0.0d0
      if (use_dewald) then
         do i = 1, ndisp
            csixi = csix(i)
            do k = 1, ndisp
               csixpr = csixpr + csixi*csix(k)
            end do
         end do
      end if
c
c     turn off the dispersion potential if not used
c
      if (ndisp .eq. 0)  use_disp = .false.
      return
      end
