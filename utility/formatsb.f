c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2009  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program formatsb  --  convert stretch-bend to new format  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "formatsb" converts stretch-bend parameters from the previous
c     Tinker format similar to MM2/MM3 to the more general format
c     used in Tinker 5 and later, derived from the MMFF formulation
c
c     the old format has a single parameter for each atom class of
c     the central atom in a bond angle, and three possible force
c     constants based on the number of hydrogens among the two
c     terminal atom classes of the angle
c
c     the new format has a parameters for each triple of atom classes
c     corresponding to a bond angle, IA-IB-IC, and two separate force
c     constants used to separately couple the IA-IB and IB-IC bond
c     length deviations to the bond angle deviation
c
c     note executables must be built against a Tinker 4 or earlier
c     source tree that implements the old stretch-bend format
c
c
      program formatsb
      use fields
      use iounit
      use kangs
      use katoms
      use kstbnd
      use sizes
      implicit none
      integer i,j,size
      integer ia,ib,ic
      integer ja,jc
      integer na,na5,na4
      integer na3,naf,nh
      integer number
      integer clsnum(maxclass)
      real*8 force1,force2
      logical mm2stbn
      character*12 blank
      external number
c
c
c     read in the force field parameters to convert
c
      call initial
      call field
c
c     set some constants and parameters for the conversion
c
      size = 4
      blank = '            '
      mm2stbn = .false.
      if (forcefield(1:3) .eq. 'MM2')  mm2stbn = .true.
c
c     write out the name of the forcefield being converted
c
      write (iout,10)  forcefield
   10 format (/,' Converting Stretch-Bend Parameters for :  ',a,/)
c
c     determine the total number of angle bending parameters
c
      na = maxna
      na5 = maxna5
      na4 = maxna4
      na3 = maxna3
      naf = maxnaf
      do i = maxna, 1, -1
         if (ka(i) .eq. blank)  na = i - 1
      end do
      do i = maxna5, 1, -1
         if (ka5(i) .eq. blank)  na5 = i - 1
      end do
      do i = maxna4, 1, -1
         if (ka4(i) .eq. blank)  na4 = i - 1
      end do
      do i = maxna3, 1, -1
         if (ka3(i) .eq. blank)  na3 = i - 1
      end do
      do i = maxnaf, 1, -1
         if (kaf(i) .eq. blank)  naf = i - 1
      end do
c
c     get the atomic number for each atom class
c
      do i = 1, maxclass
         clsnum(i) = 0
      end do
      do i = 1, maxtyp
         j = atmcls(i)
         if (j .ne. 0)  clsnum(j) = atmnum(i)
      end do
c
c     convert stretch-bend parameters for regular angles
c
      do i = 1, na
         ia = number (ka(i)(1:4))
         ib = number (ka(i)(5:8))
         ic = number (ka(i)(9:12))
         ja = clsnum(ia)
         jc = clsnum(ic)
         nh = 1
         if (ja .le. 1)  nh = nh + 1
         if (jc .le. 1)  nh = nh + 1
         if (stbn(nh,ib) .ne. 0.0d0) then
            force1 = stbn(nh,ib)
            force2 = stbn(nh,ib)
            if (mm2stbn) then
               if (ja .le. 1)  force1 = 0.0d0
               if (jc .le. 1)  force2 = 0.0d0
            end if
            if (force1.ne.0.0d0 .or. force2.ne.0.0d0) then
               write (iout,20)  ia,ib,ic,force1,force2
   20          format ('strbnd',4x,3i5,2x,2f10.3)
            end if
         end if
      end do
c
c     convert stretch-bend parameters for 5-ring angles
c
      do i = 1, na5
         ia = number (ka5(i)(1:4))
         ib = number (ka5(i)(5:8))
         ic = number (ka5(i)(9:12))
         ja = clsnum(ia)
         jc = clsnum(ic)
         nh = 1
         if (ja .le. 1)  nh = nh + 1
         if (jc .le. 1)  nh = nh + 1
         if (stbn(nh,ib) .ne. 0.0d0) then
            force1 = stbn(nh,ib)
            force2 = stbn(nh,ib)
            if (mm2stbn) then
               if (ja .le. 1)  force1 = 0.0d0
               if (jc .le. 1)  force2 = 0.0d0
            end if
            if (force1.ne.0.0d0 .or. force2.ne.0.0d0) then
               write (iout,30)  ia,ib,ic,force1,force2
   30          format ('strbnd5',3x,3i5,2x,2f10.3)
            end if
         end if
      end do
c
c     convert stretch-bend parameters for 4-ring angles
c
      do i = 1, na4
         ia = number (ka4(i)(1:4))
         ib = number (ka4(i)(5:8))
         ic = number (ka4(i)(9:12))
         ja = clsnum(ia)
         jc = clsnum(ic)
         nh = 1
         if (ja .le. 1)  nh = nh + 1
         if (jc .le. 1)  nh = nh + 1
         if (stbn(nh,ib) .ne. 0.0d0) then
            force1 = stbn(nh,ib)
            force2 = stbn(nh,ib)
            if (mm2stbn) then
               if (ja .le. 1)  force1 = 0.0d0
               if (jc .le. 1)  force2 = 0.0d0
            end if
            if (force1.ne.0.0d0 .or. force2.ne.0.0d0) then
               write (iout,40)  ia,ib,ic,force1,force2
   40          format ('strbnd4',3x,3i5,2x,2f10.3)
            end if
         end if
      end do
c
c     convert stretch-bend parameters for 3-ring angles
c
      do i = 1, na3
         ia = number (ka3(i)(1:4))
         ib = number (ka3(i)(5:8))
         ic = number (ka3(i)(9:12))
         ja = clsnum(ia)
         jc = clsnum(ic)
         nh = 1
         if (ja .le. 1)  nh = nh + 1
         if (jc .le. 1)  nh = nh + 1
         if (stbn(nh,ib) .ne. 0.0d0) then
            force1 = stbn(nh,ib)
            force2 = stbn(nh,ib)
            if (mm2stbn) then
               if (ja .le. 1)  force1 = 0.0d0
               if (jc .le. 1)  force2 = 0.0d0
            end if
            if (force1.ne.0.0d0 .or. force2.ne.0.0d0) then
               write (iout,50)  ia,ib,ic,force1,force2
   50          format ('strbnd3',3x,3i5,2x,2f10.3)
            end if
         end if
      end do
c
c     convert stretch-bend parameters for Fourier angles
c
      do i = 1, naf
         ia = number (kaf(i)(1:4))
         ib = number (kaf(i)(5:8))
         ic = number (kaf(i)(9:12))
         ja = clsnum(ia)
         jc = clsnum(ic)
         nh = 1
         if (ja .le. 1)  nh = nh + 1
         if (jc .le. 1)  nh = nh + 1
         if (stbn(nh,ib) .ne. 0.0d0) then
            force1 = stbn(nh,ib)
            force2 = stbn(nh,ib)
            if (mm2stbn) then
               if (ja .le. 1)  force1 = 0.0d0
               if (jc .le. 1)  force2 = 0.0d0
            end if
            if (force1.ne.0.0d0 .or. force2.ne.0.0d0) then
               write (iout,60)  ia,ib,ic,force1,force2
   60          format ('strbndf',3x,3i5,2x,2f10.3)
            end if
         end if
      end do
      end
