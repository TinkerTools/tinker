c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine katom  --  atom type parameter assignment  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "katom" assigns an atom type definitions to each atom in
c     the structure and processes any new or changed values
c
c     literature reference:
c
c     K. A. Feenstra, B. Hess and H. J. C. Berendsen, "Improving
c     Efficiency of Large Time-Scale Molecular Dynamics Simulations
c     of Hydrogen-Rich Systems", Journal of Computational Chemistry,
c     8, 786-798 (1999)
c
c     C. W. Hopkins, S. Le Grand, R. C. Walker and A. E. Roitberg,
c     "Long-Time-Step Molecular Dynamics through Hydrogen Mass
c     Repartitioning", Journal of Chemical Theory and Computation,
c     11, 1864-1874 (2015)
c
c
      subroutine katom
      use atomid
      use atoms
      use couple
      use inform
      use iounit
      use katoms
      use keys
      implicit none
      integer i,j,k
      integer next,nh
      integer cls,atn,lig
      real*8 wght,sum
      real*8 hmax,hmass
      real*8 dmin,dmass
      logical header,heavy
      character*3 symb
      character*20 keyword
      character*24 notice
      character*240 record
      character*240 string
c
c
c     process keywords containing atom type parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            cls = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            if (cls .eq. 0)  cls = k
            atmcls(k) = cls
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:240)
            read (string,*,err=40,end=40)  atn,wght,lig
            if (k.ge.1 .and. k.le.maxtyp) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,10)
   10             format (/,' Additional Atom Definition Parameters :',
     &                    //,5x,'Type  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               symbol(k) = symb
               describe(k) = notice
               atmnum(k) = atn
               weight(k) = wght
               ligand(k) = lig
               if (.not. silent) then
                  write (iout,20)  k,cls,symb,notice,atn,wght,lig
   20             format (1x,i8,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            else if (k .ge. maxtyp) then
               write (iout,30)
   30          format (/,' KATOM   --  Too many Atom Types;',
     &                    ' Increase MAXTYP')
               abort = .true.
            end if
   40       continue
         end if
      end do
c
c     transfer atom type values to individual atoms
c
      do i = 1, n
         k = type(i)
         if (k .eq. 0) then
            class(i) = 0
            atomic(i) = 0
            mass(i) = 0.0d0
            valence(i) = 0
            story(i) = 'Undefined Atom Type     '
         else
            if (symbol(k) .ne. '   ')  name(i) = symbol(k)
            class(i) = atmcls(k)
            atomic(i) = atmnum(k)
            mass(i) = weight(k)
            valence(i) = ligand(k)
            story(i) = describe(k)
         end if
      end do
c
c     repartition hydrogen masses to use "heavy" hydrogens
c
      heavy = .false.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:240)
         if (keyword(1:15) .eq. 'HEAVY-HYDROGEN ') then
            heavy = .true.
            hmax = 4.0d0
            read (string,*,err=50,end=50)  hmax
         end if
   50    continue
      end do
      if (heavy) then
         do i = 1, n
            nh = 0
            sum = mass(i)
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  nh = nh + 1
                  sum = sum + mass(k)
               end if
            end do
            hmass = min(hmax,sum/dble(nh+1))
            do j = 1, n12(i)
               k = i12(j,i)
               if (atomic(k) .eq. 1) then
                  dmass = hmass - mass(k)
                  mass(k) = hmass
                  mass(i) = mass(i) - dmass
               end if
            end do
         end do
         do i = 1, n
            if (mass(i) .lt. hmax) then
               dmass = hmax - mass(i)
               dmin = hmax + dmass
               do j = 1, n12(i)
                  k = i12(j,i)
                  if (mass(k) .gt. dmin) then
                     mass(k) = mass(k) - dmass
                     mass(i) = hmax
                     goto 60
                  end if
               end do
               do j = 1, n13(i)
                  k = i13(j,i)
                  if (mass(k) .gt. dmin) then
                     mass(k) = mass(k) - dmass
                     mass(i) = hmax
                     goto 60
                  end if
               end do
   60          continue
            end if
         end do
      end if
c
c     process keywords containing atom types for specific atoms
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'ATOM ') then
            k = 0
            symb = ' '
            notice = ' '
            atn = 0
            wght = 0.0d0
            lig = 0
            call getnumb (record,k,next)
            call getnumb (record,cls,next)
            call gettext (record,symb,next)
            call getstring (record,notice,next)
            string = record(next:240)
            read (string,*,err=90,end=90)  atn,wght,lig
            if (k.lt.0 .and. k.ge.-n) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,70)
   70             format (/,' Additional Atom Definitions for',
     &                       ' Specific Atoms :',
     &                    //,5x,'Atom  Class  Symbol  Description',
     &                       15x,'Atomic',4x,'Mass',3x,'Valence',/)
               end if
               k = -k
               if (cls .eq. 0)  cls = k
               class(k) = cls
               name(k) = symb
               story(k) = notice
               atomic(k) = atn
               mass(k) = wght
               valence(k) = lig
               if (.not. silent) then
                  write (iout,80)  k,cls,symb,notice,atn,wght,lig
   80             format (1x,i8,i6,5x,a3,3x,a24,i6,f11.3,i6)
               end if
            end if
   90       continue
         end if
      end do
c
c     check for presence of undefined atom types or classes
c
      header = .true.
      do i = 1, n
         k = type(i)
         cls = class(i)
         if (k.lt.1 .or. k.gt.maxtyp
     &          .or. cls.lt.1 .or. cls.gt.maxclass) then
            abort = .true.
            if (header) then
               header = .false.
               write (iout,100)
  100          format (/,' Undefined Atom Types or Classes :',
     &                 //,' Type',10x,'Atom Number',5x,'Atom Type',
     &                    5x,'Atom Class',/)
            end if
            write (iout,110)  i,k,cls
  110       format (' Atom',9x,i8,10x,i5,10x,i5)
         end if
      end do
c
c     check the number of atoms attached to each atom
c
      header = .true.
      do i = 1, n
         if (n12(i) .ne. valence(i)) then
            if (header) then
               header = .false.
               write (iout,120)
  120          format (/,' Atoms with an Unusual Number of Attached',
     &                    ' Atoms :',
     &                 //,' Type',11x,'Atom Name',6x,'Atom Type',7x,
     &                    'Expected',4x,'Found',/)
            end if
            write (iout,130)  i,name(i),type(i),valence(i),n12(i)
  130       format (' Valence',4x,i8,'-',a3,8x,i5,10x,i5,5x,i5)
         end if
      end do
      return
      end
