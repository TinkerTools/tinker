c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1994  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program noesetup  --  put QUANTA NOE's into a keyfile  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "noesetup" takes PDB and Tinker coordinate files and a QUANTA
c     NOE file, and generates a list of distance restraints suitable
c     for input to the Tinker distgeom program; modified to handle
c     the special characters "+" and "*" in PDB atom names
c
c
      program noesetup
      use atmtyp
      use atoms
      use files
      use iounit
      use pdb
      implicit none
      integer i,j
      integer inoe,ikey
      integer ires,next
      integer bound1,bound2
      integer nlist1,list1(20)
      integer nlist2,list2(20)
      integer res1,res2,freeunit
      integer mult1,mult2
      integer ni(maxres)
      integer cai(maxres)
      integer ci(maxres)
      real*8 xx,yy,zz
      real*8 distsq,dummy
      real*8 factor,upper,lower
      character*1 letter
      character*3 resnam1,resnam2
      character*3 amino(maxres)
      character*4 atom1,atom2
      character*240 noefile,keyfile
      character*240 record,string
c
c
c     read the PDB and Tinker files for the structure
c
      call initial
      call getpdb
      call getxyz
c
c     store the residue name of each amino acid in sequence
c
      do i = 1, npdb
         amino(resnum(i)) = resnam(i)
      end do
c
c     find the backbone atoms in the coordinate file
c
      ires = 0
      do i = 1, n-3
         if (name(i)(1:1).eq.'N' .and. name(i+1)(1:1).eq.'C' .and.
     &       name(i+2)(1:1).eq.'C' .and. name(i+3)(1:1).eq.'O') then
            ires = ires + 1
            ni(ires) = i
            cai(ires) = i + 1
            ci(ires) = i + 2
         end if
      end do
c
c     open the Tinker keyfile to which output will be written
c
      ikey = freeunit ()
      keyfile = filename(1:leng)
      call suffix (keyfile,'key')
      call version (keyfile,'new')
      open (unit=ikey,file=keyfile,status='new')
c
c     write the backbone torsional constraints in keyfile format
c
      write (ikey,10)
   10 format ('#',/,'# Backbone Torsional Constraints',/,'#',/)
      do i = 1, ires
         if (i .ne. 1) then
            write (ikey,20)  ci(i-1),ni(i),cai(i),ci(i),
     &                       -80.0,-40.0,1.0,i
   20       format (/,'restrain-dihedral',2x,4i6,1x,3f7.1,
     &                 2x,'[',i4,' Phi  ]')
         end if
         if (i .ne. ires) then
            write (ikey,30)  ni(i),cai(i),ci(i),ni(i+1),
     &                       -60.0,-20.0,1.0,i
   30       format ('restrain-dihedral',2x,4i6,1x,3f7.1,
     &                 2x,'[',i4,' Psi  ]')
            write (ikey,40)  cai(i),ci(i),ni(i+1),cai(i+1),
     &                       175.0,185.0,1.0,i
   40       format ('restrain-dihedral',2x,4i6,1x,3f7.1,
     &                 2x,'[',i4,' Omega]')
         end if
      end do
c
c     open the file containing the QUANTA-formatted NOE data
c
      inoe = freeunit ()
      noefile = filename(1:leng)
      call suffix (noefile,'noe')
      call version (noefile,'old')
      open (unit=inoe,file=noefile,status='old')
c
c     write a header for the NOE distance constraints
c
      write (ikey,50)
   50 format (/,'#',/,'# NOE InterProton Distance Constraints',/,'#',/)
c
c     get the first atom involved in the NOE constraint
c
      dowhile (.true.)
         read (inoe,60,err=190,end=190)  record
   60    format (a240)
         next = index(record,':') + 1
         if (next .eq. 0)  goto 180
         call getnumb (record,res1,next)
         resnam1 = amino(res1)
         atom1 = '    '
         do i = 1, 5
            letter = record(next:next)
            if (letter .eq. '''') then
               goto 70
            else if (letter .eq. '_') then
               atom1(i:i) = ' '
            else if (i .eq. 5) then
               atom1(1:1) = letter
            else
               atom1(i:i) = letter
            end if
            next = next + 1
         end do
   70    continue
         do i = 4, 1, -1
            if (atom1(i:i) .ne. ' ') then
               if ((atom1(i:i).eq.'+' .or. atom1(i:i).eq.'*')
     &                  .and. atom1(1:1).eq.' ') then
                  if ((resnam1.ne.'PHE' .and. resnam1.ne.'TYR')
     &                     .or. atom1(3:3).eq.'B') then
                     atom1(1:1) = atom1(i:i)
                     atom1(i:i) = ' '
                  end if
               end if
               goto 80
            end if
         end do
   80    continue
c
c     now, the second atom involved in the NOE constraint
c
         record = record(next:80)
         next = index(record,':') + 1
         if (next .eq. 0)  goto 180
         call getnumb (record,res2,next)
         resnam2 = amino(res2)
         atom2 = '    '
         do i = 1, 5
            letter = record(next:next)
            if (letter .eq. '''') then
               goto 90
            else if (letter .eq. '_') then
               atom2(i:i) = ' '
            else if (i .eq. 5) then
               atom2(1:1) = letter
            else
               atom2(i:i) = letter
            end if
            next = next + 1
         end do
   90    continue
         do i = 4, 1, -1
            if (atom2(i:i) .ne. ' ') then
               if ((atom2(i:i).eq.'+' .or. atom2(i:i).eq.'*')
     &                  .and. atom2(1:1).eq.' ') then
                  if ((resnam2.ne.'PHE' .and. resnam2.ne.'TYR')
     &                     .or. atom2(3:3).eq.'B') then
                     atom2(1:1) = atom2(i:i)
                     atom2(i:i) = ' '
                  end if
               end if
               goto 100
            end if
         end do
  100    continue
c
c     locate the values of the NOE distance bounds
c
         next = next + 8
         string = record(next:240)
         read (string,*)  dummy,upper,lower
c        upper = 10.0d0
c        lower = 2.0d0
c
c     clean up QUANTA's sloppy hydrogen naming convention (!&*$#!)
c
         call fixpdb (amino(res1),atom1)
         call fixpdb (amino(res2),atom2)
c
c     find all matches of the first atom type in the PDB file
c
         nlist1 = 0
         do i = 1, npdb
            if (resnum(i) .eq. res1) then
               do j = 1, 4
                  if (atom1(j:j).ne.atmnam(i)(j:j) .and.
     &                atom1(j:j).ne.'+' .and. atom1(j:j).ne.'*')
     &               goto 110
               end do
               nlist1 = nlist1 + 1
               list1(nlist1) = i
            end if
  110       continue
         end do
c
c     search for the first atom types in the Tinker file
c
         do j = 1, nlist1
            xx = xpdb(list1(j))
            yy = ypdb(list1(j))
            zz = zpdb(list1(j))
            list1(j) = 0
            do i = 1, n
               distsq = (x(i)-xx)**2 + (y(i)-yy)**2 + (z(i)-zz)**2
               if (distsq .lt. 0.01d0) then
                  list1(j) = i
                  goto 120
               end if
            end do
  120       continue
         end do
c
c     find all matches of the second atom type in the PDB file
c
         nlist2 = 0
         do i = 1, npdb
            if (resnum(i) .eq. res2) then
               do j = 1, 4
                  if (atom2(j:j).ne.atmnam(i)(j:j) .and.
     &                atom2(j:j).ne.'+' .and. atom2(j:j).ne.'*')
     &               goto 130
               end do
               nlist2 = nlist2 + 1
               list2(nlist2) = i
            end if
  130       continue
         end do
c
c     search for the second atom types in the Tinker file
c
         do j = 1, nlist2
            xx = xpdb(list2(j))
            yy = ypdb(list2(j))
            zz = zpdb(list2(j))
            list2(j) = 0
            do i = 1, n
               distsq = (x(i)-xx)**2 + (y(i)-yy)**2 + (z(i)-zz)**2
               if (distsq .lt. 0.01d0) then
                  list2(j) = i
                  goto 140
               end if
            end do
  140       continue
         end do
c
c     check for constraints that do not have matching atoms
c
         if (nlist1.eq.0 .or. nlist2.eq.0) then
            write (iout,150)
  150       format (/,' WARNING  --  NOESETUP did not Find any',
     &                 ' Matches for the Next Constraint',/)
         end if
c
c     write the bound values in a keyfile format
c
         mult1 = nlist1
         mult2 = nlist2
         do i = 1, 4
            if (atom1(i:i) .eq. '*')  mult1 = 1
            if (atom2(i:i) .eq. '*')  mult2 = 1
         end do
         factor = 100.0d0 / dble(mult1*mult2)
         write (iout,160)  res1,atom1,res2,atom2,lower,upper,
     &                     mult1*mult2
  160    format (3x,i6,3x,a4,8x,i6,3x,a4,8x,2f8.2,i8)
         do i = 1, nlist1
            bound1 = list1(i)
            do j = 1, nlist2
               bound2 = list2(j)
               write (ikey,170)  bound1,bound2,lower,upper,
     &                           factor,res1,atom1,res2,atom2
  170          format ('restrain-distance',2x,2i6,1x,2f7.2,2x,f6.2,
     &                    2x,'['i4,1x,a4,1x,i4,1x,a4,']')
            end do
         end do
  180    continue
      end do
c
c     perform any final tasks before program exit
c
  190 continue
      call final
      end
