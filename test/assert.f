c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module assert  --  unit-test assertions and refs            ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert" provides the lightweight assertion primitives,
c     fixture loader and reference reader used by the Tinker unit test
c     binary; it is the analog of the tinker-gpu COMPARE_REALS /
c     COMPARE_GRADIENT / TestReference facilities, with no external
c     dependency beyond the Tinker library itself
c
c     npass       running count of checks that have passed
c     nfail       running count of checks that have failed
c     ta_root     absolute path of the test directory tree
c     ta_filter   optional tag/name filter selecting which tests run
c     ta_saved    working directory saved across a pushdir/popdir pair
c     ta_detail   flag to also print per-check detail lines (-v option)
c
c
      module assert
      implicit none
      integer npass
      integer nfail
      character*240 ta_root
      character*240 ta_filter
      character*240 ta_saved
      logical ta_detail
      data npass / 0 /
      data nfail / 0 /
      data ta_detail / .false. /
      save
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ta_init  --  set test root and filter            ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ta_init" records the test directory root and the optional test
c     selection filter passed on the command line
c
c
      subroutine ta_init (root,filter,detail)
      use assert
      implicit none
      character*(*) root,filter
      logical detail
c
c
      ta_root = root
      ta_filter = filter
      ta_detail = detail
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  function skiptest  --  filter tests by name or tag          ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "skiptest" returns true when the active filter, if any, matches
c     neither the test name nor its space-separated tag string, so the
c     caller can return early
c
c
      logical function skiptest (name,tags)
      use assert
      implicit none
      character*(*) name,tags
c
c
      skiptest = .false.
      if (len_trim(ta_filter) .eq. 0)  return
      if (index(name,trim(ta_filter)) .gt. 0)  return
      if (index(tags,trim(ta_filter)) .gt. 0)  return
      skiptest = .true.
      write (*,10)  trim(name)
   10 format (1x,'skip  ',a)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine pushdir  --  enter a fixture directory           ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "pushdir" saves the current working directory and changes into a
c     subdirectory of the test root, so fixture-relative paths in a key
c     file resolve correctly
c
c
      subroutine pushdir (subdir)
      use assert
      implicit none
      character*(*) subdir
c
c
      call getcwd (ta_saved)
      call chdir (trim(ta_root)//'/'//subdir)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine popdir  --  return to prior directory            ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "popdir" restores the working directory saved by the matching
c     call to "pushdir"
c
c
      subroutine popdir ()
      use assert
      implicit none
c
c
      call chdir (trim(ta_saved))
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine loadfix  --  read a fixture and set up FF        ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "loadfix" performs the same work as "getxyz" for a caller-supplied
c     base name, reading the coordinates and key file and building the
c     force field, so a single binary can set up many structures
c
c
      subroutine loadfix (base)
      implicit none
      character*(*) base
      character*240 xyzfile
      integer ixyz,freeunit
c
c
      call initial
      xyzfile = base
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      ixyz = freeunit ()
      open (unit=ixyz,file=xyzfile,status='old')
      call readxyz (ixyz)
      close (unit=ixyz)
      call mechanic
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine refpath  --  build a reference file path         ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "refpath" returns the path to a reference file, given its case
c     subdirectory and file name, under the test-tree "ref" directory
c
c
      subroutine refpath (subdir,file,path)
      use assert
      implicit none
      character*(*) subdir,file,path
c
c
      path = trim(ta_root)//'/ref/'//trim(subdir)//'/'//trim(file)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine assert_real  --  compare scalar to reference     ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert_real" checks that a computed scalar agrees with a stored
c     reference value to within an absolute tolerance, recording the
c     outcome and printing one summary line
c
c
      subroutine assert_real (got,ref,eps,label)
      use assert
      implicit none
      real*8 got,ref,eps,diff
      character*(*) label
c
c
      diff = abs(got-ref)
      if (diff .le. eps) then
         npass = npass + 1
         write (*,10)  'pass ',trim(label)
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
      end if
      if (ta_detail .or. diff.gt.eps)  write (*,20)  got,ref,diff
      if (diff .gt. eps)  call assert_summary ()
   10 format (1x,a5,1x,a)
   20 format (7x,'calc=',g16.8,' ref=',g16.8,' d=',g11.3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine assert_int  --  compare integer to reference     ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert_int" checks that a computed integer exactly agrees with a
c     stored reference integer, recording the outcome and printing one
c     summary line
c
c
      subroutine assert_int (got,ref,label)
      use assert
      implicit none
      integer got,ref
      character*(*) label
c
c
      if (got .eq. ref) then
         npass = npass + 1
         write (*,10)  'pass ',trim(label)
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
      end if
      if (ta_detail .or. got.ne.ref)  write (*,20)  got,ref
      if (got .ne. ref)  call assert_summary ()
   10 format (1x,a5,1x,a)
   20 format (7x,'calc=',i0,' ref=',i0)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine assert_grad  --  compare gradient to reference   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert_grad" checks that every Cartesian gradient component of an
c     n-atom system agrees with a stored reference to within an absolute
c     tolerance, reporting the largest deviation found
c
c
      subroutine assert_grad (derivs,refg,n,eps,label)
      use assert
      implicit none
      integer n,i,j,iw,jw
      real*8 derivs(3,*),refg(3,*),eps,diff,worst
      character*(*) label
c
c
      worst = 0.0d0
      iw = 1
      jw = 1
      do i = 1, n
         do j = 1, 3
            diff = abs(derivs(j,i)-refg(j,i))
            if (diff .gt. worst) then
               worst = diff
               iw = i
               jw = j
            end if
         end do
      end do
      if (worst .le. eps) then
         npass = npass + 1
         write (*,10)  'pass ',trim(label)
         if (ta_detail)  write (*,20)  worst
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,30)  iw,jw,derivs(jw,iw),iw,jw,refg(jw,iw),worst
         call assert_summary ()
      end if
   10 format (1x,a5,1x,a)
   20 format (7x,'maxdiff=',g12.4)
   30 format (7x,'calc[',i0,',',i0,']=',g16.8,'  ref[',i0,',',i0,
     &        ']=',g16.8,'  d=',g11.3)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine load_ref  --  read energy/gradient reference     ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "load_ref" reads a Tinker reference output file and extracts the
c     total potential energy, intermolecular energy, internal virial and
c     per-atom analytical gradient, as written by "analyze", "testvir"
c     and "testgrad"; the gradient is taken from the integer-indexed
c     "Anlyt" rows and stored by atom index
c
c
      subroutine load_ref (fname,mxn,ref_e,ref_ei,refv,refg,nat)
      use assert
      implicit none
      integer mxn,nat,i,j,atom,ios,ip,next
      real*8 ref_e,ref_ei,refv(3,3),refg(3,mxn)
      character*(*) fname
      character*256 line
c
c
      ref_e = 0.0d0
      ref_ei = 0.0d0
      nat = 0
      do i = 1, 3
         do j = 1, 3
            refv(j,i) = 0.0d0
         end do
      end do
      do i = 1, mxn
         do j = 1, 3
            refg(j,i) = 0.0d0
         end do
      end do
      open (unit=15,file=fname,status='old',iostat=ios)
      if (ios .ne. 0) then
         nfail = nfail + 1
         write (*,*) 'load_ref: cannot open ',trim(fname)
         call assert_summary ()
      end if
  100 continue
      read (15,'(a)',iostat=ios) line
      if (ios .ne. 0)  go to 190
c
c     the total potential energy follows its label and a colon
c
      if (index(line,'Total Potential Energy') .gt. 0) then
         next = index(line,':') + 1
         call getfloat (line,ref_e,next)
         go to 100
      end if
c
c     the intermolecular energy is present for multi-molecule cases
c
      if (index(line,'Intermolecular Energy') .gt. 0) then
         next = index(line,':') + 1
         call getfloat (line,ref_ei,next)
         go to 100
      end if
c
c     the internal virial is written as one labeled row and two
c     continuation rows, each containing three tensor components
c
      if (index(line,'Internal Virial Tensor') .gt. 0) then
         next = index(line,':') + 1
         call getfloat (line,refv(1,1),next)
         call getfloat (line,refv(2,1),next)
         call getfloat (line,refv(3,1),next)
         do i = 2, 3
            read (15,'(a)',iostat=ios) line
            if (ios .ne. 0)  go to 190
            next = 1
            call getfloat (line,refv(1,i),next)
            call getfloat (line,refv(2,i),next)
            call getfloat (line,refv(3,i),next)
         end do
         go to 100
      end if
c
c     per-atom gradient rows read as "Anlyt <atom> dx dy dz norm",
c     using the Tinker getnumb/getfloat string parsers from libtinker
c
      ip = index(line,'Anlyt')
      if (ip .gt. 0) then
         next = ip + 5
         call getnumb (line,atom,next)
         if (atom.ge.1 .and. atom.le.mxn) then
            call getfloat (line,refg(1,atom),next)
            call getfloat (line,refg(2,atom),next)
            call getfloat (line,refg(3,atom),next)
            if (atom .gt. nat)  nat = atom
         end if
      end if
      go to 100
  190 continue
      close (unit=15)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine load_engcnt  --  read named energy and count     ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "load_engcnt" reads a Tinker reference output file and extracts
c     the energy value and interaction count from an "EngCnt" row whose
c     component name matches the requested string
c
c
      subroutine load_engcnt (fname,name,eng,cnt)
      use assert
      implicit none
      integer cnt,ios,ip,next
      real*8 eng,rcnt
      character*(*) fname,name
      character*256 line
c
c
      eng = 0.0d0
      cnt = 0
      open (unit=15,file=fname,status='old',iostat=ios)
      if (ios .ne. 0) then
         nfail = nfail + 1
         write (*,*) 'load_engcnt: cannot open ',trim(fname)
         call assert_summary ()
      end if
  100 continue
      read (15,'(a)',iostat=ios) line
      if (ios .ne. 0)  go to 190
      if (index(line,'EngCnt') .eq. 0)  go to 100
      ip = index(line,trim(name))
      if (ip .eq. 0)  go to 100
      next = ip + len_trim(name)
      call getfloat (line,eng,next)
      call getfloat (line,rcnt,next)
      cnt = nint(rcnt)
      close (unit=15)
      return
  190 continue
      close (unit=15)
      nfail = nfail + 1
      write (*,*) 'load_engcnt: missing ',trim(name),' in ',trim(fname)
      call assert_summary ()
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine run_prog  --  run a built Tinker program         ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "run_prog" runs a built Tinker program as a subprocess, capturing
c     its output; the executable is taken from the source tree (name.x)
c     if present, else from the bin directory, located absolutely via
c     the test root so it works regardless of the current directory
c
c
      subroutine run_prog (name,args,outfile,ist)
      use assert
      implicit none
      character*(*) name,args,outfile
      integer ist
      character*512 exe,cmd
      logical exist
c
c
c     prefer the built source executable, else the installed bin
c
      exe = trim(ta_root)//'/../source/'//trim(name)//'.x'
      inquire (file=exe,exist=exist)
      if (.not. exist) then
         exe = trim(ta_root)//'/../bin/'//trim(name)
         inquire (file=exe,exist=exist)
      end if
c
c     if the program was not found, report and fail this test; the
c     caller sees ist = -1 and should skip its output comparison
c
      if (.not. exist) then
         write (*,*)  trim(name), ' program wasn''t found'
         nfail = nfail + 1
         ist = -1
         return
      end if
      cmd = trim(exe)//' '//trim(args)//' < /dev/null > '//
     &         trim(outfile)//' 2>&1'
      call execute_command_line (cmd,exitstat=ist)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine assert_files  --  tolerant file comparison       ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert_files" checks that two Tinker output files agree: banner
c     (Tinker header) lines and blank lines are ignored, corresponding
c     data lines must have matching tokens, where a token pair passes if
c     the strings are identical or if both parse as reals within eps;
c     thus text and integers must match exactly while floats may differ
c     up to the tolerance
c
c
      subroutine assert_files (reffile,outfile,eps,label)
      use assert
      implicit none
      character*(*) reffile,outfile,label
      real*8 eps,rv,ov
      character*512 rline,oline
      character*64 rtok,otok
      integer ur,uo,ln,pr,po,ir,io,freeunit,fkind
      logical rdone,odone
c
c
      ur = freeunit ()
      open (unit=ur,file=reffile,status='old',iostat=ir)
      uo = freeunit ()
      open (unit=uo,file=outfile,status='old',iostat=io)
      if (ir.ne.0 .or. io.ne.0) then
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,20)  'file open failed'
         call assert_summary ()
      end if
      ln = 0
      fkind = 0
c
c     step through the non-banner, non-blank lines in parallel
c
  100 continue
      call next_content (ur,rline,rdone)
      call next_content (uo,oline,odone)
      if (rdone .and. odone)  go to 200
      ln = ln + 1
      if (rdone .neqv. odone) then
         fkind = 1
         go to 200
      end if
c
c     compare the two lines token by token
c
      pr = 1
      po = 1
  110 continue
      call gettext (rline,rtok,pr)
      call gettext (oline,otok,po)
      if (len_trim(rtok).eq.0 .and. len_trim(otok).eq.0)  go to 100
      if ((len_trim(rtok).eq.0) .neqv. (len_trim(otok).eq.0)) then
         fkind = 2
      else if (rtok .eq. otok) then
         go to 110
      else
         read (rtok,*,iostat=ir) rv
         read (otok,*,iostat=io) ov
         if (ir.eq.0 .and. io.eq.0 .and. abs(rv-ov).le.eps)  go to 110
         fkind = 2
      end if
      go to 200
c
c     report status first, then any detail lines indented beneath
c
  200 continue
      close (unit=ur)
      close (unit=uo)
      if (fkind .eq. 0) then
         npass = npass + 1
         write (*,10)  'pass ',trim(label)
         if (ta_detail)  write (*,20)  'file match'
      else if (fkind .eq. 1) then
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,30)  ln
         call assert_summary ()
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,40)  trim(rtok),trim(otok)
         write (*,20)  'ref line:  '//trim(rline)
         write (*,20)  'out line:  '//trim(oline)
         call assert_summary ()
      end if
   10 format (1x,a5,1x,a)
   20 format (7x,a)
   30 format (7x,'file length differs near data line ',i0)
   40 format (7x,'mismatch for :  ref="',a,'"  out="',a,'"')
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine next_content  --  next non-banner data line      ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "next_content" returns the next line of a file that is neither
c     blank nor part of a Tinker banner (first non-blank character is a
c     '#'); "done" is set true at end of file
c
c
      subroutine next_content (unit,line,done)
      implicit none
      integer unit,ios,i,n
      character*(*) line
      logical done
c
c
  100 continue
      read (unit,'(a)',iostat=ios) line
      if (ios .ne. 0) then
         done = .true.
         line = ' '
         return
      end if
      if (len_trim(line) .eq. 0)  go to 100
      n = len_trim(line)
      do i = 1, n
         if (line(i:i) .ne. ' ') then
            if (line(i:i) .eq. '#')  go to 100
            done = .false.
            return
         end if
      end do
      go to 100
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine assert_summary  --  report totals and exit       ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "assert_summary" prints the pass/fail totals and terminates with a
c     nonzero exit status if any check has failed, so the driver returns
c     failure to the caller
c
c
      subroutine assert_summary ()
      use assert
      implicit none
c
c
      write (*,10)  npass,npass+nfail
   10 format (/,1x,i0,' of ',i0,' checks passed')
      if (nfail .gt. 0)  call exit (1)
      call exit (0)
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine check_engcnt  --  compare energy and count       ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "check_engcnt" loads a named EngCnt reference row and compares
c     both the energy value and the interaction count; the full energy
c     term name is included in failure labels for easier diagnosis
c
c
      subroutine check_engcnt (rpath,name,eng,cnt,eps,label)
      implicit none
      integer cnt,refcnt
      real*8 eng,refeng,eps
      character*(*) rpath,name,label
c
c
      call load_engcnt (rpath,name,refeng,refcnt)
      call assert_real (eng,refeng,eps,
     &                  label//' '//trim(name)//' energy')
      call assert_int (cnt,refcnt,label//' '//trim(name)//' count')
      return
      end
