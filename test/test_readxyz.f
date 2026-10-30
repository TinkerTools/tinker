c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_readxyz  --  XYZ file header lines  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_readxyz" checks how "readxyz" reads the lines between the
c     atom count and the first atom: the comment line of a simple XYZ
c     file is kept as the title and never read as an atom or a box,
c     simple files without a comment line still load, frame by frame
c     in an archive, and Tinker XYZ files with or without a periodic
c     box line are unchanged
c
c
      subroutine test_readxyz
      use atomid
      use atoms
      use bound
      use boxes
      implicit none
      integer iunit
      logical skiptest
      character*(*) tname
      parameter (tname='test_readxyz')
c
c
      if (skiptest(tname,'readxyz'))  return
      call pushdir ('file/readxyz')
c
c     simple XYZ files with a comment line
c
      call test_readxyz_open ('comtext.xyz',iunit)
      call test_readxyz_water (tname//' comtext',iunit,
     &                         'water from somewhere',0.0d0)
      close (unit=iunit)
      call test_readxyz_open ('comnum.xyz',iunit)
      call test_readxyz_water (tname//' comnum',iunit,
     &                         '300 K snapshot',0.0d0)
      call assert_logical (use_bounds,.false.,tname//' comnum no box')
      close (unit=iunit)
      call test_readxyz_open ('comatom.xyz',iunit)
      call test_readxyz_water (tname//' comatom',iunit,
     &                         'frame 1.0 2.0 3.0',0.0d0)
      close (unit=iunit)
      call test_readxyz_open ('comblank.xyz',iunit)
      call test_readxyz_water (tname//' comblank',iunit,'',0.0d0)
      close (unit=iunit)
c
c     simple XYZ files without a comment line
c
      call test_readxyz_open ('nocom.xyz',iunit)
      call test_readxyz_water (tname//' nocom',iunit,'',0.0d0)
      close (unit=iunit)
      call test_readxyz_open ('nocom1.xyz',iunit)
      call readxyz (iunit)
      call assert_int (n,1,tname//' nocom1 atoms')
      call assert_logical (name(1).eq.'Na',.true.,tname//' nocom1 name')
      call assert_real (x(1),1.0d0,1.0d-12,tname//' nocom1 x')
      call assert_real (y(1),2.0d0,1.0d-12,tname//' nocom1 y')
      call assert_real (z(1),3.0d0,1.0d-12,tname//' nocom1 z')
      close (unit=iunit)
c
c     two-frame simple XYZ archives, with and without comments
c
      call test_readxyz_open ('nocomarc.xyz',iunit)
      call test_readxyz_water (tname//' nocomarc 1',iunit,'',0.0d0)
      call test_readxyz_water (tname//' nocomarc 2',iunit,'',1.0d0)
      close (unit=iunit)
      call test_readxyz_open ('nocom1arc.xyz',iunit)
      call readxyz (iunit)
      call assert_int (n,1,tname//' nocom1arc 1 atoms')
      call assert_real (z(1),3.0d0,1.0d-12,tname//' nocom1arc 1 z')
      call readxyz (iunit)
      call assert_int (n,1,tname//' nocom1arc 2 atoms')
      call assert_real (z(1),4.0d0,1.0d-12,tname//' nocom1arc 2 z')
      close (unit=iunit)
      call test_readxyz_open ('comarc.xyz',iunit)
      call test_readxyz_water (tname//' comarc 1',iunit,
     &                         'frame 1.0 2.0 3.0',0.0d0)
      call test_readxyz_water (tname//' comarc 2',iunit,
     &                         'frame 2.0 3.0 4.0',1.0d0)
      close (unit=iunit)
c
c     Tinker XYZ files with a full box line, none, and a single value
c
      call test_readxyz_open ('tinker.xyz',iunit)
      call test_readxyz_water (tname//' tinker',iunit,'Water Box',0.0d0)
      call test_readxyz_tinker (tname//' tinker')
      call assert_logical (use_bounds,.true.,tname//' tinker box')
      call assert_real (zbox,20.0d0,1.0d-12,tname//' tinker zbox')
      close (unit=iunit)
      call test_readxyz_open ('tinkernobox.xyz',iunit)
      call test_readxyz_water (tname//' tinkernobox',iunit,'',0.0d0)
      call test_readxyz_tinker (tname//' tinkernobox')
      call assert_logical (use_bounds,.false.,tname//' tinkernobox box')
      close (unit=iunit)
      call test_readxyz_open ('tinkerbox1.xyz',iunit)
      call test_readxyz_water (tname//' tinkerbox1',iunit,'Water',0.0d0)
      call test_readxyz_tinker (tname//' tinkerbox1')
      call assert_logical (use_bounds,.true.,tname//' tinkerbox1 box')
      call assert_real (ybox,20.0d0,1.0d-12,tname//' tinkerbox1 ybox')
      call assert_real (zbox,20.0d0,1.0d-12,tname//' tinkerbox1 zbox')
      close (unit=iunit)
      call popdir
      call final
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_readxyz_open  --  open a test file  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine test_readxyz_open (fname,iunit)
      implicit none
      integer iunit,freeunit
      character*(*) fname
c
c
      call initial
      iunit = freeunit ()
      open (unit=iunit,file=fname,status='old')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_readxyz_water  --  read and check water  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_readxyz_water" reads the next frame from the open unit
c     and checks the title and a water molecule shifted by dz along z
c
c
      subroutine test_readxyz_water (label,iunit,texp,dz)
      use atomid
      use atoms
      use titles
      implicit none
      integer iunit
      real*8 dz,eps
      logical match
      character*(*) label,texp
c
c
      call readxyz (iunit)
      eps = 1.0d-12
      call assert_int (n,3,label//' atoms')
      match = (ltitle .eq. len(texp))
      if (match .and. ltitle.gt.0)  match = (title(1:ltitle) .eq. texp)
      call assert_logical (match,.true.,label//' title')
      call assert_logical (name(1).eq.'O' .and. name(2).eq.'H' .and.
     &                     name(3).eq.'H',.true.,label//' names')
      call assert_real (z(1),0.117d0+dz,eps,label//' z1')
      call assert_real (y(2),0.757d0,eps,label//' y2')
      call assert_real (y(3),-0.757d0,eps,label//' y3')
      call assert_real (z(3),-0.469d0+dz,eps,label//' z3')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_readxyz_tinker  --  types and bonds read  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_readxyz_tinker (label)
      use atoms
      use couple
      implicit none
      character*(*) label
c
c
      call assert_int (type(1),1,label//' type1')
      call assert_int (type(3),2,label//' type3')
      call assert_int (n12(1),2,label//' n12 1')
      call assert_int (i12(2,1),3,label//' i12 2,1')
      call assert_int (i12(1,2),1,label//' i12 1,2')
      return
      end
