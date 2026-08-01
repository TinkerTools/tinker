c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_tinkernist  --  Tinker NIST tests  ##
c     ##                                                     ##
c     #########################################################
c
c
c     "test_tinkernist" runs the dynamic.x regression tests mirrored
c     from the tinker-gpu tinkernist.cpp test cases
c
c
      subroutine test_tinkernist
      implicit none
c
c
      call test_tinkernist_nocoord
      call test_tinkernist_dcd
      call test_tinkernist_save1
      call test_tinkernist_save2
      call test_tinkernist_save_only
      call test_tinkernist_save_system
      call test_tinkernist_exc_moment
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine test_tinkernist_nocoord  --  no coords  ##
c     ##                                                     ##
c     #########################################################
c
c
      subroutine test_tinkernist_nocoord
      implicit none
      integer ist
      logical skiptest
c
c
      if (skiptest('test_tinkernist_nocoord','amoeba'))  return
      call tnist_prep ('tinkernist_nocoord','water30',.false.)
      call pushdir ('file/tinkernist_nocoord')
      call tnist_append ('water30.key','nocoord')
      call run_prog ('dynamic','water30 1 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_tinkernist_nocoord status')
         call tnist_assert_file ('water30.arc',.false.,
     &                           'test_tinkernist_nocoord arc')
         call tnist_assert_file ('water30.dyn',.true.,
     &                           'test_tinkernist_nocoord dyn')
      end if
      call popdir
      call tnist_clean ('tinkernist_nocoord')
      return
      end
c
c
c     ###############################################
c     ##                                           ##
c     ##  subroutine test_tinkernist_dcd  --  DCD  ##
c     ##                                           ##
c     ###############################################
c
c
      subroutine test_tinkernist_dcd
      implicit none
      integer ist
      logical skiptest
c
c
      if (skiptest('test_tinkernist_dcd','amoeba'))  return
      call tnist_prep ('tinkernist_dcd','water30',.true.)
      call pushdir ('file/tinkernist_dcd')
      call tnist_append_saves ('water30.key')
      call tnist_append ('water30.key','dcd-archive')
      call run_prog ('dynamic','water30 2 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_tinkernist_dcd status')
         call tnist_assert_file ('water30.dcd',.true.,
     &                           'test_tinkernist_dcd dcd')
         call tnist_assert_file ('water30.dcdv',.true.,
     &                           'test_tinkernist_dcd dcdv')
         call tnist_assert_file ('water30.dcduc',.true.,
     &                           'test_tinkernist_dcd dcduc')
         call tnist_assert_file ('water30.dcdus',.true.,
     &                           'test_tinkernist_dcd dcdus')
         call tnist_assert_file ('water30.dcdui',.true.,
     &                           'test_tinkernist_dcd dcdui')
         call tnist_assert_file ('water30.dcdud',.true.,
     &                           'test_tinkernist_dcd dcdud')
         call tnist_assert_file ('water30.dcdde',.true.,
     &                           'test_tinkernist_dcd dcdde')
         call tnist_assert_file ('water30.dcdte',.true.,
     &                           'test_tinkernist_dcd dcdte')
      end if
      call popdir
      call tnist_clean ('tinkernist_dcd')
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine test_tinkernist_save1  --  AMOEBA  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine test_tinkernist_save1
      implicit none
      logical skiptest
c
c
      if (skiptest('test_tinkernist_save1','amoeba'))  return
      call tnist_save_case ('tinkernist_save1','water30',
     &                      'water30',20)
      return
      end
c
c
c     ###################################################
c     ##                                               ##
c     ##  subroutine test_tinkernist_save2  --  HIPPO  ##
c     ##                                               ##
c     ###################################################
c
c
      subroutine test_tinkernist_save2
      implicit none
      logical skiptest
c
c
      if (skiptest('test_tinkernist_save2','amoeba'))  return
      call tnist_save_case ('tinkernist_save2','water30_2',
     &                      'water30_2',20)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_tinkernist_save_only  --  SAVE-ONLY  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_tinkernist_save_only
      implicit none
      integer ist
      logical skiptest
c
c
      if (skiptest('test_tinkernist_save_only','amoeba'))  return
      call tnist_prep ('tinkernist_save_only','water30',.true.)
      call pushdir ('file/tinkernist_save_only')
      call tnist_append ('water30.key','save-only -4 102')
      call tnist_append_saves ('water30.key')
      call run_prog ('dynamic','water30 2 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_tinkernist_save_only status')
         call tnist_assert_count ('water30.arc',99,2,
     &                            'test_tinkernist_save_only arc')
         call tnist_assert_count ('water30.vel',99,2,
     &                            'test_tinkernist_save_only vel')
         call tnist_assert_count ('water30.uchg',99,2,
     &                            'test_tinkernist_save_only uchg')
         call tnist_assert_count ('water30.ustc',99,2,
     &                            'test_tinkernist_save_only ustc')
         call tnist_assert_count ('water30.uind',99,2,
     &                            'test_tinkernist_save_only uind')
         call tnist_assert_count ('water30.udir',99,2,
     &                            'test_tinkernist_save_only udir')
         call tnist_assert_count ('water30.def',99,2,
     &                            'test_tinkernist_save_only def')
         call tnist_assert_count ('water30.tef',99,2,
     &                            'test_tinkernist_save_only tef')
      end if
      call popdir
      call tnist_clean ('tinkernist_save_only')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_tinkernist_save_system  --  moments  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_tinkernist_save_system
      implicit none
      integer ist
      logical skiptest
c
c
      if (skiptest('test_tinkernist_save_system','amoeba'))  return
      call tnist_prep ('tinkernist_save_system','water30',.true.)
      call pushdir ('file/tinkernist_save_system')
      call tnist_append ('water30.key','save-usystem')
      call tnist_append ('water30.key','save-vsystem')
      call tnist_append ('water30.key','nocoord')
      call run_prog ('dynamic','water30 2 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_tinkernist_save_system status')
         call tnist_check_save_system ('out.txt')
      end if
      call popdir
      call tnist_clean ('tinkernist_save_system')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_tinkernist_exc_moment  --  excluded  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_tinkernist_exc_moment
      implicit none
      integer ist
      logical skiptest
c
c
      if (skiptest('test_tinkernist_exc_moment','amoeba'))  return
      call tnist_prep ('tinkernist_exc_moment','water30',.true.)
      call pushdir ('file/tinkernist_exc_moment')
      call tnist_append ('water30.key','exc-moment -2677 2684')
      call tnist_append ('water30.key','save-usystem')
      call tnist_append ('water30.key','nocoord')
      call run_prog ('dynamic','water30 2 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_tinkernist_exc_moment status')
         call tnist_check_exc_moment ('out.txt')
      end if
      call popdir
      call tnist_clean ('tinkernist_exc_moment')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine tnist_save_case  --  run save regression  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine tnist_save_case (work,base,refbase,nstep)
      implicit none
      integer ist,nstep
      character*(*) work,base,refbase
      character*12 tnist_itoa
      character*240 rpath
c
c
      call tnist_prep (work,base,.true.)
      call pushdir ('file/'//trim(work))
      call tnist_append_saves (trim(base)//'.key')
      call run_prog ('dynamic',trim(base)//' '//trim(tnist_itoa(nstep))
     &              //' 1 0.001 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_'//trim(work)//' status')
         call refpath ('tinkernist',trim(refbase)//'.arc',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.arc',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' arc')
         call refpath ('tinkernist',trim(refbase)//'.vel',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.vel',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' vel')
         call refpath ('tinkernist',trim(refbase)//'.uchg',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.uchg',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' uchg')
         call refpath ('tinkernist',trim(refbase)//'.ustc',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.ustc',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' ustc')
         call refpath ('tinkernist',trim(refbase)//'.uind',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.uind',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' uind')
         call refpath ('tinkernist',trim(refbase)//'.udir',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.udir',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' udir')
         call refpath ('tinkernist',trim(refbase)//'.def',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.def',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' def')
         call refpath ('tinkernist',trim(refbase)//'.tef',rpath)
         call tnist_compare_vecfile (rpath,trim(base)//'.tef',nstep,
     &                               1.0d-4,'test_'//trim(work)//
     &                               ' tef')
      end if
      call popdir
      call tnist_clean (work)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine tnist_prep  --  create scratch fixture  ##
c     ##                                                     ##
c     #########################################################
c
c
      subroutine tnist_prep (work,base,usedyn)
      implicit none
      character*(*) work,base
      logical usedyn
      character*512 cmd
c
c
      call pushdir ('file/tinkernist')
      cmd = 'rm -rf ../'//trim(work)//' ; mkdir -p ../'//trim(work)//
     &      ' ; cp '//trim(base)//'.xyz '//trim(base)//'.key ../'//
     &      trim(work)//'/'
      call execute_command_line (cmd)
      if (usedyn) then
         cmd = 'cp '//trim(base)//'.dyn ../'//trim(work)//'/'
         call execute_command_line (cmd)
      end if
      call popdir
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine tnist_clean  --  remove scratch files  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine tnist_clean (work)
      implicit none
      character*(*) work
      character*512 cmd
c
c
      call pushdir ('file')
      cmd = 'rm -rf '//trim(work)
      call execute_command_line (cmd)
      call popdir
      return
      end
c
c
c     ####################################################
c     ##                                                ##
c     ##  subroutine tnist_append  --  append key line  ##
c     ##                                                ##
c     ####################################################
c
c
      subroutine tnist_append (key,line)
      implicit none
      integer unit,freeunit
      character*(*) key,line
c
c
      unit = freeunit ()
      open (unit=unit,file=key,status='old',position='append')
      write (unit,10)  trim(line)
   10 format (a)
      close (unit=unit)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine tnist_append_saves  --  save keywords  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine tnist_append_saves (key)
      implicit none
      character*(*) key
c
c
      call tnist_append (key,'save-ucharge')
      call tnist_append (key,'save-ustatic')
      call tnist_append (key,'save-uinduce')
      call tnist_append (key,'save-udirect')
      call tnist_append (key,'save-defield')
      call tnist_append (key,'save-tefield')
      call tnist_append (key,'save-velocity')
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine tnist_assert_file  --  file existence  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine tnist_assert_file (file,want,label)
      implicit none
      integer got,ref
      logical exist,want
      character*(*) file,label
c
c
      inquire (file=file,exist=exist)
      got = 0
      ref = 0
      if (exist)  got = 1
      if (want)  ref = 1
      call assert_int (got,ref,label)
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine tnist_compare_vecfile  --  vector compare  ##
c     ##                                                        ##
c     ############################################################
c
c
      subroutine tnist_compare_vecfile (reffile,outfile,nframe,eps,
     &                                  label)
      use assert
      implicit none
      integer nframe,ir,io,iosr,ioso,fr,nr,no,i,j
      integer ar,ao,freeunit,iw,jw,fw
      real*8 eps,rv(3),ov(3),diff,worst
      character*(*) reffile,outfile,label
      character*512 lr,lo
      character*8 name
c
c
      worst = 0.0d0
      fw = 0
      iw = 0
      jw = 0
      ir = freeunit ()
      open (unit=ir,file=reffile,status='old',iostat=iosr)
      io = freeunit ()
      open (unit=io,file=outfile,status='old',iostat=ioso)
      if (iosr.ne.0 .or. ioso.ne.0) then
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,*) 'tnist_compare_vecfile: open failed'
         call assert_summary ()
      end if
      do fr = 1, nframe
         read (ir,'(a)',iostat=iosr) lr
         read (io,'(a)',iostat=ioso) lo
         if (iosr.ne.0 .or. ioso.ne.0)  go to 90
         read (lr,*)  nr
         read (lo,*)  no
         if (nr .ne. no)  go to 90
         read (ir,'(a)',iostat=iosr) lr
         read (io,'(a)',iostat=ioso) lo
         if (iosr.ne.0 .or. ioso.ne.0)  go to 90
         do i = 1, nr
            read (ir,'(a)',iostat=iosr) lr
            read (io,'(a)',iostat=ioso) lo
            if (iosr.ne.0 .or. ioso.ne.0)  go to 90
            if (i.le.2 .or. i.gt.nr-2) then
               read (lr,*)  ar,name,rv(1),rv(2),rv(3)
               read (lo,*)  ao,name,ov(1),ov(2),ov(3)
               do j = 1, 3
                  diff = abs(ov(j)-rv(j))
                  if (diff .gt. worst) then
                     worst = diff
                     fw = fr
                     iw = i
                     jw = j
                  end if
               end do
            end if
         end do
      end do
      close (unit=ir)
      close (unit=io)
      if (worst .le. eps) then
         npass = npass + 1
         write (*,10)  'pass ',trim(label)
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,20)  fw,iw,jw,worst
         call assert_summary ()
      end if
      return
   90 continue
      close (unit=ir)
      close (unit=io)
      nfail = nfail + 1
      write (*,10)  'FAIL ',trim(label)
      write (*,*) 'tnist_compare_vecfile: malformed file'
      call assert_summary ()
   10 format (1x,a5,1x,a)
   20 format (7x,'frame=',i0,' row=',i0,' comp=',i0,' d=',g11.3)
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine tnist_assert_count  --  saved site count  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine tnist_assert_count (file,want,nframe,label)
      implicit none
      integer want,nframe,fr,nat,row,unit,ios,freeunit
      character*(*) file,label
      character*512 line
c
c
      unit = freeunit ()
      open (unit=unit,file=file,status='old',iostat=ios)
      if (ios .ne. 0) then
         call assert_int (0,want,label)
         return
      end if
      do fr = 1, nframe
         read (unit,'(a)',iostat=ios) line
         if (ios .ne. 0) then
            nat = -1
            go to 90
         end if
         read (line,*)  nat
         if (nat .ne. want)  go to 90
         read (unit,'(a)',iostat=ios) line
         if (ios .ne. 0)  go to 90
         do row = 1, nat
            read (unit,'(a)',iostat=ios) line
            if (ios .ne. 0)  go to 90
         end do
      end do
   90 continue
      close (unit=unit)
      call assert_int (nat,want,label)
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine tnist_check_save_system  --  stdout refs  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine tnist_check_save_system (file)
      implicit none
      character*(*) file
      real*8 eps
      real*8 got1(3,2),got2(3,8)
      real*8 ref1(3,2),ref2(3,8)
      integer n
c
c
      eps = 1.0d-4
      call tnist_load_log1 (file,'System Charge Dipole',got1,2,n)
      ref1(1,1)=6.720926d0
      ref1(2,1)=-27.910876d0
      ref1(3,1)=75.185866d0
      ref1(1,2)=6.531183d0
      ref1(2,2)=-28.030485d0
      ref1(3,2)=75.139570d0
      call tnist_compare_vecs (got1,ref1,2,eps,
     &                         'test_tinkernist_save_system uchg')
      call tnist_load_log1 (file,'System Static Dipole',got1,2,n)
      ref1(1,1)=-7.383698d0
      ref1(2,1)=-7.883884d0
      ref1(3,1)=6.462269d0
      ref1(1,2)=-7.381055d0
      ref1(2,2)=-7.911815d0
      ref1(3,2)=6.441521d0
      call tnist_compare_vecs (got1,ref1,2,eps,
     &                         'test_tinkernist_save_system ustc')
      call tnist_load_log1 (file,'System Induced Dipole',got1,2,n)
      ref1(1,1)=-26.916924d0
      ref1(2,1)=-16.364259d0
      ref1(3,1)=4.850315d0
      ref1(1,2)=-27.126473d0
      ref1(2,2)=-16.509707d0
      ref1(3,2)=4.743480d0
      call tnist_compare_vecs (got1,ref1,2,eps,
     &                         'test_tinkernist_save_system uind')
c
      call tnist_load_log2 (file,'Charge Dipole by Atom Type:',
     &                      got2,8,n)
      call tnist_ref_save_system_auchg (ref2)
      call tnist_compare_vecs (got2,ref2,8,eps,
     &                         'test_tinkernist_save_system auchg')
      call tnist_load_log2 (file,'Static Dipole by Atom Type:',
     &                      got2,8,n)
      call tnist_ref_save_system_austc (ref2)
      call tnist_compare_vecs (got2,ref2,8,eps,
     &                         'test_tinkernist_save_system austc')
      call tnist_load_log2 (file,'Induced Dipole by Atom Type:',
     &                      got2,8,n)
      call tnist_ref_save_system_auind (ref2)
      call tnist_compare_vecs (got2,ref2,8,eps,
     &                         'test_tinkernist_save_system auind')
      call tnist_load_log2 (file,'Velocity by Atom Type:',
     &                      got2,8,n)
      call tnist_ref_save_system_avelo (ref2)
      call tnist_compare_vecs (got2,ref2,8,eps,
     &                         'test_tinkernist_save_system avelo')
      return
      end
c
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine tnist_check_exc_moment  --  stdout refs  ##
c     ##                                                      ##
c     ##########################################################
c
c
      subroutine tnist_check_exc_moment (file)
      implicit none
      character*(*) file
      real*8 eps
      real*8 got(3,2),ref(3,2)
      integer n
c
c
      eps = 1.0d-4
      call tnist_load_log1 (file,'System Charge Dipole',got,2,n)
      ref(1,1)=-28.410462d0
      ref(2,1)=-29.982911d0
      ref(3,1)=24.644113d0
      ref(1,2)=-28.554250d0
      ref(2,2)=-30.125104d0
      ref(3,2)=24.570208d0
      call tnist_compare_vecs (got,ref,2,eps,
     &                         'test_tinkernist_exc_moment uchg')
      call tnist_load_log1 (file,'System Static Dipole',got,2,n)
      ref(1,1)=-7.383698d0
      ref(2,1)=-7.883884d0
      ref(3,1)=6.462269d0
      ref(1,2)=-7.381055d0
      ref(2,2)=-7.911815d0
      ref(3,2)=6.441521d0
      call tnist_compare_vecs (got,ref,2,eps,
     &                         'test_tinkernist_exc_moment ustc')
      call tnist_load_log1 (file,'System Induced Dipole',got,2,n)
      ref(1,1)=-25.738935d0
      ref(2,1)=-19.017606d0
      ref(3,1)=6.832111d0
      ref(1,2)=-25.932641d0
      ref(2,2)=-19.158320d0
      ref(3,2)=6.678709d0
      call tnist_compare_vecs (got,ref,2,eps,
     &                         'test_tinkernist_exc_moment uind')
      return
      end
c
c
c     ######################################################
c     ##                                                  ##
c     ##  subroutine tnist_load_log1  --  system moments  ##
c     ##                                                  ##
c     ######################################################
c
c
      subroutine tnist_load_log1 (file,label,vec,mx,n)
      implicit none
      integer mx,n,unit,ios,next,freeunit
      real*8 vec(3,mx)
      character*(*) file,label
      character*512 line
c
c
      n = 0
      unit = freeunit ()
      open (unit=unit,file=file,status='old',iostat=ios)
      if (ios .ne. 0)  return
   10 continue
      read (unit,'(a)',iostat=ios) line
      if (ios .ne. 0)  go to 20
      if (index(line,trim(label)) .gt. 0) then
         if (n .lt. mx) then
            n = n + 1
            next = index(line,trim(label)) + len_trim(label)
            call getfloat (line,vec(1,n),next)
            call getfloat (line,vec(2,n),next)
            call getfloat (line,vec(3,n),next)
         end if
      end if
      go to 10
   20 continue
      close (unit=unit)
      return
      end
c
c
c     #########################################################
c     ##                                                     ##
c     ##  subroutine tnist_load_log2  --  atom type moments  ##
c     ##                                                     ##
c     #########################################################
c
c
      subroutine tnist_load_log2 (file,label,vec,mx,n)
      implicit none
      integer mx,n,unit,ios,k,itype,freeunit
      real*8 vec(3,mx)
      character*(*) file,label
      character*512 line
c
c
      n = 0
      unit = freeunit ()
      open (unit=unit,file=file,status='old',iostat=ios)
      if (ios .ne. 0)  return
   10 continue
      read (unit,'(a)',iostat=ios) line
      if (ios .ne. 0)  go to 30
      if (index(line,trim(label)) .gt. 0) then
         read (unit,'(a)',iostat=ios) line
         do k = 1, 4
            read (unit,'(a)',iostat=ios) line
            if (ios .ne. 0)  go to 30
            if (n .lt. mx) then
               n = n + 1
               read (line,*)  itype,vec(1,n),vec(2,n),vec(3,n)
            end if
         end do
      end if
      go to 10
   30 continue
      close (unit=unit)
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  subroutine tnist_compare_vecs  --  vector arrays  ##
c     ##                                                    ##
c     ########################################################
c
c
      subroutine tnist_compare_vecs (got,ref,n,eps,label)
      use assert
      implicit none
      integer n,i,j,iw,jw
      real*8 got(3,n),ref(3,n),eps,diff,worst
      character*(*) label
c
c
      worst = 0.0d0
      iw = 0
      jw = 0
      do i = 1, n
         do j = 1, 3
            diff = abs(got(j,i)-ref(j,i))
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
      else
         nfail = nfail + 1
         write (*,10)  'FAIL ',trim(label)
         write (*,20)  iw,jw,got(jw,iw),ref(jw,iw),worst
         call assert_summary ()
      end if
   10 format (1x,a5,1x,a)
   20 format (7x,'row=',i0,' comp=',i0,' calc=',g16.8,
     &        ' ref=',g16.8,' d=',g11.3)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine tnist_ref_save_system_auchg  --  charge refs  ##
c     ##                                                           ##
c     ###############################################################
c
c
      subroutine tnist_ref_save_system_auchg (ref)
      implicit none
      real*8 ref(3,8)
c
c
      ref(1,1)=303.681720d0
      ref(2,1)=230.022755d0
      ref(3,1)=292.223414d0
      ref(1,2)=-332.092182d0
      ref(2,2)=-260.005666d0
      ref(3,2)=-267.579302d0
      ref(1,3)=203.333735d0
      ref(2,3)=139.690855d0
      ref(3,3)=202.363501d0
      ref(1,4)=-168.202346d0
      ref(2,4)=-137.618821d0
      ref(3,4)=-151.821748d0
      ref(1,5)=303.676734d0
      ref(2,5)=230.010554d0
      ref(3,5)=292.205929d0
      ref(1,6)=-332.230984d0
      ref(2,6)=-260.135658d0
      ref(3,6)=-267.635721d0
      ref(1,7)=203.312444d0
      ref(2,7)=139.706756d0
      ref(3,7)=202.374783d0
      ref(1,8)=-168.227011d0
      ref(2,8)=-137.612138d0
      ref(3,8)=-151.805421d0
      return
      end
c
c
      subroutine tnist_ref_save_system_austc (ref)
      implicit none
      real*8 ref(3,8)
c
c
      ref(1,1)=-6.783194d0
      ref(2,1)=-7.362902d0
      ref(3,1)=6.004808d0
      ref(1,2)=-0.600503d0
      ref(2,2)=-0.520982d0
      ref(3,2)=0.457461d0
      ref(1,3)=0.0d0
      ref(2,3)=0.0d0
      ref(3,3)=0.0d0
      ref(1,4)=0.0d0
      ref(2,4)=0.0d0
      ref(3,4)=0.0d0
      ref(1,5)=-6.766523d0
      ref(2,5)=-7.383731d0
      ref(3,5)=5.971497d0
      ref(1,6)=-0.614533d0
      ref(2,6)=-0.528084d0
      ref(3,6)=0.470024d0
      ref(1,7)=0.0d0
      ref(2,7)=0.0d0
      ref(3,7)=0.0d0
      ref(1,8)=0.0d0
      ref(2,8)=0.0d0
      ref(3,8)=0.0d0
      return
      end
c
c
      subroutine tnist_ref_save_system_auind (ref)
      implicit none
      real*8 ref(3,8)
c
c
      ref(1,1)=-13.726947d0
      ref(2,1)=-10.959873d0
      ref(3,1)=4.394567d0
      ref(1,2)=-12.011988d0
      ref(2,2)=-8.057732d0
      ref(3,2)=2.437544d0
      ref(1,3)=0.049399d0
      ref(2,3)=-0.011118d0
      ref(3,3)=-0.005801d0
      ref(1,4)=-1.227388d0
      ref(2,4)=2.664465d0
      ref(3,4)=-1.975995d0
      ref(1,5)=-13.830041d0
      ref(2,5)=-11.025101d0
      ref(3,5)=4.293926d0
      ref(1,6)=-12.102600d0
      ref(2,6)=-8.133219d0
      ref(3,6)=2.384783d0
      ref(1,7)=0.049148d0
      ref(2,7)=-0.010589d0
      ref(3,7)=-0.005383d0
      ref(1,8)=-1.242980d0
      ref(2,8)=2.659202d0
      ref(3,8)=-1.929846d0
      return
      end
c
c
      subroutine tnist_ref_save_system_avelo (ref)
      implicit none
      real*8 ref(3,8)
c
c
      ref(1,1)=6.427442d0
      ref(2,1)=9.137643d0
      ref(3,1)=8.614916d0
      ref(1,2)=-115.797508d0
      ref(2,2)=-113.198323d0
      ref(3,2)=-15.953530d0
      ref(1,3)=-4.357217d0
      ref(2,3)=3.195626d0
      ref(3,3)=2.390313d0
      ref(1,4)=5.174714d0
      ref(2,4)=-1.274009d0
      ref(3,4)=-3.400850d0
      ref(1,5)=-0.121014d0
      ref(2,5)=3.595936d0
      ref(3,5)=8.494728d0
      ref(1,6)=-75.282660d0
      ref(2,6)=-82.537639d0
      ref(3,6)=-68.066174d0
      ref(1,7)=-4.484956d0
      ref(2,7)=3.433754d0
      ref(3,7)=2.314456d0
      ref(1,8)=5.111564d0
      ref(2,8)=-1.492540d0
      ref(3,8)=-3.387783d0
      return
      end
c
c
c     ########################################################
c     ##                                                    ##
c     ##  function tnist_itoa  --  integer to text         ##
c     ##                                                    ##
c     ########################################################
c
c
      character*12 function tnist_itoa (ival)
      implicit none
      integer ival
c
c
      write (tnist_itoa,10)  ival
   10 format (i0)
      return
      end
