c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_analyzelmda  --  lambda free energy tests  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_analyzelmda" runs "analyzelmda" on a short orthogonal
c     space tempering and adaptive biasing force history, comparing
c     every file written against its stored reference; on a hand made
c     metadynamics history with one gaussian, whose free energy and
c     profile have a closed form; on the history written by a short
c     metadynamics run; on hand made TI files whose totals are exact,
c     whose errors follow from the spread of their block averages, or
c     whose totals cannot be found; and on history files written by
c     another method than the one chosen; the bin count and theta map
c     keyword rules shared with dynamics are checked directly, and the
c     FEP mode is checked by "test_bar" on its FEP files
c
c
      subroutine test_analyzelmda
      implicit none
      logical skiptest
c
c
      if (skiptest('test_analyzelmda','lambda'))  return
      call test_lmda_ost
      call test_lmda_abf
      call test_lmda_meta
      call test_lmda_metadyn
      call test_lmda_ti
      call test_lmda_wrong
      call test_lmda_keys
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_lmda_ost  --  OST output regression  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_lmda_ost" analyzes a short orthogonal space tempering
c     history on the coarser grid set by its keyfile, and compares the
c     time series, free energy and g kernel files to their references
c
c
      subroutine test_lmda_ost
      implicit none
      integer ist
c
c
      call pushdir ('file/analyzelmda')
      call lmdawipe ('wost')
      call run_prog ('analyzelmda','1 wost','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_analyzelmda ost status')
         call lmdacmp ('wost.series','test_analyzelmda ost series')
         call lmdacmp ('wost.freeenergy',
     &                 'test_analyzelmda ost free energy')
         call lmdacmp ('wost.gkernel','test_analyzelmda ost gkernel')
      end if
      call lmdawipe ('wost')
      call popdir
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_lmda_abf  --  ABF output regression  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_lmda_abf" analyzes a short adaptive biasing force history
c     on the finer grid set by its keyfile, compares the time series
c     and free energy files to their references, and checks that no
c     g kernel file is written
c
c
      subroutine test_lmda_abf
      implicit none
      integer ist
      logical exist
c
c
      call pushdir ('file/analyzelmda')
      call lmdawipe ('wabf')
      call run_prog ('analyzelmda','3 wabf','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_analyzelmda abf status')
         call lmdacmp ('wabf.series','test_analyzelmda abf series')
         call lmdacmp ('wabf.freeenergy',
     &                 'test_analyzelmda abf free energy')
         inquire (file='wabf.gkernel',exist=exist)
         call assert_logical (exist,.false.,
     &                        'test_analyzelmda abf no gkernel')
      end if
      call lmdawipe ('wabf')
      call popdir
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_lmda_meta  --  one gaussian metadynamics  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_lmda_meta" analyzes a metadynamics history holding one
c     gaussian of height 0.5 and width 0.1 at lambda 0.3, cut off at
c     four widths; with its image at -0.3 the bias gives a free energy
c     change of Vbias(0) - Vbias(1) = 0.044318484119 and a bias part
c     of G(0.5) of Vbias(0) - Vbias(0.5) = -0.225636348447 Kcal/mol,
c     to which the profile adds kT ln(J(0.5)/J(0)) for the theta map,
c     with the default smoothed triangle map and with the sine squared
c     map; the profile runs from zero to the free energy change, and
c     a keyfile bin count below three rebuilds it on three bins
c
c
      subroutine test_lmda_meta
      use dlmda
      use units
      implicit none
      integer ist
      real*8 eps
      real*8 value,err
      real*8 dg,g0,g1
      real*8 sig,cut,rt
      real*8 gbias,gmid
      real*8 lmdathjac
      logical found
c
c
      eps = 1.0d-9
      sig = 0.1d0
      cut = 4.0d0 * sig
      rt = gasconst * 298.0d0
      gbias = -0.225636348447d0
      call pushdir ('file/analyzelmda')
      call lmdawipe ('meta1')
      call run_prog ('analyzelmda','2 meta1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_analyzelmda meta status')
         call lmdacmp ('meta1.series','test_analyzelmda meta series')
         call barfind ('meta1.freeenergy','Delta G',found,dg,err)
         call assert_logical (found,.true.,
     &                        'test_analyzelmda meta delta g found')
         call assert_real (dg,0.044318484119d0,eps,
     &                     'test_analyzelmda meta delta g')
         call lmdaprof ('meta1.freeenergy',0.5d0,found,value)
         call assert_logical (found,.true.,
     &                        'test_analyzelmda meta profile found')
         lmdathmap = 'TRI'
         lmdathalpha = 0.999999999d0
         gmid = gbias + rt*log(lmdathjac(0.5d0,sig,cut)
     &                         /lmdathjac(0.0d0,sig,cut))
         call assert_real (value,gmid,eps,
     &                     'test_analyzelmda meta profile middle')
         call lmdaprof ('meta1.freeenergy',0.0d0,found,g0)
         call assert_real (g0,0.0d0,eps,
     &                     'test_analyzelmda meta profile start')
         call lmdaprof ('meta1.freeenergy',1.0d0,found,g1)
         call assert_real (g1,dg,eps,
     &                     'test_analyzelmda meta profile end')
      end if
      call lmdawipe ('meta1')
c
c     the sine squared theta map named in the keyfile changes the
c     profile but not the free energy change
c
      call lmdawipe ('meta1sin')
      call run_prog ('analyzelmda','2 meta1sin','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_analyzelmda meta sin status')
         call barfind ('meta1sin.freeenergy','Delta G',found,value,err)
         call assert_real (value,dg,eps,
     &                     'test_analyzelmda meta sin delta g')
         call lmdaprof ('meta1sin.freeenergy',0.5d0,found,value)
         lmdathmap = 'SIN'
         gmid = gbias + rt*log(lmdathjac(0.5d0,sig,cut)
     &                         /lmdathjac(0.0d0,sig,cut))
         call assert_real (value,gmid,eps,
     &                     'test_analyzelmda meta sin profile middle')
         call lmdaprof ('meta1sin.freeenergy',1.0d0,found,g1)
         call assert_real (g1,dg,eps,
     &                     'test_analyzelmda meta sin profile end')
      end if
      call lmdawipe ('meta1sin')
      lmdathmap = 'TRI'
c
c     a bin count below three rebuilds the history on three bins, as
c     in dynamics, leaving the free energy and the profile unchanged
c
      call lmdawipe ('meta1nb1')
      call run_prog ('analyzelmda','2 meta1nb1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,'test_analyzelmda meta nbin 1 status')
         call barfind ('meta1nb1.freeenergy','Delta G',found,value,err)
         call assert_real (value,dg,eps,
     &                     'test_analyzelmda meta nbin 1 delta g')
         call lmdaprof ('meta1nb1.freeenergy',0.5d0,found,value)
         call assert_logical (found,.true.,
     &                        'test_analyzelmda meta nbin 1 middle bin')
         lmdathmap = 'TRI'
         lmdathalpha = 0.999999999d0
         gmid = gbias + rt*log(lmdathjac(0.5d0,sig,cut)
     &                         /lmdathjac(0.0d0,sig,cut))
         call assert_real (value,gmid,eps,
     &                     'test_analyzelmda meta nbin 1 middle')
         call lmdaprof ('meta1nb1.freeenergy',0.1d0,found,value)
         call assert_logical (found,.false.,
     &                        'test_analyzelmda meta nbin 1 three bins')
      end if
      call lmdawipe ('meta1nb1')
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_lmda_metadyn  --  metadynamics round trip  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_lmda_metadyn" runs a short metadynamics simulation that
c     writes its history through the shared lambda bias header, then
c     checks that the history is read back with the number of gaussians
c     and the free energy given in the header by dynamics; the wide
c     gaussians reach the end points, and the run starts off center
c     with taller gaussians, so the free energy compared is not zero
c
c
      subroutine test_lmda_metadyn
      implicit none
      integer i,ist,ios
      integer nhead(4)
      integer unit,freeunit
      real*8 value,err
      real*8 head(7)
      logical found
      character*240 line
      character*512 cmd
c
c
c     copy the structure and keyfile into a scratch directory, and
c     switch the keyfile to metadynamics with frequent, wide and tall
c     deposits that start away from the middle of the lambda range
c
      call pushdir ('file/mutate')
      cmd = 'rm -rf ../lmda_metadyn ; mkdir -p ../lmda_metadyn ;'//
     &      ' cp water2.xyz ../lmda_metadyn/ ;'//
     &      ' cp vsoft_default_ost.key ../lmda_metadyn/water2.key'
      call execute_command_line (cmd)
      call popdir
      call pushdir ('file/lmda_metadyn')
      call tnist_append ('water2.key','lambda-mode META')
      call tnist_append ('water2.key','lambda-interval 10')
      call tnist_append ('water2.key','wlhist 0.1')
      call tnist_append ('water2.key','lambda 0.3')
      call tnist_append ('water2.key','hbias 0.01')
      call run_prog ('dynamic','water2 40 1.0 0.01 2 298','out.txt',
     &               ist)
      if (ist .eq. -1)  goto 30
      call assert_int (ist,0,'test_analyzelmda metadyn dynamics')
c
c     read the gaussian count and free energy from the history header
c
      found = .false.
      unit = freeunit ()
      open (unit=unit,file='water2.meta',status='old',iostat=ios)
      if (ios .eq. 0) then
         read (unit,'(a)',iostat=ios)  line
         found = (index(line,'Metadynamics History') .ne. 0)
         read (unit,'(a)',iostat=ios)  line
         read (unit,*,iostat=ios)  (nhead(i),i=1,4)
         read (unit,'(a)',iostat=ios)  line
         read (unit,'(a)',iostat=ios)  line
         read (unit,'(a)',iostat=ios)  line
         read (unit,*,iostat=ios)  (head(i),i=1,7)
         if (ios .ne. 0)  found = .false.
         close (unit=unit)
      end if
      call assert_logical (found,.true.,
     &                     'test_analyzelmda metadyn header')
      if (.not. found)  goto 30
      call assert_logical (abs(head(7)).gt.1.0d-6,.true.,
     &                     'test_analyzelmda metadyn nonzero delta g')
c
c     analyze the history and compare with the header values
c
      call run_prog ('analyzelmda','2 water2','out2.txt',ist)
      if (ist .eq. -1)  goto 30
      call assert_int (ist,0,'test_analyzelmda metadyn status')
      call barfind ('water2.freeenergy','Number of Gaussians',found,
     &              value,err)
      call assert_real (value,dble(nhead(3)),0.0d0,
     &                  'test_analyzelmda metadyn gaussians')
      call barfind ('water2.freeenergy','Delta G',found,value,err)
      call assert_real (value,head(7),1.0d-9,
     &                  'test_analyzelmda metadyn delta g')
   30 continue
      call popdir
      call tnist_clean ('lmda_metadyn')
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_lmda_ti  --  free energy from TI files  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_lmda_ti" analyzes hand made TI files; a derivative linear
c     in lambda and constant within each window, written from lambda
c     1 down to 0, gives a trapezoid total of exactly 4 Kcal/mol with
c     no error; two blocks per window at 0, 0.5 and 1 with means 2, 4
c     and 5 and standard errors 1, 2 and 0 give 3.75 +/- sqrt(1.0625);
c     windows short of lambda 1, or fewer than the file header gives,
c     leave the totals out, and a single window stops with a message
c
c
      subroutine test_lmda_ti
      implicit none
      integer ist
      real*8 eps
      real*8 value,err
      logical found,exist
c
c
      eps = 1.0d-4
      call pushdir ('file/analyzelmda')
c
c     linear derivative, with the default and a given trial count
c
      call lmdawipe ('ti_linear')
      call run_prog ('analyzelmda','5 ti_linear','out.txt',ist)
      if (ist .eq. -1)  goto 10
      call assert_int (ist,0,'test_analyzelmda ti linear status')
      call barfind ('ti_linear.freeenergy','TI Full Sample',found,
     &              value,err)
      call assert_real (value,4.0d0,eps,
     &                  'test_analyzelmda ti linear total')
      call assert_real (err,0.0d0,eps,
     &                  'test_analyzelmda ti linear error')
      call barfind ('ti_linear.freeenergy','TI Bootstrap',found,
     &              value,err)
      call assert_real (value,4.0d0,eps,
     &                  'test_analyzelmda ti linear bootstrap')
      call assert_real (err,0.0d0,eps,
     &                  'test_analyzelmda ti linear bootstrap error')
      call lmdawipe ('ti_linear')
      call run_prog ('analyzelmda','5 ti_linear 50','out.txt',ist)
      call barfind ('ti_linear.freeenergy','Bootstrap Trials :',found,
     &              value,err)
      call assert_real (value,50.0d0,0.0d0,
     &                  'test_analyzelmda ti trials given')
      call lmdawipe ('ti_linear')
c
c     standard errors from the spread of the block averages
c
      call lmdawipe ('ti_spread')
      call run_prog ('analyzelmda','5 ti_spread','out.txt',ist)
      call barfind ('ti_spread.freeenergy','TI Full Sample',found,
     &              value,err)
      call assert_real (value,3.75d0,eps,
     &                  'test_analyzelmda ti spread total')
      call assert_real (err,sqrt(1.0625d0),eps,
     &                  'test_analyzelmda ti spread error')
      call barfind ('ti_spread.freeenergy','TI Bootstrap',found,
     &              value,err)
      call assert_logical (found,.true.,
     &                     'test_analyzelmda ti spread bootstrap')
      call lmdawipe ('ti_spread')
c
c     windows short of lambda 1, and fewer than in the header
c
      call lmdawipe ('ti_partial')
      call run_prog ('analyzelmda','5 ti_partial','out.txt',ist)
      call barfind ('ti_partial.freeenergy','do not Span 0 to 1',
     &              found,value,err)
      call assert_logical (found,.true.,
     &                     'test_analyzelmda ti partial note')
      call barfind ('ti_partial.freeenergy','TI Full Sample',found,
     &              value,err)
      call assert_logical (found,.false.,
     &                     'test_analyzelmda ti partial no total')
      call lmdawipe ('ti_partial')
      call lmdawipe ('ti_header')
      call run_prog ('analyzelmda','5 ti_header','out.txt',ist)
      call barfind ('ti_header.freeenergy','Differ from the File',
     &              found,value,err)
      call assert_logical (found,.true.,
     &                     'test_analyzelmda ti header note')
      call barfind ('ti_header.freeenergy','Windows in Header :',
     &              found,value,err)
      call assert_real (value,5.0d0,0.0d0,
     &                  'test_analyzelmda ti header count')
      call lmdawipe ('ti_header')
c
c     a single window stops before any file is written
c
      call lmdawipe ('ti_one')
      call run_prog ('analyzelmda','5 ti_one','out.txt',ist)
      call barfind ('out.txt','Fewer than Two Lambda Windows',found,
     &              value,err)
      call assert_logical (found,.true.,
     &                     'test_analyzelmda ti one window note')
      inquire (file='ti_one.freeenergy',exist=exist)
      call assert_logical (exist,.false.,
     &                     'test_analyzelmda ti one window no file')
      call lmdawipe ('ti_one')
   10 continue
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_lmda_wrong  --  history of another method  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_lmda_wrong" chooses a method other than the one that wrote
c     each history file, and checks that the run stops with a message
c     before any output file is written
c
c
      subroutine test_lmda_wrong
      implicit none
      integer ist
      real*8 value,err
      logical found,exist
c
c
      call pushdir ('file/analyzelmda')
      call lmdawipe ('wabf')
      call run_prog ('analyzelmda','1 wabf.abf','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','OST History not Found',found,
     &                 value,err)
         call assert_logical (found,.true.,
     &                        'test_analyzelmda ost on abf file')
         inquire (file='wabf.series',exist=exist)
         call assert_logical (exist,.false.,
     &                        'test_analyzelmda ost on abf no file')
      end if
      call lmdawipe ('wost')
      call run_prog ('analyzelmda','2 wost.ost','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Metadynamics History not Found',
     &                 found,value,err)
         call assert_logical (found,.true.,
     &                        'test_analyzelmda meta on ost file')
      end if
      call lmdawipe ('wabf')
      call lmdawipe ('wost')
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine lmdawipe  --  remove the analysis outputs  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "lmdawipe" removes every version of the files "analyzelmda"
c     writes for the given base name, so each run starts clean
c
c
      subroutine lmdawipe (base)
      implicit none
      character*(*) base
c
c
      call execute_command_line ('rm -f '//trim(base)//'.series* '//
     &                           trim(base)//'.freeenergy* '//
     &                           trim(base)//'.gkernel*')
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine lmdacmp  --  compare one analysis output  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "lmdacmp" compares an output file of "analyzelmda" with the
c     reference of the same name
c
c
      subroutine lmdacmp (fname,label)
      implicit none
      real*8 eps
      character*(*) fname,label
      character*240 rpath
c
c
      eps = 1.0d-6
      call refpath ('analyzelmda',fname,rpath)
      call assert_files (rpath,fname,eps,label)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine lmdaprof  --  free energy profile at lambda  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "lmdaprof" finds the value of the free energy profile at the
c     given lambda in a metadynamics free energy file
c
c
      subroutine lmdaprof (fname,lambda,found,value)
      implicit none
      integer ios,k
      integer unit,freeunit
      real*8 lambda,value
      real*8 x,g
      logical found,table
      character*(*) fname
      character*240 line
c
c
      found = .false.
      table = .false.
      value = 0.0d0
      unit = freeunit ()
      open (unit=unit,file=fname,status='old',iostat=ios)
      do while (ios .eq. 0)
         read (unit,'(a)',iostat=ios)  line
         if (ios .ne. 0)  goto 10
         if (index(line,'Free Energy Profile') .ne. 0) then
            table = .true.
         else if (table) then
            read (line,*,iostat=k)  x,g
            if (k.eq.0 .and. abs(x-lambda).lt.1.0d-6) then
               found = .true.
               value = g
               goto 10
            end if
         end if
      end do
   10 continue
      close (unit=unit)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_lmda_keys  --  shared lambda keywords  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_lmda_keys" checks the bin count rule shared by dynamics
c     and the analysis keys, which keeps an odd count of at least
c     three, and the theta map keywords with their defaults and the
c     replacement of an unknown map or a sharpness outside (0,1)
c
c
      subroutine test_lmda_keys
      use dlmda
      use keys
      implicit none
      integer lmdanbin
      real*8 eps
c
c
      eps = 1.0d-12
      call assert_int (lmdanbin(-4),3,'lmdanbin negative count')
      call assert_int (lmdanbin(0),3,'lmdanbin zero count')
      call assert_int (lmdanbin(1),3,'lmdanbin one bin')
      call assert_int (lmdanbin(2),3,'lmdanbin two bins')
      call assert_int (lmdanbin(4),5,'lmdanbin even count')
      call assert_int (lmdanbin(11),11,'lmdanbin odd count')
c
c     no keywords give the smoothed triangle map and its sharpness
c
      if (allocated(keyline))  deallocate (keyline)
      allocate (keyline(2))
      nkey = 0
      lmdathmap = 'SIN'
      lmdathalpha = 0.5d0
      call lmdathkey
      call assert_logical (lmdathmap.eq.'TRI',.true.,
     &                     'lmdathkey default map')
      call assert_real (lmdathalpha,0.999999999d0,eps,
     &                  'lmdathkey default sharpness')
c
c     a named map and a valid sharpness are kept
c
      nkey = 2
      keyline(1) = 'lambda-theta-map sin'
      keyline(2) = 'LAMBDA-THETA-ALPHA 0.5'
      call lmdathkey
      call assert_logical (lmdathmap.eq.'SIN',.true.,
     &                     'lmdathkey sine squared map')
      call assert_real (lmdathalpha,0.5d0,eps,
     &                  'lmdathkey given sharpness')
c
c     an unknown map and a sharpness outside (0,1) take the defaults
c
      keyline(1) = 'LAMBDA-THETA-MAP COS'
      keyline(2) = 'LAMBDA-THETA-ALPHA 1.5'
      call lmdathkey
      call assert_logical (lmdathmap.eq.'TRI',.true.,
     &                     'lmdathkey unknown map')
      call assert_real (lmdathalpha,0.999999999d0,eps,
     &                  'lmdathkey sharpness out of range')
c
c     restore an empty keyword list and the default theta map
c
      nkey = 0
      deallocate (keyline)
      call lmdathkey
      return
      end
