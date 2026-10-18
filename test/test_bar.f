c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_bar  --  BAR perturbed energy file tests  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_bar" runs "bar" in mode 1 on two pairs of water trajectory
c     archives, creating a Tinker BAR file of perturbed potential
c     energies for each, and checks each generated file against its
c     stored reference with the tolerant "assert_files" comparison;
c     it then checks the shared FEP and BAR calculation directly,
c     including that frame volumes do not change the free energy,
c     the free energy found in mode 3 from FEP files built so that
c     the results do not depend on the random bootstrap samples, the
c     frame selection of mode 2, the PV term in the enthalpy at a
c     given pressure, bootstrap errors on data with constant energy
c     differences, the input errors that stop mode 2, the volumes
c     and pressure read in mode 3, the BAR enthalpy for trajectories
c     with different numbers of frames and the shift of Bennett, the
c     bootstrap versions of the estimators, the errors of the mode 3
c     totals, mode 2 on the water trajectories, exponents that would
c     overflow, bootstrap trials that fail, and mode 1 energies taken
c     from trajectory logs
c
c
      subroutine test_bar
      implicit none
      logical skiptest
c
c
      if (skiptest('test_bar','amoeba'))  return
      call test_bar_case ('water-00-0000-0950','water-00-0000-0900')
      call test_bar_case ('water-00-0950-1000','water-00-0900-1000')
      call test_bar_vol
      call test_bar_overlap
      call test_bar_fep
      call test_bar_cross
      call test_bar_frames
      call test_bar_pv
      call test_bar_err
      call test_bar_input
      call test_bar_fepvol
      call test_bar_hbar
      call test_bar_bennett
      call test_bar_boot
      call test_bar_chain
      call test_bar_ref
      call test_bar_overflow
      call test_bar_nonconv
      call test_bar_logs
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_bar_case  --  one BAR file regression  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "test_bar_case" runs "bar 1 <A>.arc 298 <B>.arc 298 N" from the
c     fixture directory, which writes "<A>.bar", and compares it to the
c     reference of the same name
c
c
      subroutine test_bar_case (basea,baseb)
      implicit none
      real*8 eps
      integer ist
      character*(*) basea,baseb
      character*240 rpath
      character*512 args
c
c
      eps = 1.0d-3
      call pushdir ('file/bar')
      call execute_command_line
     &   ('rm -f '//trim(basea)//'.bar out.txt')
      args = '1 '//trim(basea)//'.arc 298 '//trim(baseb)//'.arc 298 N'
      call run_prog ('bar',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('bar',trim(basea)//'.bar',rpath)
      call assert_files (rpath,trim(basea)//'.bar',eps,
     &                   'test_bar '//trim(basea))
      call execute_command_line
     &   ('rm -f '//trim(basea)//'.bar out.txt')
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_bar_vol  --  volumes in FEP and BAR sums  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_bar_vol" checks that frame volumes, which vary from frame
c     to frame, do not enter the FEP and BAR free energy sums; with a
c     constant energy difference every estimate is exact and every
c     bootstrap error is zero only if no volume term is applied
c
c
      subroutine test_bar_vol
      implicit none
      integer i,n
      integer nbst,istat
      integer nboot(2)
      real*8 temp
      real*8 ua0(20),ua1(20)
      real*8 ub0(20),ub1(20)
      real*8 vola(20),volb(20)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
c
c
c     state 1 lies 0.8 Kcal/mol above state 0 in every frame of
c     both trajectories, while the frame volumes vary by 20 percent
c
      n = 20
      nbst = 200
      temp = 298.0d0
      do i = 1, n
         vola(i) = 1000.0d0 + 200.0d0*sin(dble(i))
         volb(i) = 1000.0d0 + 200.0d0*cos(dble(i))
         ua0(i) = 0.0d0
         ua1(i) = 0.8d0
         ub1(i) = 0.0d0
         ub0(i) = -0.8d0
      end do
      call barpair (n,n,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,0,'test_bar vol status')
      call assert_int (nboot(1),nbst,'test_bar vol converged')
      call assert_real (dgfep(1),0.8d0,1.0d-6,'test_bar vol fore FEP')
      call assert_real (dgfep(2),0.0d0,1.0d-6,'test_bar vol fore boot')
      call assert_real (dgfep(3),0.8d0,1.0d-6,'test_bar vol back FEP')
      call assert_real (dgfep(4),0.0d0,1.0d-6,'test_bar vol back boot')
      call assert_real (dgbar(1),0.8d0,1.0d-6,'test_bar vol BAR')
      call assert_real (dgbar(3),0.8d0,1.0d-6,'test_bar vol BAR mean')
      call assert_real (dgbar(4),0.0d0,1.0d-6,'test_bar vol BAR boot')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_bar_overlap  --  states without overlap  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_bar_overlap" checks that a pair of states with no overlap
c     ends the BAR iteration with the status for no overlap, while
c     the FEP sums still give the exact energy differences, that mode
c     3 marks such a pair and computes no totals, and that mode 2 stops
c     with a message about the overlap and prints no NaN value
c
c
      subroutine test_bar_overlap
      implicit none
      integer i,n
      integer nbst,istat,ist
      integer nboot(2)
      real*8 temp,value,err
      real*8 ua0(10),ua1(10)
      real*8 ub0(10),ub1(10)
      real*8 vola(10),volb(10)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      logical found
c
c
c     each state lies 2000 Kcal/mol above the other in its own frames
c
      n = 10
      nbst = 50
      temp = 298.0d0
      do i = 1, n
         ua0(i) = 0.0d0
         ua1(i) = 2000.0d0
         ub0(i) = 2000.0d0
         ub1(i) = 0.0d0
         vola(i) = 0.0d0
         volb(i) = 0.0d0
      end do
      call barpair (n,n,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,2,'test_bar overlap status')
      call assert_real (dgfep(1),2000.0d0,1.0d-6,
     &                  'test_bar overlap fore FEP')
      call assert_real (dgfep(3),-2000.0d0,1.0d-6,
     &                  'test_bar overlap back FEP')
c
c     mode 3 marks the pair and mode 2 stops with a message
c
      call barrun3 ('overlap.fep 0 0 0','overlap.out',
     &              'test_bar overlap mode 3')
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','2 overlap.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','the States have no Overlap',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar overlap mode 2')
         call barfind ('out.txt','NaN',found,value,err)
         call assert_logical (found,.false.,'test_bar overlap no NaN')
      end if
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_bar_fep  --  free energy from FEP files  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_bar_fep" runs "bar" in mode 3 on FEP files with a constant
c     energy within each window and a constant energy difference over
c     each pair, so that every estimate is exact and every error is
c     zero; also checks the sample selection, a window without samples
c     and a file header that disagrees with the windows found
c
c
      subroutine test_bar_fep
      implicit none
c
c
      call barrun3 ('analytic.fep 0 0 0','analytic.out',
     &              'test_bar fep analytic')
      call barrun3 ('analytic.fep 2 4 2 500','analytic-sel.out',
     &              'test_bar fep selection')
      call barrun3 ('gap.fep 0 0 0','gap.out',
     &              'test_bar fep missing window')
      call barrun3 ('gap6.fep 0 0 0','gap6.out',
     &              'test_bar fep header mismatch')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_bar_cross  --  mode 2 agrees with mode 3  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_bar_cross" runs the same pair of states through mode 2
c     from a BAR file and through mode 3 from an FEP file, and checks
c     that both give the same BAR iteration value and analytical error,
c     which mode 3 gives for each pair, as its total error comes from
c     the bootstrap
c
c
      subroutine test_bar_cross
      implicit none
      integer i,k,ist,ios
      integer unit,freeunit
      real*8 eps
      real*8 val2,err2
      real*8 val3,err3
      real*8 v(5)
      logical found
      character*240 line
c
c
      eps = 1.0d-4
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out2.txt out3.txt')
      call run_prog ('bar','2 cross.bar 0 0 0 0 0 0','out2.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call run_prog ('bar','3 cross.fep 0 0 0','out3.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call barfind ('out2.txt','Free Energy via BAR Estimate',found,
     &              val2,err2)
      if (.not. found) then
         call barfind ('out2.txt','BAR Free Energy Estimate',found,
     &                 val2,err2)
      end if
      call assert_logical (found,.true.,'test_bar cross mode 2')
      call barfind ('out3.txt','Free Energy via BAR Iteration',found,
     &              val3,err3)
      call assert_logical (found,.true.,'test_bar cross mode 3')
      call assert_real (val3,val2,eps,'test_bar cross BAR value')
c
c     read the analytical error from the row of the only pair
c
      err3 = -1.0d0
      found = .false.
      unit = freeunit ()
      open (unit=unit,file='out3.txt',status='old',iostat=ios)
      do while (ios .eq. 0)
         read (unit,'(a)',iostat=ios)  line
         if (ios .ne. 0)  goto 10
         if (index(line,'Free Energy by Window Pair') .ne. 0) then
            found = .true.
         else if (found) then
            read (line,*,iostat=k)  (v(i),i=1,5)
            if (k .eq. 0) then
               err3 = v(5)
               goto 10
            end if
         end if
      end do
   10 continue
      close (unit=unit)
      call assert_real (err3,err2,eps,'test_bar cross BAR error')
      call execute_command_line ('rm -f out2.txt out3.txt')
      call popdir
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_bar_frames  --  mode 2 frame selection  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_bar_frames" runs mode 2 on a BAR file with ten frames per
c     trajectory, selecting frames 1 to 10 by 3 from trajectory A and
c     2 to 10 by 4 from trajectory B, and checks that the last frame
c     of each range is included, giving four and three frames
c
c
      subroutine test_bar_frames
      implicit none
      integer ist
      integer nfrm(2)
c
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','2 frames.bar 1 10 3 2 10 4','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call barframes ('out.txt',nfrm)
      call assert_int (nfrm(1),4,'test_bar frames trajectory A')
      call assert_int (nfrm(2),3,'test_bar frames trajectory B')
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_bar_pv  --  PV term in the enthalpy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "test_bar_pv" runs mode 2 at a pressure of 1000 atm given on
c     the command line, on an exactly solvable NPT system, an ideal
c     gas of two particles where state 1 adds energy a*V; the volume
c     density is V**2*exp(-beta*P*V) in state 0 and
c     V**2*exp(-beta*(P+a)*V) in state 1, so the free energy change
c     is 3kT*ln((P+a)/P), here 1 Kcal/mol, and the enthalpy change
c     with the PV term is zero, while enthalpies without the PV term,
c     or with it at the wrong pressure, are about 0.76 Kcal/mol; the
c     reported values do not depend on the bootstrap, and are checked
c     against full-sample values found independently for these frames
c
c
      subroutine test_bar_pv
      use units
      implicit none
      integer i,k,m,n,ist
      integer unit,freeunit
      real*8 temp,beta,pres
      real*8 ratio,a,u
      real*8 value,err
      real*8 alpha(3)
      real*8 expect(9)
      real*8 vola(2000),volb(2000)
      logical found
      character*32 label(9)
      data label  / 'Free Energy via Forward FEP',
     &              'Free Energy via Backward FEP',
     &              'Free Energy via BAR Estimate',
     &              'Enthalpy via Direct Estimate',
     &              'Enthalpy via Forward FEP',
     &              'Enthalpy via Backward FEP',
     &              'Enthalpy via BAR Estimate',
     &              'BAR Estimate -T*dS Value',
     &              'PdV Work Term at' /
      data expect  / 0.9975d0, 0.9959d0, 0.9985d0, 0.0d0,
     &               -0.0043d0, 0.0350d0, -0.0063d0, 1.0048d0,
     &               1000.0d0 /
c
c
c     choose a so that the free energy change is 1 Kcal/mol, and
c     draw the volumes from a three-dimensional Kronecker sequence
c
      n = 2000
      temp = 298.0d0
      beta = 1.0d0 / (gasconst*temp)
      pres = 1000.0d0 / prescon
      ratio = exp(beta/3.0d0)
      a = pres * (ratio-1.0d0)
      alpha(1) = sqrt(2.0d0)
      alpha(2) = sqrt(3.0d0)
      alpha(3) = sqrt(5.0d0)
      do i = 1, n
         u = 0.0d0
         do k = 1, 3
            u = u - log(mod(dble(i)*alpha(k),1.0d0))
         end do
         vola(i) = u / (beta*pres)
         volb(i) = u / (beta*(pres+a))
      end do
c
c     write the BAR file with state 0 and state 1 trajectories
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f pv.bar out.txt')
      unit = freeunit ()
      open (unit=unit,file='pv.bar',status='new')
      write (unit,10)  n,temp,0
   10 format (i8,f10.2,2x,'PV Toy State',i2)
      do i = 1, n
         write (unit,20)  i,0.0d0,a*vola(i),vola(i)
   20    format (i8,3f20.8)
      end do
      write (unit,10)  n,temp,1
      do i = 1, n
         write (unit,20)  i,0.0d0,a*volb(i),volb(i)
      end do
      close (unit=unit)
      call run_prog ('bar','2 pv.bar 0 0 0 0 0 0 1000','out.txt',ist)
      if (ist .eq. -1) then
         call execute_command_line ('rm -f pv.bar out.txt')
         call popdir
         return
      end if
c
c     check each value against its full-sample reference
c
      do m = 1, 9
         call barfind ('out.txt',trim(label(m)),found,value,err)
         if (m.eq.3 .and. .not.found) then
            call barfind ('out.txt','BAR Free Energy Estimate',found,
     &                    value,err)
         end if
         if (.not. found)  value = 1.0d10
         call assert_real (value,expect(m),2.0d-4,
     &                     'test_bar pv '//trim(label(m)))
      end do
      call execute_command_line ('rm -f pv.bar out.txt')
      call popdir
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_bar_err  --  bootstrap error estimates  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_bar_err" runs mode 2 on a BAR file whose energies and
c     volumes vary from frame to frame while the energy difference
c     is a constant 0.8 Kcal/mol, with no pressure given so that the
c     default of 1 atm is used; every free energy and enthalpy is
c     then exactly 0.8 in every bootstrap resample, so each error
c     must be zero if the resampled values come from the same frames
c
c
      subroutine test_bar_err
      implicit none
      integer m,ist
      real*8 value,err
      real*8 expect(7)
      logical found,haserr(7)
      character*32 label(7)
      data label  / 'Free Energy via Forward FEP',
     &              'Free Energy via Backward FEP',
     &              'Enthalpy via Forward FEP',
     &              'Enthalpy via Backward FEP',
     &              'Enthalpy via BAR Estimate',
     &              'BAR Estimate -T*dS Value',
     &              'PdV Work Term at' /
      data expect  / 0.8d0, 0.8d0, 0.8d0, 0.8d0, 0.8d0,
     &               0.0d0, 1.0d0 /
      data haserr  / .true., .true., .true., .true., .true.,
     &               .false., .false. /
c
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','2 err.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      do m = 1, 7
         call barfind ('out.txt',trim(label(m)),found,value,err)
         if (.not. found) then
            value = 1.0d10
            err = 1.0d10
         end if
         call assert_real (value,expect(m),1.0d-4,
     &                     'test_bar err '//trim(label(m)))
         if (haserr(m)) then
            call assert_real (err,0.0d0,1.0d-4,
     &                        'test_bar err SD '//trim(label(m)))
         end if
      end do
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_bar_input  --  input that stops mode 2  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "test_bar_input" checks that mode 2 stops with a message when a
c     frame selection is empty, including a first frame past the end
c     of the file, when the two trajectories are at different
c     temperatures, when a header has no temperature, when frame
c     volumes are given for only one trajectory or only some frames,
c     and when a row cannot be read; it also checks input that mode 2
c     accepts, an explicit zero volume, negative frame values, headers
c     with no title and a file that ends early
c
c
      subroutine test_bar_input
      implicit none
      integer i,ist
      integer unit,freeunit
      integer nfrm(2)
      real*8 value,err,u0
      logical found
c
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','2 frames.bar 11 20 1 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','No Frames Selected from Trajectory A',
     &                 found,value,err)
         call assert_logical (found,.true.,'test_bar input no frames')
      end if
      call run_prog ('bar','2 temps.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Require a Single Temperature',
     &                 found,value,err)
         call assert_logical (found,.true.,'test_bar input temperature')
      end if
      call run_prog ('bar','2 frames.bar 1000 0 0 1 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','No Frames Selected from Trajectory A',
     &                 found,value,err)
         call assert_logical (found,.true.,'test_bar input past end A')
      end if
      call run_prog ('bar','2 frames.bar 1 0 0 1000 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','No Frames Selected from Trajectory B',
     &                 found,value,err)
         call assert_logical (found,.true.,'test_bar input past end B')
      end if
      call execute_command_line ('printf "       2\n'//
     &   '       1  -10.0  -9.5\n       2  -10.1  -9.6\n       2\n'//
     &   '       1  -10.5  -10.0\n       2  -10.6  -10.1\n"'//
     &   ' > notemp.bar')
      call run_prog ('bar','2 notemp.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','No Temperature Found for Trajectory',
     &                 found,value,err)
         call assert_logical (found,.true.,'test_bar input no temp')
      end if
c
c     volumes for only one trajectory, or for only some frames of a
c     trajectory, stop mode 2 before any PdV term is found
c
      call execute_command_line ('printf "       2    298.00\n'//
     &   '       1  -10.0  -9.5  1000.0\n'//
     &   '       2  -10.1  -9.6  1010.0\n'//
     &   '       2    298.00\n       1  -10.5  -10.0\n'//
     &   '       2  -10.6  -10.1\n" > onevol.bar')
      call run_prog ('bar','2 onevol.bar 0 0 0 0 0 0 1','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Volumes Missing for       2 of'//
     &                 '       2 Frames of Trajectory B',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar input one volume')
         call barfind ('out.txt','PdV',found,value,err)
         call assert_logical (found,.false.,'test_bar input no PdV')
      end if
      call execute_command_line ('printf "       2    298.00\n'//
     &   '       1  -10.0  -9.5  1000.0\n'//
     &   '       2  -10.1  -9.6\n'//
     &   '       2    298.00\n       1  -10.5  -10.0  1000.0\n'//
     &   '       2  -10.6  -10.1  1010.0\n" > partvol.bar')
      call run_prog ('bar','2 partvol.bar 0 0 0 0 0 0 1','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Volumes Missing for       1 of'//
     &                 '       2 Frames of Trajectory A',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar input part volume')
      end if
c
c     a row that cannot be read stops mode 2, whether an energy is
c     unreadable, the last row is cut short, or a volume is unreadable
c
      call execute_command_line ('printf "       2    298.00\n'//
     &   '       1  -10.0  -9.5\n       2  ******  -9.6\n'//
     &   '       2    298.00\n       1  -10.5  -10.0\n'//
     &   '       2  -10.6  -10.1\n" > badrow1.bar')
      call execute_command_line ('printf "       2    298.00\n'//
     &   '       1  -10.0  -9.5\n       2  -10.1  -9.6\n'//
     &   '       2    298.00\n       1  -10.5  -10.0\n'//
     &   '       2  -10.6\n" > badrow2.bar')
      call execute_command_line ('printf "       2    298.00\n'//
     &   '       1  -10.0  -9.5  1000.0\n'//
     &   '       2  -10.1  -9.6  ******\n'//
     &   '       2    298.00\n       1  -10.5  -10.0  1000.0\n'//
     &   '       2  -10.6  -10.1  1010.0\n" > badrow3.bar')
      call run_prog ('bar','2 badrow1.bar 0 0 0 0 0 0 1','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Unreadable Line in BAR File',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar input bad energy')
      end if
      call run_prog ('bar','2 badrow2.bar 0 0 0 0 0 0 1','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Unreadable Line in BAR File',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar input cut row')
      end if
      call run_prog ('bar','2 badrow3.bar 0 0 0 0 0 0 1','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Unreadable Line in BAR File',found,
     &                 value,err)
         call assert_logical (found,.true.,'test_bar input bad volume')
      end if
c
c     every frame has an energy difference of 0.5 Kcal/mol in both
c     directions, and an explicit zero volume means no volume, so
c     there is no PdV term and no pressure is asked for
c
      unit = freeunit ()
      open (unit=unit,file='zerovol.bar',status='replace')
      write (unit,10)  10
   10 format (i8,'    298.00')
      do i = 1, 10
         u0 = -100.0d0 + 0.37d0*dble(i-1)
         write (unit,20)  i,u0,u0+0.5d0,0.0d0
   20    format (i8,3f18.4)
      end do
      write (unit,10)  10
      do i = 1, 10
         u0 = -90.0d0 + 0.21d0*dble(i-1)
         write (unit,20)  i,u0,u0+0.5d0,0.0d0
      end do
      close (unit=unit)
      call run_prog ('bar','2 zerovol.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barfind ('out.txt','Free Energy via BAR Estimate',found,
     &                 value,err)
         if (.not. found) then
            call barfind ('out.txt','BAR Free Energy Estimate',found,
     &                    value,err)
         end if
         call assert_real (value,0.5d0,1.0d-4,'test_bar input zero vol')
         call assert_real (err,0.0d0,1.0d-6,'test_bar input zero err')
         call barfind ('out.txt','PdV',found,value,err)
         call assert_logical (found,.false.,'test_bar input zero PdV')
         call barfind ('out.txt','Pressure',found,value,err)
         call assert_logical (found,.false.,'test_bar input zero pres')
      end if
c
c     negative frame values select the first and last frames and a
c     step of one, with no prompt
c
      call run_prog ('bar','2 frames.bar -5 -1 -2 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barframes ('out.txt',nfrm)
         call assert_int (nfrm(1),10,'test_bar input negative A')
         call assert_int (nfrm(2),10,'test_bar input negative B')
         call barfind ('out.txt','First & Last Frame',found,value,err)
         call assert_logical (found,.false.,'test_bar input no prompt')
      end if
c
c     headers with no title are read, and a file that ends before the
c     frame count of its last block keeps only the frames it has
c
      call execute_command_line ('printf "       3    298.00\n'//
     &   '       1  -10.0  -9.5\n       2  -10.1  -9.6\n'//
     &   '       3  -10.2  -9.7\n       4    298.00\n'//
     &   '       1  -10.5  -10.0\n       2  -10.6  -10.1\n"'//
     &   ' > short.bar')
      call run_prog ('bar','2 short.bar 0 0 0 0 0 0','out.txt',ist)
      if (ist .ne. -1) then
         call barframes ('out.txt',nfrm)
         call assert_int (nfrm(1),3,'test_bar input short A')
         call assert_int (nfrm(2),2,'test_bar input short B')
         call barfind ('out.txt','Free Energy via BAR Estimate',found,
     &                 value,err)
         if (.not. found) then
            call barfind ('out.txt','BAR Free Energy Estimate',found,
     &                    value,err)
         end if
         call assert_real (value,0.5d0,1.0d-4,'test_bar input short')
      end if
      call execute_command_line ('rm -f notemp.bar onevol.bar '//
     &                           'partvol.bar badrow1.bar '//
     &                           'badrow2.bar badrow3.bar '//
     &                           'zerovol.bar short.bar out.txt')
      call popdir
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_bar_fepvol  --  volumes in mode 3 files  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_bar_fepvol" runs mode 3 on a two-window FEP file holding
c     the same frames and volumes as "test_bar_pv", so the totals must
c     match the mode 2 values found there; the file is then written
c     again without the pressure line, so the volumes are ignored and
c     the enthalpies are those without the PV term; with the pressure
c     line but only zero volumes, the PV terms are dropped with a note
c     and the values match those without the pressure line, and with
c     one volume missing among the others, mode 3 stops with a message
c
c
      subroutine test_bar_fepvol
      use units
      implicit none
      integer i,k,m,n
      integer ist,icase,iexp
      integer unit,freeunit
      real*8 temp,beta,patm
      real*8 pres,ratio,a,u
      real*8 value,err,v
      real*8 alpha(3)
      real*8 expect(6,2)
      real*8 vola(2000),volb(2000)
      logical found
      character*16 tag(4)
      character*32 label(6)
      data label  / 'Free Energy via Forward FEP',
     &              'Free Energy via Backward FEP',
     &              'Free Energy via BAR Iteration',
     &              'Enthalpy via Direct Estimate',
     &              'Enthalpy via BAR Estimate',
     &              'BAR Estimate -T*dS Value' /
      data expect  / 0.9975d0, 0.9959d0, 0.9985d0,
     &               0.0d0, -0.0063d0, 1.0048d0,
     &               0.9975d0, 0.9959d0, 0.9985d0,
     &               0.7639d0, 0.7619d0, 0.2365d0 /
      data tag  / 'pressure', 'no pressure', 'zero volume',
     &            'missing volume' /
c
c
c     build the volumes of the two windows as in "test_bar_pv"
c
      n = 2000
      temp = 298.0d0
      patm = 1000.0d0
      beta = 1.0d0 / (gasconst*temp)
      pres = patm / prescon
      ratio = exp(beta/3.0d0)
      a = pres * (ratio-1.0d0)
      alpha(1) = sqrt(2.0d0)
      alpha(2) = sqrt(3.0d0)
      alpha(3) = sqrt(5.0d0)
      do i = 1, n
         u = 0.0d0
         do k = 1, 3
            u = u - log(mod(dble(i)*alpha(k),1.0d0))
         end do
         vola(i) = u / (beta*pres)
         volb(i) = u / (beta*(pres+a))
      end do
c
c     write the FEP file with and then without the pressure line,
c     then with the pressure line and zero or one missing volume,
c     where window 1 at lambda 1 has energy 0 and window 2 at lambda
c     0 has energy a*V
c
      call pushdir ('file/bar')
      do icase = 1, 4
         call execute_command_line ('rm -f fepvol.fep out.txt')
         unit = freeunit ()
         open (unit=unit,file='fepvol.fep',status='new')
         write (unit,10)  2,temp
   10    format (/,' Free Energy Perturbation Parameters :',
     &           //,' Number of Lambda Windows',9x,i10,
     &           /,' Simulation Temperature',13x,f8.2)
         if (icase .ne. 2)  write (unit,20)  patm
   20    format (' Simulation Pressure',16x,f8.2)
         write (unit,30)
   30    format (/,' Potential Energies at Trial Lambda Values :',
     &           //,7x,'Step',4x,'Sim Lambda',2x,'Trial Lambda',
     &              14x,'Energy',12x,'Volume',/)
         do i = 1, n
            v = vola(i)
            if (icase.eq.3 .or. (icase.eq.4 .and. i.eq.1))  v = 0.0d0
            write (unit,40)  i,1.0d0,1.0d0,0.0d0,v
            write (unit,40)  i,1.0d0,0.0d0,a*vola(i),v
   40       format (i11,2f14.8,f20.8,f18.4)
         end do
         do i = 1, n
            v = volb(i)
            if (icase .eq. 3)  v = 0.0d0
            write (unit,40)  n+i,0.0d0,0.0d0,a*volb(i),v
            write (unit,40)  n+i,0.0d0,1.0d0,0.0d0,v
         end do
         close (unit=unit)
         call run_prog ('bar','3 fepvol.fep 0 0 0','out.txt',ist)
         if (ist .eq. -1)  goto 50
c
c     zero volumes give a note, and a missing volume stops the run
c
         if (icase .eq. 3) then
            call barfind ('out.txt','No Box Volumes',found,value,err)
            call assert_logical (found,.true.,
     &                           'test_bar fepvol zero volume note')
         end if
         if (icase .eq. 4) then
            call barfind ('out.txt','Box Volume Missing at Step',
     &                    found,value,err)
            call assert_logical (found,.true.,
     &                           'test_bar fepvol missing volume stop')
            call assert_real (value,1.0d0,0.0d0,
     &                        'test_bar fepvol missing volume step')
            goto 50
         end if
c
c     check the pressure echo and each total against its reference
c
         call barfind ('out.txt','Simulation Pressure :',found,
     &                 value,err)
         call assert_logical (found,icase.eq.1,
     &                        'test_bar fepvol '//trim(tag(icase))//
     &                        ' echo')
         if (icase .eq. 1) then
            call assert_real (value,patm,1.0d-6,
     &                        'test_bar fepvol pressure value')
         end if
         iexp = min(icase,2)
         do m = 1, 6
            call barfind ('out.txt',trim(label(m)),found,value,err)
            if (.not. found)  value = 1.0d10
            call assert_real (value,expect(m,iexp),2.0d-4,
     &                        'test_bar fepvol '//trim(tag(icase))//
     &                        ' '//trim(label(m)))
         end do
      end do
   50 continue
      call execute_command_line ('rm -f fepvol.fep out.txt')
      call popdir
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_bar_hbar  --  BAR enthalpy, unequal frames  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_bar_hbar" checks the BAR free energy and enthalpy when the
c     trajectories have different numbers of frames; the frames are the
c     exact Boltzmann populations of two configurations, with weights
c     of 3:1 in state 0 over 4 frames and 1:5 in state 1 over 12 frames,
c     so both values must equal the exact ones, while sums over frames
c     that are not averaged per trajectory give an enthalpy of -0.25
c     instead of 0.19 Kcal/mol
c
c
      subroutine test_bar_hbar
      use units
      implicit none
      integer i,k,na,nb
      integer nbst,istat
      integer nboot(2)
      real*8 temp,rt,c
      real*8 dg,dh
      real*8 u0(2),u1(2)
      real*8 ua0(4),ua1(4)
      real*8 ub0(12),ub1(12)
      real*8 vola(4),volb(12)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
c
c
c     set the configuration energies and the exact values
c
      temp = 298.0d0
      rt = gasconst * temp
      c = 0.5d0
      u0(1) = -rt * log(3.0d0)
      u0(2) = 0.0d0
      u1(1) = c
      u1(2) = -rt*log(5.0d0) + c
      dg = -rt*log(6.0d0/4.0d0) + c
      dh = (u1(1)+5.0d0*u1(2))/6.0d0 - (3.0d0*u0(1)+u0(2))/4.0d0
c
c     fill each trajectory with its exact populations
c
      na = 4
      nb = 12
      do i = 1, na
         k = 2
         if (i .le. 3)  k = 1
         ua0(i) = u0(k)
         ua1(i) = u1(k)
         vola(i) = 0.0d0
      end do
      do i = 1, nb
         k = 2
         if (i .le. 2)  k = 1
         ub0(i) = u0(k)
         ub1(i) = u1(k)
         volb(i) = 0.0d0
      end do
      nbst = 200
      call barpair (na,nb,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,0,'test_bar hbar status')
      call assert_real (dgbar(1),dg,2.0d-4,'test_bar hbar BAR')
      call assert_real (hbar(1),dh,2.0d-4,'test_bar hbar enthalpy')
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_bar_bennett  --  BAR with unequal sizes  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_bar_bennett" checks that the BAR free energy solves the
c     equation of Bennett with the optimal shift of ln(nb/na) when the
c     trajectories have 5 and 15 frames; the root is found here by
c     bisection, and differs by 0.1 Kcal/mol from the solution without
c     the shift
c
c
      subroutine test_bar_bennett
      use units
      implicit none
      integer i,k,na,nb
      integer nbst,istat
      integer nboot(2)
      real*8 temp,rt,s,x
      real*8 lnr,glo,ghi,g
      real*8 res,reslo
      real*8 ua0(5),ua1(5)
      real*8 ub0(15),ub1(15)
      real*8 vola(5),volb(15)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
      real*8 xa(5)
      data xa  / -0.9d0, -0.4d0, 0.0d0, 0.3d0, 0.8d0 /
c
c
c     harmonic states displaced from each other, sampled unevenly
c
      na = 5
      nb = 15
      temp = 298.0d0
      rt = gasconst * temp
      s = 1.2d0
      do i = 1, na
         ua0(i) = 0.5d0 * xa(i)**2
         ua1(i) = 0.5d0*(xa(i)-s)**2 + 0.5d0
         vola(i) = 0.0d0
      end do
      do i = 1, nb
         x = s - 0.9d0 + 1.8d0*dble(i-1)/dble(nb-1)
         ub0(i) = 0.5d0 * x**2
         ub1(i) = 0.5d0*(x-s)**2 + 0.5d0
         volb(i) = 0.0d0
      end do
c
c     solve the equation of Bennett by bisection
c
      lnr = log(dble(nb)/dble(na))
      glo = -10.0d0
      ghi = 10.0d0
      do k = 1, 200
         g = 0.5d0 * (glo+ghi)
         res = 0.0d0
         reslo = 0.0d0
         do i = 1, na
            res = res + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-g)/rt-lnr))
            reslo = reslo + 1.0d0/(1.0d0+exp((ua1(i)-ua0(i)-glo)/rt
     &                                      -lnr))
         end do
         do i = 1, nb
            res = res - 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+g)/rt+lnr))
            reslo = reslo - 1.0d0/(1.0d0+exp((ub0(i)-ub1(i)+glo)/rt
     &                                      +lnr))
         end do
         if (reslo*res .le. 0.0d0) then
            ghi = g
         else
            glo = g
         end if
      end do
      g = 0.5d0 * (glo+ghi)
      nbst = 200
      call barpair (na,nb,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,0,'test_bar bennett status')
      call assert_real (dgbar(1),g,2.0d-4,'test_bar bennett BAR')
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_bar_boot  --  bootstrap kernels are exact  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_bar_boot" checks that the bootstrap versions of the FEP,
c     BAR and BAR enthalpy estimators, which use the exponentials found
c     once by "bootprep", match the direct ones on a selection of frames
c     with repeats, at a reference free energy away from the solution,
c     and checks the default number of bootstrap trials
c
c
      subroutine test_bar_boot
      use units
      implicit none
      integer i,na,nb,niter
      integer barnbst
      integer ia(5),ib(15)
      real*8 temp,rt,pres,s,x
      real*8 cref,xmaxa,xmaxb
      real*8 dg1,dh1,u1,v1
      real*8 dg2,dh2,u2,v2
      real*8 c1,c2,sd,h1,h2
      real*8 ua0(5),ua1(5)
      real*8 ub0(15),ub1(15)
      real*8 vola(5),volb(15)
      real*8 wa(5),hwa(5),ea(5)
      real*8 wb(15),hwb(15),eb(15)
      real*8 hist(0:100)
      real*8 xa(5)
      logical conv1,conv2
      data xa  / -0.9d0, -0.4d0, 0.0d0, 0.3d0, 0.8d0 /
      data ia  / 1, 1, 2, 4, 5 /
c
c
c     harmonic states with frame volumes, at a pressure of 1 atm
c
      na = 5
      nb = 15
      temp = 298.0d0
      rt = gasconst * temp
      pres = 1.0d0 / prescon
      s = 1.2d0
      do i = 1, na
         ua0(i) = 0.5d0 * xa(i)**2
         ua1(i) = 0.5d0*(xa(i)-s)**2 + 0.5d0
         vola(i) = 1000.0d0 + 10.0d0*dble(i)
      end do
      do i = 1, nb
         x = s - 0.9d0 + 1.8d0*dble(i-1)/dble(nb-1)
         ub0(i) = 0.5d0 * x**2
         ub1(i) = 0.5d0*(x-s)**2 + 0.5d0
         volb(i) = 1000.0d0 - 5.0d0*dble(i)
         ib(i) = mod(4*i,nb) + 1
      end do
      cref = 0.3d0
      call bootprep (na,nb,ua0,ua1,ub0,ub1,vola,volb,rt,rt,pres,cref,
     &               xmaxa,xmaxb,wa,hwa,wb,hwb,ea,eb)
c
c     compare the forward and backward FEP values
c
      call barfep (na,ia,ua0,ua1,vola,rt,pres,dg1,dh1,u1,v1)
      call bootfep (na,ia,ua0,ua1,vola,wa,hwa,xmaxa,rt,pres,dg2,dh2,
     &              u2,v2)
      call assert_real (dg2,dg1,1.0d-8,'test_bar boot fore FEP')
      call assert_real (dh2,dh1,1.0d-8,'test_bar boot fore enthalpy')
      call assert_real (u2,u1,1.0d-8,'test_bar boot fore energy')
      call assert_real (v2,v1,1.0d-8,'test_bar boot fore volume')
      call barfep (nb,ib,ub1,ub0,volb,rt,pres,dg1,dh1,u1,v1)
      call bootfep (nb,ib,ub1,ub0,volb,wb,hwb,xmaxb,rt,pres,dg2,dh2,
     &              u2,v2)
      call assert_real (dg2,dg1,1.0d-8,'test_bar boot back FEP')
      call assert_real (dh2,dh1,1.0d-8,'test_bar boot back enthalpy')
c
c     compare the BAR free energy and enthalpy
c
      call bariter (na,ia,nb,ib,ua0,ua1,ub0,ub1,rt,rt,0.0d0,100,c1,sd,
     &              niter,hist,conv1)
      call bootiter (na,ia,nb,ib,ea,eb,cref,rt,rt,0.0d0,100,c2,conv2)
      call assert_logical (conv1,.true.,'test_bar boot BAR converged')
      call assert_logical (conv2,.true.,'test_bar boot trial converged')
      call assert_real (c2,c1,1.0d-8,'test_bar boot BAR')
      call barhbar (na,ia,nb,ib,ua0,ua1,ub0,ub1,vola,volb,rt,rt,pres,
     &              c1,h1)
      call boothbar (na,ia,nb,ib,ua0,ua1,ub0,ub1,vola,volb,ea,eb,cref,
     &               rt,rt,pres,c1,h2)
      call assert_real (h2,h1,1.0d-8,'test_bar boot BAR enthalpy')
c
c     the default number of trials is about 1.0d8 frame evaluations,
c     but at least 2000 trials for errors to about 1.5 percent
c
      call assert_int (barnbst(1000),10000,'test_bar boot trials small')
      call assert_int (barnbst(40000),2500,'test_bar boot trials mid')
      call assert_int (barnbst(1000000),2000,'test_bar boot trials big')
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine test_bar_chain  --  mode 3 errors of the totals  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "test_bar_chain" runs mode 3 on five harmonic windows with 100
c     samples each, U = (x-lambda*s)**2/2 + 2*lambda with beta*s**2
c     equal to 9, so the total free energy change is 2 Kcal/mol; the
c     samples come from a normal generator with seed 4, and pairs
c     sharing a window are correlated, so the error of the total from
c     the joint bootstrap is about 1.3 times the sum in quadrature of
c     the errors of the pairs, which would understate it
c
c
      subroutine test_bar_chain
      implicit none
      integer i,k,m,ist,ios
      integer unit,freeunit
      real*8 value,err,quad
      real*8 v(8)
      logical found
      character*240 line
c
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','3 chain.fep 0 0 0','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
c
c     sum the bootstrap errors of the four pairs in quadrature
c
      quad = 0.0d0
      m = 0
      unit = freeunit ()
      open (unit=unit,file='out.txt',status='old',iostat=ios)
      found = .false.
      do while (ios .eq. 0)
         read (unit,'(a)',iostat=ios)  line
         if (ios .ne. 0)  goto 10
         if (index(line,'Free Energy by Window Pair') .ne. 0) then
            found = .true.
         else if (found .and. m.lt.4) then
            read (line,*,iostat=k)  (v(i),i=1,8)
            if (k .eq. 0) then
               m = m + 1
               quad = quad + v(7)**2
            end if
         end if
      end do
   10 continue
      close (unit=unit)
      quad = sqrt(quad)
      call assert_int (m,4,'test_bar chain pairs')
c
c     the total is within twice its error of the exact value, and
c     its error exceeds the sum in quadrature over the pairs
c
      call barfind ('out.txt','Free Energy via BAR Iteration',found,
     &              value,err)
      call assert_logical (found,.true.,'test_bar chain total')
      call assert_logical (abs(value-2.0d0).lt.2.0d0*err,.true.,
     &                     'test_bar chain total value')
      call assert_logical (err.gt.1.1d0*quad,.true.,
     &                     'test_bar chain total error')
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_bar_ref  --  mode 2 on water trajectories  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_bar_ref" runs mode 2 at 1 atm on the two BAR files of the
c     water trajectories that "test_bar_case" builds in mode 1, and
c     checks every value found from all of the frames against those
c     of Tinker9, which keeps the same values in its own test
c
c
      subroutine test_bar_ref
      implicit none
      integer j,m,ist
      real*8 value,err
      real*8 sdbar(2)
      real*8 expect(9,2)
      logical found
      character*18 base(2)
      character*32 label(9)
      character*240 rpath
      data base  / 'water-00-0000-0950', 'water-00-0950-1000' /
      data label  / 'Free Energy via Forward FEP',
     &              'Free Energy via Backward FEP',
     &              'Free Energy via BAR Estimate',
     &              'Atm',
     &              'Enthalpy via Direct Estimate',
     &              'Enthalpy via Forward FEP',
     &              'Enthalpy via Backward FEP',
     &              'Enthalpy via BAR Estimate',
     &              'BAR Estimate -T*dS Value' /
      data expect  / -0.0282d0, -0.2117d0, -0.2072d0, -0.0009d0,
     &               -18.5438d0, -161.2045d0, -348.4886d0,
     &               -342.1090d0, 341.9018d0,
     &               1.1255d0, 1.5060d0, 1.3773d0, 0.0002d0,
     &               -31.5853d0, 268.6882d0, 290.2645d0,
     &               398.8347d0, -397.4574d0 /
      data sdbar  / 0.2593d0, 0.2533d0 /
c
c
      call pushdir ('file/bar')
      do j = 1, 2
         call refpath ('bar',base(j)//'.bar',rpath)
         call execute_command_line ('cp '//trim(rpath)//' ref.bar')
         call run_prog ('bar','2 ref.bar 0 0 0 0 0 0 1','out.txt',ist)
         if (ist .eq. -1)  goto 10
         do m = 1, 9
c
c     the PdV value follows the pressure on its line
c
            if (m .eq. 4) then
               call barfind ('out.txt','Atm       ',found,value,err)
            else
               call barfind ('out.txt',trim(label(m)),found,value,err)
            end if
            if (m.eq.3 .and. .not.found) then
               call barfind ('out.txt','BAR Free Energy Estimate',found,
     &                       value,err)
            end if
            if (.not. found)  value = 1.0d10
            call assert_real (value,expect(m,j),1.0d-4,'test_bar ref '//
     &                        base(j)//' '//trim(label(m)))
            if (m .eq. 3) then
               call assert_real (err,sdbar(j),1.0d-4,'test_bar ref '//
     &                           base(j)//' BAR error')
            end if
         end do
      end do
   10 continue
      call execute_command_line ('rm -f ref.bar out.txt')
      call popdir
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_bar_overflow  --  exponents that overflow  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_bar_overflow" checks a pair of states 1000 Kcal/mol apart
c     in every frame, so that exp(-beta*dU) overflows in the forward
c     FEP sum and underflows in the backward one unless the exponents
c     are shifted; every value must then be -1000 with zero error
c
c
      subroutine test_bar_overflow
      implicit none
      integer i,n
      integer nbst,istat
      integer nboot(2)
      real*8 temp
      real*8 ua0(10),ua1(10)
      real*8 ub0(10),ub1(10)
      real*8 vola(10),volb(10)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
c
c
      n = 10
      nbst = 200
      temp = 298.0d0
      do i = 1, n
         ua0(i) = -100.0d0 + 0.37d0*dble(i-1)
         ua1(i) = ua0(i) - 1000.0d0
         ub1(i) = -1090.0d0 + 0.21d0*dble(i-1)
         ub0(i) = ub1(i) + 1000.0d0
         vola(i) = 0.0d0
         volb(i) = 0.0d0
      end do
      call barpair (n,n,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,0,'test_bar overflow status')
      call assert_real (dgfep(1),-1000.0d0,1.0d-4,
     &                  'test_bar overflow fore FEP')
      call assert_real (dgfep(2),0.0d0,1.0d-4,
     &                  'test_bar overflow fore error')
      call assert_real (dgfep(3),-1000.0d0,1.0d-4,
     &                  'test_bar overflow back FEP')
      call assert_real (dgfep(4),0.0d0,1.0d-4,
     &                  'test_bar overflow back error')
      call assert_real (dgbar(1),-1000.0d0,1.0d-4,
     &                  'test_bar overflow BAR')
      call assert_real (dgbar(3),-1000.0d0,1.0d-4,
     &                  'test_bar overflow BAR mean')
      call assert_real (hbar(1),-1000.0d0,1.0d-4,
     &                  'test_bar overflow BAR enthalpy')
      call assert_real (hbar(2),0.0d0,1.0d-4,
     &                  'test_bar overflow BAR enthalpy error')
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_bar_nonconv  --  failed bootstrap trials  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "test_bar_nonconv" checks a pair where nine of the ten state 0
c     frames have no overlap with state 1, so a bootstrap trial that
c     misses the remaining frame cannot converge; about a third of
c     the trials fail, and they must be skipped rather than loop
c     forever or enter the bootstrap mean
c
c
      subroutine test_bar_nonconv
      implicit none
      integer i,n
      integer nbst,istat
      integer nboot(2)
      real*8 temp
      real*8 ua0(10),ua1(10)
      real*8 ub0(10),ub1(10)
      real*8 vola(10),volb(10)
      real*8 dgfep(4),dgbar(4)
      real*8 sdu(2),uave(4)
      real*8 hbar(3)
c
c
      n = 10
      nbst = 2000
      temp = 298.0d0
      do i = 1, n
         ua0(i) = -100.0d0 - 0.1d0*dble(i-1)
         ua1(i) = ua0(i) + 1000.0d0
         if (i .eq. 1)  ua1(i) = ua0(i)
         ub1(i) = -100.0d0 + 0.1d0*dble(i-1)
         ub0(i) = ub1(i) + 0.3d0*dble(mod(i-1,3)) - 0.3d0
         vola(i) = 0.0d0
         volb(i) = 0.0d0
      end do
      call barpair (n,n,ua0,ua1,ub0,ub1,vola,volb,temp,temp,1.0d0,
     &              nbst,.false.,dgfep,dgbar,nboot,sdu,uave,hbar,
     &              istat)
      call assert_int (istat,0,'test_bar nonconv status')
      call assert_logical (nboot(2).gt.0,.true.,
     &                     'test_bar nonconv some trials fail')
      call assert_logical (nboot(1).ge.2,.true.,
     &                     'test_bar nonconv some trials converge')
      call assert_logical (abs(dgbar(3)).lt.1.0d100,.true.,
     &                     'test_bar nonconv mean finite')
      call assert_logical (abs(dgbar(4)).lt.1.0d100,.true.,
     &                     'test_bar nonconv error finite')
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine barframes  --  frame counts from mode 2 output  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "barframes" reads the number of frames that mode 2 reports for
c     each of the two trajectories, or zero if not found
c
c
      subroutine barframes (fname,nfrm)
      implicit none
      integer k,m,ios
      integer unit,freeunit
      integer nfrm(2)
      character*(*) fname
      character*20 label
      character*240 line
c
c
      label = 'Number of Frames :'
      nfrm(1) = 0
      nfrm(2) = 0
      m = 0
      unit = freeunit ()
      open (unit=unit,file=fname,status='old',iostat=ios)
      do while (ios.eq.0 .and. m.lt.2)
         read (unit,'(a)',iostat=ios)  line
         k = index(line,trim(label))
         if (ios.eq.0 .and. k.gt.0) then
            m = m + 1
            read (line(k+len_trim(label):),*,iostat=ios)  nfrm(m)
         end if
      end do
      close (unit=unit)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine test_bar_logs  --  mode 1 with trajectory logs  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "test_bar_logs" runs mode 1 with energies taken from trajectory
c     logs, on two-frame archives of water with no periodic box, base
c     names that contain a dot and a final blank line; the logs must
c     be found, their energies written to "<A>.bar" with two frames
c     per trajectory and no volume column, a log with no energies must
c     give recomputed ones, and a log with fewer energies than the
c     archive has frames must stop without leaving a BAR file
c
c
      subroutine test_bar_logs
      implicit none
      integer i,k,ist,nlong
      integer unit,freeunit
      integer nfrm(2)
      real*8 u0(4),u1(4)
      real*8 value,err
      logical exist,found
      character*240 args
      character*240 line
c
c
      call pushdir ('file/bar')
      call execute_command_line ('rm -f bar.np-a.* bar.np-b.* out.txt')
      call execute_command_line ('cat ../mutate/water2.xyz '//
     &   '../mutate/water2.xyz > bar.np-a.arc')
      call execute_command_line ('echo >> bar.np-a.arc')
      call execute_command_line ('cp bar.np-a.arc bar.np-b.arc')
      call execute_command_line ('echo "parameters '//
     &   '../commit_6fe8e913/water03" > bar.np-a.key')
      call execute_command_line ('cp bar.np-a.key bar.np-b.key')
      args = '1 bar.np-a.arc 298 bar.np-b.arc 298 Y'
c
c     both logs have one energy per frame
c
      call execute_command_line ('printf " Current Potential'//
     &   '        123.0000\n Current Potential        124.0000\n"'//
     &   ' > bar.np-a.log')
      call execute_command_line ('printf " Current Potential'//
     &   '        223.0000\n Current Potential        224.0000\n"'//
     &   ' > bar.np-b.log')
      call run_prog ('bar',trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call execute_command_line ('rm -f bar.np-a.* bar.np-b.* '//
     &                              'out.txt')
         call popdir
         return
      end if
      do i = 1, 4
         u0(i) = 0.0d0
         u1(i) = 0.0d0
      end do
      unit = freeunit ()
      open (unit=unit,file='bar.np-a.bar',status='old',err=10)
      read (unit,*,err=10,end=10)
      read (unit,*,err=10,end=10)  k,u0(1),u1(1)
      read (unit,*,err=10,end=10)  k,u0(2),u1(2)
      read (unit,*,err=10,end=10)
      read (unit,*,err=10,end=10)  k,u0(3),u1(3)
      read (unit,*,err=10,end=10)  k,u0(4),u1(4)
   10 continue
      close (unit=unit)
      call assert_real (u0(1),123.0d0,1.0d-6,'test_bar logs A frame 1')
      call assert_real (u0(2),124.0d0,1.0d-6,'test_bar logs A frame 2')
      call assert_real (u1(3),223.0d0,1.0d-6,'test_bar logs B frame 1')
      call assert_real (u1(4),224.0d0,1.0d-6,'test_bar logs B frame 2')
c
c     the blank line that ends each archive adds no frame, and the
c     rows of a system with no periodic box have no volume column
c
      do i = 1, 2
         nfrm(i) = 0
      end do
      nlong = 0
      unit = freeunit ()
      open (unit=unit,file='bar.np-a.bar',status='old',err=30)
      do i = 1, 6
         read (unit,'(a)',err=30,end=30)  line
         if (i.eq.1 .or. i.eq.4) then
            read (line,*,err=30,end=30)  nfrm((i+2)/3)
         else if (len_trim(line) .gt. 46) then
            nlong = nlong + 1
         end if
      end do
   30 continue
      close (unit=unit)
      call assert_int (nfrm(1),2,'test_bar logs A frames')
      call assert_int (nfrm(2),2,'test_bar logs B frames')
      call assert_int (nlong,0,'test_bar logs no volume column')
c
c     a log with no energies gives energies recomputed in state 0,
c     the same as those of the same frame in state 1
c
      call execute_command_line ('rm -f bar.np-a.bar')
      call execute_command_line ('echo "no energies" > bar.np-a.log')
      call run_prog ('bar',trim(args),'out.txt',ist)
      u0(1) = 0.0d0
      u1(1) = 1.0d0
      unit = freeunit ()
      open (unit=unit,file='bar.np-a.bar',status='old',err=20)
      read (unit,*,err=20,end=20)
      read (unit,*,err=20,end=20)  k,u0(1),u1(1)
   20 continue
      close (unit=unit)
      call assert_real (u0(1),u1(1),1.0d-4,'test_bar logs recompute')
c
c     a log with fewer energies than frames stops mode 1, and the
c     BAR file already begun for trajectory A is removed
c
      call execute_command_line ('rm -f bar.np-a.bar')
      call execute_command_line ('printf " Current Potential'//
     &   '        223.0000\n" > bar.np-b.log')
      call run_prog ('bar',trim(args),'out.txt',ist)
      call barfind ('out.txt','Energies but Archive has',found,
     &              value,err)
      call assert_logical (found,.true.,'test_bar logs mismatch')
      inquire (file='bar.np-a.bar',exist=exist)
      call assert_logical (exist,.false.,'test_bar logs no file')
      call execute_command_line ('rm -f bar.np-a.* bar.np-b.* out.txt')
      call popdir
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine barrun3  --  one mode 3 output regression  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "barrun3" runs "bar 3" with the given arguments from the fixture
c     directory and compares its output to the named reference
c
c
      subroutine barrun3 (args,refname,label)
      implicit none
      integer ist
      real*8 eps
      character*(*) args,refname,label
      character*240 rpath
c
c
      eps = 1.0d-4
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','3 '//trim(args),'out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
      call refpath ('bar',refname,rpath)
      call assert_files (rpath,'out.txt',eps,label)
      call execute_command_line ('rm -f out.txt')
      call popdir
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine barfind  --  value on a labeled output line  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "barfind" finds the first line of a file containing the label,
c     and returns the first number following the label and the first
c     number following any "+/-" on that line
c
c
      subroutine barfind (fname,label,found,value,err)
      implicit none
      integer k,ios
      integer unit,freeunit
      real*8 value,err
      logical found
      character*(*) fname,label
      character*240 line
c
c
      found = .false.
      value = 0.0d0
      err = 0.0d0
      unit = freeunit ()
      open (unit=unit,file=fname,status='old',iostat=ios)
      if (ios .ne. 0)  return
   10 continue
      read (unit,'(a)',iostat=ios)  line
      if (ios .ne. 0)  goto 20
      k = index(line,label)
      if (k .eq. 0)  goto 10
      found = .true.
      read (line(k+len(label):),*,iostat=ios)  value
      k = index(line,'+/-')
      if (k .gt. 0)  read (line(k+3:),*,iostat=ios)  err
   20 continue
      close (unit=unit)
      return
      end
