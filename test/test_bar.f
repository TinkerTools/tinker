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
c     differences, the input errors that stop mode 2, and the volumes
c     and pressure read in mode 3
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
c     ends the BAR iteration with a failure status while the FEP sums
c     still give the exact energy differences, that mode 3 marks
c     such a pair and computes no totals, and that mode 2 stops with
c     a message rather than iterating forever
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
      call assert_int (istat,1,'test_bar overlap status')
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
         call barfind ('out.txt','not Converged',found,value,err)
         call assert_logical (found,.true.,'test_bar overlap mode 2')
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
c     that both give the same BAR iteration value and analytical error
c
c
      subroutine test_bar_cross
      implicit none
      integer ist
      real*8 eps
      real*8 val2,err2
      real*8 val3,err3
      logical found
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
      integer k,m,ist,ios
      integer unit,freeunit
      integer nfrm(2)
      character*20 label
      character*240 line
c
c
      label = 'Number of Frames :'
      nfrm(1) = 0
      nfrm(2) = 0
      call pushdir ('file/bar')
      call execute_command_line ('rm -f out.txt')
      call run_prog ('bar','2 frames.bar 1 10 3 2 10 4','out.txt',ist)
      if (ist .eq. -1) then
         call popdir
         return
      end if
c
c     read the frame count reported for each trajectory
c
      m = 0
      unit = freeunit ()
      open (unit=unit,file='out.txt',status='old',iostat=ios)
      do while (ios.eq.0 .and. m.lt.2)
         read (unit,'(a)',iostat=ios)  line
         k = index(line,trim(label))
         if (ios.eq.0 .and. k.gt.0) then
            m = m + 1
            read (line(k+len_trim(label):),*,iostat=ios)  nfrm(m)
         end if
      end do
      close (unit=unit)
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
c     frame selection is empty and when the two trajectories are at
c     different temperatures
c
c
      subroutine test_bar_input
      implicit none
      integer ist
      real*8 value,err
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
      call execute_command_line ('rm -f out.txt')
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
c     the enthalpies are those without the PV term
c
c
      subroutine test_bar_fepvol
      use units
      implicit none
      integer i,k,m,n
      integer ist,icase
      integer unit,freeunit
      real*8 temp,beta,patm
      real*8 pres,ratio,a,u
      real*8 value,err
      real*8 alpha(3)
      real*8 expect(6,2)
      real*8 vola(2000),volb(2000)
      logical found
      character*12 tag(2)
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
      data tag  / 'pressure', 'no pressure' /
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
c     where window 1 at lambda 1 has energy 0 and window 2 at lambda
c     0 has energy a*V
c
      call pushdir ('file/bar')
      do icase = 1, 2
         call execute_command_line ('rm -f fepvol.fep out.txt')
         unit = freeunit ()
         open (unit=unit,file='fepvol.fep',status='new')
         write (unit,10)  2,temp
   10    format (/,' Free Energy Perturbation Parameters :',
     &           //,' Number of Lambda Windows',9x,i10,
     &           /,' Simulation Temperature',13x,f8.2)
         if (icase .eq. 1)  write (unit,20)  patm
   20    format (' Simulation Pressure',16x,f8.2)
         write (unit,30)
   30    format (/,' Potential Energies at Trial Lambda Values :',
     &           //,7x,'Step',4x,'Sim Lambda',2x,'Trial Lambda',
     &              14x,'Energy',12x,'Volume',/)
         do i = 1, n
            write (unit,40)  i,1.0d0,1.0d0,0.0d0,vola(i)
            write (unit,40)  i,1.0d0,0.0d0,a*vola(i),vola(i)
   40       format (i11,2f14.8,f20.8,f18.4)
         end do
         do i = 1, n
            write (unit,40)  n+i,0.0d0,0.0d0,a*volb(i),volb(i)
            write (unit,40)  n+i,0.0d0,1.0d0,0.0d0,volb(i)
         end do
         close (unit=unit)
         call run_prog ('bar','3 fepvol.fep 0 0 0','out.txt',ist)
         if (ist .eq. -1)  goto 50
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
         do m = 1, 6
            call barfind ('out.txt',trim(label(m)),found,value,err)
            if (.not. found)  value = 1.0d10
            call assert_real (value,expect(m,icase),2.0d-4,
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
