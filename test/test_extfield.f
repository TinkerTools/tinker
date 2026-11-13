c
c
c     ##################################################################
c     ##  COPYRIGHT (C) 2026 by  Moses K. J. Chung and Jay W. Ponder  ##
c     ##                     All Rights Reserved                      ##
c     ##################################################################
c
c     ##########################################################
c     ##                                                      ##
c     ##  subroutine test_extfield  --  external field tests  ##
c     ##                                                      ##
c     ##########################################################
c
c
c     "test_extfield" checks the external field cases exercised by
c     the tinker-gpu extfield.cpp tests
c
c
      subroutine test_extfield
      implicit none
c
c
      call test_extfield_mpole
      call test_extfield_polar
      call test_extfield_mpolar
      call test_extfield_vdwpchg
      call test_extfield_dynfreq
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine test_extfield_mpole  --  multipole field  ##
c     ##                                                       ##
c     ###########################################################
c
c
      subroutine test_extfield_mpole
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpole.key','extfield.1.txt',
     &   'extfield_mpole','Atomic Multipoles','amoeba',.true.)
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpole.key','extfield.1.txt',
     &   'extfield_mpole','Atomic Multipoles','amoeba',.false.)
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine test_extfield_polar  --  polarization field  ##
c     ##                                                          ##
c     ##############################################################
c
c
      subroutine test_extfield_polar
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_polar.key','extfield.2.txt',
     &   'extfield_polar','Polarization','amoeba',.true.)
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_polar.key','extfield.2.txt',
     &   'extfield_polar','Polarization','amoeba',.false.)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_extfield_mpolar  --  mpole + polar field  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_extfield_mpolar
      implicit none
c
c
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpolar.key','extfield.3.txt',
     &   'extfield_mpolar',' ','amoeba',.true.)
      call test_extfield_case ('water4_Na_shifted',
     &   'water4_Na_shifted_mpolar.key','extfield.3.txt',
     &   'extfield_mpolar',' ','amoeba',.false.)
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine test_extfield_vdwpchg  --  vdw + charge field  ##
c     ##                                                            ##
c     ################################################################
c
c
      subroutine test_extfield_vdwpchg
      implicit none
c
c
      call test_extfield_case ('amber','amber.key','extfield.4.txt',
     &   'extfield_vdwpchg',' ','amber',.true.)
      call test_extfield_case ('amber','amber.key','extfield.4.txt',
     &   'extfield_vdwpchg',' ','amber',.false.)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine test_extfield_dynfreq  --  oscillating field  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "test_extfield_dynfreq" runs ten Verlet steps of "dynamic" in an
c     external field oscillating at EXFLD-FREQ and checks positions and
c     velocities of the final restart file; the virial is left out, as
c     "exfield1" has an open analytical virial bug
c
c     "dynamic" starts from the accelerations stored in water4.dyn,
c     so they are those of the full field, as tinker-gpu computes them;
c     the result matches the tinker-gpu extfield.5 references
c
c
      subroutine test_extfield_dynfreq
      implicit none
      integer ist
      logical skiptest
      character*240 rpath
      character*(*) tname
      parameter (tname='test_extfield_dynfreq')
c
c
      if (skiptest(tname,'amoeba',4))  return
      call mdprep ('extfield','extfield_dynfreq','water4')
      call pushdir ('file/extfield_dynfreq')
      call tnist_append ('water4.key','external-field  150  -300  450')
      call tnist_append ('water4.key','exfld-freq  20.5')
      call tnist_append ('water4.key','integrator verlet')
      call run_prog ('dynamic','water4 10 2 0.02 1','out.txt',ist)
      if (ist .ne. -1) then
         call assert_int (ist,0,tname//' status')
         call refpath ('extfield','water4_dynfreq.dyn',rpath)
         call assert_dynfile (rpath,'water4.dyn',1.0d-6,1.0d-5,tname)
      end if
      call popdir
      call tnist_clean ('extfield_dynfreq')
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine test_extfield_case  --  run external field  ##
c     ##                                                         ##
c     #############################################################
c
c
      subroutine test_extfield_case (base,key,reffile,tname,term,fftag,
     &                               uselist)
      use action
      use atoms
      use energi
      use extfld
      use units
      use virial
      implicit none
      integer nat,refcnt,i
      real*8 energy,e,ref_e,ref_ei,refeng
      real*8 eps_e,eps_g,eps_v
      real*8, allocatable :: derivs(:,:)
      real*8 refv(3,3)
      real*8, allocatable :: refg(:,:)
      logical skiptest,hasterm,uselist
      character*(*) base,key,reffile,tname,term,fftag
      character*48 cname
      character*240 rpath
      character*(*) tpre
      parameter (tpre='test_')
c
c
      if (uselist) then
         cname = tpre//trim(tname)//'_list'
      else
         cname = tpre//trim(tname)//'_nolist'
      end if
      if (skiptest(cname,fftag,1))  return
      call pushdir ('file/extfield')
      if (uselist) then
         call loadfix_keyadd (base,key,'neighbor-list')
      else
         call loadfix (base,key)
      end if
      allocate (derivs(3,n))
      allocate (refg(3,n))
      call refpath ('extfield',reffile,rpath)
      call load_ref (rpath,n,ref_e,ref_ei,refv,refg,nat)
      hasterm = len_trim(term) .ne. 0
      if (hasterm) then
         call load_engcnt (rpath,term,refeng,refcnt)
         ref_e = refeng
      end if
      eps_e = 1.0d-4
      eps_g = 1.0d-4
      eps_v = 1.0d-3
c
c     level 0  --  external field energy
c
      e = energy ()
      call assert_real (esum,ref_e,eps_e,trim(cname)//' energy (v0)')
      if (trim(term) .eq. 'Atomic Multipoles')
     &   call assert_real (em,refeng,eps_e,
     &                     trim(cname)//' mpole (v0)')
      if (trim(term) .eq. 'Polarization')
     &   call assert_real (ep,refeng,eps_e,
     &                     trim(cname)//' polar (v0)')
c
c     level 1  --  energy, Cartesian gradient and internal virial
c
      call gradient (e,derivs)
      call assert_real (esum,ref_e,eps_e,trim(cname)//' grad-e (v1)')
      call assert_grad (derivs,refg,n,eps_g,trim(cname)//' grad (v1)')
      call assert_grad (vir,refv,3,eps_v,trim(cname)//' virial (v1)')
c
c     level 3  --  total energy and optional per-term analysis
c
      call analysis (e)
      call assert_real (esum,ref_e,eps_e,trim(cname)//' analysis (v3)')
      if (trim(term) .eq. 'Atomic Multipoles')
     &   call check_engcnt (rpath,'Atomic Multipoles',em,nem,
     &                      eps_e,trim(cname)//' mpole (v3)')
      if (trim(term) .eq. 'Polarization')
     &   call check_engcnt (rpath,'Polarization',ep,nep,
     &                      eps_e,trim(cname)//' polar (v3)')
      deallocate (derivs)
      deallocate (refg)
      call popdir
      call final
      return
      end
