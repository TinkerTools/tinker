c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2022  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine setprm  --  allocate force field parameters  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "setprm" counts and allocates memory space for force field
c     parameter values involving multiple atom types or classes as
c     found in the parameter file and keyfile
c
c
      subroutine setprm
      use kangs
      use kantor
      use kbonds
      use kcflux
      use kdipol
      use keys
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdwpr
      use params
      use restrn
      implicit none
      integer i,next
      character*20 keyword
      character*240 record
c
c
c     zero out the count of each force field parameter type
c
      maxnb = 0
      maxnb5 = 0
      maxnb4 = 0
      maxnb3 = 0
      maxnel = 0
      maxna = 0
      maxna5 = 0
      maxna4 = 0
      maxna3 = 0
      maxnap = 0
      maxnaf = 0
      maxnsb = 0
      maxnu = 0
      maxnopb = 0
      maxnopd = 0
      maxndi = 0
      maxnti = 0
      maxnt = 0
      maxnt5 = 0
      maxnt4 = 0
      maxnpt = 0
      maxnbt = 0
      maxnat = 0
      maxntt = 0
      maxnvp = 0
      maxnhb = 0
      maxnd = 0
      maxnd5 = 0
      maxnd4 = 0
      maxnd3 = 0
      maxnmp = 0
      maxncfb = 0
      maxncfa = 0
      maxnpi = 0
      maxnpi5 = 0
      maxnpi4 = 0
      maxfix = 0
c
c     find any parameter values found in the parameter file
c
      do i = 1, nprm
         record = prmline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'BOND ')  maxnb = maxnb + 1
         if (keyword(1:6) .eq. 'BOND5 ')  maxnb5 = maxnb5 + 1
         if (keyword(1:6) .eq. 'BOND4 ')  maxnb4 = maxnb4 + 1
         if (keyword(1:6) .eq. 'BOND3 ')  maxnb3 = maxnb3 + 1
         if (keyword(1:9) .eq. 'ELECTNEG ')  maxnel = maxnel + 1
         if (keyword(1:6) .eq. 'ANGLE ')  maxna = maxna + 1
         if (keyword(1:7) .eq. 'ANGLE5 ')  maxna5 = maxna5 + 1
         if (keyword(1:7) .eq. 'ANGLE4 ')  maxna4 = maxna4 + 1
         if (keyword(1:7) .eq. 'ANGLE3 ')  maxna3 = maxna3 + 1
         if (keyword(1:7) .eq. 'ANGLEP ')  maxnap = maxnap + 1
         if (keyword(1:7) .eq. 'ANGLEF ')  maxnaf = maxnaf + 1
         if (keyword(1:7) .eq. 'STRBND ')  maxnsb = maxnsb + 1
         if (keyword(1:9) .eq. 'UREYBRAD ')  maxnu = maxnu + 1
         if (keyword(1:7) .eq. 'OPBEND ')  maxnopb = maxnopb + 1
         if (keyword(1:7) .eq. 'OPDIST ')  maxnopd = maxnopd + 1
         if (keyword(1:9) .eq. 'IMPROPER ')  maxndi = maxndi + 1
         if (keyword(1:8) .eq. 'IMPTORS ')  maxnti = maxnti + 1
         if (keyword(1:8) .eq. 'TORSION ')  maxnt = maxnt + 1
         if (keyword(1:9) .eq. 'TORSION5 ')  maxnt5 = maxnt5 + 1
         if (keyword(1:9) .eq. 'TORSION4 ')  maxnt4 = maxnt4 + 1
         if (keyword(1:7) .eq. 'PITORS ')  maxnpt = maxnpt + 1
         if (keyword(1:8) .eq. 'STRTORS ')  maxnbt = maxnbt + 1
         if (keyword(1:8) .eq. 'ANGTORS ')  maxnat = maxnat + 1
         if (keyword(1:8) .eq. 'TORTORS ')  maxntt = maxntt + 1
         if (keyword(1:6) .eq. 'VDWPR ')  maxnvp = maxnvp + 1
         if (keyword(1:6) .eq. 'HBOND ')  maxnhb = maxnhb + 1
         if (keyword(1:7) .eq. 'DIPOLE ')  maxnd = maxnd + 1
         if (keyword(1:8) .eq. 'DIPOLE5 ')  maxnd5 = maxnd5 + 1
         if (keyword(1:8) .eq. 'DIPOLE4 ')  maxnd4 = maxnd4 + 1
         if (keyword(1:8) .eq. 'DIPOLE3 ')  maxnd3 = maxnd3 + 1
         if (keyword(1:10) .eq. 'MULTIPOLE ')  maxnmp = maxnmp + 1
         if (keyword(1:9) .eq. 'BNDCFLUX ')  maxncfb = maxncfb + 1
         if (keyword(1:9) .eq. 'ANGCFLUX ')  maxncfa = maxncfa + 1
         if (keyword(1:7) .eq. 'PIBOND ')  maxnpi = maxnpi + 1
         if (keyword(1:8) .eq. 'PIBOND5 ')  maxnpi5 = maxnpi5 + 1
         if (keyword(1:8) .eq. 'PIBOND4 ')  maxnpi4 = maxnpi4 + 1
         if (keyword(1:9) .eq. 'RESTRAIN-')  maxfix = maxfix + 1
      end do
c
c     find additional parameter values found in the keyfile
c
      do i = 1, nkey
         record = keyline(i)
         next = 1
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:5) .eq. 'BOND ')  maxnb = maxnb + 1
         if (keyword(1:6) .eq. 'BOND5 ')  maxnb5 = maxnb5 + 1
         if (keyword(1:6) .eq. 'BOND4 ')  maxnb4 = maxnb4 + 1
         if (keyword(1:6) .eq. 'BOND3 ')  maxnb3 = maxnb3 + 1
         if (keyword(1:9) .eq. 'ELECTNEG ')  maxnel = maxnel + 1
         if (keyword(1:6) .eq. 'ANGLE ')  maxna = maxna + 1
         if (keyword(1:7) .eq. 'ANGLE5 ')  maxna5 = maxna5 + 1
         if (keyword(1:7) .eq. 'ANGLE4 ')  maxna4 = maxna4 + 1
         if (keyword(1:7) .eq. 'ANGLE3 ')  maxna3 = maxna3 + 1
         if (keyword(1:7) .eq. 'ANGLEP ')  maxnap = maxnap + 1
         if (keyword(1:7) .eq. 'ANGLEF ')  maxnaf = maxnaf + 1
         if (keyword(1:7) .eq. 'STRBND ')  maxnsb = maxnsb + 1
         if (keyword(1:9) .eq. 'UREYBRAD ')  maxnu = maxnu + 1
         if (keyword(1:7) .eq. 'OPBEND ')  maxnopb = maxnopb + 1
         if (keyword(1:7) .eq. 'OPDIST ')  maxnopd = maxnopd + 1
         if (keyword(1:9) .eq. 'IMPROPER ')  maxndi = maxndi + 1
         if (keyword(1:8) .eq. 'IMPTORS ')  maxnti = maxnti + 1
         if (keyword(1:8) .eq. 'TORSION ')  maxnt = maxnt + 1
         if (keyword(1:9) .eq. 'TORSION5 ')  maxnt5 = maxnt5 + 1
         if (keyword(1:9) .eq. 'TORSION4 ')  maxnt4 = maxnt4 + 1
         if (keyword(1:7) .eq. 'PITORS ')  maxnpt = maxnpt + 1
         if (keyword(1:8) .eq. 'STRTORS ')  maxnbt = maxnbt + 1
         if (keyword(1:8) .eq. 'ANGTORS ')  maxnat = maxnat + 1
         if (keyword(1:8) .eq. 'TORTORS ')  maxntt = maxntt + 1
         if (keyword(1:6) .eq. 'VDWPR ')  maxnvp = maxnvp + 1
         if (keyword(1:6) .eq. 'HBOND ')  maxnhb = maxnhb + 1
         if (keyword(1:7) .eq. 'DIPOLE ')  maxnd = maxnd + 1
         if (keyword(1:8) .eq. 'DIPOLE5 ')  maxnd5 = maxnd5 + 1
         if (keyword(1:8) .eq. 'DIPOLE4 ')  maxnd4 = maxnd4 + 1
         if (keyword(1:8) .eq. 'DIPOLE3 ')  maxnd3 = maxnd3 + 1
         if (keyword(1:10) .eq. 'MULTIPOLE ')  maxnmp = maxnmp + 1
         if (keyword(1:9) .eq. 'BNDCFLUX ')  maxncfb = maxncfb + 1
         if (keyword(1:9) .eq. 'ANGCFLUX ')  maxncfa = maxncfa + 1
         if (keyword(1:7) .eq. 'PIBOND ')  maxnpi = maxnpi + 1
         if (keyword(1:8) .eq. 'PIBOND5 ')  maxnpi5 = maxnpi5 + 1
         if (keyword(1:8) .eq. 'PIBOND4 ')  maxnpi4 = maxnpi4 + 1
         if (keyword(1:9) .eq. 'RESTRAIN-')  maxfix = maxfix + 1
      end do
c
c     set the allocated memory for each parameter type
c
      maxnb = max(500,maxnb+100)
      maxnb5 = max(500,maxnb5+100)
      maxnb4 = max(500,maxnb4+100)
      maxnb3 = max(500,maxnb3+100)
      maxnel = max(500,maxnel+100)
      maxna = max(500,maxna+100)
      maxna5 = max(500,maxna5+100)
      maxna4 = max(500,maxna4+100)
      maxna3 = max(500,maxna3+100)
      maxnap = max(500,maxnap+100)
      maxnaf = max(500,maxnaf+100)
      maxnsb = max(500,maxnsb+100)
      maxnu = max(500,maxnu+100)
      maxnopb = max(500,maxnopb+100)
      maxnopd = max(500,maxnopd+100)
      maxndi = max(500,maxndi+100)
      maxnti = max(500,maxnti+100)
      maxnt = max(500,maxnt+100)
      maxnt5 = max(500,maxnt5+100)
      maxnt4 = max(500,maxnt4+100)
      maxnpt = max(500,maxnpt+100)
      maxnbt = max(500,maxnbt+100)
      maxnat = max(500,maxnat+100)
      maxntt = max(50,maxntt+10)
      maxnvp = max(500,maxnvp+100)
      maxnhb = max(500,maxnhb+100)
      maxnd = max(500,maxnd+100)
      maxnd5 = max(500,maxnd5+100)
      maxnd4 = max(500,maxnd4+100)
      maxnd3 = max(500,maxnd3+100)
      maxnmp = max(500,maxnmp+100)
      maxncfb = max(500,maxncfb+100)
      maxncfa = max(500,maxncfa+100)
      maxnpi = max(500,maxnpi+100)
      maxnpi5 = max(500,maxnpi5+100)
      maxnpi4 = max(500,maxnpi4+100)
      maxfix = max(500,maxfix+100)
c
c     allocate bond stretching forcefield parameters
c
      if (allocated(bcon))  deallocate (bcon)
      allocate (bcon(maxnb))
      if (allocated(bcon5))  deallocate (bcon5)
      allocate (bcon5(maxnb5))
      if (allocated(bcon4))  deallocate (bcon4)
      allocate (bcon4(maxnb4))
      if (allocated(bcon3))  deallocate (bcon3)
      allocate (bcon3(maxnb3))
      if (allocated(blen))  deallocate (blen)
      allocate (blen(maxnb))
      if (allocated(blen5))  deallocate (blen5)
      allocate (blen5(maxnb5))
      if (allocated(blen4))  deallocate (blen4)
      allocate (blen4(maxnb4))
      if (allocated(blen3))  deallocate (blen3)
      allocate (blen3(maxnb3))
      if (allocated(dlen))  deallocate (dlen)
      allocate (dlen(maxnel))
      if (allocated(kb))  deallocate (kb)
      allocate (kb(maxnb))
      if (allocated(kb5))  deallocate (kb5)
      allocate (kb5(maxnb5))
      if (allocated(kb4))  deallocate (kb4)
      allocate (kb4(maxnb4))
      if (allocated(kb3))  deallocate (kb3)
      allocate (kb3(maxnb3))
      if (allocated(kel))  deallocate (kel)
      allocate (kel(maxnel))
c
c     allocate bond angle bend forcefield parameters
c
      if (allocated(acon))  deallocate (acon)
      allocate (acon(maxna))
      if (allocated(acon5))  deallocate (acon5)
      allocate (acon5(maxna5))
      if (allocated(acon4))  deallocate (acon4)
      allocate (acon4(maxna4))
      if (allocated(acon3))  deallocate (acon3)
      allocate (acon3(maxna3))
      if (allocated(aconp))  deallocate (aconp)
      allocate (aconp(maxnap))
      if (allocated(aconf))  deallocate (aconf)
      allocate (aconf(maxnaf))
      if (allocated(ang))  deallocate (ang)
      allocate (ang(3,maxna))
      if (allocated(ang5))  deallocate (ang5)
      allocate (ang5(3,maxna5))
      if (allocated(ang4))  deallocate (ang4)
      allocate (ang4(3,maxna4))
      if (allocated(ang3))  deallocate (ang3)
      allocate (ang3(3,maxna3))
      if (allocated(angp))  deallocate (angp)
      allocate (angp(2,maxnap))
      if (allocated(angf))  deallocate (angf)
      allocate (angf(2,maxnaf))
      if (allocated(ka))  deallocate (ka)
      allocate (ka(maxna))
      if (allocated(ka5))  deallocate (ka5)
      allocate (ka5(maxna5))
      if (allocated(ka4))  deallocate (ka4)
      allocate (ka4(maxna4))
      if (allocated(ka3))  deallocate (ka3)
      allocate (ka3(maxna3))
      if (allocated(kap))  deallocate (kap)
      allocate (kap(maxnap))
      if (allocated(kaf))  deallocate (kaf)
      allocate (kaf(maxnaf))
c
c     allocate stretch-bend forcefield parameters
c
      if (allocated(stbn))  deallocate (stbn)
      allocate (stbn(2,maxnsb))
      if (allocated(ksb))  deallocate (ksb)
      allocate (ksb(maxnsb))
c
c     allocate Urey-Bradley term forcefield parameters
c
      if (allocated(ucon))  deallocate (ucon)
      allocate (ucon(maxnu))
      if (allocated(dst13))  deallocate (dst13)
      allocate (dst13(maxnu))
      if (allocated(ku))  deallocate (ku)
      allocate (ku(maxnu))
c
c     allocate out-of-plane bend forcefield parameters
c
      if (allocated(opbn))  deallocate (opbn)
      allocate (opbn(maxnopb))
      if (allocated(kopb))  deallocate (kopb)
      allocate (kopb(maxnopb))
c
c     allocate out-of-plane distance forcefield parameters
c
      if (allocated(opds))  deallocate (opds)
      allocate (opds(maxnopd))
      if (allocated(kopd))  deallocate (kopd)
      allocate (kopd(maxnopd))
c
c     allocate improper dihedral forcefield parameters
c
      if (allocated(dcon))  deallocate (dcon)
      allocate (dcon(maxndi))
      if (allocated(tdi))  deallocate (tdi)
      allocate (tdi(maxndi))
      if (allocated(kdi))  deallocate (kdi)
      allocate (kdi(maxndi))
c
c     allocate improper torsion forcefield parameters
c
      if (allocated(ti1))  deallocate (ti1)
      allocate (ti1(2,maxnti))
      if (allocated(ti2))  deallocate (ti2)
      allocate (ti2(2,maxnti))
      if (allocated(ti3))  deallocate (ti3)
      allocate (ti3(2,maxnti))
      if (allocated(kti))  deallocate (kti)
      allocate (kti(maxnti))
c
c     allocate torsion angle forcefield parameters
c
      if (allocated(t1))  deallocate (t1)
      allocate (t1(2,maxnt))
      if (allocated(t2))  deallocate (t2)
      allocate (t2(2,maxnt))
      if (allocated(t3))  deallocate (t3)
      allocate (t3(2,maxnt))
      if (allocated(t4))  deallocate (t4)
      allocate (t4(2,maxnt))
      if (allocated(t5))  deallocate (t5)
      allocate (t5(2,maxnt))
      if (allocated(t6))  deallocate (t6)
      allocate (t6(2,maxnt))
      if (allocated(t15))  deallocate (t15)
      allocate (t15(2,maxnt5))
      if (allocated(t25))  deallocate (t25)
      allocate (t25(2,maxnt5))
      if (allocated(t35))  deallocate (t35)
      allocate (t35(2,maxnt5))
      if (allocated(t45))  deallocate (t45)
      allocate (t45(2,maxnt5))
      if (allocated(t55))  deallocate (t55)
      allocate (t55(2,maxnt5))
      if (allocated(t65))  deallocate (t65)
      allocate (t65(2,maxnt5))
      if (allocated(t14))  deallocate (t14)
      allocate (t14(2,maxnt4))
      if (allocated(t24))  deallocate (t24)
      allocate (t24(2,maxnt4))
      if (allocated(t34))  deallocate (t34)
      allocate (t34(2,maxnt4))
      if (allocated(t44))  deallocate (t44)
      allocate (t44(2,maxnt4))
      if (allocated(t54))  deallocate (t54)
      allocate (t54(2,maxnt4))
      if (allocated(t64))  deallocate (t64)
      allocate (t64(2,maxnt4))
      if (allocated(kt))  deallocate (kt)
      allocate (kt(maxnt))
      if (allocated(kt5))  deallocate (kt5)
      allocate (kt5(maxnt5))
      if (allocated(kt4))  deallocate (kt4)
      allocate (kt4(maxnt4))
c
c     allocate pi-system torsion forcefield parameters
c
      if (allocated(ptcon))  deallocate (ptcon)
      allocate (ptcon(maxnpt))
      if (allocated(kpt))  deallocate (kpt)
      allocate (kpt(maxnpt))
c
c     allocate stretch-torsion forcefield parameters
c
      if (allocated(btcon))  deallocate (btcon)
      allocate (btcon(9,maxnbt))
      if (allocated(kbt))  deallocate (kbt)
      allocate (kbt(maxnbt))
c
c     allocate angle-torsion forcefield parameters
c
      if (allocated(atcon))  deallocate (atcon)
      allocate (atcon(6,maxnat))
      if (allocated(kat))  deallocate (kat)
      allocate (kat(maxnat))
c
c     allocate torsion-torsion forcefield parameters
c
      if (allocated(tnx))  deallocate (tnx)
      allocate (tnx(maxntt))
      if (allocated(tny))  deallocate (tny)
      allocate (tny(maxntt))
      if (allocated(ttx))  deallocate (ttx)
      allocate (ttx(maxtgrd,maxntt))
      if (allocated(tty))  deallocate (tty)
      allocate (tty(maxtgrd,maxntt))
      if (allocated(tbf))  deallocate (tbf)
      allocate (tbf(maxtgrd2,maxntt))
      if (allocated(tbx))  deallocate (tbx)
      allocate (tbx(maxtgrd2,maxntt))
      if (allocated(tby))  deallocate (tby)
      allocate (tby(maxtgrd2,maxntt))
      if (allocated(tbxy))  deallocate (tbxy)
      allocate (tbxy(maxtgrd2,maxntt))
      if (allocated(ktt))  deallocate (ktt)
      allocate (ktt(maxntt))
c
c     allocate special vdw term forcefield parameters
c
      if (allocated(radpr))  deallocate (radpr)
      allocate (radpr(maxnvp))
      if (allocated(epspr))  deallocate (epspr)
      allocate (epspr(maxnvp))
      if (allocated(kvpr))  deallocate (kvpr)
      allocate (kvpr(maxnvp))
c
c     allocate H-bonding term forcefield parameters
c
      if (allocated(radhb))  deallocate (radhb)
      allocate (radhb(maxnhb))
      if (allocated(epshb))  deallocate (epshb)
      allocate (epshb(maxnhb))
      if (allocated(khb))  deallocate (khb)
      allocate (khb(maxnhb))
c
c     allocate bond dipole forcefield parameters
c
      if (allocated(dpl))  deallocate (dpl)
      allocate (dpl(maxnd))
      if (allocated(dpl5))  deallocate (dpl5)
      allocate (dpl5(maxnd5))
      if (allocated(dpl4))  deallocate (dpl4)
      allocate (dpl4(maxnd4))
      if (allocated(dpl3))  deallocate (dpl3)
      allocate (dpl3(maxnd3))
      if (allocated(pos))  deallocate (pos)
      allocate (pos(maxnd))
      if (allocated(pos5))  deallocate (pos5)
      allocate (pos5(maxnd5))
      if (allocated(pos4))  deallocate (pos4)
      allocate (pos4(maxnd4))
      if (allocated(pos3))  deallocate (pos3)
      allocate (pos3(maxnd3))
      if (allocated(kd))  deallocate (kd)
      allocate (kd(maxnd))
      if (allocated(kd5))  deallocate (kd5)
      allocate (kd5(maxnd5))
      if (allocated(kd4))  deallocate (kd4)
      allocate (kd4(maxnd4))
      if (allocated(kd3))  deallocate (kd3)
      allocate (kd3(maxnd3))
c
c     allocate atomic multipole forcefield parameters
c
      if (allocated(multip))  deallocate (multip)
      allocate (multip(13,maxnmp))
      if (allocated(mpaxis))  deallocate (mpaxis)
      allocate (mpaxis(maxnmp))
      if (allocated(kmp))  deallocate (kmp)
      allocate (kmp(maxnmp))
c
c     allocate charge flux term forcefield parameters
c
      if (allocated(cflb))  deallocate (cflb)
      allocate (cflb(maxncfb))
      if (allocated(cfla))  deallocate (cfla)
      allocate (cfla(2,maxncfa))
      if (allocated(cflab))  deallocate (cflab)
      allocate (cflab(2,maxncfa))
      if (allocated(kcfb))  deallocate (kcfb)
      allocate (kcfb(maxncfb))
      if (allocated(kcfa))  deallocate (kcfa)
      allocate (kcfa(maxncfa))
c
c     allocate pisystem orbital forcefield parameters
c
      if (allocated(sslope))  deallocate (sslope)
      allocate (sslope(maxnpi))
      if (allocated(sslope5))  deallocate (sslope5)
      allocate (sslope5(maxnpi5))
      if (allocated(sslope4))  deallocate (sslope4)
      allocate (sslope4(maxnpi4))
      if (allocated(tslope))  deallocate (tslope)
      allocate (tslope(maxnpi))
      if (allocated(tslope5))  deallocate (tslope5)
      allocate (tslope5(maxnpi5))
      if (allocated(tslope4))  deallocate (tslope4)
      allocate (tslope4(maxnpi4))
      if (allocated(kpi))  deallocate (kpi)
      allocate (kpi(maxnpi))
      if (allocated(kpi5))  deallocate (kpi5)
      allocate (kpi5(maxnpi5))
      if (allocated(kpi4))  deallocate (kpi4)
      allocate (kpi4(maxnpi4))
      return
      end
