c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1997  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine final  --  final actions before program exit  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "final" performs any final program actions such as deallocation
c     of global memory, prints a status message, and then pauses if
c     necessary to avoid closing the execution window
c
c
      subroutine final
      use align
      use analyz
      use angang
      use angbnd
      use angtor
      use atmlst
      use bitor
      use bndstr
      use cell
      use cflux
      use charge
      use chgpen
      use chgtrn
      use chunks
      use couple
      use deriv
      use dipole
      use disgeo
      use domega
      use expol
      use faces
      use fft
      use fields
      use fracs
      use freeze
      use group
      use hessn
      use hpmf
      use ielscf
      use improp
      use imptor
      use inform
      use iounit
      use kanang
      use kangs
      use kantor
      use katoms
      use kbonds
      use kcflux
      use kchrge
      use kcpen
      use kctrn
      use kdipol
      use kdsp
      use kexpl
      use khbond
      use kiprop
      use kitors
      use kmulti
      use kopbnd
      use kopdst
      use korbs
      use kpitor
      use kpolpr
      use kpolr
      use krepl
      use ksolut
      use kstbnd
      use ksttor
      use ktorsn
      use ktrtor
      use kurybr
      use kvdwpr
      use kvdws
      use kxrepl
      use light
      use limits
      use merck
      use molcul
      use moldyn
      use mpole
      use mrecip
      use mutant
      use neigh
      use nonpol
      use omega
      use opbend
      use opdist
      use orbits
      use paths
      use pbstuf
      use pdb
      use piorbs
      use pistuf
      use pitors
      use pme
      use polar
      use polgrp
      use polopt
      use polpcg
      use poltcg
      use potfit
      use qmstuf
      use refer
      use repel
      use restrn
      use rgddyn
      use rigid
      use ring
      use rotbnd
      use socket
      use solpot
      use solute
      use stodyn
      use strbnd
      use strtor
      use syntrn
      use tarray
      use tors
      use tortor
      use uprior
      use urey
      use usage
      use vdw
      use vibs
      use warp
      use xrepel
      implicit none
c
c
c     free memory used by the APBS Poisson-Boltzmann solver
c
      if (solvtyp(1:2) .eq. 'PB') then
         call apbsfinal
      end if
c
c     free memory used by the Fourier transform routines
c
      if (use_ewald .or. use_dewald) then
         call fftclose
      end if
c
c     close any open socket used for external communication
c
      if (use_socket) then
         call sktkill
      end if
c
c     print a final status message before exiting Tinker
c
      if (debug) then
         write (iout,10)
   10    format (/,' Tinker is Exiting following Normal Termination')
      end if
c
c     ensure any output is written to the storage device
c
      flush (iout)
c
c     deallocation of global arrays from module align
c
      if (allocated(ifit))  deallocate (ifit)
      if (allocated(wfit))  deallocate (wfit)
c
c     deallocation of global arrays from module analyz
c
      if (allocated(aesum))  deallocate (aesum)
      if (allocated(aeb))  deallocate (aeb)
      if (allocated(aea))  deallocate (aea)
      if (allocated(aeba))  deallocate (aeba)
      if (allocated(aeub))  deallocate (aeub)
      if (allocated(aeaa))  deallocate (aeaa)
      if (allocated(aeopb))  deallocate (aeopb)
      if (allocated(aeopd))  deallocate (aeopd)
      if (allocated(aeid))  deallocate (aeid)
      if (allocated(aeit))  deallocate (aeit)
      if (allocated(aet))  deallocate (aet)
      if (allocated(aept))  deallocate (aept)
      if (allocated(aebt))  deallocate (aebt)
      if (allocated(aeat))  deallocate (aeat)
      if (allocated(aett))  deallocate (aett)
      if (allocated(aev))  deallocate (aev)
      if (allocated(aer))  deallocate (aer)
      if (allocated(aedsp))  deallocate (aedsp)
      if (allocated(aec))  deallocate (aec)
      if (allocated(aecd))  deallocate (aecd)
      if (allocated(aed))  deallocate (aed)
      if (allocated(aem))  deallocate (aem)
      if (allocated(aep))  deallocate (aep)
      if (allocated(aect))  deallocate (aect)
      if (allocated(aerxf))  deallocate (aerxf)
      if (allocated(aes))  deallocate (aes)
      if (allocated(aelf))  deallocate (aelf)
      if (allocated(aeg))  deallocate (aeg)
      if (allocated(aex))  deallocate (aex)
c
c     deallocation of global arrays from module angang
c
      if (allocated(iaa))  deallocate (iaa)
      if (allocated(kaa))  deallocate (kaa)
c
c     deallocation of global arrays from module angbnd
c
      if (allocated(iang))  deallocate (iang)
      if (allocated(ak))  deallocate (ak)
      if (allocated(anat))  deallocate (anat)
      if (allocated(afld))  deallocate (afld)
c
c     deallocation of global arrays from module angtor
c
      if (allocated(iat))  deallocate (iat)
      if (allocated(kant))  deallocate (kant)
c
c     deallocation of global arrays from module atmlst
c
      if (allocated(bndlist))  deallocate (bndlist)
      if (allocated(anglist))  deallocate (anglist)
      if (allocated(balist))  deallocate (balist)
c
c     deallocation of global arrays from module bitor
c
      if (allocated(ibitor))  deallocate (ibitor)
c
c     deallocation of global arrays from module bndstr
c
      if (allocated(ibnd))  deallocate (ibnd)
      if (allocated(bk))  deallocate (bk)
      if (allocated(bl))  deallocate (bl)
c
c     deallocation of global arrays from module cell
c
      if (allocated(icell))  deallocate (icell)
c
c     deallocation of global arrays from module cflux
c
      if (allocated(bflx))  deallocate (bflx)
      if (allocated(aflx))  deallocate (aflx)
      if (allocated(abflx))  deallocate (abflx)
c
c     deallocation of global arrays from module charge
c
      if (allocated(iion))  deallocate (iion)
      if (allocated(jion))  deallocate (jion)
      if (allocated(kion))  deallocate (kion)
      if (allocated(pchg))  deallocate (pchg)
      if (allocated(pchg0))  deallocate (pchg0)
c
c     deallocation of global arrays from module chgpen
c
      if (allocated(pcore))  deallocate (pcore)
      if (allocated(pval))  deallocate (pval)
      if (allocated(pval0))  deallocate (pval0)
      if (allocated(palpha))  deallocate (palpha)
c
c     deallocation of global arrays from module chgtrn
c
      if (allocated(chgct))  deallocate (chgct)
      if (allocated(dmpct))  deallocate (dmpct)
c
c     deallocation of global arrays from module chunks
c
      if (allocated(pmetable))  deallocate (pmetable)
c
c     deallocation of global arrays from module couple
c
      if (allocated(n13))  deallocate (n13)
      if (allocated(n14))  deallocate (n14)
      if (allocated(n15))  deallocate (n15)
      if (allocated(i13))  deallocate (i13)
      if (allocated(i14))  deallocate (i14)
      if (allocated(i15))  deallocate (i15)
c
c     deallocation of global arrays from module deriv
c
      if (allocated(desum))  deallocate (desum)
      if (allocated(deb))  deallocate (deb)
      if (allocated(dea))  deallocate (dea)
      if (allocated(deba))  deallocate (deba)
      if (allocated(deub))  deallocate (deub)
      if (allocated(deaa))  deallocate (deaa)
      if (allocated(deopb))  deallocate (deopb)
      if (allocated(deopd))  deallocate (deopd)
      if (allocated(deid))  deallocate (deid)
      if (allocated(deit))  deallocate (deit)
      if (allocated(det))  deallocate (det)
      if (allocated(dept))  deallocate (dept)
      if (allocated(debt))  deallocate (debt)
      if (allocated(deat))  deallocate (deat)
      if (allocated(dett))  deallocate (dett)
      if (allocated(dev))  deallocate (dev)
      if (allocated(der))  deallocate (der)
      if (allocated(dedsp))  deallocate (dedsp)
      if (allocated(dec))  deallocate (dec)
      if (allocated(decd))  deallocate (decd)
      if (allocated(ded))  deallocate (ded)
      if (allocated(dem))  deallocate (dem)
      if (allocated(dep))  deallocate (dep)
      if (allocated(dect))  deallocate (dect)
      if (allocated(derxf))  deallocate (derxf)
      if (allocated(des))  deallocate (des)
      if (allocated(delf))  deallocate (delf)
      if (allocated(deg))  deallocate (deg)
      if (allocated(dex))  deallocate (dex)
c
c     deallocation of global arrays from module dipole
c
      if (allocated(idpl))  deallocate (idpl)
      if (allocated(bdpl))  deallocate (bdpl)
      if (allocated(sdpl))  deallocate (sdpl)
c
c     deallocation of global arrays from module disgeo
c
      if (allocated(dbnd))  deallocate (dbnd)
      if (allocated(georad))  deallocate (georad)
c
c     deallocation of global arrays from module domega
c
      if (allocated(tesum))  deallocate (tesum)
      if (allocated(teb))  deallocate (teb)
      if (allocated(tea))  deallocate (tea)
      if (allocated(teba))  deallocate (teba)
      if (allocated(teub))  deallocate (teub)
      if (allocated(teaa))  deallocate (teaa)
      if (allocated(teopb))  deallocate (teopb)
      if (allocated(teopd))  deallocate (teopd)
      if (allocated(teid))  deallocate (teid)
      if (allocated(teit))  deallocate (teit)
      if (allocated(tet))  deallocate (tet)
      if (allocated(tept))  deallocate (tept)
      if (allocated(tebt))  deallocate (tebt)
      if (allocated(teat))  deallocate (teat)
      if (allocated(tett))  deallocate (tett)
      if (allocated(tev))  deallocate (tev)
      if (allocated(ter))  deallocate (ter)
      if (allocated(tedsp))  deallocate (tedsp)
      if (allocated(tec))  deallocate (tec)
      if (allocated(tecd))  deallocate (tecd)
      if (allocated(ted))  deallocate (ted)
      if (allocated(tem))  deallocate (tem)
      if (allocated(tep))  deallocate (tep)
      if (allocated(tect))  deallocate (tect)
      if (allocated(terxf))  deallocate (terxf)
      if (allocated(tes))  deallocate (tes)
      if (allocated(telf))  deallocate (telf)
      if (allocated(teg))  deallocate (teg)
      if (allocated(tex))  deallocate (tex)
c
c     deallocation of global arrays from module expol
c
      if (allocated(kpep))  deallocate (kpep)
      if (allocated(prepep))  deallocate (prepep)
      if (allocated(dmppep))  deallocate (dmppep)
      if (allocated(lpep))  deallocate (lpep)
      if (allocated(polscale))  deallocate (polscale)
      if (allocated(polinv))  deallocate (polinv)
c
c     deallocation of global arrays from module faces
c
      if (allocated(ar))  deallocate (ar)
      if (allocated(axyz))  deallocate (axyz)
      if (allocated(skip))  deallocate (skip)
      if (allocated(nosurf))  deallocate (nosurf)
      if (allocated(afree))  deallocate (afree)
      if (allocated(abur))  deallocate (abur)
      if (allocated(cls))  deallocate (cls)
      if (allocated(clst))  deallocate (clst)
      if (allocated(acls))  deallocate (acls)
      if (allocated(ttfe))  deallocate (ttfe)
      if (allocated(ttle))  deallocate (ttle)
      if (allocated(enext))  deallocate (enext)
      if (allocated(tta))  deallocate (tta)
      if (allocated(ttbur))  deallocate (ttbur)
      if (allocated(ttfree))  deallocate (ttfree)
      if (allocated(tfe))  deallocate (tfe)
      if (allocated(ta))  deallocate (ta)
      if (allocated(tr))  deallocate (tr)
      if (allocated(t))  deallocate (t)
      if (allocated(tax))  deallocate (tax)
      if (allocated(tfree))  deallocate (tfree)
      if (allocated(pa))  deallocate (pa)
      if (allocated(p))  deallocate (p)
      if (allocated(va))  deallocate (va)
      if (allocated(vp))  deallocate (vp)
      if (allocated(vxyz))  deallocate (vxyz)
      if (allocated(env))  deallocate (env)
      if (allocated(fnen))  deallocate (fnen)
      if (allocated(ca))  deallocate (ca)
      if (allocated(ct))  deallocate (ct)
      if (allocated(cr))  deallocate (cr)
      if (allocated(c))  deallocate (c)
      if (allocated(eqc))  deallocate (eqc)
      if (allocated(eqv))  deallocate (eqv)
      if (allocated(afe))  deallocate (afe)
      if (allocated(ale))  deallocate (ale)
      if (allocated(eqnext))  deallocate (eqnext)
      if (allocated(fsen))  deallocate (fsen)
      if (allocated(fseq))  deallocate (fseq)
      if (allocated(cyneq))  deallocate (cyneq)
      if (allocated(cyeq))  deallocate (cyeq)
      if (allocated(fqa))  deallocate (fqa)
      if (allocated(fqncy))  deallocate (fqncy)
      if (allocated(fqcy))  deallocate (fqcy)
c
c     deallocation of global arrays from module fft
c
      if (allocated(ffttable))  deallocate (ffttable)
c
c     deallocation of global arrays from module fields
c
      if (allocated(biotyp))  deallocate (biotyp)
c
c     deallocation of global arrays from module fracs
c
      if (allocated(xfrac))  deallocate (xfrac)
      if (allocated(yfrac))  deallocate (yfrac)
      if (allocated(zfrac))  deallocate (zfrac)
c
c     deallocation of global arrays from module freeze
c
      if (allocated(iratx))  deallocate (iratx)
      if (allocated(kratx))  deallocate (kratx)
      if (allocated(irat))  deallocate (irat)
      if (allocated(krat))  deallocate (krat)
      if (allocated(ratimage))  deallocate (ratimage)
c
c     deallocation of global arrays from module group
c
      if (allocated(kgrp))  deallocate (kgrp)
      if (allocated(grplist))  deallocate (grplist)
      if (allocated(igrp))  deallocate (igrp)
      if (allocated(grpmass))  deallocate (grpmass)
      if (allocated(wgrp))  deallocate (wgrp)
c
c     deallocation of global arrays from module hessn
c
      if (allocated(hessx))  deallocate (hessx)
      if (allocated(hessy))  deallocate (hessy)
      if (allocated(hessz))  deallocate (hessz)
c
c     deallocation of global arrays from module hpmf
c
      if (allocated(ipmf))  deallocate (ipmf)
      if (allocated(rpmf))  deallocate (rpmf)
      if (allocated(acsa))  deallocate (acsa)
c
c     deallocation of global arrays from module ielscf
c
      if (allocated(uaux))  deallocate (uaux)
      if (allocated(upaux))  deallocate (upaux)
      if (allocated(vaux))  deallocate (vaux)
      if (allocated(vpaux))  deallocate (vpaux)
      if (allocated(aaux))  deallocate (aaux)
      if (allocated(apaux))  deallocate (apaux)
c
c     deallocation of global arrays from module improp
c
      if (allocated(iiprop))  deallocate (iiprop)
      if (allocated(kprop))  deallocate (kprop)
      if (allocated(vprop))  deallocate (vprop)
c
c     deallocation of global arrays from module imptor
c
      if (allocated(iitors))  deallocate (iitors)
      if (allocated(itors1))  deallocate (itors1)
      if (allocated(itors2))  deallocate (itors2)
      if (allocated(itors3))  deallocate (itors3)
c
c     deallocation of global arrays from module kanang
c
      if (allocated(anan))  deallocate (anan)
c
c     deallocation of global arrays from module kangs
c
      if (allocated(acon))  deallocate (acon)
      if (allocated(acon5))  deallocate (acon5)
      if (allocated(acon4))  deallocate (acon4)
      if (allocated(acon3))  deallocate (acon3)
      if (allocated(aconp))  deallocate (aconp)
      if (allocated(aconf))  deallocate (aconf)
      if (allocated(ang))  deallocate (ang)
      if (allocated(ang5))  deallocate (ang5)
      if (allocated(ang4))  deallocate (ang4)
      if (allocated(ang3))  deallocate (ang3)
      if (allocated(angp))  deallocate (angp)
      if (allocated(angf))  deallocate (angf)
      if (allocated(ka))  deallocate (ka)
      if (allocated(ka5))  deallocate (ka5)
      if (allocated(ka4))  deallocate (ka4)
      if (allocated(ka3))  deallocate (ka3)
      if (allocated(kap))  deallocate (kap)
      if (allocated(kaf))  deallocate (kaf)
c
c     deallocation of global arrays from module kantor
c
      if (allocated(atcon))  deallocate (atcon)
      if (allocated(kat))  deallocate (kat)
c
c     deallocation of global arrays from module katoms
c
      if (allocated(atmcls))  deallocate (atmcls)
      if (allocated(atmnum))  deallocate (atmnum)
      if (allocated(ligand))  deallocate (ligand)
      if (allocated(weight))  deallocate (weight)
      if (allocated(symbol))  deallocate (symbol)
      if (allocated(describe))  deallocate (describe)
c
c     deallocation of global arrays from module kbonds
c
      if (allocated(bcon))  deallocate (bcon)
      if (allocated(bcon5))  deallocate (bcon5)
      if (allocated(bcon4))  deallocate (bcon4)
      if (allocated(bcon3))  deallocate (bcon3)
      if (allocated(blen))  deallocate (blen)
      if (allocated(blen5))  deallocate (blen5)
      if (allocated(blen4))  deallocate (blen4)
      if (allocated(blen3))  deallocate (blen3)
      if (allocated(dlen))  deallocate (dlen)
      if (allocated(kb))  deallocate (kb)
      if (allocated(kb5))  deallocate (kb5)
      if (allocated(kb4))  deallocate (kb4)
      if (allocated(kb3))  deallocate (kb3)
      if (allocated(kel))  deallocate (kel)
c
c     deallocation of global arrays from module kcflux
c
      if (allocated(cflb))  deallocate (cflb)
      if (allocated(cfla))  deallocate (cfla)
      if (allocated(cflab))  deallocate (cflab)
      if (allocated(kcfb))  deallocate (kcfb)
      if (allocated(kcfa))  deallocate (kcfa)
c
c     deallocation of global arrays from module kchrge
c
      if (allocated(chg))  deallocate (chg)
c
c     deallocation of global arrays from module kcpen
c
      if (allocated(cpele))  deallocate (cpele)
      if (allocated(cpalp))  deallocate (cpalp)
c
c     deallocation of global arrays from module kctrn
c
      if (allocated(ctchg))  deallocate (ctchg)
      if (allocated(ctdmp))  deallocate (ctdmp)
c
c     deallocation of global arrays from module kdipol
c
      if (allocated(dpl))  deallocate (dpl)
      if (allocated(dpl5))  deallocate (dpl5)
      if (allocated(dpl4))  deallocate (dpl4)
      if (allocated(dpl3))  deallocate (dpl3)
      if (allocated(pos))  deallocate (pos)
      if (allocated(pos5))  deallocate (pos5)
      if (allocated(pos4))  deallocate (pos4)
      if (allocated(pos3))  deallocate (pos3)
      if (allocated(kd))  deallocate (kd)
      if (allocated(kd5))  deallocate (kd5)
      if (allocated(kd4))  deallocate (kd4)
      if (allocated(kd3))  deallocate (kd3)
c
c     deallocation of global arrays from module kdsp
c
      if (allocated(dspsix))  deallocate (dspsix)
      if (allocated(dspdmp))  deallocate (dspdmp)
c
c     deallocation of global arrays from module kexpl
c
      if (allocated(pepk))  deallocate (pepk)
      if (allocated(peppre))  deallocate (peppre)
      if (allocated(pepdmp))  deallocate (pepdmp)
      if (allocated(pepl))  deallocate (pepl)
c
c     deallocation of global arrays from module khbond
c
      if (allocated(radhb))  deallocate (radhb)
      if (allocated(epshb))  deallocate (epshb)
      if (allocated(khb))  deallocate (khb)
c
c     deallocation of global arrays from module kiprop
c
      if (allocated(dcon))  deallocate (dcon)
      if (allocated(tdi))  deallocate (tdi)
      if (allocated(kdi))  deallocate (kdi)
c
c     deallocation of global arrays from module kitors
c
      if (allocated(ti1))  deallocate (ti1)
      if (allocated(ti2))  deallocate (ti2)
      if (allocated(ti3))  deallocate (ti3)
      if (allocated(kti))  deallocate (kti)
c
c     deallocation of global arrays from module kmulti
c
      if (allocated(multip))  deallocate (multip)
      if (allocated(mpaxis))  deallocate (mpaxis)
      if (allocated(kmp))  deallocate (kmp)
c
c     deallocation of global arrays from module kopbnd
c
      if (allocated(opbn))  deallocate (opbn)
      if (allocated(kopb))  deallocate (kopb)
c
c     deallocation of global arrays from module kopdst
c
      if (allocated(opds))  deallocate (opds)
      if (allocated(kopd))  deallocate (kopd)
c
c     deallocation of global arrays from module korbs
c
      if (allocated(electron))  deallocate (electron)
      if (allocated(ionize))  deallocate (ionize)
      if (allocated(repulse))  deallocate (repulse)
      if (allocated(sslope))  deallocate (sslope)
      if (allocated(sslope5))  deallocate (sslope5)
      if (allocated(sslope4))  deallocate (sslope4)
      if (allocated(tslope))  deallocate (tslope)
      if (allocated(tslope5))  deallocate (tslope5)
      if (allocated(tslope4))  deallocate (tslope4)
      if (allocated(kpi))  deallocate (kpi)
      if (allocated(kpi5))  deallocate (kpi5)
      if (allocated(kpi4))  deallocate (kpi4)
c
c     deallocation of global arrays from module kpitor
c
      if (allocated(ptcon))  deallocate (ptcon)
      if (allocated(kpt))  deallocate (kpt)
c
c     deallocation of global arrays from module kpolpr
c
      if (allocated(thlpr))  deallocate (thlpr)
      if (allocated(thdpr))  deallocate (thdpr)
      if (allocated(kppr))  deallocate (kppr)
c
c     deallocation of global arrays from module kpolr
c
      if (allocated(pgrp))  deallocate (pgrp)
      if (allocated(polr))  deallocate (polr)
      if (allocated(athl))  deallocate (athl)
      if (allocated(dthl))  deallocate (dthl)
c
c     deallocation of global arrays from module krepl
c
      if (allocated(prsiz))  deallocate (prsiz)
      if (allocated(prdmp))  deallocate (prdmp)
      if (allocated(prele))  deallocate (prele)
c
c     deallocation of global arrays from module ksolut
c
      if (allocated(pbr))  deallocate (pbr)
      if (allocated(csr))  deallocate (csr)
      if (allocated(gkr))  deallocate (gkr)
      if (allocated(snk))  deallocate (snk)
c
c     deallocation of global arrays from module kstbnd
c
      if (allocated(stbn))  deallocate (stbn)
      if (allocated(ksb))  deallocate (ksb)
c
c     deallocation of global arrays from module ksttor
c
      if (allocated(btcon))  deallocate (btcon)
      if (allocated(kbt))  deallocate (kbt)
c
c     deallocation of global arrays from module ktorsn
c
      if (allocated(t1))  deallocate (t1)
      if (allocated(t2))  deallocate (t2)
      if (allocated(t3))  deallocate (t3)
      if (allocated(t4))  deallocate (t4)
      if (allocated(t5))  deallocate (t5)
      if (allocated(t6))  deallocate (t6)
      if (allocated(t15))  deallocate (t15)
      if (allocated(t25))  deallocate (t25)
      if (allocated(t35))  deallocate (t35)
      if (allocated(t45))  deallocate (t45)
      if (allocated(t55))  deallocate (t55)
      if (allocated(t65))  deallocate (t65)
      if (allocated(t14))  deallocate (t14)
      if (allocated(t24))  deallocate (t24)
      if (allocated(t34))  deallocate (t34)
      if (allocated(t44))  deallocate (t44)
      if (allocated(t54))  deallocate (t54)
      if (allocated(t64))  deallocate (t64)
      if (allocated(kt))  deallocate (kt)
      if (allocated(kt5))  deallocate (kt5)
      if (allocated(kt4))  deallocate (kt4)
c
c     deallocation of global arrays from module ktrtor
c
      if (allocated(tnx))  deallocate (tnx)
      if (allocated(tny))  deallocate (tny)
      if (allocated(ttx))  deallocate (ttx)
      if (allocated(tty))  deallocate (tty)
      if (allocated(tbf))  deallocate (tbf)
      if (allocated(tbx))  deallocate (tbx)
      if (allocated(tby))  deallocate (tby)
      if (allocated(tbxy))  deallocate (tbxy)
      if (allocated(ktt))  deallocate (ktt)
c
c     deallocation of global arrays from module kurybr
c
      if (allocated(ucon))  deallocate (ucon)
      if (allocated(dst13))  deallocate (dst13)
      if (allocated(ku))  deallocate (ku)
c
c     deallocation of global arrays from module kvdwpr
c
      if (allocated(radpr))  deallocate (radpr)
      if (allocated(epspr))  deallocate (epspr)
      if (allocated(kvpr))  deallocate (kvpr)
c
c     deallocation of global arrays from module kvdws
c
      if (allocated(rad))  deallocate (rad)
      if (allocated(eps))  deallocate (eps)
      if (allocated(rad4))  deallocate (rad4)
      if (allocated(eps4))  deallocate (eps4)
      if (allocated(reduct))  deallocate (reduct)
c
c     deallocation of global arrays from module kxrepl
c
      if (allocated(pxrz))  deallocate (pxrz)
      if (allocated(pxrdmp))  deallocate (pxrdmp)
      if (allocated(pxrcr))  deallocate (pxrcr)
c
c     deallocation of global arrays from module light
c
      if (allocated(kbx))  deallocate (kbx)
      if (allocated(kby))  deallocate (kby)
      if (allocated(kbz))  deallocate (kbz)
      if (allocated(kex))  deallocate (kex)
      if (allocated(key))  deallocate (key)
      if (allocated(kez))  deallocate (kez)
      if (allocated(locx))  deallocate (locx)
      if (allocated(locy))  deallocate (locy)
      if (allocated(locz))  deallocate (locz)
      if (allocated(rgx))  deallocate (rgx)
      if (allocated(rgy))  deallocate (rgy)
      if (allocated(rgz))  deallocate (rgz)
c
c     deallocation of global arrays from module merck
c
      if (allocated(mmff_ka))  deallocate (mmff_ka)
      if (allocated(mmff_ka1))  deallocate (mmff_ka1)
      if (allocated(mmff_ka2))  deallocate (mmff_ka2)
      if (allocated(mmff_ka3))  deallocate (mmff_ka3)
      if (allocated(mmff_ka4))  deallocate (mmff_ka4)
      if (allocated(mmff_ka5))  deallocate (mmff_ka5)
      if (allocated(mmff_ka6))  deallocate (mmff_ka6)
      if (allocated(mmff_ka7))  deallocate (mmff_ka7)
      if (allocated(mmff_ka8))  deallocate (mmff_ka8)
      if (allocated(mmff_ang0))  deallocate (mmff_ang0)
      if (allocated(mmff_ang1))  deallocate (mmff_ang1)
      if (allocated(mmff_ang2))  deallocate (mmff_ang2)
      if (allocated(mmff_ang3))  deallocate (mmff_ang3)
      if (allocated(mmff_ang4))  deallocate (mmff_ang4)
      if (allocated(mmff_ang5))  deallocate (mmff_ang5)
      if (allocated(mmff_ang6))  deallocate (mmff_ang6)
      if (allocated(mmff_ang7))  deallocate (mmff_ang7)
      if (allocated(mmff_ang8))  deallocate (mmff_ang8)
      if (allocated(stbn_abc))  deallocate (stbn_abc)
      if (allocated(stbn_cba))  deallocate (stbn_cba)
      if (allocated(stbn_abc1))  deallocate (stbn_abc1)
      if (allocated(stbn_cba1))  deallocate (stbn_cba1)
      if (allocated(stbn_abc2))  deallocate (stbn_abc2)
      if (allocated(stbn_cba2))  deallocate (stbn_cba2)
      if (allocated(stbn_abc3))  deallocate (stbn_abc3)
      if (allocated(stbn_cba3))  deallocate (stbn_cba3)
      if (allocated(stbn_abc4))  deallocate (stbn_abc4)
      if (allocated(stbn_cba4))  deallocate (stbn_cba4)
      if (allocated(stbn_abc5))  deallocate (stbn_abc5)
      if (allocated(stbn_cba5))  deallocate (stbn_cba5)
      if (allocated(stbn_abc6))  deallocate (stbn_abc6)
      if (allocated(stbn_cba6))  deallocate (stbn_cba6)
      if (allocated(stbn_abc7))  deallocate (stbn_abc7)
      if (allocated(stbn_cba7))  deallocate (stbn_cba7)
      if (allocated(stbn_abc8))  deallocate (stbn_abc8)
      if (allocated(stbn_cba8))  deallocate (stbn_cba8)
      if (allocated(stbn_abc9))  deallocate (stbn_abc9)
      if (allocated(stbn_cba9))  deallocate (stbn_cba9)
      if (allocated(stbn_abc10))  deallocate (stbn_abc10)
      if (allocated(stbn_cba10))  deallocate (stbn_cba10)
      if (allocated(stbn_abc11))  deallocate (stbn_abc11)
      if (allocated(stbn_cba11))  deallocate (stbn_cba11)
c
c     deallocation of global arrays from module molcul
c
      if (allocated(imol))  deallocate (imol)
      if (allocated(kmol))  deallocate (kmol)
      if (allocated(molcule))  deallocate (molcule)
      if (allocated(molmass))  deallocate (molmass)
c
c     deallocation of global arrays from module moldyn
c
      if (allocated(v))  deallocate (v)
      if (allocated(a))  deallocate (a)
      if (allocated(aalt))  deallocate (aalt)
c
c     deallocation of global arrays from module mpole
c
      if (allocated(ipole))  deallocate (ipole)
      if (allocated(polsiz))  deallocate (polsiz)
      if (allocated(pollist))  deallocate (pollist)
      if (allocated(zaxis))  deallocate (zaxis)
      if (allocated(xaxis))  deallocate (xaxis)
      if (allocated(yaxis))  deallocate (yaxis)
      if (allocated(pole))  deallocate (pole)
      if (allocated(rpole))  deallocate (rpole)
      if (allocated(mono0))  deallocate (mono0)
      if (allocated(polaxe))  deallocate (polaxe)
c
c     deallocation of global arrays from module mrecip
c
      if (allocated(cmp))  deallocate (cmp)
      if (allocated(fmp))  deallocate (fmp)
      if (allocated(cphi))  deallocate (cphi)
      if (allocated(fphi))  deallocate (fphi)
c
c     deallocation of global arrays from module mutant
c
      if (allocated(imut))  deallocate (imut)
      if (allocated(type0))  deallocate (type0)
      if (allocated(class0))  deallocate (class0)
      if (allocated(type1))  deallocate (type1)
      if (allocated(class1))  deallocate (class1)
      if (allocated(mut))  deallocate (mut)
c
c     deallocation of global arrays from module neigh
c
      if (allocated(nvlst))  deallocate (nvlst)
      if (allocated(vlst))  deallocate (vlst)
      if (allocated(nelst))  deallocate (nelst)
      if (allocated(elst))  deallocate (elst)
      if (allocated(nulst))  deallocate (nulst)
      if (allocated(ulst))  deallocate (ulst)
      if (allocated(xvold))  deallocate (xvold)
      if (allocated(yvold))  deallocate (yvold)
      if (allocated(zvold))  deallocate (zvold)
      if (allocated(xeold))  deallocate (xeold)
      if (allocated(yeold))  deallocate (yeold)
      if (allocated(zeold))  deallocate (zeold)
      if (allocated(xuold))  deallocate (xuold)
      if (allocated(yuold))  deallocate (yuold)
      if (allocated(zuold))  deallocate (zuold)
c
c     deallocation of global arrays from module nonpol
c
      if (allocated(radcav))  deallocate (radcav)
      if (allocated(raddsp))  deallocate (raddsp)
      if (allocated(epsdsp))  deallocate (epsdsp)
      if (allocated(cdsp))  deallocate (cdsp)
c
c     deallocation of global arrays from module omega
c
      if (allocated(iomega))  deallocate (iomega)
      if (allocated(zline))  deallocate (zline)
      if (allocated(dihed))  deallocate (dihed)
c
c     deallocation of global arrays from module opbend
c
      if (allocated(iopb))  deallocate (iopb)
      if (allocated(opbk))  deallocate (opbk)
c
c     deallocation of global arrays from module opdist
c
      if (allocated(iopd))  deallocate (iopd)
      if (allocated(opdk))  deallocate (opdk)
c
c     deallocation of global arrays from module orbits
c
      if (allocated(qorb))  deallocate (qorb)
      if (allocated(worb))  deallocate (worb)
      if (allocated(emorb))  deallocate (emorb)
c
c     deallocation of global arrays from module paths
c
      if (allocated(pc0))  deallocate (pc0)
      if (allocated(pc1))  deallocate (pc1)
      if (allocated(pvect))  deallocate (pvect)
      if (allocated(pstep))  deallocate (pstep)
      if (allocated(pzet))  deallocate (pzet)
      if (allocated(gc))  deallocate (gc)
c
c     deallocation of global arrays from module pbstuf
c
      if (allocated(apbe))  deallocate (apbe)
      if (allocated(pbep))  deallocate (pbep)
      if (allocated(pbfp))  deallocate (pbfp)
      if (allocated(pbtp))  deallocate (pbtp)
      if (allocated(pbeuind))  deallocate (pbeuind)
      if (allocated(pbeuinp))  deallocate (pbeuinp)
c
c     deallocation of global arrays from module pdb
c
      if (allocated(resnum))  deallocate (resnum)
      if (allocated(resatm))  deallocate (resatm)
      if (allocated(npdb12))  deallocate (npdb12)
      if (allocated(ipdb12))  deallocate (ipdb12)
      if (allocated(pdblist))  deallocate (pdblist)
      if (allocated(xpdb))  deallocate (xpdb)
      if (allocated(ypdb))  deallocate (ypdb)
      if (allocated(zpdb))  deallocate (zpdb)
      if (allocated(pdbres))  deallocate (pdbres)
      if (allocated(pdbsym))  deallocate (pdbsym)
      if (allocated(pdbatm))  deallocate (pdbatm)
      if (allocated(pdbtyp))  deallocate (pdbtyp)
c
c     deallocation of global arrays from module piorbs
c
      if (allocated(iorbit))  deallocate (iorbit)
      if (allocated(iconj))  deallocate (iconj)
      if (allocated(kconj))  deallocate (kconj)
      if (allocated(piperp))  deallocate (piperp)
      if (allocated(ibpi))  deallocate (ibpi)
      if (allocated(itpi))  deallocate (itpi)
      if (allocated(pbpl))  deallocate (pbpl)
      if (allocated(pnpl))  deallocate (pnpl)
      if (allocated(listpi))  deallocate (listpi)
c
c     deallocation of global arrays from module pistuf
c
      if (allocated(bkpi))  deallocate (bkpi)
      if (allocated(blpi))  deallocate (blpi)
      if (allocated(kslope))  deallocate (kslope)
      if (allocated(lslope))  deallocate (lslope)
      if (allocated(torsp2))  deallocate (torsp2)
c
c     deallocation of global arrays from module pitors
c
      if (allocated(ipit))  deallocate (ipit)
      if (allocated(kpit))  deallocate (kpit)
c
c     deallocation of global arrays from module pme
c
      if (allocated(igrid))  deallocate (igrid)
      if (allocated(thetai1))  deallocate (thetai1)
      if (allocated(thetai2))  deallocate (thetai2)
      if (allocated(thetai3))  deallocate (thetai3)
      if (allocated(qgrid))  deallocate (qgrid)
      if (allocated(qfac))  deallocate (qfac)
c
c     deallocation of global arrays from module polar
c
      if (allocated(ipolar))  deallocate (ipolar)
      if (allocated(polarity))  deallocate (polarity)
      if (allocated(thole))  deallocate (thole)
      if (allocated(tholed))  deallocate (tholed)
      if (allocated(pdamp))  deallocate (pdamp)
      if (allocated(thlval))  deallocate (thlval)
      if (allocated(thdval))  deallocate (thdval)
      if (allocated(udir))  deallocate (udir)
      if (allocated(udirp))  deallocate (udirp)
      if (allocated(udirs))  deallocate (udirs)
      if (allocated(udirps))  deallocate (udirps)
      if (allocated(uind))  deallocate (uind)
      if (allocated(uinp))  deallocate (uinp)
      if (allocated(uinds))  deallocate (uinds)
      if (allocated(uinps))  deallocate (uinps)
      if (allocated(uexact))  deallocate (uexact)
      if (allocated(douind))  deallocate (douind)
c
c     deallocation of global arrays from module polgrp
c
      if (allocated(ip11))  deallocate (ip11)
      if (allocated(ip12))  deallocate (ip12)
      if (allocated(ip13))  deallocate (ip13)
      if (allocated(ip14))  deallocate (ip14)
c
c     deallocation of global arrays from module polopt
c
      if (allocated(copt))  deallocate (copt)
      if (allocated(copm))  deallocate (copm)
      if (allocated(uopt))  deallocate (uopt)
      if (allocated(uoptp))  deallocate (uoptp)
      if (allocated(uopts))  deallocate (uopts)
      if (allocated(uoptps))  deallocate (uoptps)
c
c     deallocation of global arrays from module polpcg
c
      if (allocated(mindex))  deallocate (mindex)
      if (allocated(minv))  deallocate (minv)
c
c     deallocation of global arrays from module poltcg
c
      if (allocated(uad))  deallocate (uad)
      if (allocated(uap))  deallocate (uap)
      if (allocated(ubd))  deallocate (ubd)
      if (allocated(ubp))  deallocate (ubp)
c
c     deallocation of global arrays from module potfit
c
      if (allocated(ipgrid))  deallocate (ipgrid)
      if (allocated(fit0))  deallocate (fit0)
      if (allocated(fchg))  deallocate (fchg)
      if (allocated(fpol))  deallocate (fpol)
      if (allocated(pgrid))  deallocate (pgrid)
      if (allocated(epot))  deallocate (epot)
      if (allocated(fitchg))  deallocate (fitchg)
      if (allocated(fitpol))  deallocate (fitpol)
      if (allocated(fitpol))  deallocate (fitcpen)
      if (allocated(gatm))  deallocate (gatm)
      if (allocated(fatm))  deallocate (fatm)
      if (allocated(fxdpl))  deallocate (fxdpl)
      if (allocated(fydpl))  deallocate (fydpl)
      if (allocated(fzdpl))  deallocate (fzdpl)
      if (allocated(vchg))  deallocate (vchg)
      if (allocated(vpol))  deallocate (vpol)
      if (allocated(vcpen))  deallocate (vcpen)
      if (allocated(varpot))  deallocate (varpot)
c
c     deallocation of global arrays from module qmstuf
c
      if (allocated(gx))  deallocate (gx)
      if (allocated(gy))  deallocate (gy)
      if (allocated(gz))  deallocate (gz)
      if (allocated(gfreq))  deallocate (gfreq)
      if (allocated(gforce))  deallocate (gforce)
      if (allocated(gh))  deallocate (gh)
c
c     deallocation of global arrays from module refer
c
      if (allocated(reftyp))  deallocate (reftyp)
      if (allocated(n12ref))  deallocate (n12ref)
      if (allocated(i12ref))  deallocate (i12ref)
      if (allocated(xref))  deallocate (xref)
      if (allocated(yref))  deallocate (yref)
      if (allocated(zref))  deallocate (zref)
      if (allocated(refnam))  deallocate (refnam)
c
c     deallocation of global arrays from module repel
c
      if (allocated(irep))  deallocate (irep)
      if (allocated(replist))  deallocate (replist)
      if (allocated(sizpr))  deallocate (sizpr)
      if (allocated(dmppr))  deallocate (dmppr)
      if (allocated(elepr))  deallocate (elepr)
      if (allocated(repole))  deallocate (repole)
      if (allocated(rrepole))  deallocate (rrepole)
c
c     deallocation of global arrays from module xrepel
c
      if (allocated(ixrep))  deallocate (ixrep)
      if (allocated(xreplist))  deallocate (xreplist)
      if (allocated(zpxr))  deallocate (zpxr)
      if (allocated(dmppxr))  deallocate (dmppxr)
      if (allocated(crpxr))  deallocate (crpxr)
      if (allocated(cpxr))  deallocate (cpxr)
      if (allocated(rcpxr))  deallocate (rcpxr)
      if (allocated(xrepole))  deallocate (xrepole)
c
c     deallocation of global arrays from module restrn
c
      if (allocated(ipfix))  deallocate (ipfix)
      if (allocated(kpfix))  deallocate (kpfix)
      if (allocated(idfix))  deallocate (idfix)
      if (allocated(iafix))  deallocate (iafix)
      if (allocated(itfix))  deallocate (itfix)
      if (allocated(igfix))  deallocate (igfix)
      if (allocated(ichir))  deallocate (ichir)
      if (allocated(xpfix))  deallocate (xpfix)
      if (allocated(ypfix))  deallocate (ypfix)
      if (allocated(zpfix))  deallocate (zpfix)
      if (allocated(pfix))  deallocate (pfix)
      if (allocated(dfix))  deallocate (dfix)
      if (allocated(afix))  deallocate (afix)
      if (allocated(tfix))  deallocate (tfix)
      if (allocated(gfix))  deallocate (gfix)
      if (allocated(chir))  deallocate (chir)
c
c     deallocation of global arrays from module rgddyn
c
      if (allocated(xcmo))  deallocate (xcmo)
      if (allocated(ycmo))  deallocate (ycmo)
      if (allocated(zcmo))  deallocate (zcmo)
      if (allocated(vcm))  deallocate (vcm)
      if (allocated(wcm))  deallocate (wcm)
      if (allocated(lm))  deallocate (lm)
      if (allocated(vc))  deallocate (vc)
      if (allocated(wc))  deallocate (wc)
      if (allocated(linear))  deallocate (linear)
c
c     deallocation of global arrays from module rigid
c
      if (allocated(xrb))  deallocate (xrb)
      if (allocated(yrb))  deallocate (yrb)
      if (allocated(zrb))  deallocate (zrb)
      if (allocated(rbc))  deallocate (rbc)
c
c     deallocation of global arrays from module ring
c
      if (allocated(iring3))  deallocate (iring3)
      if (allocated(iring4))  deallocate (iring4)
      if (allocated(iring5))  deallocate (iring5)
      if (allocated(iring6))  deallocate (iring6)
      if (allocated(iring7))  deallocate (iring7)
c
c     deallocation of global arrays from module rotbnd
c
      if (allocated(rot))  deallocate (rot)
c
c     deallocation of global arrays from module solute
c
      if (allocated(rsolv))  deallocate (rsolv)
      if (allocated(asolv))  deallocate (asolv)
      if (allocated(rborn))  deallocate (rborn)
      if (allocated(drb))  deallocate (drb)
      if (allocated(drbp))  deallocate (drbp)
      if (allocated(drobc))  deallocate (drobc)
      if (allocated(gpol))  deallocate (gpol)
      if (allocated(shct))  deallocate (shct)
      if (allocated(aobc))  deallocate (aobc)
      if (allocated(bobc))  deallocate (bobc)
      if (allocated(gobc))  deallocate (gobc)
      if (allocated(vsolv))  deallocate (vsolv)
      if (allocated(wace))  deallocate (wace)
      if (allocated(s2ace))  deallocate (s2ace)
      if (allocated(uace))  deallocate (uace)
      if (allocated(sneck))  deallocate (sneck)
      if (allocated(bornint))  deallocate (bornint)
c
c     deallocation of global arrays from module stodyn
c
      if (allocated(fgamma))  deallocate (fgamma)
c
c     deallocation of global arrays from module strbnd
c
      if (allocated(isb))  deallocate (isb)
      if (allocated(sbk))  deallocate (sbk)
c
c     deallocation of global arrays from module strtor
c
      if (allocated(ist))  deallocate (ist)
      if (allocated(kst))  deallocate (kst)
c
c     deallocation of global arrays from module syntrn
c
      if (allocated(xmin1))  deallocate (xmin1)
      if (allocated(xmin2))  deallocate (xmin2)
      if (allocated(xm))  deallocate (xm)
c
c     deallocation of global arrays from module tarray
c
      if (allocated(tindex))  deallocate (tindex)
      if (allocated(tdipdip))  deallocate (tdipdip)
c
c     deallocation of global arrays from module tors
c
      if (allocated(itors))  deallocate (itors)
      if (allocated(tors1))  deallocate (tors1)
      if (allocated(tors2))  deallocate (tors2)
      if (allocated(tors3))  deallocate (tors3)
      if (allocated(tors4))  deallocate (tors4)
      if (allocated(tors5))  deallocate (tors5)
      if (allocated(tors6))  deallocate (tors6)
c
c     deallocation of global arrays from module tortor
c
      if (allocated(itt))  deallocate (itt)
c
c     deallocation of global arrays from module uprior
c
      if (allocated(udalt))  deallocate (udalt)
      if (allocated(upalt))  deallocate (upalt)
      if (allocated(usalt))  deallocate (usalt)
      if (allocated(upsalt))  deallocate (upsalt)
c
c     deallocation of global arrays from module urey
c
      if (allocated(iury))  deallocate (iury)
      if (allocated(uk))  deallocate (uk)
      if (allocated(ul))  deallocate (ul)
c
c     deallocation of global arrays from module usage
c
      if (allocated(iuse))  deallocate (iuse)
      if (allocated(use))  deallocate (use)
c
c     deallocation of global arrays from module vdw
c
      if (allocated(ivdw))  deallocate (ivdw)
      if (allocated(jvdw))  deallocate (jvdw)
      if (allocated(mvdw))  deallocate (mvdw)
      if (allocated(ired))  deallocate (ired)
      if (allocated(kred))  deallocate (kred)
      if (allocated(radmin))  deallocate (radmin)
      if (allocated(epsilon))  deallocate (epsilon)
      if (allocated(radmin4))  deallocate (radmin4)
      if (allocated(epsilon4))  deallocate (epsilon4)
      if (allocated(radhbnd))  deallocate (radhbnd)
      if (allocated(epshbnd))  deallocate (epshbnd)
c
c     deallocation of global arrays from module vibs
c
      if (allocated(rho))  deallocate (rho)
      if (allocated(rhok))  deallocate (rhok)
      if (allocated(rwork))  deallocate (rwork)
c
c     deallocation of global arrays from module warp
c
      if (allocated(m2))  deallocate (m2)
c
c     may need a pause to avoid closing the execution window
c
      if (holdup) then
         read (input,20)
   20    format ()
      end if
      return
      end
