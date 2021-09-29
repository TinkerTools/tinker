c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2017  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine ommdata  --  transfer Tinker data to OpenMM  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "ommdata" uses calls to the OpenMM interface to copy data
c     from Tinker modules to the corresponding OpenMM structures
c
c
      subroutine ommdata ()
      use angbnd
      use angpot
      use angtor
      use atomid
      use atoms
      use bath
      use bitor
      use bndpot
      use bndstr
      use bound
      use boxes
      use cell
      use charge
      use chgpot
      use couple
      use deriv
      use energi
      use ewald
      use freeze
      use group
      use imptor
      use inform
      use ktrtor
      use kvdwpr
      use kvdws
      use limits
      use mdstuf
      use molcul
      use moldyn
      use mplpot
      use mpole
      use mutant
      use nonpol
      use opbend
      use openmm
      use pitors
      use pme
      use polar
      use polgrp
      use polopt
      use polpot
      use potent
      use restrn
      use sizes
      use solpot
      use solute
      use stodyn
      use strbnd
      use strtor
      use torpot
      use tors
      use tortor
      use units
      use urey
      use urypot
      use usage
      use vdw
      use vdwpot
      implicit none
c
c
c     use C++ interface calls to map Tinker variables to OpenMM
c
      call set_angbnd_data (nangle,iang,ak,anat,afld)
      call set_angpot_data (angunit,stbnunit,aaunit,opbunit,opdunit,
     &                      cang,qang,pang,sang,copb,qopb,popb,sopb,
     &                      copd,qopd,popd,sopd,opbtyp,angtyp)
      call set_angtor_data (nangtor,iat,kant)
      call set_atomid_data (tag,class,atomic,valence,mass,name,story)
      call set_atoms_data (n,type,x,y,z)
      call set_bath_data (maxnose,voltrial,kelvin,atmsph,tautemp,
     &                    taupres,compress,collide,eta,volmove,vbar,
     &                    qbar,gbar,vnh,qnh,gnh,isothermal,isobaric,
     &                    anisotrop,volscale,barostat,thermostat)
      call set_bitor_data (nbitor,ibitor)
      call set_bndpot_data (cbnd,qbnd,bndunit,bndtyp)
      call set_bndstr_data (nbond,ibnd,bk,bl)
      call set_bound_data (polycut,polycut2,use_bounds,use_replica,
     &                     use_polymer)
      call set_boxes_data (xbox,ybox,zbox,alpha,beta,gamma,xbox2,
     &                     ybox2,zbox2,box23,volbox,alpha_sin,alpha_cos,
     &                     beta_sin,beta_cos,gamma_sin,gamma_cos,
     &                     beta_term,gamma_term,lvec,recip,orthogonal,
     &                     monoclinic,triclinic,octahedron,dodecadron,
     &                     nonprism,spacegrp)
      call set_cell_data (ncell,icell,xcell,ycell,zcell,
     &                    xcell2,ycell2,zcell2)
      call set_charge_data (nion,iion,jion,kion,pchg)
      call set_chgpot_data (electric,dielec,ebuffer,c2scale,c3scale,
     &                      c4scale,c5scale,neutnbr,neutcut)
      call set_couple_data (n12,n13,n14,n15,i12,i13,i14,i15)
      call set_deriv_data (desum,deb,dea,deba,deub,deaa,deopb,deopd,
     &                     deid,deit,det,dept,debt,deat,dett,dev,der,
     &                     dedsp,dec,decd,ded,dem,dep,dect,derxf,
     &                     des,delf,deg,dex)
      call set_energi_data (esum,eb,ea,eba,eub,eaa,eopb,eopd,eid,eit,
     &                      et,ept,ebt,eat,ett,ev,er,edsp,ec,ecd,ed,
     &                      em,ep,ect,erxf,es,elf,eg,ex)
      call set_ewald_data (aewald,aeewald,apewald,adewald,boundary)
      call set_freeze_data (nrat,nratx,iratx,kratx,irat,rateps,
     &                      krat,use_rattle,ratimage)
      call set_group_data (ngrp,kgrp,grplist,igrp,grpmass,wgrp,
     &                     use_group,use_intra,use_inter)
      call set_imptor_data (nitors,iitors,itors1,itors2,itors3)
      call set_inform_data (maxask,digits,iprint,iwrite,isend,
     &                      silent,verbose,debug,holdup,abort)
      call set_ktrtor_data (maxntt,maxtgrd,maxtgrd2,tnx,tny,
     &                      ttx,tty,tbf,tbx,tby,tbxy,ktt)
      call set_kvdwpr_data (maxnvp,radpr,epspr,kvpr)
      call set_kvdws_data (rad,eps,rad4,eps4,reduct)
      call set_limits_data (vdwcut,repcut,dispcut,chgcut,dplcut,
     &                      mpolecut,ctrncut,vdwtaper,reptaper,
     &                      disptaper,chgtaper,dpltaper,mpoletaper,
     &                      ctrntaper,ewaldcut,dewaldcut,usolvcut,
     &                      use_ewald,use_dewald,use_lights,use_list,
     &                      use_vlist,use_dlist,use_clist,use_mlist,
     &                      use_ulist)
      call set_mdstuf_data (nfree,irest,bmnmix,arespa,dorest,integrate)
      call set_molcul_data (nmol,imol,kmol,molcule,totmass,molmass)
      call set_moldyn_data (v,a,aalt)
      call set_mplpot_data (m2scale,m3scale,m4scale,m5scale,use_chgpen)
      call set_mpole_data (maxpole,npole,ipole,polsiz,pollist,
     &                     zaxis,xaxis,yaxis,pole,rpole,spole,
     &                     srpole,polaxe)
      call set_mutant_data (nmut,vcouple,imut,type0,class0,type1,
     &                      class1,lambda,tlambda,vlambda,elambda,
     &                      scexp,scalpha,mut)
      call set_nonpol_data (epso,epsh,rmino,rminh,awater,slevy,
     &                      solvprs,surften,spcut,spoff,stcut,
     &                      stoff,rcav,rdisp,cdisp)
      call set_opbend_data (nopbend,iopb,opbk)
      call set_openmm_data (ommHandle,cudaPrecision,
     &                      ommPlatform,cudaDevice)
      call set_pitors_data (npitors,ipit,kpit)
      call set_pme_data (nfft1,nfft2,nfft3,nefft1,nefft2,nefft3,ndfft1,
     &                   ndfft2,ndfft3,bsorder,bseorder,bsporder,
     7                   bsdorder,igrid,bsmod1,bsmod2,bsmod3,bsbuild,
     &                   thetai1,thetai2,thetai3,qgrid,qfac)
      call set_polar_data (npolar,ipolar,polarity,thole,dirdamp,pdamp,
     &                     udir,udirp,udirs,udirps,uind,uinp,uinds,
     &                     uinps,uexact,douind)
      call set_polgrp_data (maxp11,maxp12,maxp13,maxp14,np11,
     &                      np12,np13,np14,ip11,ip12,ip13,ip14)
      call set_polopt_data (maxopt,optorder,optlevel,copt,copm,
     &                      uopt,uoptp,uopts,uoptps,fopt,foptp)
      call set_polpot_data (politer,poleps,p2scale,p3scale,p4scale,
     &                      p5scale,p2iscale,p3iscale,p4iscale,p5iscale,
     &                      d1scale,d2scale,d3scale,d4scale,u1scale,
     &                      u2scale,u3scale,u4scale,w2scale,w3scale,
     &                      w4scale,w5scale,udiag,polprt,dpequal,
     &                      use_thole,use_dirdamp,poltyp)
      call set_potent_data (use_bond,use_angle,use_strbnd,use_urey,
     &                      use_angang,use_opbend,use_opdist,use_improp,
     &                      use_imptor,use_tors,use_pitors,use_strtor,
     &                      use_angtor,use_tortor,use_vdw,use_repuls,
     &                      use_disp,use_charge,use_chgdpl,use_dipole,
     &                      use_mpole,use_polar,use_chgtrn,use_chgflx,
     &                      use_rxnfld,use_solv,use_metal,use_geom,
     &                      use_extra,use_born,use_orbit)
      call set_restrn_data (npfix,ndfix,nafix,ntfix,ngfix,nchir,ipfix,
     &                      kpfix,idfix,iafix,itfix,igfix,ichir,depth,
     &                      width,rwall,xpfix,ypfix,zpfix,pfix,dfix,
     &                      afix,tfix,gfix,chir,use_basin,use_wall)
      call set_sizes_data (maxatm,maxtyp,maxclass,maxval,maxref,
     &                     maxgrp,maxres,maxfix)
      call set_solpot_data (solvtyp,borntyp)
      call set_solute_data (doffset,p1,p2,p3,p4,p5,rsolv,asolv,rborn,
     &                      drb,drbp,drobc,gpol,shct,aobc,bobc,gobc,
     &                      vsolv,wace,s2ace,uace)
      call set_stodyn_data (friction,fgamma,use_sdarea)
      call set_strbnd_data (nstrbnd,isb,sbk)
      call set_strtor_data (nstrtor,ist,kst)
      call set_torpot_data (idihunit,itorunit,torsunit,ptorunit,
     &                      storunit,atorunit,ttorunit)
      call set_tors_data (ntors,itors,tors1,tors2,tors3,
     &                    tors4,tors5,tors6)
      call set_tortor_data (ntortor,itt)
      call set_units_data (avogadro,lightspd,boltzmann,gasconst,elemchg,
     &                     vacperm,emass,planck,joule,ekcal,bohr,
     &                     hartree,evolt,efreq,coulomb,debye,prescon)
      call set_urey_data (nurey,iury,uk,ul)
      call set_urypot_data (cury,qury,ureyunit)
      call set_usage_data (nuse,iuse,use)
      call set_vdw_data (nvdw,ivdw,jvdw,ired,kred,xred,yred,zred,radmin,
     &                   epsilon,radmin4,epsilon4,radhbnd,epshbnd)
      call set_vdwpot_data (maxgauss,ngauss,igauss,abuck,bbuck,cbuck,
     &                      ghal,dhal,v2scale,v3scale,v4scale,v5scale,
     &                      use_vcorr,vdwindex,radtyp,radsiz,gausstyp,
     &                      radrule,epsrule,vdwtyp)
      return
      end
