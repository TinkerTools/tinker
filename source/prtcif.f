c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2025  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine prtcif  --  output of PDB in PDBx/mmCIF format  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "prtcif" writes out a set of PDB coordinates in PDBx/mmCIF
c     format to an external file
c
c
      subroutine prtcif (icif,imodel)
      use ascii
      use bound
      use boxes
      use files
      use pdb
      use sequen
      use titles
      implicit none
      integer i,k
      integer icif,imodel
      integer start,stop
      integer resmax,resnumb
      integer, allocatable :: resid(:)
      real*8 crdmin,crdmax
      real*8 occupy,bfac
      logical opened
      logical rename
      logical reformat
      character*1 letter
      character*1 chnname,entity
      character*1 insert,formal
      character*1, allocatable :: chain(:)
      character*2 atmc,resc
      character*3 resname
      character*4 atmname
      character*6 crdc
      character*240 fstr
      character*240 ciffile
c
c
c     set flags for residue naming and large value formatting
c
      rename = .false.
      reformat = .true.
c
c     open the output unit if not already done
c
      inquire (unit=icif,opened=opened)
      if (.not. opened) then
         ciffile = filename(1:leng)//'.cif'
         call version (ciffile,'new')
         open (unit=icif,file=ciffile,status='new')
      end if
c
c     write out the structure title as the initial line
c
      if (ltitle .eq. 0) then
         fstr = '(''_struct.title '')'
         write (icif,fstr(1:18))
      else
         letter = char(apostrophe)
         fstr = '(''_struct.title'',10x,a1,a,a1)'
         write (icif,fstr(1:29))  letter,title(1:ltitle),letter
      end if
      fstr = '(''#'')'
      write (icif,fstr(1:5))
c
c     include any lattice parameters in the header info
c
      if (use_bounds) then
         fstr = '(''_cell.length_a'',8x,f9.3)'
         write (icif,fstr(1:26))  xbox
         fstr = '(''_cell.length_b'',8x,f9.3)'
         write (icif,fstr(1:26))  ybox
         fstr = '(''_cell.length_c'',8x,f9.3)'
         write (icif,fstr(1:26))  zbox
         fstr = '(''_cell.angle_alpha'',5x,f9.3)'
         write (icif,fstr(1:29))  alpha
         fstr = '(''_cell.angle_beta'',6x,f9.3)'
         write (icif,fstr(1:28))  beta
         fstr = '(''_cell.angle_gamma'',5x,f9.3)'
         write (icif,fstr(1:29))  gamma
         fstr = '(''#'')'
         write (icif,fstr(1:5))
      end if
c
c     perform dynamic allocation of some local arrays
c
      allocate (resid(maxres))
      allocate (chain(maxres))
c
c     find the chain name and chain position for each residue
c
      do i = 1, nchain
         start = ichain(1,i)
         stop = ichain(2,i)
         do k = start, stop
            resid(k) = k - start + 1
            chain(k) = chnnam(i)
         end do
      end do
c
c     change some Tinker residue names to match PDB standards
c
      if (rename) then
         do i = 1, npdb
            if (pdbres(i) .eq. 'CYX')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'CYD')  pdbres(i) = 'CYS'
            if (pdbres(i) .eq. 'TYD')  pdbres(i) = 'TYR'
            if (pdbres(i) .eq. 'HID')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIE')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'HIP')  pdbres(i) = 'HIS'
            if (pdbres(i) .eq. 'ASH')  pdbres(i) = 'ASP'
            if (pdbres(i) .eq. 'GLH')  pdbres(i) = 'GLU'
            if (pdbres(i) .eq. 'LYD')  pdbres(i) = 'LYS'
         end do
      end if
c
c     set formatting to match the PDB fixed format standard
c
      atmc = 'i3'
      resc = 'i3'
      crdc = '3f7.3 '
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         resmax = 0
         crdmin = 0.0d0
         crdmax = 0.0d0
         do i = 1, npdb
            if (pdbrec(i) .eq. 'ATOM  ') then
               resmax = max(resmax,resid(resnum(i)))
            else
               resmax = max(resmax,resnum(i))
            end if
            crdmin = min(crdmin,xpdb(i),ypdb(i),zpdb(i))
            crdmax = max(crdmax,xpdb(i),ypdb(i),zpdb(i))
         end do
         if (npdb .ge. 1000)  atmc = 'i4'
         if (npdb .ge. 10000)  atmc = 'i5'
         if (npdb .ge. 100000)  atmc = 'i6'
         if (npdb .ge. 1000000)  atmc = 'i7'
         if (resmax .ge. 1000)  resc = 'i4'
         if (resmax .ge. 10000)  resc = 'i5'
         if (resmax .ge. 100000)  resc = 'i6'
         if (resmax .ge. 1000000)  resc = 'i7'
         if (crdmin .le. -10.0d0)  crdc = '3f8.3 '
         if (crdmin .le. 100.0d0)  crdc = '3f8.3 '
         if (crdmin .le. -100.0d0)  crdc = '3f9.3 '
         if (crdmax .ge. 1000.0d0)  crdc = '3f9.3 '
         if (crdmin .le. -1000.0d0)  crdc = '3f10.3'
         if (crdmax .ge. 10000.0d0)  crdc = '3f10.3'
      end if
c
c     write the loop structure for the coordinates section
c
      fstr = '(''loop_'')'
      write (icif,fstr(1:9))
      fstr = '(''_atom_site.group_PDB '')'
      write (icif,fstr(1:25))
      fstr = '(''_atom_site.id '')'
      write (icif,fstr(1:25))
      fstr = '(''_atom_site.type_symbol '')'
      write (icif,fstr(1:27))
      fstr = '(''_atom_site.label_atom_id '')'
      write (icif,fstr(1:29))
      fstr = '(''_atom_site.label_alt_id '')'
      write (icif,fstr(1:28))
      fstr = '(''_atom_site.label_comp_id '')'
      write (icif,fstr(1:29))
      fstr = '(''_atom_site.label_asym_id '')'
      write (icif,fstr(1:29))
      fstr = '(''_atom_site.label_entity_id '')'
      write (icif,fstr(1:31))
      fstr = '(''_atom_site.label_seq_id '')'
      write (icif,fstr(1:28))
      fstr = '(''_atom_site.pdbx_PDB_ins_code '')'
      write (icif,fstr(1:33))
      fstr = '(''_atom_site.Cartn_x '')'
      write (icif,fstr(1:23))
      fstr = '(''_atom_site.Cartn_y '')'
      write (icif,fstr(1:23))
      fstr = '(''_atom_site.Cartn_z '')'
      write (icif,fstr(1:23))
      fstr = '(''_atom_site.occupancy '')'
      write (icif,fstr(1:25))
      fstr = '(''_atom_site.B_iso_or_equiv '')'
      write (icif,fstr(1:30))
      fstr = '(''_atom_site.pdbx_formal_charge '')'
      write (icif,fstr(1:34))
      fstr = '(''_atom_site.auth_seq_id '')'
      write (icif,fstr(1:27))
      fstr = '(''_atom_site.auth_comp_id '')'
      write (icif,fstr(1:28))
      fstr = '(''_atom_site.auth_asym_id '')'
      write (icif,fstr(1:28))
      fstr = '(''_atom_site.auth_atom_id '')'
      write (icif,fstr(1:28))
      fstr = '(''_atom_site.pdbx_PDB_model_num '')'
      write (icif,fstr(1:34))
c
c     write information and coordinates for each PDB atom
c
      fstr = '(a6,'//atmc//',1x,a3,1x,a4,1x,a1,1x,a3,1x,a1,1x,a1,'
     &          //resc//',1x,a1,1x,'//crdc//',2f5.2,1x,a1,'
     &          //resc//',1x,a3,1x,a1,1x,a4,i3)'
      altsym = '.'
      entity = '1'
      insert = '?'
      occupy = 1.0d0
      bfac = 0.0d0
      formal = '?'
      imodel = 1
      do i = 1, npdb
         atmname = pdbatm(i)
         if (atmname(1:1) .eq. ' ')  atmname = atmname(2:4)//' '
         resname = pdbres(i)
         if (resname(2:3) .eq. '  ')  resname = '  '//resname(1:1)
         if (resname(3:3) .eq. ' ')  resname = ' '//resname(1:2)
         if (pdbrec(i) .eq. 'ATOM  ') then
            resnumb = resid(resnum(i))
            chnname = chain(resnum(i))
         else
            resnumb = resnum(i)
            chnname = ' '
         end if
         write (icif,fstr)  pdbrec(i),i,pdbsym(i),atmname,altsym,
     &                      resname,chnname,entity,resnumb,insert,
     &                      xpdb(i),ypdb(i),zpdb(i),occupy,bfac,formal,
     &                      resnumb,resname,chnname,atmname,imodel
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (resid)
      deallocate (chain)
c
c     check for large values requiring extended formatting
c
      if (reformat) then
         if (npdb .ge. 10000)  atmc = 'i6'
         if (npdb .ge. 100000)  atmc = 'i7'
         if (npdb .ge. 1000000)  atmc = 'i8'
      end if
c
c     close the output unit if opened by this routine
c
c     if (.not. opened)  close (unit=icif)
      return
      end
