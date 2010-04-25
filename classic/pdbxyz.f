c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program pdbxyz  --  Protein Data Bank to XYZ coordinates  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "pdbxyz" takes as input a Protein Data Bank file and then
c     converts to and writes out a Cartesian coordinates file and,
c     for biopolymers, a sequence file
c
c
      program pdbxyz
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'files.i'
      include 'inform.i'
      include 'katoms.i'
      include 'pdb.i'
      include 'resdue.i'
      integer i,j,it,next
      integer ipdb,ixyz,iseq
      integer last,freeunit
      integer size(maxatm)
      real*8 xi,yi,zi,rij
      real*8 rcut,rmax(0:9)
      logical peptide
      logical nucacid
      logical clash
      character*1 letter
      character*3 resname
      character*3 reslast
      character*120 pdbfile
      character*120 xyzfile
      character*120 seqfile
c
c
c     get the Protein Data Bank file and a parameter set
c
      call initial
      call getpdb
      call field
c
c     decide whether the system contains only polypeptides
c
      peptide = .false.
      reslast = '***'
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resname = resnam(i)
            if (resname .ne. reslast) then
               reslast = resname
               do j = 1, maxamino
                  if (resname .eq. amino(j)) then
                     peptide = .true.
                     goto 10
                  end if
               end do
               peptide = .false.
               goto 20
   10          continue
            end if
         end if
      end do
   20 continue
c
c     decide whether the system contains only nucleic acids
c
      nucacid = .false.
      reslast = '***'
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            resname = resnam(i)
            if (resname .ne. reslast) then
               reslast = resname
               do j = 1, maxnuc
                  if (resname .eq. nuclz(j)) then
                     nucacid = .true.
                     goto 30
                  end if
               end do
               nucacid = .false.
               goto 40
   30          continue
            end if
         end if
      end do
   40 continue
c
c     open the TINKER coordinates file to be used for output
c
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
c
c     reopen the PDB file and read the first coordinate set
c
      ipdb = freeunit ()
      pdbfile = filename
      call suffix (pdbfile,'pdb')
      call version (pdbfile,'old')
      open (unit=ipdb,file=pdbfile,status ='old')
      rewind (unit=ipdb)
      call readpdb (ipdb)
c
c     use special translation mechanisms used for biopolymers
c
      dowhile (.not. abort)
         if (peptide .or. nucacid) then
            if (peptide)  call ribosome
            if (nucacid)  call ligase
            call hetatom
            last = n
            do i = last, 1, -1
               if (type(i) .eq. 0)  call delete (i)
            end do
c
c     get general atom properties for distance-based connectivity
c
         else
            n = npdb
            do i = 1, n
               x(i) = xpdb(i)
               y(i) = ypdb(i)
               z(i) = zpdb(i)
               name(i) = atmnam(i)(2:4)
               n12(i) = 0
               next = 1
               call getnumb (resnam(i),type(i),next)
               it = type(i)
               if (it .eq. 0) then
                  letter = name(i)(1:1)
                  call upcase (letter)
                  if (letter .eq. 'H') then
                     size(i) = 1
                  else if (letter .eq. 'C') then
                     size(i) = 2
                  else if (letter .eq. 'N') then
                     size(i) = 2
                  else if (letter .eq. 'O') then
                     size(i) = 2
                  else if (letter .eq. 'P') then
                     size(i) = 3
                  else if (letter .eq. 'S') then
                     size(i) = 3
                  else
                     size(i) = 0
                  end if
               else if (ligand(it) .eq. 0) then
                  size(i) = 0
               else if (atmnum(it) .le. 2) then
                  size(i) = 1
               else if (atmnum(it) .le. 10) then
                  size(i) = 2
               else
                  size(i) = 3
               end if
            end do
c
c     set the maximum bonded distance between atom type pairs
c
            rmax(0) = -1.0d0
            rmax(1) = -1.0d0
            rmax(2) = 1.3d0
            rmax(3) = 1.55d0
            rmax(4) = 1.75d0
            rmax(6) = 2.0d0
            rmax(9) = 2.2d0
c
c     find and connect atom pairs within bonding distance
c
            do i = 1, n-1
               xi = x(i)
               yi = y(i)
               zi = z(i)
               do j = i+1, n
                  rcut = rmax(size(i)*size(j))
                  rij = sqrt((xi-x(j))**2 + (yi-y(j))**2
     &                           + (zi-z(j))**2)
                  if (rij .le. rcut) then
                     n12(i) = n12(i) + 1
                     i12(n12(i),i) = j
                     n12(j) = n12(j) + 1
                     i12(n12(j),j) = i
                  end if
               end do
            end do
         end if
c
c     sort the attached atom lists into ascending order
c
         do i = 1, n
            call sort (n12(i),i12(1,i))
         end do
c
c     check for atom pairs with identical coordinates
c
         clash = .false.
         call chkxyz (clash)
c
c     write the TINKER coordinates and reset the connectivities
c
         call prtxyz (ixyz)
         do i = 1, n
            n12(i) = 0
         end do
c
c     read the next coordinate set from Protein Data Bank file
c
         call readpdb (ipdb)
      end do
c
c     write a sequence file for proteins and nucleic acids
c
      if (peptide .or. nucacid) then
         iseq = freeunit ()
         seqfile = filename(1:leng)//'.seq'
         call version (seqfile,'new')
         open (unit=iseq,file=seqfile,status='new')
         call prtseq (iseq)
         close (unit=iseq)
      end if
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      call final
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine oldatm  --  transfer coordinates from PDB  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "oldatm" get the Cartesian coordinates for an atom from
c     the Protein Data Bank file, then assigns the atom type
c     and atomic connectivities
c
c
      subroutine oldatm (i,bionum,i1,ires)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'fields.i'
      include 'iounit.i'
      include 'katoms.i'
      include 'sequen.i'
      include 'pdb.i'
      integer i,bionum
      integer i1,ires
c
c
c     get coordinates, assign atom type, and update connectivities
c
      if (bionum .ne. 0) then
         if (i .ne. 0) then
            type(n) = biotyp(bionum)
            if (type(n) .ne. 0) then
               name(n) = symbol(type(n))
            else
               name(n) = '   '
            end if
            x(n) = xpdb(i)
            y(n) = ypdb(i)
            z(n) = zpdb(i)
            if (i1 .ne. 0)  call addbond (n,i1)
            n = n + 1
         else
            write (iout,10)  bionum,ires,seq(ires)
   10       format (/,' OLDATM  --  A PDB Atom of Biotype',i5,
     &                 ' is Missing in Residue',i5,'-',a3)
            call fatal
         end if
      end if
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine newatm  --  create and define a new atom  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "newatm" creates and defines an atom needed for the
c     Cartesian coordinates file, but which may not present
c     in the original Protein Data Bank file
c
c
      subroutine newatm (i,bionum,ia,bond,ib,angle1,ic,angle2,chiral)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'fields.i'
      include 'katoms.i'
      include 'pdb.i'
      integer i,bionum
      integer ia,ib,ic
      integer chiral
      real*8 bond
      real*8 angle1
      real*8 angle2
c
c
c     set the atom type, compute coordinates, assign
c     connectivities and increment the atom counter
c
      if (bionum .ne. 0) then
         type(n) = biotyp(bionum)
         if (type(n) .ne. 0) then
            name(n) = symbol(type(n))
         else
            name(n) = '   '
         end if
         if (i .eq. 0) then
            call xyzatm (n,ia,bond,ib,angle1,ic,angle2,chiral)
         else
            x(n) = xpdb(i)
            y(n) = ypdb(i)
            z(n) = zpdb(i)
         end if
         call addbond (n,ia)
         n = n + 1
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine addbond  --  add a bond between two atoms  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "addbond" adds entries to the attached atoms list in
c     order to generate a direct connection between two atoms
c
c
      subroutine addbond (i,j)
      implicit none
      include 'sizes.i'
      include 'couple.i'
      integer i,j
c
c
c     add connectivity between the two atoms
c
      if (i.ne.0 .and. j.ne.0) then
         n12(i) = n12(i) + 1
         i12(n12(i),i) = j
         n12(j) = n12(j) + 1
         i12(n12(j),j) = i
      end if
      return
      end
c
c
c     ############################################################
c     ##                                                        ##
c     ##  subroutine findatm  --  locate PDB atom in a residue  ##
c     ##                                                        ##
c     ############################################################
c
c
c     "findatm" locates a specific PDB atom name type within a
c     range of atoms from the PDB file, returns zero if the name
c     type was not found
c
c
      subroutine findatm (name,start,stop,ipdb)
      implicit none
      include 'sizes.i'
      include 'pdb.i'
      integer i,ipdb
      integer start,stop
      character*4 name
c
c
c     search for the specified atom within the residue
c
      ipdb = 0
      do i = start, stop
         if (atmnam(i) .eq. name) then
            ipdb = i
            goto 10
         end if
      end do
   10 continue
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine ribosome  --  coordinates from PDB polypeptide  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "ribosome" translates a polypeptide structure in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine ribosome
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'resdue.i'
      include 'sequen.i'
      integer i,j,k
      integer ityp,nres
      integer jres,kres
      integer start,stop
      integer cyxtyp
      integer ncys,ndisulf
      integer resatm(2,maxres)
      integer ni(maxres),cai(maxres)
      integer ci(maxres),oi(maxres)
      integer icyclic(maxres)
      integer icys(maxres),idisulf(2,maxres)
      integer ntyp(maxamino),catyp(maxamino)
      integer ctyp(maxamino),hntyp(maxamino)
      integer otyp(maxamino),hatyp(maxamino)
      integer nntyp(maxamino),cantyp(maxamino)
      integer cntyp(maxamino),hnntyp(maxamino)
      integer ontyp(maxamino),hantyp(maxamino)
      integer nctyp(maxamino),cactyp(maxamino)
      integer cctyp(maxamino),hnctyp(maxamino)
      integer octyp(maxamino),hactyp(maxamino)
      real*8 rik,xcys(maxres)
      real*8 ycys(maxres)
      real*8 zcys(maxres)
      logical newchain
      logical midchain
      logical endchain
      logical cyclic
      logical header
      character*3 resname
      character*4 atmname
c
c     biopolymer atom types for amino acid backbone atoms
c
      data ntyp    /   1,   7,  15,  27,  41,  55,  65,  77,  87,
     &                96, 107, 122, 138, 161, 178, 194, 210, 220,
     &               232, 244, 258, 271, 287, 304, 318, 325,   0,
     &                 0,   0,   0,   1 /
      data catyp   /   2,   8,  16,  28,  42,  56,  66,  78,  88,
     &                97, 108, 123, 139, 162, 179, 195, 211, 221,
     &               233, 245, 259, 272, 288, 305, 319, 326,   0,
     &                 0,   0,   0,   2 /
      data ctyp    /   3,   9,  17,  29,  43,  57,  67,  79,  89,
     &                98, 109, 124, 140, 163, 180, 196, 212, 222,
     &               234, 246, 260, 273, 289, 306, 320, 327,   0,
     &                 0,   0,   0,   3 /
      data hntyp   /   4,  10,  18,  30,  44,  58,  68,  80,  90,
     &                 0, 110, 125, 141, 164, 181, 197, 213, 223,
     &               235, 247, 261, 274, 290, 307, 321, 328,   0,
     &                 0,   0,   0,   4 /
      data otyp    /   5,  11,  19,  31,  45,  59,  69,  81,  91,
     &                99, 111, 126, 142, 165, 182, 198, 214, 224,
     &               236, 248, 262, 275, 291, 308, 322, 329,   0,
     &                 0,   0,   0,   5 /
      data hatyp   /   6,  12,  20,  32,  46,  60,  70,  82,  92,
     &               100, 112, 127, 143, 166, 183, 199, 215, 225,
     &               237, 249, 263, 276, 292, 309,   0, 330,   0,
     &                 0,   0,   0,   6 /
c
c     biopolymer atom types for N-terminal backbone atoms
c
      data nntyp   / 350, 356, 362, 368, 374, 380, 386, 392, 398,
     &               404, 412, 418, 424, 430, 436, 442, 448, 454,
     &               460, 466, 472, 478, 484, 490, 496, 325,   0,
     &                 0,   0,   0, 350 /
      data cantyp  / 351, 357, 363, 369, 375, 381, 387, 393, 399,
     &               405, 413, 419, 425, 431, 437, 443, 449, 455,
     &               461, 467, 473, 479, 485, 491, 497, 326,   0,
     &               340,   0,   0, 351 /
      data cntyp   / 352, 358, 364, 370, 376, 382, 388, 394, 400,
     &               406, 414, 420, 426, 432, 438, 444, 450, 456,
     &               462, 468, 474, 480, 486, 492, 498, 327, 337,
     &               342,   0,   0, 352 /
      data hnntyp  / 353, 359, 365, 371, 377, 383, 389, 395, 401,
     &               407, 415, 421, 427, 433, 439, 445, 451, 457,
     &               463, 469, 475, 481, 487, 493, 499, 328,   0,
     &                 0,   0,   0, 353 /
      data ontyp   / 354, 360, 366, 372, 378, 384, 390, 396, 402,
     &               408, 416, 422, 428, 434, 440, 446, 452, 458,
     &               464, 470, 476, 482, 488, 494, 500, 329, 339,
     &               343,   0,   0, 354 /
      data hantyp  / 355, 361, 367, 373, 379, 385, 391, 397, 403,
     &               409, 417, 423, 429, 435, 441, 447, 453, 459,
     &               465, 471, 477, 483, 489, 495,   0, 330, 338,
     &               341,   0,   0, 355 /
c
c     biopolymer atom types for C-terminal backbone atoms
c
      data nctyp   / 501, 507, 513, 519, 525, 531, 537, 543, 549,
     &               555, 560, 566, 572, 578, 584, 590, 596, 602,
     &               608, 614, 620, 626, 632, 638, 644,   0,   0,
     &                 0, 344, 346, 501 /
      data cactyp  / 502, 508, 514, 520, 526, 532, 538, 544, 550,
     &               556, 561, 567, 573, 579, 585, 591, 597, 603,
     &               609, 615, 621, 627, 633, 639, 645,   0,   0,
     &                 0,   0, 348, 502 /
      data cctyp   / 503, 509, 515, 521, 527, 533, 539, 545, 551,
     &               557, 562, 568, 574, 580, 586, 592, 598, 604,
     &               610, 616, 622, 628, 634, 640, 646,   0,   0,
     &                 0,   0,   0, 503 /
      data hnctyp  / 504, 510, 516, 522, 528, 534, 540, 546, 552,
     &                 0, 563, 569, 575, 581, 587, 593, 599, 605,
     &               611, 617, 623, 629, 635, 641, 647,   0,   0,
     &                 0, 345, 347, 504 /
      data octyp   / 505, 511, 517, 523, 529, 535, 541, 547, 553,
     &               558, 564, 570, 576, 582, 588, 594, 600, 606,
     &               612, 618, 624, 630, 636, 642, 648,   0,   0,
     &                 0,   0,   0, 505 /
      data hactyp  / 506, 512, 518, 524, 530, 536, 542, 548, 554,
     &               559, 565, 571, 577, 583, 589, 595, 601, 607,
     &               613, 619, 625, 631, 637, 643,   0,   0,   0,
     &                 0,   0, 349, 506 /
c
c
c     set a pointer to the first and last atom of each residue
c
      nres = 0
      k = 0
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            if (resnum(i) .ne. k) then
               k = resnum(i)
               if (nres .ne. 0)  resatm(2,nres) = i - 1
               nres = nres + 1
               resatm(1,nres) = i
            end if
         end if
      end do
      if (nres .ne. 0)  resatm(2,nres) = npdb
c
c     get the three-letter sequence and code for each residue
c
      nseq = nres
      do i = 1, nres
         start = resatm(1,i)
         resname = resnam(start)
         seq(i) = 'UNK'
         seqtyp(i) = maxamino
         do k = 1, maxamino
            if (resname .eq. amino(k)) then
               seq(i) = amino(k)
               seqtyp(i) = k
               goto 10
            end if
         end do
   10    continue
      end do
c
c     search for the presence of cyclic polypeptide chains
c
      do i = 1, nres
         icyclic(i) = 0
      end do
      do i = 1, nchain
         jres = ichain(1,i)
         kres = ichain(2,i)
         start = resatm(1,jres)
         stop = resatm(2,jres)
         call findatm (' N  ',start,stop,j)
         start = resatm(1,kres)
         stop = resatm(2,kres)
         call findatm (' C  ',start,stop,k)
         if (j.ne.0 .and. k.ne.0) then
            rik = sqrt((xpdb(k)-xpdb(j))**2 + (ypdb(k)-ypdb(j))**2
     &                         + (zpdb(k)-zpdb(j))**2)
            if (rik .le. 3.0d0) then
               ni(jres) = j
               ci(kres) = k
               icyclic(jres) = kres
               icyclic(kres) = jres
            end if
         end if
      end do
c
c     search for any potential cystine disulfide bonds
c
      do i = 1, maxamino
         if (amino(i) .eq. 'CYX')  cyxtyp = i
      end do
      ncys = 0
      do i = 1, nres
         start = resatm(1,i)
         resname = resnam(start)
         if (resname.eq.'CYS' .or. resname.eq.'CYX') then
            stop = resatm(2,i)
            call findatm (' SG ',start,stop,k)
            ncys = ncys + 1
            icys(ncys) = i
            xcys(ncys) = xpdb(k)
            ycys(ncys) = ypdb(k)
            zcys(ncys) = zpdb(k)
         end if
      end do
      ndisulf = 0
      do i = 1, ncys-1
         do k = i+1, ncys
            rik = sqrt((xcys(i)-xcys(k))**2 + (ycys(i)-ycys(k))**2
     &                         + (zcys(i)-zcys(k))**2)
            if (rik .le. 3.0d0) then
               ndisulf = ndisulf + 1
               idisulf(1,ndisulf) = min(icys(i),icys(k))
               idisulf(2,ndisulf) = max(icys(i),icys(k))
            end if
         end do
      end do
      do i = 1, nres
         icys(i) = 0
      end do
      do i = 1, ndisulf
         icys(idisulf(1,i)) = 1
         icys(idisulf(2,i)) = 1
      end do
c
c     set the current atom to be the first atom
c
      n = 1
      newchain = .true.
c
c     locate and assign the atoms that make up each residue
c
      do i = 1, nres
         ityp = seqtyp(i)
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = resnam(start)
c
c     check that the maximum allowed atoms is not exceeded
c
         if (n+25 .gt. maxatm) then
            write (iout,20)  maxatm
   20       format (/,' RIBOSOME  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     check for a cysteine that should be changed to cystine
c
         if (icys(i) .eq. 1) then
            resname = 'CYX'
            seqtyp(i) = cyxtyp
            ityp = cyxtyp
         end if
c
c     test for the final residue of a peptide chain
c
         endchain = .false.
         do k = 1, nchain
            if (i .eq. ichain(2,k))  endchain = .true.
         end do
c
c     residue not at start or end of chain is in the middle
c
         if (newchain .or. endchain) then
            midchain = .false.
         else
            midchain = .true.
         end if
c
c     check to see if the current chain is cyclic
c
         cyclic = .false.
         if (newchain .or. endchain) then
            if (icyclic(i) .ne. 0)  cyclic = .true.
         end if
c
c     build the amide nitrogen of the current residue
c
         call findatm (' N  ',start,stop,k)
         if (k .ne. 0)  ni(i) = n
         if (midchain) then
            j = ntyp(ityp)
            call oldatm (k,j,ci(i-1),i)
         else if (newchain) then
            if (cyclic) then
               j = ntyp(ityp)
            else
               j = nntyp(ityp)
            end if
            call oldatm (k,j,0,i)
         else if (endchain) then
            if (cyclic) then
               j = ntyp(ityp)
            else
               j = nctyp(ityp)
            end if
            call oldatm (k,j,ci(i-1),i)
         end if
c
c     build the alpha carbon of the current residue
c
         atmname = ' CA '
         if (resname .eq. 'ACE')  atmname = ' CH3'
         if (resname .eq. 'NME')  atmname = ' CH3'
         call findatm (atmname,start,stop,k)
         if (k .ne. 0)  cai(i) = n
         if (midchain .or. cyclic .or. nres.eq.1) then
            j = catyp(ityp)
            call oldatm (k,j,ni(i),i)
         else if (newchain) then
            j = cantyp(ityp)
            call oldatm (k,j,ni(i),i)
         else if (endchain) then
            j = cactyp(ityp)
            call oldatm (k,j,ni(i),i)
         end if
c
c     build the carbonyl carbon of the current residue
c
         call findatm (' C  ',start,stop,k)
         if (k .ne. 0)  ci(i) = n
         if (midchain .or. cyclic) then
            j = ctyp(ityp)
            call oldatm (k,j,cai(i),i)
         else if (newchain) then
            j = cntyp(ityp)
            call oldatm (k,j,cai(i),i)
         else if (endchain) then
            j = cctyp(ityp)
            call oldatm (k,j,cai(i),i)
         end if
c
c     build the carbonyl oxygen of the current residue
c
         call findatm (' O  ',start,stop,k)
         if (k .ne. 0)  oi(i) = n
         if (midchain .or. cyclic) then
            j = otyp(ityp)
            call oldatm (k,j,ci(i),i)
         else if (newchain) then
            j = ontyp(ityp)
            call oldatm (k,j,ci(i),i)
         else if (endchain) then
            j = octyp(ityp)
            call oldatm (k,j,ci(i),i)
         end if
c
c     build the amide hydrogens of the current residue
c
         if (midchain .or. (endchain.and.cyclic)) then
            j = hntyp(ityp)
            call findatm (' H  ',start,stop,k)
            call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                      cai(i),119.0d0,1)
         else if (newchain .and. cyclic) then
            j = hntyp(ityp)
            call findatm (' H  ',start,stop,k)
            call newatm (k,j,ni(i),1.01d0,ci(icyclic(i)),119.0d0,
     &                      cai(i),119.0d0,1)
         else if (newchain) then
            j = hnntyp(ityp)
            if (resname .eq. 'PRO') then
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),0.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-120.0d0,0)
            else if (resname .eq. 'PCA') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-60.0d0,0)
            else
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,cai(i),109.5d0,
     &                         ci(i),-60.0d0,0)
            end if
         else if (endchain) then
            j = hnctyp(ityp)
            if (resname .eq. 'NH2') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),120.9d0,
     &                         cai(i-1),0.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),120.3d0,
     &                         cai(i-1),180.0d0,0)
            else if (resname .eq. 'NME') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                         cai(i),119.0d0,1)
            else
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ni(i),1.01d0,ci(i-1),119.0d0,
     &                         cai(i),119.0d0,1)
            end if
         end if
c
c     build the alpha hydrogen of the current residue
c
         if (resname .eq. 'GLY') then
            call findatm (' HA2',start,stop,k)
         else
            call findatm (' HA ',start,stop,k)
         end if
         if (midchain .or. cyclic) then
            j = hatyp(ityp)
            call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                      ci(i),109.5d0,-1)
         else if (newchain) then
            j = hantyp(ityp)
            if (resname .eq. 'FOR') then
               call findatm (' H  ',start,stop,k)
               call newatm (k,j,ci(i),1.12d0,oi(i),0.0d0,0,0.0d0,0)
            else if (resname .eq. 'ACE') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ci(i),109.5d0,
     &                         oi(i),-60.0d0,0)
            else
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i),109.5d0,-1)
            end if
         else if (endchain) then
            j = hactyp(ityp)
            if (resname .eq. 'NME') then
               call findatm (' H1 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),180.0d0,0)
               call findatm (' H2 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),60.0d0,0)
               call findatm (' H3 ',start,stop,k)
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i-1),-60.0d0,0)
            else
               call newatm (k,j,cai(i),1.10d0,ni(i),109.5d0,
     &                         ci(i),109.5d0,-1)
            end if
         end if
c
c     build the side chain atoms of the current residue
c
         call addside (resname,i,start,stop,cai(i),ni(i),ci(i),icys(i))
c
c     build the terminal oxygen at the end of a peptide chain
c
         if (endchain .and. .not.cyclic) then
            call findatm (' OXT',start,stop,k)
            if (k .eq. 0)  call findatm (' OT2',start,stop,k)
            j = octyp(ityp)
            call newatm (k,j,ci(i),1.25d0,cai(i),117.0d0,
     &                      oi(i),126.0d0,1)
         end if
c
c     next residue starts a new chain if current one just ended
c
         if (endchain) then
            newchain = .true.
         else
            newchain = .false.
         end if
      end do
c
c     connect the terminal residues of a cyclic polypeptide
c
      header = .true.
      do i = 1, nchain
         j = ichain(1,i)
         k = icyclic(j)
         if (k .ne. 0) then
            call addbond (ni(j),ci(k))
            if (verbose) then
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format ()
               end if
               write (iout,40)  j,k
   40          format (' Peptide Cyclization between Residues :  ',2i5)
            end if
         end if
      end do
c
c     connect the sulfur atoms involved in disulfide bonds
c
      header = .true.
      do i = 1, ndisulf
         call addbond (icys(idisulf(1,i)),icys(idisulf(2,i)))
         if (verbose) then
            if (header) then
               header = .false.
               write (iout,50)
   50          format ()
            end if
            write (iout,60)  idisulf(1,i),idisulf(2,i)
   60       format (' Disulfide Bond between Residues :  ',2i5)
         end if
      end do
c
c     total number of atoms is one less than the current atom
c
      n = n - 1
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine addside  --  build the amino acid side chains  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "addside" builds the Cartesian coordinates for a single amino
c     acid side chain; coordinates are read from the Protein Data
c     Bank file or found from internal coordinates, then atom types
c     are assigned and connectivity data generated
c
c
      subroutine addside (resname,ires,start,stop,cai,ni,ci,cystype)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'sequen.i'
      integer i,ires
      integer start,stop
      integer cai,ni,ci
      integer cystype
      character*3 resname
c
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         call findatm (' HA3',start,stop,i)
         if (ires .eq. 1) then
            call newatm (i,355,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         else if (ires .eq. nseq) then
            call newatm (i,506,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         else
            call newatm (i,6,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         end if
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,13,cai,ires)
         call findatm (' HB1',start,stop,i)
         call newatm (i,14,n-1,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm (' HB2',start,stop,i)
         call newatm (i,14,n-2,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm (' HB3',start,stop,i)
         call newatm (i,14,n-3,1.10d0,cai,110.2d0,ni,-60.0d0,0)
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,21,cai,ires)
         call findatm (' CG1',start,stop,i)
         call oldatm (i,23,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,25,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,22,n-3,1.10d0,cai,107.0d0,n-2,108.2d0,1)
         call findatm ('HG11',start,stop,i)
         call newatm (i,24,n-3,1.10d0,n-4,111.6d0,cai,180.0d0,0)
         call findatm ('HG12',start,stop,i)
         call newatm (i,24,n-4,1.10d0,n-5,111.6d0,cai,60.0d0,0)
         call findatm ('HG13',start,stop,i)
         call newatm (i,24,n-5,1.10d0,n-6,111.6d0,cai,-60.0d0,0)
         call findatm ('HG21',start,stop,i)
         call newatm (i,26,n-5,1.10d0,n-7,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,26,n-6,1.10d0,n-8,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,26,n-7,1.10d0,n-9,111.6d0,cai,-60.0d0,0)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,33,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,35,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,37,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,39,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,34,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,34,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,36,n-5,1.10d0,n-6,107.0d0,n-4,108.2d0,1)
         call findatm ('HD11',start,stop,i)
         call newatm (i,38,n-5,1.10d0,n-6,111.6d0,n-7,180.0d0,0)
         call findatm ('HD12',start,stop,i)
         call newatm (i,38,n-6,1.10d0,n-7,111.6d0,n-8,60.0d0,0)
         call findatm ('HD13',start,stop,i)
         call newatm (i,38,n-7,1.10d0,n-8,111.6d0,n-9,-60.0d0,0)
         call findatm ('HD21',start,stop,i)
         call newatm (i,40,n-7,1.10d0,n-9,111.6d0,n-10,180.0d0,0)
         call findatm ('HD22',start,stop,i)
         call newatm (i,40,n-8,1.10d0,n-10,111.6d0,n-11,60.0d0,0)
         call findatm ('HD23',start,stop,i)
         call newatm (i,40,n-9,1.10d0,n-11,111.6d0,n-12,-60.0d0,0)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,47,cai,ires)
         call findatm (' CG1',start,stop,i)
         call oldatm (i,49,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,51,n-2,ires)
         call findatm (' CD1',start,stop,i)
         if (i .eq. 0)  call findatm (' CD ',start,stop,i)
         call oldatm (i,53,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,48,n-4,1.10d0,cai,107.0d0,n-3,108.2d0,-1)
         call findatm ('HG12',start,stop,i)
         call newatm (i,50,n-4,1.10d0,n-5,109.5d0,n-2,109.5d0,1)
         call findatm ('HG13',start,stop,i)
         call newatm (i,50,n-5,1.10d0,n-6,109.5d0,n-3,109.5d0,-1)
         call findatm ('HG21',start,stop,i)
         call newatm (i,52,n-5,1.10d0,n-7,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,52,n-6,1.10d0,n-8,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,52,n-7,1.10d0,n-9,111.6d0,cai,-60.0d0,0)
         call findatm ('HD11',start,stop,i)
         call newatm (i,54,n-7,1.10d0,n-9,111.6d0,n-10,180.0d0,0)
         call findatm ('HD12',start,stop,i)
         call newatm (i,54,n-8,1.10d0,n-10,111.6d0,n-11,60.0d0,0)
         call findatm ('HD13',start,stop,i)
         call newatm (i,54,n-9,1.10d0,n-11,111.6d0,n-12,-60.0d0,0)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,61,cai,ires)
         call findatm (' OG ',start,stop,i)
         call oldatm (i,63,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,62,n-2,1.10d0,cai,109.2d0,n-1,109.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,62,n-3,1.10d0,cai,109.2d0,n-2,109.5d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,64,n-3,0.94d0,n-4,106.9d0,cai,180.0d0,0)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,71,cai,ires)
         call findatm (' OG1',start,stop,i)
         call oldatm (i,73,n-1,ires)
         call findatm (' CG2',start,stop,i)
         call oldatm (i,75,n-2,ires)
         call findatm (' HB ',start,stop,i)
         call newatm (i,72,n-3,1.10d0,cai,107.0d0,n-2,108.2d0,-1)
         call findatm (' HG1',start,stop,i)
         call newatm (i,74,n-3,0.94d0,n-4,106.9d0,cai,180.0d0,0)
         call findatm ('HG21',start,stop,i)
         call newatm (i,76,n-3,1.10d0,n-5,111.6d0,cai,180.0d0,0)
         call findatm ('HG22',start,stop,i)
         call newatm (i,76,n-4,1.10d0,n-6,111.6d0,cai,60.0d0,0)
         call findatm ('HG23',start,stop,i)
         call newatm (i,76,n-5,1.10d0,n-7,111.6d0,cai,-60.0d0,0)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,83,cai,ires)
         call findatm (' SG ',start,stop,i)
         call oldatm (i,85,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,84,n-2,1.10d0,cai,109.5d0,n-1,107.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,84,n-3,1.10d0,cai,109.5d0,n-2,107.5d0,-1)
         call findatm (' HG ',start,stop,i)
         call newatm (i,86,n-3,1.34d0,n-4,96.0d0,cai,180.0d0,0)
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,93,cai,ires)
         call findatm (' SG ',start,stop,i)
         cystype = n
         call oldatm (i,95,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,94,n-2,1.10d0,cai,109.5d0,n-1,107.5d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,94,n-3,1.10d0,cai,109.5d0,n-2,107.5d0,-1)
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,101,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,103,n-1,ires)
         call findatm (' CD ',start,stop,i)
         if (ires .eq. 1) then
            call oldatm (i,410,n-1,ires)
         else
            call oldatm (i,105,n-1,ires)
         end if
         call addbond (n-1,ni)
         call findatm (' HB2',start,stop,i)
         call newatm (i,102,n-3,1.10d0,cai,111.2d0,n-2,111.2d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,102,n-4,1.10d0,cai,111.2d0,n-3,111.2d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,104,n-4,1.10d0,n-5,111.2d0,n-3,111.2d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,104,n-5,1.10d0,n-6,111.2d0,n-4,111.2d0,-1)
         if (ires .eq. 1) then
            call findatm (' HD2',start,stop,i)
            call newatm (i,411,n-5,1.10d0,n-6,111.2d0,ni,111.2d0,1)
            call findatm (' HD3',start,stop,i)
            call newatm (i,411,n-6,1.10d0,n-7,111.2d0,ni,111.2d0,-1)
         else
            call findatm (' HD2',start,stop,i)
            call newatm (i,106,n-5,1.10d0,n-6,111.2d0,ni,111.2d0,1)
            call findatm (' HD3',start,stop,i)
            call newatm (i,106,n-6,1.10d0,n-7,111.2d0,ni,111.2d0,-1)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,113,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,115,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,116,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,116,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,118,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,118,n-2,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,120,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,114,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,114,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,117,n-7,1.09d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,117,n-7,1.09d0,n-9,120.0d0,n-10,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,119,n-7,1.09d0,n-9,120.0d0,n-10,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,119,n-7,1.09d0,n-9,120.0d0,n-11,180.0d0,0)
         call findatm (' HZ ',start,stop,i)
         call newatm (i,121,n-7,1.09d0,n-8,120.0d0,n-10,180.0d0,0)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,128,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,130,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,131,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,131,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,133,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,133,n-2,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,135,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' OH ',start,stop,i)
         call oldatm (i,136,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,129,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,129,n-9,1.10d0,cai,107.9d0,n-8,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,132,n-8,1.09d0,n-9,120.0d0,n-10,0.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,132,n-8,1.09d0,n-10,120.0d0,n-11,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,134,n-8,1.09d0,n-10,120.0d0,n-11,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,134,n-8,1.09d0,n-10,120.0d0,n-12,180.0d0,0)
         call findatm (' HH ',start,stop,i)
         call newatm (i,137,n-7,0.97d0,n-8,108.0d0,n-9,0.0d0,0)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,144,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,146,n-1,ires)
         call findatm (' CD1',start,stop,i)
         call oldatm (i,147,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,149,n-2,ires)
         call findatm (' NE1',start,stop,i)
         call oldatm (i,150,n-2,ires)
         call findatm (' CE2',start,stop,i)
         call oldatm (i,152,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' CE3',start,stop,i)
         call oldatm (i,153,n-3,ires)
         call findatm (' CZ2',start,stop,i)
         call oldatm (i,155,n-2,ires)
         call findatm (' CZ3',start,stop,i)
         call oldatm (i,157,n-2,ires)
         call findatm (' CH2',start,stop,i)
         call oldatm (i,159,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,145,n-10,1.10d0,cai,107.9d0,n-9,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,145,n-11,1.10d0,cai,107.9d0,n-10,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,148,n-10,1.09d0,n-11,126.0d0,n-12,0.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,151,n-9,1.01d0,n-11,126.3d0,n-12,180.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,154,n-8,1.09d0,n-6,120.0d0,n-5,180.0d0,0)
         call findatm (' HZ2',start,stop,i)
         call newatm (i,156,n-8,1.09d0,n-6,120.0d0,n-7,180.0d0,0)
         call findatm (' HZ3',start,stop,i)
         call newatm (i,158,n-8,1.09d0,n-7,120.0d0,n-9,180.0d0,0)
         call findatm (' HH2',start,stop,i)
         call newatm (i,160,n-8,1.09d0,n-9,120.0d0,n-11,180.0d0,0)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,167,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,169,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,170,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,172,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,174,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,176,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,168,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,168,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,171,n-6,1.02d0,n-4,126.0d0,n-3,180.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,173,n-6,1.09d0,n-4,126.0d0,n-5,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,175,n-6,1.09d0,n-5,126.0d0,n-7,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,177,n-6,1.02d0,n-7,126.0d0,n-9,180.0d0,0)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,184,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,186,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,187,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,189,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,191,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,193,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,185,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,185,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD1',start,stop,i)
         call newatm (i,188,n-6,1.02d0,n-4,126.0d0,n-3,180.0d0,0)
         call findatm (' HD2',start,stop,i)
         call newatm (i,190,n-6,1.09d0,n-4,126.0d0,n-5,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,192,n-6,1.09d0,n-5,126.0d0,n-7,180.0d0,0)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,200,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,202,n-1,ires)
         call findatm (' ND1',start,stop,i)
         call oldatm (i,203,n-1,ires)
         call findatm (' CD2',start,stop,i)
         call oldatm (i,204,n-2,ires)
         call findatm (' CE1',start,stop,i)
         call oldatm (i,206,n-2,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,208,n-1,ires)
         call addbond (n-1,n-3)
         call findatm (' HB2',start,stop,i)
         call newatm (i,201,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,201,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,205,n-5,1.09d0,n-3,126.0d0,n-4,180.0d0,0)
         call findatm (' HE1',start,stop,i)
         call newatm (i,207,n-5,1.09d0,n-4,126.0d0,n-6,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,209,n-5,1.02d0,n-6,126.0d0,n-8,180.0d0,0)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,216,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,218,n-1,ires)
         call findatm (' OD1',start,stop,i)
         call oldatm (i,219,n-1,ires)
         call findatm (' OD2',start,stop,i)
         call oldatm (i,219,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,217,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,217,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,226,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,228,n-1,ires)
         call findatm (' OD1',start,stop,i)
         call oldatm (i,229,n-1,ires)
         call findatm (' ND2',start,stop,i)
         call oldatm (i,230,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,227,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,227,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm ('HD21',start,stop,i)
         call newatm (i,231,n-3,1.01d0,n-5,120.9d0,n-6,0.0d0,0)
         call findatm ('HD22',start,stop,i)
         call newatm (i,231,n-4,1.01d0,n-6,120.3d0,n-7,180.0d0,0)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,238,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,240,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,242,n-1,ires)
         call findatm (' OE1',start,stop,i)
         call oldatm (i,243,n-1,ires)
         call findatm (' OE2',start,stop,i)
         call oldatm (i,243,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,239,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,239,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,241,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,241,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,250,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,252,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,254,n-1,ires)
         call findatm (' OE1',start,stop,i)
         call oldatm (i,255,n-1,ires)
         call findatm (' NE2',start,stop,i)
         call oldatm (i,256,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,251,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,251,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,253,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,253,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm ('HE21',start,stop,i)
         call newatm (i,257,n-5,1.01d0,n-7,120.9d0,n-8,0.0d0,0)
         call findatm ('HE22',start,stop,i)
         call newatm (i,257,n-6,1.01d0,n-8,120.3d0,n-9,180.0d0,0)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,264,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,266,n-1,ires)
         call findatm (' SD ',start,stop,i)
         call oldatm (i,268,n-1,ires)
         call findatm (' CE ',start,stop,i)
         call oldatm (i,269,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,265,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,265,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,267,n-5,1.10d0,n-6,109.5d0,n-4,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,267,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,-1)
         call findatm (' HE1',start,stop,i)
         call newatm (i,270,n-5,1.10d0,n-6,110.2d0,n-7,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,270,n-6,1.10d0,n-7,110.2d0,n-8,60.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,270,n-7,1.10d0,n-8,110.2d0,n-9,-60.0d0,0)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,277,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,279,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,281,n-1,ires)
         call findatm (' CE ',start,stop,i)
         call oldatm (i,283,n-1,ires)
         call findatm (' NZ ',start,stop,i)
         call oldatm (i,285,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,278,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,278,n-6,1.10d0,cai,107.9d0,n-5,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,280,n-6,1.10d0,n-7,109.5d0,n-5,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,280,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,282,n-7,1.10d0,n-8,109.5d0,n-6,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,282,n-8,1.10d0,n-9,109.5d0,n-7,109.5d0,-1)
         call findatm (' HE2',start,stop,i)
         call newatm (i,284,n-8,1.10d0,n-9,110.9d0,n-7,107.3d0,1)
         call findatm (' HE3',start,stop,i)
         call newatm (i,284,n-9,1.10d0,n-10,110.9d0,n-8,107.3d0,-1)
         call findatm (' HZ1',start,stop,i)
         call newatm (i,286,n-9,1.04d0,n-10,110.5d0,n-11,180.0d0,0)
         call findatm (' HZ2',start,stop,i)
         call newatm (i,286,n-10,1.04d0,n-11,110.5d0,n-12,60.0d0,0)
         call findatm (' HZ3',start,stop,i)
         call newatm (i,286,n-11,1.04d0,n-12,110.5d0,n-13,-60.0d0,0)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,293,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,295,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,297,n-1,ires)
         call findatm (' NE ',start,stop,i)
         call oldatm (i,299,n-1,ires)
         call findatm (' CZ ',start,stop,i)
         call oldatm (i,301,n-1,ires)
         call findatm (' NH1',start,stop,i)
         call oldatm (i,302,n-1,ires)
         call findatm (' NH2',start,stop,i)
         call oldatm (i,302,n-2,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,294,n-7,1.10d0,cai,107.9d0,n-6,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,294,n-8,1.10d0,cai,107.9d0,n-7,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,296,n-8,1.10d0,n-9,109.5d0,n-7,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,296,n-9,1.10d0,n-10,109.5d0,n-8,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,298,n-9,1.10d0,n-10,109.5d0,n-8,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,298,n-10,1.10d0,n-11,109.5d0,n-9,109.5d0,-1)
         call findatm (' HE ',start,stop,i)
         call newatm (i,300,n-10,1.01d0,n-11,118.5d0,n-9,120.0d0,1)
         call findatm ('HH11',start,stop,i)
         call newatm (i,303,n-9,1.01d0,n-10,122.5d0,n-11,0.0d0,0)
         call findatm ('HH12',start,stop,i)
         call newatm (i,303,n-10,1.01d0,n-11,118.8d0,n-12,180.0d0,0)
         call findatm ('HH21',start,stop,i)
         call newatm (i,303,n-10,1.01d0,n-12,122.5d0,n-13,0.0d0,0)
         call findatm ('HH22',start,stop,i)
         call newatm (i,303,n-11,1.01d0,n-13,118.8d0,n-14,180.0d0,0)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,310,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,312,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,314,n-1,ires)
         call findatm (' NE ',start,stop,i)
         call oldatm (i,316,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,311,n-4,1.10d0,cai,107.9d0,n-3,110.0d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,311,n-5,1.10d0,cai,107.9d0,n-4,110.0d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,313,n-5,1.10d0,n-7,109.5d0,n-4,109.5d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,313,n-6,1.10d0,n-8,109.5d0,n-5,109.5d0,-1)
         call findatm (' HD2',start,stop,i)
         call newatm (i,315,n-6,1.10d0,n-8,109.5d0,n-5,109.5d0,1)
         call findatm (' HD3',start,stop,i)
         call newatm (i,315,n-7,1.10d0,n-9,109.5d0,n-6,109.5d0,-1)
         call findatm (' HE1',start,stop,i)
         call newatm (i,317,n-7,1.04d0,n-8,110.5d0,n-9,180.0d0,0)
         call findatm (' HE2',start,stop,i)
         call newatm (i,317,n-8,1.04d0,n-9,110.5d0,n-10,60.0d0,0)
         call findatm (' HE3',start,stop,i)
         call newatm (i,317,n-9,1.04d0,n-10,110.5d0,n-11,-60.0d0,0)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call findatm (' CB1',start,stop,i)
         call oldatm (i,323,cai,ires)
         call findatm (' CB2',start,stop,i)
         call oldatm (i,323,cai,ires)
         call findatm ('HB11',start,stop,i)
         call newatm (i,324,n-2,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm ('HB12',start,stop,i)
         call newatm (i,324,n-3,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm ('HB13',start,stop,i)
         call newatm (i,324,n-4,1.10d0,cai,110.2d0,ni,-60.0d0,0)
         call findatm ('HB21',start,stop,i)
         call newatm (i,324,n-4,1.10d0,cai,110.2d0,ni,180.0d0,0)
         call findatm ('HB22',start,stop,i)
         call newatm (i,324,n-5,1.10d0,cai,110.2d0,ni,60.0d0,0)
         call findatm ('HB23',start,stop,i)
         call newatm (i,324,n-6,1.10d0,cai,110.2d0,ni,-60.0d0,0)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call findatm (' CB ',start,stop,i)
         call oldatm (i,331,cai,ires)
         call findatm (' CG ',start,stop,i)
         call oldatm (i,333,n-1,ires)
         call findatm (' CD ',start,stop,i)
         call oldatm (i,335,n-1,ires)
         call addbond (n-1,ni)
         call findatm (' OE ',start,stop,i)
         call oldatm (i,336,n-1,ires)
         call findatm (' HB2',start,stop,i)
         call newatm (i,332,n-4,1.10d0,cai,111.2d0,n-3,111.2d0,1)
         call findatm (' HB3',start,stop,i)
         call newatm (i,332,n-5,1.10d0,cai,111.2d0,n-4,111.2d0,-1)
         call findatm (' HG2',start,stop,i)
         call newatm (i,334,n-5,1.10d0,n-6,111.2d0,n-4,111.2d0,1)
         call findatm (' HG3',start,stop,i)
         call newatm (i,334,n-6,1.10d0,n-7,111.2d0,n-5,111.2d0,-1)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         if (ires .eq. 1) then
            call newatm (0,355,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         else if (ires .eq. nseq) then
            call newatm (0,506,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         else
            call newatm (0,6,cai,1.10d0,ni,109.5d0,ci,109.5d0,1)
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine ligase  --  coordinates from PDB nucleic acid  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "ligase" translates a nucleic acid structure in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine ligase
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'pdb.i'
      include 'resdue.i'
      include 'sequen.i'
      integer i,j,k
      integer ityp,nres
      integer start,stop
      integer poi,o5i,c5i
      integer c4i,o4i,c1i
      integer c3i,c2i,o3i,o2i
      integer resatm(2,maxres)
      integer o5typ(maxnuc),c5typ(maxnuc)
      integer h51typ(maxnuc),h52typ(maxnuc)
      integer c4typ(maxnuc),h4typ(maxnuc)
      integer o4typ(maxnuc),c1typ(maxnuc)
      integer h1typ(maxnuc),c3typ(maxnuc)
      integer h3typ(maxnuc),c2typ(maxnuc)
      integer o3typ(maxnuc),o2typ(maxnuc)
      integer h21typ(maxnuc),h22typ(maxnuc)
      integer ptyp(maxnuc),optyp(maxnuc)
      integer h5ttyp(maxnuc),h3ttyp(maxnuc)
      logical newchain,endchain
      logical deoxy(maxres)
      character*3 resname
c
c     biopolymer atom types for nucleic acid backbone atoms
c
      data o5typ   / 1001, 1031, 1062, 1090, 1117, 1146, 1176, 1203,
     &                  0,    0,    0,    0 /
      data c5typ   / 1002, 1032, 1063, 1091, 1118, 1147, 1177, 1204,
     &                  0,    0,    0,    0 /
      data h51typ  / 1003, 1033, 1064, 1092, 1119, 1148, 1178, 1205,
     &                  0,    0,    0,    0 /
      data h52typ  / 1004, 1034, 1065, 1093, 1120, 1149, 1179, 1206,
     &                  0,    0,    0,    0 /
      data c4typ   / 1005, 1035, 1066, 1094, 1121, 1150, 1180, 1207,
     &                  0,    0,    0,    0 /
      data h4typ   / 1006, 1036, 1067, 1095, 1122, 1151, 1181, 1208,
     &                  0,    0,    0,    0 /
      data o4typ   / 1007, 1037, 1068, 1096, 1123, 1152, 1182, 1209,
     &                  0,    0,    0,    0 /
      data c1typ   / 1008, 1038, 1069, 1097, 1124, 1153, 1183, 1210,
     &                  0,    0,    0,    0 /
      data h1typ   / 1009, 1039, 1070, 1098, 1125, 1154, 1184, 1211,
     &                  0,    0,    0,    0 /
      data c3typ   / 1010, 1040, 1071, 1099, 1126, 1155, 1185, 1212,
     &                  0,    0,    0,    0 /
      data h3typ   / 1011, 1041, 1072, 1100, 1127, 1156, 1186, 1213,
     &                  0,    0,    0,    0 /
      data c2typ   / 1012, 1042, 1073, 1101, 1128, 1157, 1187, 1214,
     &                  0,    0,    0,    0 /
      data h21typ  / 1013, 1043, 1074, 1102, 1129, 1158, 1188, 1215,
     &                  0,    0,    0,    0 /
      data h22typ  / 1015, 1045, 1076, 1104, 1130, 1159, 1189, 1216,
     &                  0,    0,    0,    0 /
      data o3typ   / 1016, 1046, 1077, 1105, 1131, 1160, 1190, 1217,
     &                  0,    0,    0,    0 /
      data o2typ   / 1014, 1044, 1075, 1103,    0,    0,    0,    0,
     &                  0,    0,    0,    0 /
      data ptyp    / 1230, 1230, 1230, 1230, 1242, 1242, 1242, 1242,
     &                  0,    0,    0,    0 /
      data optyp   / 1231, 1231, 1231, 1231, 1243, 1243, 1243, 1243,
     &                  0,    0,    0,    0 /
      data h5ttyp  / 1233, 1233, 1233, 1233, 1245, 1245, 1245, 1245,
     &                  0,    0,    0,    0 /
      data h3ttyp  / 1238, 1238, 1238, 1238, 1250, 1250, 1250, 1250,
     &                  0,    0,    0,    0 /
c
c
c     set a pointer to the first and last atom of each residue
c
      nres = 0
      k = 0
      do i = 1, npdb
         if (pdbtyp(i) .eq. 'ATOM  ') then
            if (resnum(i) .ne. k) then
               k = resnum(i)
               if (nres .ne. 0)  resatm(2,nres) = i - 1
               nres = nres + 1
               resatm(1,nres) = i
            end if
         end if
      end do
      if (nres .ne. 0)  resatm(2,nres) = npdb
c
c     check for deoxyribose and change residue name if necessary
c
      do i = 1, nres
         deoxy(i) = .false.
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = resnam(start)
         call findatm (' O2''',start,stop,k)
         if (k .eq. 0) then
            deoxy(i) = .true.
            do j = start, stop
               if (resname .eq. 'A  ')  resnam(j) = 'DA '
               if (resname .eq. 'G  ')  resnam(j) = 'DG '
               if (resname .eq. 'C  ')  resnam(j) = 'DC '
               if (resname .eq. 'U  ')  resnam(j) = 'DU '
               if (resname .eq. 'T  ')  resnam(j) = 'DT '
            end do
         end if
      end do
c
c     get the three-letter sequence and code for each residue
c
      nseq = nres
      do i = 1, nres
         start = resatm(1,i)
         resname = resnam(start)
         seq(i) = 'UNK'
         seqtyp(i) = maxnuc
         do k = 1, maxnuc
            if (resname .eq. nuclz(k)) then
               seq(i) = nuclz(k)
               seqtyp(i) = k
               goto 10
            end if
         end do
   10    continue
      end do
c
c     set the current atom to be the first atom
c
      n = 1
c
c     locate and assign the atoms that make up each residue
c
      do i = 1, nres
         ityp = seqtyp(i)
         start = resatm(1,i)
         stop = resatm(2,i)
         resname = resnam(start)
c
c     check that the maximum allowed atoms is not exceeded
c
         if (n+25 .gt. maxatm) then
            write (iout,20)  maxatm
   20       format (/,' LIGASE  --  The Maximum of',i8,' Atoms',
     &                 ' has been Exceeded')
            call fatal
         end if
c
c     test for initial or final residue of a nucleotide chain
c
         newchain = .false.
         endchain = .false.
         do j = 1, nchain
            if (i .eq. ichain(1,j)) then
               newchain = .true.
               poi = 0
               o3i = 0
            end if
            if (i .eq. ichain(2,j))  endchain = .true.
         end do
c
c     build the phosphate atoms of the current residue
c
         if (resname .eq. 'TP ') then

         else if (resname .eq. 'DP ') then

         else if (resname .eq. 'MP ') then

         else if (.not. newchain) then
            call findatm (' P  ',start,stop,k)
            if (k .ne. 0)  poi = n
            j = ptyp(ityp)
            call oldatm (k,j,o3i,i)
            call findatm (' OP1',start,stop,k)
            j = optyp(ityp)
            call oldatm (k,j,n-1,i)
            call findatm (' OP2',start,stop,k)
            j = optyp(ityp)
            call oldatm (k,j,n-2,i)
         end if
c
c     build the ribose sugar atoms of the current residue
c
         call findatm (' O5''',start,stop,k)
         if (k .ne. 0)  o5i = n
         j = o5typ(ityp)
         if (newchain) then
            if (deoxy(i)) then
               j = 1244
            else
               j = 1232
            end if
         end if
         call oldatm (k,j,poi,i)
         call findatm (' C5''',start,stop,k)
         if (k .ne. 0)  c5i = n
         j = c5typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C4''',start,stop,k)
         if (k .ne. 0)  c4i = n
         j = c4typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' O4''',start,stop,k)
         if (k .ne. 0)  o4i = n
         j = o4typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C1''',start,stop,k)
         if (k .ne. 0)  c1i = n
         j = c1typ(ityp)
         call oldatm (k,j,n-1,i)
         call findatm (' C3''',start,stop,k)
         if (k .ne. 0)  c3i = n
         j = c3typ(ityp)
         call oldatm (k,j,n-3,i)
         call findatm (' C2''',start,stop,k)
         if (k .ne. 0)  c2i = n
         j = c2typ(ityp)
         call oldatm (k,j,n-1,i)
         call addbond (n-1,n-3)
         call findatm (' O3''',start,stop,k)
         if (k .ne. 0)  o3i = n
         j = o3typ(ityp)
         if (endchain) then
            if (deoxy(i)) then
               j = 1249
            else
               j = 1237
            end if
         end if
         call oldatm (k,j,n-2,i)
         if (.not. deoxy(i)) then
            call findatm (' O2''',start,stop,k)
            if (k .ne. 0)  o2i = n
            j = o2typ(ityp)
            call oldatm (k,j,n-2,i)
         end if
c
c     build the hydrogen atoms of the current residue
c
         if (newchain) then
            call findatm (' H5T',start,stop,k)
            j = h5ttyp(ityp)
            call newatm (k,j,o5i,1.00d0,c5i,109.5d0,c4i,180.0d0,0)
         end if
         call findatm (' H5''',start,stop,k)
         j = h51typ(ityp)
         call newatm (k,j,c5i,1.09d0,o5i,109.5d0,c4i,109.5d0,1)
         call findatm ('H5''''',start,stop,k)
         j = h52typ(ityp)
         call newatm (k,j,c5i,1.09d0,o5i,109.5d0,c4i,109.5d0,-1)
         call findatm (' H4''',start,stop,k)
         j = h4typ(ityp)
         call newatm (k,j,c4i,1.09d0,c5i,109.5d0,c3i,109.5d0,-1)
         call findatm (' H3''',start,stop,k)
         j = h3typ(ityp)
         call newatm (k,j,c3i,1.09d0,c4i,109.5d0,c2i,109.5d0,-1)
         if (deoxy(i)) then
            call findatm (' H2''',start,stop,k)
            j = h21typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,-1)
            call findatm ('H2''''',start,stop,k)
            j = h22typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,1)
         else
            call findatm (' H2''',start,stop,k)
            j = h21typ(ityp)
            call newatm (k,j,c2i,1.09d0,c3i,109.5d0,c1i,109.5d0,-1)
            call findatm ('HO2''',start,stop,k)
            j = h22typ(ityp)
            call newatm (k,j,o2i,1.00d0,c2i,109.5d0,c3i,180.0d0,0)
         end if
         call findatm (' H1''',start,stop,k)
         j = h1typ(ityp)
         call newatm (k,j,c1i,1.09d0,o4i,109.5d0,c2i,109.5d0,-1)
         if (endchain) then
            call findatm (' H3T',start,stop,k)
            j = h3ttyp(ityp)
            call newatm (k,j,o3i,1.00d0,c3i,109.5d0,c4i,180.0d0,0)
         end if
c
c     build the standard base atoms of the current residue
c
         call addbase (resname,i,start,stop,c1i)
      end do
c
c     total number of atoms is one less than the current atom
c
      n = n - 1
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine addbase  --  build a single nucleic acid base  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "addbase" builds the Cartesian coordinates for a single nucleic
c     acid base; coordinates are read from the Protein Data Bank file
c     or found from internal coordinates, then atom types are assigned
c     and connectivity data generated
c
c
      subroutine addbase (resname,ires,start,stop,c1i)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      integer i,ires
      integer start,stop
      integer c1i
      character*3 resname
c
c
c     adenine in adenosine residue  (A)
c
      if (resname .eq. 'A ') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1017,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1021,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1020,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1019,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1025,n-1,ires)
         call findatm (' N6 ',start,stop,i)
         call oldatm (i,1027,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1024,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1023,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1022,n-1,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1018,n-1,ires)
         call addbond (n-1,n-7)
         call addbond (n-1,n-10)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1030,n-9,1.08d0,n-8,123.1d0,n-7,180.0d0,0)
         call findatm (' H61',start,stop,i)
         call newatm (i,1028,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
         call findatm (' H62',start,stop,i)
         call newatm (i,1029,n-7,1.00d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' H2 ',start,stop,i)
         call newatm (i,1026,n-6,1.08d0,n-5,115.4d0,n-4,180.0d0,0)
c
c     guanine in guanosine residue  (G)
c
      else if (resname .eq. 'G  ') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1047,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1051,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1050,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1049,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1055,n-1,ires)
         call findatm (' O6 ',start,stop,i)
         call oldatm (i,1060,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1054,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1053,n-1,ires)
         call findatm (' N2 ',start,stop,i)
         call oldatm (i,1057,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1052,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1048,n-1,ires)
         call addbond (n-1,n-8)
         call addbond (n-1,n-11)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1061,n-10,1.08d0,n-9,123.0d0,n-8,180.0d0,0)
         call findatm (' H1 ',start,stop,i)
         call newatm (i,1056,n-6,1.00d0,n-8,117.4d0,n-9,180.0d0,0)
         call findatm (' H21',start,stop,i)
         call newatm (i,1058,n-5,1.00d0,n-6,120.0d0,n-7,0.0d0,0)
         call findatm (' H22',start,stop,i)
         call newatm (i,1059,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
c
c     cytosine in cytidine residue  (C)
c
      else if (resname .eq. 'C  ') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1078,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1079,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1084,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1080,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1081,n-1,ires)
         call findatm (' N4 ',start,stop,i)
         call oldatm (i,1085,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1082,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1083,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H41',start,stop,i)
         call newatm (i,1086,n-3,1.00d0,n-4,120.0d0,n-5,0.0d0,0)
         call findatm (' H42',start,stop,i)
         call newatm (i,1087,n-4,1.00d0,n-5,120.0d0,n-6,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1088,n-4,1.08d0,n-6,121.6d0,n-7,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1089,n-4,1.08d0,n-5,119.4d0,n-7,180.0d0,0)
c
c     uracil in uridine residue  (U)
c
      else if (resname .eq. 'U  ') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1106,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1107,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1112,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1108,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1109,n-1,ires)
         call findatm (' O4 ',start,stop,i)
         call oldatm (i,1114,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1110,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1111,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H3 ',start,stop,i)
         call newatm (i,1113,n-5,1.00d0,n-7,116.5d0,n-8,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1115,n-3,1.08d0,n-5,120.4d0,n-6,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1116,n-3,1.08d0,n-4,118.6d0,n-6,180.0d0,0)
c
c     adenine in deoxyadenosine residue  (DA)
c
      else if (resname .eq. 'DA ') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1132,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1136,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1135,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1134,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1140,n-1,ires)
         call findatm (' N6 ',start,stop,i)
         call oldatm (i,1142,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1139,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1138,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1137,n-1,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1133,n-1,ires)
         call addbond (n-1,n-7)
         call addbond (n-1,n-10)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1145,n-9,1.08d0,n-8,123.1d0,n-7,180.0d0,0)
         call findatm (' H61',start,stop,i)
         call newatm (i,1143,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
         call findatm (' H62',start,stop,i)
         call newatm (i,1144,n-7,1.00d0,n-8,120.0d0,n-9,0.0d0,0)
         call findatm (' H2 ',start,stop,i)
         call newatm (i,1141,n-6,1.08d0,n-5,115.4d0,n-4,180.0d0,0)
c
c     guanine in deoxyguanosine residue  (DG)
c
      else if (resname .eq. 'DG ') then
         call findatm (' N9 ',start,stop,i)
         call oldatm (i,1161,c1i,ires)
         call findatm (' C8 ',start,stop,i)
         call oldatm (i,1165,n-1,ires)
         call findatm (' N7 ',start,stop,i)
         call oldatm (i,1164,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1163,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1169,n-1,ires)
         call findatm (' O6 ',start,stop,i)
         call oldatm (i,1174,n-1,ires)
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1168,n-2,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1167,n-1,ires)
         call findatm (' N2 ',start,stop,i)
         call oldatm (i,1171,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1166,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1162,n-1,ires)
         call addbond (n-1,n-8)
         call addbond (n-1,n-11)
         call findatm (' H8 ',start,stop,i)
         call newatm (i,1175,n-10,1.08d0,n-9,123.0d0,n-8,180.0d0,0)
         call findatm (' H1 ',start,stop,i)
         call newatm (i,1170,n-6,1.00d0,n-8,117.4d0,n-9,180.0d0,0)
         call findatm (' H21',start,stop,i)
         call newatm (i,1172,n-5,1.00d0,n-6,120.0d0,n-7,0.0d0,0)
         call findatm (' H22',start,stop,i)
         call newatm (i,1173,n-6,1.00d0,n-7,120.0d0,n-8,180.0d0,0)
c
c     cytosine in deoxycytidine residue  (DC)
c
      else if (resname .eq. 'DC ') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1191,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1192,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1197,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1193,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1194,n-1,ires)
         call findatm (' N4 ',start,stop,i)
         call oldatm (i,1198,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1195,n-2,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1196,n-1,ires)
         call addbond (n-1,n-8)
         call findatm (' H41',start,stop,i)
         call newatm (i,1199,n-3,1.00d0,n-4,120.0d0,n-5,0.0d0,0)
         call findatm (' H42',start,stop,i)
         call newatm (i,1200,n-4,1.00d0,n-5,120.0d0,n-6,180.0d0,0)
         call findatm (' H5 ',start,stop,i)
         call newatm (i,1201,n-4,1.08d0,n-6,121.6d0,n-7,180.0d0,0)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1202,n-4,1.08d0,n-5,119.4d0,n-7,180.0d0,0)
c
c     thymine in deoxythymidine residue  (DT)
c
      else if (resname .eq. 'DT ') then
         call findatm (' N1 ',start,stop,i)
         call oldatm (i,1218,c1i,ires)
         call findatm (' C2 ',start,stop,i)
         call oldatm (i,1219,n-1,ires)
         call findatm (' O2 ',start,stop,i)
         call oldatm (i,1224,n-1,ires)
         call findatm (' N3 ',start,stop,i)
         call oldatm (i,1220,n-2,ires)
         call findatm (' C4 ',start,stop,i)
         call oldatm (i,1221,n-1,ires)
         call findatm (' O4 ',start,stop,i)
         call oldatm (i,1226,n-1,ires)
         call findatm (' C5 ',start,stop,i)
         call oldatm (i,1222,n-2,ires)
         call findatm (' C7 ',start,stop,i)
         call oldatm (i,1227,n-1,ires)
         call findatm (' C6 ',start,stop,i)
         call oldatm (i,1223,n-2,ires)
         call addbond (n-1,n-9)
         call findatm (' H3 ',start,stop,i)
         call newatm (i,1225,n-6,1.00d0,n-8,116.8d0,n-9,180.0d0,0)
         call findatm (' H71',start,stop,i)
         call newatm (i,1228,n-3,1.09d0,n-4,109.5d0,n-6,0.0d0,0)
         call findatm (' H72',start,stop,i)
         call newatm (i,1228,n-4,1.09d0,n-5,109.5d0,n-1,109.5d0,1)
         call findatm (' H73',start,stop,i)
         call newatm (i,1228,n-5,1.09d0,n-6,109.5d0,n-2,109.5d0,-1)
         call findatm (' H6 ',start,stop,i)
         call newatm (i,1229,n-5,1.08d0,n-7,119.4d0,n-9,180.0d0,0)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine hetatom  --  coordinates of PDB water and ions  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "hetatom" translates water molecules and ions in Protein Data
c     Bank format to a Cartesian coordinate file and sequence file
c
c
      subroutine hetatom
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'pdb.i'
      integer i
c
c
c     find water molecules and ions in PDB HETATM records
c
      n = n + 1
      i = 0
      dowhile (i .lt. npdb)
         i = i + 1
         if (pdbtyp(i) .eq. 'HETATM') then
            if (resnam(i) .eq. 'HOH') then
               if (atmnam(i) .eq. ' O  ') then
                  call oldatm (i,2001,0,0)
                  if (atmnam(i+1).eq.' H  ' .and.
     &                atmnam(i+2).eq.' H  ') then
                     call oldatm (i+1,2002,n-1,0)
                     call oldatm (i+2,2002,n-2,0)
                     i = i + 2
                  else
                     call newatm (0,2002,n-1,0.96d0,n-2,109.5d0,
     &                               n-3,120.0d0,0)
                     call newatm (0,2002,n-2,0.96d0,n-1,109.5d0,
     &                               n-3,120.0d0,0)
                  end if
               end if
            else if (resnam(i) .eq. 'NA ') then
               call oldatm (i,2003,0,0)
            else if (resnam(i) .eq. 'K  ') then
               call oldatm (i,2004,0,0)
            else if (resnam(i) .eq. 'MG ') then
               call oldatm (i,2005,0,0)
            else if (resnam(i) .eq. 'CA ') then
               call oldatm (i,2006,0,0)
            else if (resnam(i) .eq. 'CL ') then
               call oldatm (i,2007,0,0)
            end if
         end if
      end do
      n = n - 1
      return
      end
