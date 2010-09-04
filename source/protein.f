c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  program protein  --  build a polypeptide from sequence  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "protein" builds the internal and Cartesian coordinates
c     of a polypeptide from amino acid sequence and torsional
c     angle values for the peptide backbone and side chains
c
c
      program protein
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'iounit.i'
      include 'sequen.i'
      include 'titles.i'
      integer i,izmt
      integer ixyz,iseq
      integer natom,mode
      integer freeunit,trimtext
      logical exist,clash
      character*120 seqfile
      character*120 intfile
      character*120 xyzfile
c
c
c     get the name to use for the output structure files
c
      call initial
      call nextarg (filename,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' Enter Name to be used for Output Files :  ',$)
         read (input,20)  filename
   20    format (a120)
      end if
      call basefile (filename)
c
c     get the title line for the output files
c
      write (iout,30)
   30 format (/,' Enter Title :  ',$)
      read (input,40)  title
   40 format (a120)
      ltitle = trimtext (title)
c
c     read the keyfile and force field parameter file
c
      call getkey
      call field
c
c     get the sequence, build a Z-matrix, convert to Cartesians
c
      call getseq
      call prochain
      call connect
      call attach
      call molecule
      call makexyz
c
c     perform a packing calculation for multiple chains
c
      if (nchain .gt. 1) then
         call pauling
         call inertia (2)
      end if
c
c     remove any dummy atoms from Cartesian coordinates
c
      natom = n
      do i = natom, 1, -1
         if (type(i) .eq. 0)  call delete (i)
      end do
c
c     convert to internal and Cartesian coordinates
c
      mode = 0
      call makeint (mode)
      call makexyz
c
c     check for atom pairs with identical coordinates
c
      clash = .false.
      call chkxyz (clash)
c
c     write out a amino acid sequence file
c
      iseq = freeunit ()
      seqfile = filename(1:leng)//'.seq'
      call version (seqfile,'new')
      open (unit=iseq,file=seqfile,status='new')
      call prtseq (iseq)
      close (unit=iseq)
c
c     write out an internal coordinates file
c
      izmt = freeunit ()
      intfile = filename(1:leng)//'.int'
      call version (intfile,'new')
      open (unit=izmt,file=intfile,status='new')
      call prtint (izmt)
      close (unit=izmt)
c
c     write out a Cartesian coordinates file
c
      ixyz = freeunit ()
      xyzfile = filename(1:leng)//'.xyz'
      call version (xyzfile,'new')
      open (unit=ixyz,file=xyzfile,status='new')
      call prtxyz (ixyz)
      close (unit=ixyz)
c
c     perform any final tasks before program exit
c
      call final
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getseq  --  amino acid sequence and angles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getseq" asks the user for the amino acid sequence
c     and torsional angle values needed to define a peptide
c
c
      subroutine getseq
      implicit none
      include 'sizes.i'
      include 'iounit.i'
      include 'phipsi.i'
      include 'resdue.i'
      include 'sequen.i'
      integer i,j,next
      integer length,trimtext
      logical done
      character*1 ucase(26)
      character*1 chir(maxres)
      character*3 name
      character*120 record
      character*120 string
      data ucase  / 'A','B','C','D','E','F','G','H','I','J','K','L',
     &              'M','N','O','P','Q','R','S','T','U','V','W','X',
     &              'Y','Z' /
c
c
c     provide a header to explain the method of sequence input
c
      write (iout,10)
   10 format (/,' Enter One Residue per Line, Free Format: ',
     &           ' 1 or 3 Letter Code, then',
     &        /,' Phi/Psi/Omega (3F), Chi Angles (4F), then',
     &           ' Disulfide Partner if a',
     &        /,' CYS (I), and D/L Chirality as desired (A1)',
     &        //,' The allowed N-Cap Residues are Acetyl=ACE or',
     &           ' Formyl=FOR, and the',
     &        /,' possible C-Cap Residues are N-MethylAmide=NME',
     &           ' or Amide=NH2',
     &        //,' Use Residue=MOL to Begin a New Strand,',
     &           ' Residue=<CR> to End Entry')
c
c     initially, assume that only a single strand is present
c
      nchain = 1
      ichain(1,1) = 1
      chnnam(1) = ' '
c
c     get the amino acid sequence data and dihedral angle values
c
      i = 0
      done = .false.
      do while (.not. done)
         i = i + 1
         phi(i) = 0.0d0
         psi(i) = 0.0d0
         omega(i) = 0.0d0
         do j = 1, 4
            chi(j,i) = 0.0d0
         end do
         disulf(i) = 0
         chir(i) = ' '
         write (iout,20)  i
   20    format (/,' Enter Residue',i4,' :  ',$)
         read (input,30)  record
   30    format (a120)
         call upcase (record)
         next = 1
         call getword (record,name,next)
         length = trimtext (name)
         string = record(next:120)
         read (string,*,err=40,end=40)  phi(i),psi(i),omega(i),
     &                                  (chi(j,i),j=1,4),disulf(i)
   40    continue
         call getword (record,chir(i),next)
c
c     handle special names used for certain amino acids
c
         if (name .eq. 'CYH')  name = 'CYS'
         if (name .eq. 'CSS')  name = 'CYX'
         if (name .eq. 'HIP')  name = 'HIS'
c
c     disulfide bridged residues are cystine instead of cysteine
c
         if (name(1:1).eq.'C' .and. disulf(i).ne.0) then
            length = 3
            name = 'CYX'
         end if
c
c     process and store the current amino acid residue type
c
         if (name .eq. 'MOL') then
            i = i - 1
            ichain(2,nchain) = i
            nchain = nchain + 1
            ichain(1,nchain) = i + 1
         else
            if (name .eq. '   ') then
               done = .true.
               nseq = i - 1
               ichain(2,nchain) = nseq
            else
               seq(i) = amino(maxamino)
               seqtyp(i) = 0
               if (length .eq. 1) then
                  do j = 1, maxamino
                     if (name(1:1) .eq. amino1(j)) then
                        seq(i) = amino(j)
                        seqtyp(i) = j
                     end if
                  end do
               else if (length .eq. 3) then
                  do j = 1, maxamino
                     if (name .eq. amino(j)) then
                        seq(i) = amino(j)
                        seqtyp(i) = j
                     end if
                  end do
               end if
               if (seqtyp(i) .eq. 0) then
                  i = i - 1
                  write (iout,50)  name
   50             format (/,' GETSEQ  --  Amino Acid Type ',a3,
     &                       ' is Not Supported')
               end if
            end if
         end if
      end do
c
c     set chain identifiers if multiple chains are present
c
      if (nchain .gt. 1) then
         do i = 1, nchain
            chnnam(i) = ucase(i)
         end do
      end if
c
c     set default values for the phi-psi-omega-chi angles;
c     use extended values if no phi-psi values were given
c
      do i = 1, nseq
         if (phi(i).eq.0.0d0 .and. psi(i).eq.0.0d0) then
            phi(i) = -135.0d0
            psi(i) = 135.0d0
         end if
         if (omega(i) .eq. 0.0d0) then
            omega(i) = 180.0d0
         end if
         if (chi(1,i) .eq. 0.0d0) then
            do j = 1, 4
               chi(j,i) = 180.0d0
               if (seq(i) .eq. 'PRO')  chi(j,i) = 0.0d0
               if (seq(i) .eq. 'PCA')  chi(j,i) = 0.0d0
            end do
            if (seq(i) .eq. 'PHE')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'TYR')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'TRP')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HIS')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HID')  chi(2,i) = 90.0d0
            if (seq(i) .eq. 'HIE')  chi(2,i) = 90.0d0
         end if
c
c     check for the presence of any disulfide bonds
c
         if (disulf(i) .ne. 0) then
            if (seq(i) .ne. 'CYX') then
               write (iout,60)  i
   60          format (' GETSEQ  --  Error in Disulfide Bond',
     &                    ' at Residue',i5)
            end if
            if (i.lt.disulf(i) .and. disulf(disulf(i)).ne.i) then
               write (iout,70)  i,disulf(i)
   70          format (' GETSEQ  --  Error in Disulfide Bond',
     &                    ' at Residue',i5,' or',i5)
            end if
         end if
c
c     check the D/L chirality of the residues
c
         if (chir(i) .eq. 'D') then
            chiral(i) = -1
         else
            chiral(i) = 1
         end if
      end do
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine prochain  --  build polypeptide backbone  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "prochain" builds up the internal coordinates for an amino
c     acid sequence from the phi, psi, omega and chi values
c
c
      subroutine prochain
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'iounit.i'
      include 'phipsi.i'
      include 'resdue.i'
      include 'sequen.i'
      integer i,k,m,next
      integer ntyp(maxamino),catyp(maxamino)
      integer ctyp(maxamino),hntyp(maxamino)
      integer otyp(maxamino),hatyp(maxamino)
      integer nntyp(maxamino),cantyp(maxamino)
      integer cntyp(maxamino),hnntyp(maxamino)
      integer ontyp(maxamino),hantyp(maxamino)
      integer nctyp(maxamino),cactyp(maxamino)
      integer cctyp(maxamino),hnctyp(maxamino)
      integer octyp(maxamino),hactyp(maxamino)
      integer ni(maxres),cai(maxres),ci(maxres)
      logical single,cyclic
      character*1 answer
      character*3 resname
      character*120 record
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
c     determine whether the peptide chain is cyclic
c
      cyclic = .false.
      write (iout,10)
   10 format (/,' Cyclize the Polypeptide Chain [N] :  ',$)
      read (input,20)  record
   20 format (a120)
      next = 1
      call gettext (record,answer,next)
      call upcase (answer)
      if (answer .eq. 'Y')  cyclic = .true.
c
c     initialize the atom counter to the first atom
c
      n = 1
c
c     set atom counter and the first residue number and type
c
      do m = 1, nchain
         single = .false.
         if (ichain(1,m) .eq. ichain(2,m))  single = .true.
         i = ichain(1,m)
         k = seqtyp(i)
         resname = amino(k)
c
c     build the first residue for a cyclic peptide
c
         if (cyclic) then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (ntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (catyp(k),1.46d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               call zatom (ctyp(k),1.51d0,110.7d0,0.0d0,
     &                     cai(i),ni(i),0,0)
            else
               ni(i) = n
               call zatom (ntyp(k),30.0d0,150.0d0,180.0d0,n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (catyp(k),1.46d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               call zatom (ctyp(k),1.51d0,110.7d0,180.0d0,
     &                     cai(i),ni(i),n-3,0)
            end if
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue as an N-terminal formyl group
c
         else if (resname .eq. 'FOR') then
            if (m .eq. 1) then
               ci(i) = n
               call zatom (cntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               ni(i) = n
               call zatom (ontyp(k),1.22d0,0.0d0,0.0d0,n-1,0,0,0)
               cai(i) = n
               call zatom (hantyp(k),1.12d0,120.0d0,0.0d0,n-2,n-1,0,0)
            else
               ci(i) = n
               call zatom (cntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               ni(i) = n
               call zatom (ontyp(k),1.22d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               cai(i) = n
               call zatom (hantyp(k),1.12d0,120.0d0,0.0d0,n-2,n-1,n-3,0)
            end if
            psi(i) = 180.0d0
c
c     build the first residue as an N-terminal acetyl group
c
         else if (resname .eq. 'ACE') then
            if (m .eq. 1) then
               cai(i) = n
               call zatom (cantyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               ci(i) = n
               call zatom (cntyp(k),1.51d0,0.0d0,0.0d0,n-1,0,0,0)
               call zatom (ontyp(k),1.22d0,122.5d0,0.0d0,n-1,n-2,0,0)
            else
               cai(i) = n
               call zatom (cantyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               ci(i) = n
               call zatom (cntyp(k),1.51d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (ontyp(k),1.22d0,122.5d0,0.0d0,
     &                     n-1,n-2,n-3,0)
            end if
            ni(i) = n
            call zatom (hantyp(k),1.11d0,107.9d0,0.0d0,n-3,n-2,n-1,0)
            call zatom (hantyp(k),1.11d0,107.9d0,109.4d0,n-4,n-3,n-1,1)
            call zatom (hantyp(k),1.11d0,107.9d0,109.4d0,n-5,n-4,n-2,-1)
            psi(i) = 180.0d0
c
c     build the first residue as a proline
c
         else if (resname .eq. 'PRO') then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,0.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hnntyp(k),1.02d0,109.5d0,-120.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue as a pyroglutamic acid
c
         else if (resname .eq. 'PCA') then
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,-60.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the first residue for all other standard amino acids
c
         else
            if (m .eq. 1) then
               ni(i) = n
               call zatom (nntyp(k),0.0d0,0.0d0,0.0d0,0,0,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,0.0d0,0.0d0,ni(i),0,0,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,0.0d0,
     &                        cai(i),ni(i),0,0)
               end if
            else
               ni(i) = n
               call zatom (nntyp(k),30.0d0,150.0d0,180.0d0,
     &                     n-1,n-2,n-3,0)
               call zatom (-2,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
               cai(i) = n
               call zatom (cantyp(k),1.50d0,150.0d0,180.0d0,
     &                     ni(i),n-2,n-3,0)
               ci(i) = n
               if (single) then
                  call zatom (cctyp(k),1.51d0,111.6d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               else
                  call zatom (cntyp(k),1.51d0,110.7d0,180.0d0,
     &                        cai(i),ni(i),n-3,0)
               end if
            end if
            if (single) then
               call zatom (octyp(k),1.25d0,117.0d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            else
               call zatom (ontyp(k),1.22d0,122.5d0,psi(1)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
            end if
            call zatom (hnntyp(k),1.02d0,109.5d0,phi(i),
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hnntyp(k),1.02d0,109.5d0,108.0d0,
     &                  ni(i),cai(i),n-1,1)
            call zatom (hnntyp(k),1.02d0,109.5d0,108.0d0,
     &                  ni(i),cai(i),n-2,-1)
            call zatom (hantyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
         end if
c
c     build atoms for residues in the middle of the chain
c
         do i = ichain(1,m)+1, ichain(2,m)-1
            k = seqtyp(i)
            resname = amino(k)
            ni(i) = n
            call zatom (ntyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            cai(i) = n
            call zatom (catyp(k),1.46d0,121.0d0,omega(i-1),
     &                  ni(i),ci(i-1),cai(i-1),0)
            ci(i) = n
            call zatom (ctyp(k),1.51d0,111.6d0,phi(i),
     &                  cai(i),ni(i),ci(i-1),0)
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
         end do
c
c     set the number and type of the last residue
c
         i = ichain(2,m)
         k = seqtyp(i)
         resname = amino(k)
c
c     build the last residue for a cyclic peptide
c
         if (cyclic) then
            ni(i) = n
            call zatom (ntyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            cai(i) = n
            call zatom (catyp(k),1.46d0,121.0d0,omega(i-1),
     &                  ni(i),ci(i-1),cai(i-1),0)
            ci(i) = n
            call zatom (ctyp(k),1.51d0,111.6d0,phi(i),
     &                  cai(i),ni(i),ci(i-1),0)
            call zatom (-1,0.0d0,0.0d0,0.0d0,ni(1),ci(i),0,0)
            call zatom (otyp(k),1.22d0,122.5d0,psi(i)-180.0d0,
     &                  ci(i),cai(i),ni(i),0)
            call zatom (hntyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                  ni(i),cai(i),ci(i),0)
            call zatom (hatyp(k),1.11d0,109.5d0,107.9d0,
     &                  cai(i),ni(i),ci(i),-chiral(i))
            call proside (resname,i,cai(i),ni(i),ci(i))
c
c     build the last residue as a C-terminal amide
c
         else if (resname .eq. 'NH2') then
            call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            call zatom (hnctyp(k),1.02d0,119.0d0,0.0d0,
     &                  n-1,ci(i-1),cai(i-1),0)
            call zatom (hnctyp(k),1.02d0,119.0d0,180.0d0,
     &                  n-2,ci(i-1),cai(i-1),0)
c
c     build the last residue as a C-terminal N-methylamide
c
         else if (resname .eq. 'NME') then
            call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                  ci(i-1),cai(i-1),ni(i-1),0)
            call zatom (cactyp(k),1.46d0,121.0d0,180.0d0,
     &                  n-1,ci(i-1),cai(i-1),0)
            call zatom (hnctyp(k),1.02d0,118.0d0,121.0d0,
     &                  n-2,ci(i-1),n-1,1)
            call zatom (hactyp(k),1.11d0,109.5d0,180.0d0,
     &                  n-2,n-3,ci(i-1),0)
            call zatom (hactyp(k),1.11d0,109.5d0,109.5d0,
     &                  n-3,n-4,n-1,1)
            call zatom (hactyp(k),1.11d0,109.5d0,109.5d0,
     &                  n-4,n-5,n-2,-1)
c
c     build the last residue for all other standard amino acids
c
         else
            if (.not. single) then
               ni(i) = n
               call zatom (nctyp(k),1.34d0,112.7d0,psi(i-1),
     &                     ci(i-1),cai(i-1),ni(i-1),0)
               cai(i) = n
               call zatom (cactyp(k),1.46d0,121.0d0,omega(i-1),
     &                     ni(i),ci(i-1),cai(i-1),0)
               ci(i) = n
               call zatom (cctyp(k),1.51d0,111.6d0,phi(i),
     &                     cai(i),ni(i),ci(i-1),0)
               call zatom (octyp(k),1.25d0,117.0d0,psi(i)-180.0d0,
     &                     ci(i),cai(i),ni(i),0)
               call zatom (hnctyp(k),1.02d0,121.0d0,phi(i)-180.0d0,
     &                     ni(i),cai(i),ci(i),0)
               call zatom (hactyp(k),1.11d0,109.5d0,107.9d0,
     &                     cai(i),ni(i),ci(i),-chiral(i))
               call proside (resname,i,cai(i),ni(i),ci(i))
            end if
            call zatom (octyp(k),1.25d0,117.0d0,psi(i),
     &                  ci(i),cai(i),ni(i),0)
         end if
      end do
c
c     finally, set the total number of atoms
c
      n = n - 1
      return
      end
c
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine proside  --  build amino acid side chain  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "proside" builds the side chain for a single amino acid
c     residue in terms of internal coordinates
c
c     resname   3-letter name of current amino acid residue
c     i         number of the current amino acid residue
c     cai       atom number of alpha carbon in residue i
c     ni        atom number of amide nitrogen in residue i
c     ci        atom number of carbonyl carbon in residue i
c
c
      subroutine proside (resname,i,cai,ni,ci)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'phipsi.i'
      include 'sequen.i'
      integer i,cai,ni,ci
      character*3 resname
c
c
c     glycine residue  (GLY)
c
      if (resname .eq. 'GLY') then
         if (i .eq. 1) then
            call zatom (355,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         else if (i .eq. nseq) then
            call zatom (506,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         else
            call zatom (6,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         end if
c
c     alanine residue  (ALA)
c
      else if (resname .eq. 'ALA') then
         call zatom (13,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (14,1.11d0,109.4d0,chi(1,i),n-1,cai,ni,0)
         call zatom (14,1.11d0,109.4d0,109.4d0,n-2,cai,n-1,1)
         call zatom (14,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,-1)
c
c     valine residue  (VAL)
c
      else if (resname .eq. 'VAL') then
         call zatom (21,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (23,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (25,1.54d0,109.5d0,109.5d0,n-2,cai,n-1,-1)
         call zatom (22,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,1)
         call zatom (24,1.11d0,109.4d0,180.0d0,n-3,n-4,cai,0)
         call zatom (24,1.11d0,109.4d0,109.4d0,n-4,n-5,n-1,1)
         call zatom (24,1.11d0,109.4d0,109.4d0,n-5,n-6,n-2,-1)
         call zatom (26,1.11d0,109.4d0,180.0d0,n-5,n-7,cai,0)
         call zatom (26,1.11d0,109.4d0,109.4d0,n-6,n-8,n-1,1)
         call zatom (26,1.11d0,109.4d0,109.4d0,n-7,n-9,n-2,-1)
c
c     leucine residue  (LEU)
c
      else if (resname .eq. 'LEU') then
         call zatom (33,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (35,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (37,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (39,1.54d0,109.5d0,109.4d0,n-2,n-3,n-1,-1)
         call zatom (34,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (34,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (36,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (38,1.11d0,109.4d0,180.0d0,n-5,n-6,n-7,0)
         call zatom (38,1.11d0,109.4d0,109.4d0,n-6,n-7,n-1,1)
         call zatom (38,1.11d0,109.4d0,109.4d0,n-7,n-8,n-2,-1)
         call zatom (40,1.11d0,109.4d0,180.0d0,n-7,n-9,n-10,0)
         call zatom (40,1.11d0,109.4d0,109.4d0,n-8,n-10,n-1,1)
         call zatom (40,1.11d0,109.4d0,109.4d0,n-9,n-11,n-2,-1)
c
c     isoleucine residue  (ILE)
c
      else if (resname .eq. 'ILE') then
         call zatom (47,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (49,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (51,1.54d0,109.5d0,109.5d0,n-2,cai,n-1,1)
         call zatom (53,1.54d0,109.5d0,chi(2,i),n-2,n-3,cai,0)
         call zatom (48,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,-1)
         call zatom (50,1.11d0,109.4d0,109.4d0,n-4,n-5,n-2,1)
         call zatom (50,1.11d0,109.4d0,109.4d0,n-5,n-6,n-3,-1)
         call zatom (52,1.11d0,110.0d0,180.0d0,n-5,n-7,n-6,0)
         call zatom (52,1.11d0,110.0d0,109.0d0,n-6,n-8,n-1,1)
         call zatom (52,1.11d0,110.0d0,109.0d0,n-7,n-9,n-2,-1)
         call zatom (54,1.11d0,110.0d0,180.0d0,n-7,n-9,n-10,0)
         call zatom (54,1.11d0,110.0d0,109.0d0,n-8,n-10,n-1,1)
         call zatom (54,1.11d0,110.0d0,109.0d0,n-9,n-11,n-2,-1)
c
c     serine residue  (SER)
c
      else if (resname .eq. 'SER') then
         call zatom (61,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (63,1.41d0,107.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (62,1.11d0,109.4d0,106.7d0,n-2,cai,n-1,1)
         call zatom (62,1.11d0,109.4d0,106.7d0,n-3,cai,n-2,-1)
         call zatom (64,0.94d0,106.9d0,chi(2,i),n-3,n-4,cai,0)
c
c     threonine residue  (THR)
c
      else if (resname .eq. 'THR') then
         call zatom (71,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (73,1.41d0,107.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (75,1.54d0,109.5d0,107.7d0,n-2,cai,n-1,1)
         call zatom (72,1.11d0,109.4d0,106.7d0,n-3,cai,n-2,-1)
         call zatom (74,0.94d0,106.9d0,chi(2,i),n-3,n-4,cai,0)
         call zatom (76,1.11d0,110.0d0,180.0d0,n-3,n-5,cai,0)
         call zatom (76,1.11d0,110.0d0,109.0d0,n-4,n-6,n-1,1)
         call zatom (76,1.11d0,110.0d0,109.0d0,n-5,n-7,n-2,-1)
c
c     cysteine residue  (CYS)
c
      else if (resname .eq. 'CYS') then
         call zatom (83,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (85,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (84,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
         call zatom (84,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
         call zatom (86,1.34d0,96.0d0,chi(2,i),n-3,n-4,cai,0)
c
c     cystine residue  (CYX)
c
      else if (resname .eq. 'CYX') then
         if (disulf(i) .gt. i) then
            call zatom (93,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
            call zatom (95,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
            call zatom (94,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
            call zatom (94,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
            disulf(i) = n - 3
         else if (disulf(i) .lt. i) then
            call zatom (93,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
            call zatom (95,1.82d0,109.0d0,chi(1,i),n-1,cai,ni,0)
            call zatom (94,1.11d0,109.4d0,112.0d0,n-2,cai,n-1,1)
            call zatom (94,1.11d0,109.4d0,112.0d0,n-3,cai,n-2,-1)
            call zatom (-1,0.0d0,0.0d0,0.0d0,disulf(disulf(i)),n-3,0,0)
         end if
c
c     proline residue  (PRO)
c
      else if (resname .eq. 'PRO') then
         call zatom (101,1.54d0,107.0d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (103,1.54d0,107.0d0,chi(1,i),n-1,cai,ni,0)
         if (i .eq. 1) then
            call zatom (410,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         else
            call zatom (105,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         end if
         call zatom (-1,0.0d0,0.0d0,0.0d0,ni,n-1,0,0)
         call zatom (102,1.11d0,109.4d0,109.4d0,n-3,cai,n-2,1)
         call zatom (102,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,-1)
         call zatom (104,1.11d0,109.4d0,109.4d0,n-4,n-5,n-3,1)
         call zatom (104,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,-1)
         if (i .eq. 1) then
            call zatom (411,1.11d0,109.4d0,109.4d0,n-5,n-6,ni,1)
            call zatom (411,1.11d0,109.4d0,109.4d0,n-6,n-7,ni,-1)
         else
            call zatom (106,1.11d0,109.4d0,109.4d0,n-5,n-6,ni,1)
            call zatom (106,1.11d0,109.4d0,109.4d0,n-6,n-7,ni,-1)
         end if
c
c     phenylalanine residue  (PHE)
c
      else if (resname .eq. 'PHE') then
         call zatom (113,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (115,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (116,1.39d0,120.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (116,1.39d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (118,1.39d0,120.0d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (118,1.39d0,120.0d0,180.0d0,n-2,n-4,n-5,0)
         call zatom (120,1.39d0,120.0d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (114,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,1)
         call zatom (114,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,-1)
         call zatom (117,1.10d0,120.0d0,120.0d0,n-7,n-8,n-5,1)
         call zatom (117,1.10d0,120.0d0,120.0d0,n-7,n-9,n-5,1)
         call zatom (119,1.10d0,120.0d0,120.0d0,n-7,n-9,n-5,1)
         call zatom (119,1.10d0,120.0d0,120.0d0,n-7,n-9,n-6,1)
         call zatom (121,1.10d0,120.0d0,120.0d0,n-7,n-9,n-8,1)
c
c     tyrosine residue  (TYR)
c
      else if (resname .eq. 'TYR') then
         call zatom (128,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (130,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (131,1.39d0,120.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (131,1.39d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (133,1.39d0,120.0d0,180.0d0,n-2,n-3,n-4,0)
         call zatom (133,1.39d0,120.0d0,180.0d0,n-2,n-4,n-5,0)
         call zatom (135,1.39d0,120.0d0,0.0d0,n-2,n-4,n-5,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (136,1.36d0,120.0d0,120.0d0,n-1,n-2,n-3,1)
         call zatom (129,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,1)
         call zatom (129,1.11d0,109.4d0,109.4d0,n-9,cai,n-8,-1)
         call zatom (132,1.10d0,120.0d0,120.0d0,n-8,n-9,n-6,1)
         call zatom (132,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (134,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (134,1.10d0,120.0d0,120.0d0,n-8,n-10,n-7,1)
         call zatom (137,0.97d0,108.0d0,0.0d0,n-7,n-8,n-9,0)
c
c     tryptophan residue  (TRP)
c
      else if (resname .eq. 'TRP') then
         call zatom (144,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (146,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (147,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (149,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (150,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (152,1.35d0,108.0d0,0.0d0,n-1,n-3,n-4,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-3,n-1,0,0)
         call zatom (153,1.35d0,120.0d0,180.0d0,n-3,n-1,n-2,0)
         call zatom (155,1.35d0,120.0d0,0.0d0,n-2,n-4,n-1,0)
         call zatom (157,1.35d0,120.0d0,0.0d0,n-2,n-5,n-3,0)
         call zatom (159,1.35d0,120.0d0,0.0d0,n-2,n-4,n-6,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,n-2,n-1,0,0)
         call zatom (145,1.11d0,109.4d0,109.4d0,n-10,cai,n-9,1)
         call zatom (145,1.11d0,109.4d0,109.4d0,n-11,cai,n-10,-1)
         call zatom (148,1.10d0,126.0d0,126.0d0,n-10,n-11,n-8,1)
         call zatom (151,1.05d0,126.0d0,126.0d0,n-9,n-11,n-8,1)
         call zatom (154,1.10d0,120.0d0,120.0d0,n-8,n-11,n-6,1)
         call zatom (156,1.10d0,120.0d0,120.0d0,n-8,n-10,n-6,1)
         call zatom (158,1.10d0,120.0d0,120.0d0,n-8,n-10,n-7,1)
         call zatom (160,1.10d0,120.0d0,120.0d0,n-8,n-10,n-9,1)
c
c     histidine (HD and HE) residue  (HIS)
c
      else if (resname .eq. 'HIS') then
         call zatom (167,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (169,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (170,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (172,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (174,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (176,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (168,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (168,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (171,1.02d0,126.0d0,0.0d0,n-6,n-7,n-8,0)
         call zatom (173,1.10d0,126.0d0,126.0d0,n-6,n-8,n-4,1)
         call zatom (175,1.10d0,126.0d0,126.0d0,n-6,n-8,n-5,1)
         call zatom (177,1.02d0,126.0d0,126.0d0,n-6,n-8,n-7,1)
c
c     histidine (HD only) residue  (HID)
c
      else if (resname .eq. 'HID') then
         call zatom (184,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (186,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (187,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (189,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (191,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (193,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (185,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (185,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (188,1.02d0,126.0d0,0.0d0,n-6,n-7,n-8,0)
         call zatom (190,1.10d0,126.0d0,126.0d0,n-6,n-8,n-4,1)
         call zatom (192,1.10d0,126.0d0,126.0d0,n-6,n-8,n-5,1)
c
c     histidine (HE only) residue  (HIE)
c
      else if (resname .eq. 'HIE') then
         call zatom (200,1.54d0,109.5d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (202,1.50d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (203,1.35d0,126.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (204,1.35d0,126.0d0,108.0d0,n-2,n-3,n-1,1)
         call zatom (206,1.35d0,108.0d0,0.0d0,n-2,n-3,n-1,0)
         call zatom (208,1.35d0,108.0d0,0.0d0,n-2,n-4,n-3,0)
         call zatom (-1,0.0d0,0.0d0,.00d0,n-2,n-1,0,0)
         call zatom (201,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,1)
         call zatom (201,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,-1)
         call zatom (205,1.10d0,126.0d0,126.0d0,n-5,n-7,n-3,1)
         call zatom (207,1.10d0,126.0d0,126.0d0,n-5,n-7,n-4,1)
         call zatom (209,1.02d0,126.0d0,126.0d0,n-5,n-7,n-6,1)
c
c     aspartic acid residue  (ASP)
c
      else if (resname .eq. 'ASP') then
         call zatom (216,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (218,1.51d0,107.8d0,chi(1,i),n-1,cai,ni,0)
         call zatom (219,1.25d0,117.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (219,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (217,1.11d0,109.4d0,107.9d0,n-4,cai,n-3,1)
         call zatom (217,1.11d0,109.4d0,107.9d0,n-5,cai,n-4,-1)
c
c     asparagine residue  (ASN)
c
      else if (resname .eq. 'ASN') then
         call zatom (226,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (228,1.51d0,107.8d0,chi(1,i),n-1,cai,ni,0)
         call zatom (229,1.22d0,122.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (230,1.34d0,112.7d0,124.0d0,n-2,n-3,n-1,1)
         call zatom (227,1.11d0,109.4d0,107.9d0,n-4,cai,n-3,1)
         call zatom (227,1.11d0,109.4d0,107.9d0,n-5,cai,n-4,-1)
         call zatom (231,1.02d0,119.0d0,0.0d0,n-3,n-5,n-6,0)
         call zatom (231,1.02d0,119.0d0,120.0d0,n-4,n-6,n-1,1)
c
c     glutamic acid residue  (GLU)
c
      else if (resname .eq. 'GLU') then
         call zatom (238,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (240,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (242,1.51d0,107.8d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (243,1.25d0,117.0d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (243,1.25d0,117.0d0,126.0d0,n-2,n-3,n-1,1)
         call zatom (239,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (239,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (241,1.11d0,109.4d0,107.9d0,n-6,n-7,n-5,1)
         call zatom (241,1.11d0,109.4d0,107.9d0,n-7,n-8,n-6,-1)
c
c     glutamine residue  (GLN)
c
      else if (resname .eq. 'GLN') then
         call zatom (250,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (252,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (254,1.51d0,107.8d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (255,1.22d0,122.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (256,1.34d0,112.7d0,124.0d0,n-2,n-3,n-1,1)
         call zatom (251,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (251,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (253,1.11d0,109.4d0,107.9d0,n-6,n-7,n-5,1)
         call zatom (253,1.11d0,109.4d0,107.9d0,n-7,n-8,n-6,-1)
         call zatom (257,1.02d0,119.0d0,0.0d0,n-5,n-7,n-8,0)
         call zatom (257,1.02d0,119.0d0,120.0d0,n-6,n-8,n-1,1)
c
c     methionine residue  (MET)
c
      else if (resname .eq. 'MET') then
         call zatom (264,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (266,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (268,1.82d0,109.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (269,1.82d0,96.3d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (265,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (265,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (267,1.11d0,109.4d0,112.0d0,n-5,n-6,n-4,1)
         call zatom (267,1.11d0,109.4d0,112.0d0,n-6,n-7,n-5,-1)
         call zatom (270,1.11d0,112.0d0,180.0d0,n-5,n-6,n-7,0)
         call zatom (270,1.11d0,112.0d0,109.4d0,n-6,n-7,n-1,1)
         call zatom (270,1.11d0,112.0d0,109.4d0,n-7,n-8,n-2,-1)
c
c     lysine residue  (LYS)
c
      else if (resname .eq. 'LYS') then
         call zatom (277,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (279,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (281,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (283,1.54d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (285,1.50d0,109.5d0,chi(4,i),n-1,n-2,n-3,0)
         call zatom (278,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,1)
         call zatom (278,1.11d0,109.4d0,109.4d0,n-6,cai,n-5,-1)
         call zatom (280,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,1)
         call zatom (280,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,-1)
         call zatom (282,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,1)
         call zatom (282,1.11d0,109.4d0,109.4d0,n-8,n-9,n-7,-1)
         call zatom (284,1.11d0,109.4d0,108.8d0,n-8,n-9,n-7,1)
         call zatom (284,1.11d0,109.4d0,108.8d0,n-9,n-10,n-8,-1)
         call zatom (286,1.02d0,109.5d0,180.0d0,n-9,n-10,n-11,0)
         call zatom (286,1.02d0,109.5d0,109.5d0,n-10,n-11,n-1,1)
         call zatom (286,1.02d0,109.5d0,109.5d0,n-11,n-12,n-2,-1)
c
c     arginine residue  (ARG)
c
      else if (resname .eq. 'ARG') then
         call zatom (293,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (295,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (297,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (299,1.45d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (301,1.35d0,120.0d0,chi(4,i),n-1,n-2,n-3,0)
         call zatom (302,1.35d0,120.0d0,180.0d0,n-1,n-2,n-3,0)
         call zatom (302,1.35d0,120.0d0,120.0d0,n-2,n-3,n-1,1)
         call zatom (294,1.11d0,109.4d0,109.4d0,n-7,cai,n-6,1)
         call zatom (294,1.11d0,109.4d0,109.4d0,n-8,cai,n-7,-1)
         call zatom (296,1.11d0,109.4d0,109.4d0,n-8,n-9,n-7,1)
         call zatom (296,1.11d0,109.4d0,109.4d0,n-9,n-10,n-8,-1)
         call zatom (298,1.11d0,109.4d0,109.4d0,n-9,n-10,n-8,1)
         call zatom (298,1.11d0,109.4d0,109.4d0,n-10,n-11,n-9,-1)
         call zatom (300,1.02d0,120.0d0,120.0d0,n-10,n-11,n-9,1)
         call zatom (303,1.02d0,120.0d0,180.0d0,n-9,n-10,n-11,0)
         call zatom (303,1.02d0,120.0d0,120.0d0,n-10,n-11,n-1,1)
         call zatom (303,1.02d0,120.0d0,180.0d0,n-10,n-12,n-13,0)
         call zatom (303,1.02d0,120.0d0,120.0d0,n-11,n-13,n-1,1)
c
c     ornithine residue  (ORN)
c
      else if (resname .eq. 'ORN') then
         call zatom (310,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (312,1.54d0,109.5d0,chi(1,i),n-1,cai,ni,0)
         call zatom (314,1.54d0,109.5d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (316,1.50d0,109.5d0,chi(3,i),n-1,n-2,n-3,0)
         call zatom (311,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (311,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (313,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (313,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,-1)
         call zatom (315,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,1)
         call zatom (315,1.11d0,109.4d0,109.4d0,n-7,n-8,n-6,-1)
         call zatom (317,1.02d0,109.5d0,180.0d0,n-7,n-8,n-9,0)
         call zatom (317,1.02d0,109.5d0,109.5d0,n-8,n-9,n-1,1)
         call zatom (317,1.02d0,109.5d0,109.5d0,n-9,n-10,n-2,-1)
c
c     methylalanine residue  (AIB)
c
      else if (resname .eq. 'AIB') then
         call zatom (323,1.54d0,109.5d0,107.8d0,cai,ni,ci,-chiral(i))
         call zatom (323,1.54d0,109.5d0,107.8d0,cai,ni,ci,chiral(i))
         call zatom (324,1.11d0,109.4d0,chi(1,i),n-2,cai,ni,0)
         call zatom (324,1.11d0,109.4d0,109.4d0,n-3,cai,n-1,1)
         call zatom (324,1.11d0,109.4d0,109.4d0,n-4,cai,n-2,-1)
         call zatom (324,1.11d0,109.4d0,chi(1,i),n-4,cai,ni,0)
         call zatom (324,1.11d0,109.4d0,109.4d0,n-5,cai,n-1,1)
         call zatom (324,1.11d0,109.4d0,109.4d0,n-6,cai,n-2,-1)
c
c     pyroglutamic acid residue  (PCA)
c
      else if (resname .eq. 'PCA') then
         call zatom (331,1.54d0,107.0d0,109.5d0,cai,ni,ci,chiral(i))
         call zatom (333,1.54d0,107.0d0,chi(1,i),n-1,cai,ni,0)
         call zatom (335,1.54d0,107.0d0,chi(2,i),n-1,n-2,cai,0)
         call zatom (-1,0.0d0,0.0d0,0.0d0,ni,n-1,0,0)
         call zatom (336,1.22d0,126.0d0,126.0d0,n-1,ni,n-2,1)
         call zatom (332,1.11d0,109.4d0,109.4d0,n-4,cai,n-3,1)
         call zatom (332,1.11d0,109.4d0,109.4d0,n-5,cai,n-4,-1)
         call zatom (334,1.11d0,109.4d0,109.4d0,n-5,n-6,n-4,1)
         call zatom (334,1.11d0,109.4d0,109.4d0,n-6,n-7,n-5,-1)
c
c     unknown residue  (UNK)
c
      else if (resname .eq. 'UNK') then
         if (i .eq. 1) then
            call zatom (355,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         else if (i .eq. nseq) then
            call zatom (506,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         else
            call zatom (6,1.11d0,109.5d0,107.9d0,cai,ni,ci,chiral(i))
         end if
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine pauling  --  pack multiple polypeptide chains  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "pauling" uses a rigid body optimization to approximately
c     pack multiple polypeptide chains
c
c
      subroutine pauling
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'katoms.i'
      include 'kgeoms.i'
      include 'molcul.i'
      include 'output.i'
      include 'potent.i'
      include 'rigid.i'
      include 'usage.i'
      integer i,j,k,nvar
      real*8 minimum,grdmin
      real*8 pauling1
      real*8 xx(maxopt)
      external pauling1,optsave
c
c
c     set all atoms to be active during energy evaluations
c
      nuse = n
      do i = 1, n
         use(i) = .true.
      end do
c
c     only geometric restraints will by used in optimization
c
      call potoff
      use_geom = .true.
c
c     set the default values for the restraint variables
c
      npfix = 0
      ndfix = 0
      ntfix = 0
      ngfix = 0
      nchir = 0
      use_basin = .true.
      depth = 3.0d0
      width = 1.5d0
      use_wall = .false.
c
c     assign each chain to a separate molecule-based group
c
      use_group = .true.
      ngrp = nmol
      do i = 1, ngrp
         igrp(1,i) = imol(1,i)
         igrp(2,i) = imol(2,i)
         do j = igrp(1,i), igrp(2,i)
            kgrp(j) = kmol(j)
            grplist(kgrp(j)) = i
         end do
      end do
      do i = 0, ngrp
         do j = 0, ngrp
            wgrp(j,i) = 1.0d0
         end do
         wgrp(i,i) = 1.0d0
      end do
c
c     assume unit mass for each atom and set group masses
c
      do i = 1, n
         mass(i) = 1.0d0
      end do
      do i = 1, ngrp
         grpmass(i) = dble(igrp(2,i)-igrp(1,i)+1)
      end do
c
c     set pairwise restraints between the centers of chains
c
      do i = 1, ngrp-1
         do j = i+1, ngrp
            ngfix = ngfix + 1
            igfix(1,ngfix) = i
            igfix(2,ngfix) = j
            gfix(1,ngfix) = 1.0d0
            gfix(2,ngfix) = 11.0d0 * dble(j-i)
            gfix(3,ngfix) = 11.0d0 * dble(j-i)
         end do
      end do
c
c     set position restraints on alpha carbons of each chain
c
      do i = 1, n
         if (atmnum(type(i)) .eq. 6) then
            do j = 1, n12(i)
               if (atmnum(type(i12(j,i))) .eq. 7) then
                  do k = 1, n13(i)
                     if (atmnum(type(i13(k,i))) .eq. 8) then
                        npfix = npfix + 1
                        ipfix(npfix) = i
                        kpfix(1,npfix) = 1
                        kpfix(2,npfix) = 1
                        kpfix(3,npfix) = 0
                        xpfix(npfix) = 11.0d0 * dble(grplist(i)-1)
                        ypfix(npfix) = 0.0d0
                        zpfix(npfix) = 0.0d0
                        pfix(1,npfix) = 1.0d0
                        pfix(2,npfix) = 0.0d0
                        goto 10
                     end if
                  end do
               end if
            end do
         end if
   10    continue
      end do
c
c     get rigid body reference coordinates for each chain
c
      call orient
c
c     transfer rigid body coordinates to optimization parameters
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            xx(nvar) = rbc(j,i)
         end do
      end do
c
c     make the call to the optimization routine
c
      iprint = 0
      iwrite = 0
      grdmin = 0.1d0
      coordtype = 'NONE'
      call ocvm (nvar,xx,minimum,grdmin,pauling1,optsave)
c
c     transfer optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     convert from rigid body to Cartesian coordinates
c
      call rigidxyz
      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  function pauling1  --  energy and gradient for pauling  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "pauling1" is a service routine that computes the energy
c     and gradient for optimally conditioned variable metric
c     optimization of rigid bodies
c
c
      function pauling1 (xx,g)
      implicit none
      include 'sizes.i'
      include 'group.i'
      include 'math.i'
      include 'rigid.i'
      integer i,j,nvar
      real*8 pauling1,e
      real*8 xx(maxopt)
      real*8 g(maxopt)
      real*8 derivs(6,maxgrp)
c
c
c     translate optimization parameters to rigid body coordinates
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            rbc(j,i) = xx(nvar)
         end do
      end do
c
c     compute and store the energy and gradient
c
      call rigidxyz
      call gradrgd (e,derivs)
      pauling1 = e
c
c     translate rigid body gradient to optimization gradient
c
      nvar = 0
      do i = 1, ngrp
         do j = 1, 6
            nvar = nvar + 1
            g(nvar) = derivs(j,i)
         end do
      end do
      return
      end
