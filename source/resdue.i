c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  resdue.i  --  biopolymer residue names and biotype codes  ##
c     ##                                                            ##
c     ################################################################
c
c
c     ntyp     biotypes for mid-chain peptide backbone N atoms
c     catyp    biotypes for mid-chain peptide backbone CA atoms
c     ctyp     biotypes for mid-chain peptide backbone C atoms
c     hntyp    biotypes for mid-chain peptide backbone HN atoms
c     otyp     biotypes for mid-chain peptide backbone O atoms
c     hatyp    biotypes for mid-chain peptide backbone HA atoms
c     cbtyp    biotypes for mid-chain peptide backbone CB atoms
c     nntyp    biotypes for N-terminal peptide backbone N atoms
c     cantyp   biotypes for N-terminal peptide backbone CA atoms
c     cntyp    biotypes for N-terminal peptide backbone C atoms
c     hnntyp   biotypes for N-terminal peptide backbone HN atoms
c     ontyp    biotypes for N-terminal peptide backbone O atoms
c     hantyp   biotypes for N-terminal peptide backbone HA atoms
c     nctyp    biotypes for C-terminal peptide backbone N atoms
c     cactyp   biotypes for C-terminal peptide backbone CA atoms
c     cctyp    biotypes for C-terminal peptide backbone C atoms
c     hnctyp   biotypes for C-terminal peptide backbone HN atoms
c     octyp    biotypes for C-terminal peptide backbone O atoms
c     hactyp   biotypes for C-terminal peptide backbone HA atoms
c     o5typ    biotypes for nucleotide backbone and sugar O5' atoms
c     c5typ    biotypes for nucleotide backbone and sugar C5' atoms
c     h51typ   biotypes for nucleotide backbone and sugar H5' atoms
c     h52typ   biotypes for nucleotide backbone and sugar H5'' atoms
c     c4typ    biotypes for nucleotide backbone and sugar C4' atoms
c     h4typ    biotypes for nucleotide backbone and sugar H4' atoms
c     o4typ    biotypes for nucleotide backbone and sugar O4' atoms
c     c1typ    biotypes for nucleotide backbone and sugar C1' atoms
c     h1typ    biotypes for nucleotide backbone and sugar H1' atoms
c     c3typ    biotypes for nucleotide backbone and sugar C3' atoms
c     h3typ    biotypes for nucleotide backbone and sugar H3' atoms
c     c2typ    biotypes for nucleotide backbone and sugar C2' atoms
c     h21typ   biotypes for nucleotide backbone and sugar H2' atoms
c     o2typ    biotypes for nucleotide backbone and sugar O2' atoms
c     h22typ   biotypes for nucleotide backbone and sugar H2'' atoms
c     o3typ    biotypes for nucleotide backbone and sugar O3' atoms
c     ptyp     biotypes for nucleotide backbone and sugar P atoms
c     optyp    biotypes for nucleotide backbone and sugar OP atoms
c     h5ttyp   biotypes for nucleotide backbone and sugar H5T atoms
c     h3ttyp   biotypes for nucleotide backbone and sugar H3T atoms
c     amino    three-letter abbreviations for amino acids types
c     nuclz    three-letter abbreviations for nucleic acids types
c     amino1   one-letter abbreviations for amino acids types
c     nuclz1   one-letter abbreviations for nucleic acids types
c
c
      integer ntyp,catyp
      integer ctyp,hntyp
      integer otyp,hatyp,cbtyp
      integer nntyp,cantyp,cntyp
      integer hnntyp,ontyp,hantyp
      integer nctyp,cactyp,cctyp
      integer hnctyp,octyp,hactyp
      integer o5typ,c5typ,h51typ
      integer h52typ,c4typ,h4typ
      integer o4typ,c1typ,h1typ
      integer c3typ,h3typ,c2typ
      integer h21typ,o2typ,h22typ
      integer o3typ,ptyp,optyp
      integer h5ttyp,h3ttyp
      character*1 amino1,nuclz1
      character*3 amino,nuclz
      common /resdue/ ntyp(maxamino),catyp(maxamino),ctyp(maxamino),
     &                hntyp(maxamino),otyp(maxamino),hatyp(maxamino),
     &                cbtyp(maxamino),nntyp(maxamino),cantyp(maxamino),
     &                cntyp(maxamino),hnntyp(maxamino),ontyp(maxamino),
     &                hantyp(maxamino),nctyp(maxamino),cactyp(maxamino),
     &                cctyp(maxamino),hnctyp(maxamino),octyp(maxamino),
     &                hactyp(maxamino),o5typ(maxnuc),c5typ(maxnuc),
     &                h51typ(maxnuc),h52typ(maxnuc),c4typ(maxnuc),
     &                h4typ(maxnuc),o4typ(maxnuc),c1typ(maxnuc),
     &                h1typ(maxnuc),c3typ(maxnuc),h3typ(maxnuc),
     &                c2typ(maxnuc),h21typ(maxnuc),o2typ(maxnuc),
     &                h22typ(maxnuc),o3typ(maxnuc),ptyp(maxnuc),
     &                optyp(maxnuc),h5ttyp(maxnuc),h3ttyp(maxnuc),
     &                amino(maxamino),nuclz(maxnuc),amino1(maxamino),
     &                nuclz1(maxnuc)
