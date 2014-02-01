c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  sequen.i  --  sequence information for a biopolymer  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     nseq     total number of residues in biopolymer sequences
c     nchain   number of separate biopolymer sequence chains
c     ichain   first and last residue in each biopolymer chain
c     seqtyp   residue type for each residue in the sequence
c     seq      three-letter code for each residue in the sequence
c     chnnam   one-letter identifier for each sequence chain
c     chntyp   contents of each chain (GENERIC, PEPTIDE or NUCLEIC)
c
c
      integer nseq,seqtyp
      integer nchain,ichain
      character*1 chnnam
      character*3 seq
      character*7 chntyp
      common /sequen/ nseq,nchain,ichain(2,maxres),seqtyp(maxres),
     &                seq(maxres),chnnam(maxres),chntyp(maxres)
