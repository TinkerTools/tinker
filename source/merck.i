c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 2007 by Nicolas Staelens & Jay W. Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  merck.i  --  parameters specific for MMFF force field  ##
c     ##                                                         ##
c     #############################################################
c
c
c     bt_1       atom pairs having Bond Type 1
c     nlignes    number of atom pairs having Bond Type 1
c     eqclass    table of atom class equivalencies used to find
c                default parameters if explicit values are missing
c                (see J. Comput. Chem., 17, 490-519, '95, Table IV)
c     crd        number of attached neighbors    |
c     val        valency value                   |  see T. A. Halgren,
c     pilp       if 0, no lone pair              |  J. Comput. Chem.,
c                if 1, one or more lone pair(s)  |  17, 616-645 (1995)
c     mltb       multibond indicator             |
c     arom       aromaticity indicator           |
c     lin        linearity indicator             |
c     sbmb       single- vs multiple-bond flag   |
c     mmffarom   aromatic rings parameters
c     mmffaromc  cationic aromatic rings parameters
c     mmffaroma  anionic aromatic rings parameters
c
c
      integer bt_1
      integer nlignes
      integer eqclass
      integer crd,val,pilp,mltb
      integer arom,lin,sbmb
      integer mmffarom
      integer mmffaromc
      integer mmffaroma
      common /merck1/ bt_1(500,2),nlignes,eqclass(500,5),crd(100),
     &                val(100),pilp(100),mltb(100),arom(100),lin(100),
     &                sbmb(100),mmffarom(maxtyp,6),
     &                mmffaromc(maxtyp,6),mmffaroma(maxtyp,6)
c
c
c     mmff_kb   bond force constant for pairs of atom classes
c     mmff_kb1  bond force constant for class pairs with Bond Type 1
c     mmff_b0   bond length value for pairs of atom classes
c     mmff_b01  bond length value for class pairs with Bond Type 1
c     rad0      covalent atomic radius for empirical bond rules
c     paulel    Pauling electronegativities for empirical bond rules
c     r0ref     reference bond length for empirical bond rules
c     kbref     reference force constant for empirical bond rules
c
c
      real*8 mmff_kb,mmff_kb1
      real*8 mmff_b0,mmff_b01
      real*8 rad0,r0ref
      real*8 kbref,paulel
      common /merck2/ mmff_kb(100,100),mmff_kb1(100,100),
     &                mmff_b0(100,100),mmff_b01(100,100),
     &                rad0(100),r0ref(100,100),kbref(100,100),
     &                paulel(100)
c
c
c     mmff_ka     angle force constant for triples of atom classes
c     mmff_ka1    angle force constant with one bond of Type 1
c     mmff_ka2    angle force constant with both bonds of Type 1
c     mmff_ka3    angle force constant for 3-membered ring
c     mmff_ka4    angle force constant for 4-membered ring
c     mmff_ka5    angle force constant for 3-ring and one Bond Type 1
c     mmff_ka6    angle force constant for 3-ring and both Bond Type 1
c     mmff_ka7    angle force constant for 4-ring and one Bond Type 1
c     mmff_ka8    angle force constant for 4-ring and both Bond Type 1
c     mmff_ang0   ideal bond angle for triples of atom classes
c     mmff_ang01  ideal bond angle with one bond of Type 1
c     mmff_ang02  ideal bond angle with both bonds of Type 1
c     mmff_ang03  ideal bond angle for 3-membered ring
c     mmff_ang04  ideal bond angle for 4-membered ring
c     mmff_ang05  ideal bond angle for 3-ring and one Bond Type 1
c     mmff_ang06  ideal bond angle for 3-ring and both Bond Type 1
c     mmff_ang07  ideal bond angle for 4-ring and one Bond Type 1
c     mmff_ang08  ideal bond angle for 4-ring and both Bond Type 1
c
c
      real*8 mmff_ka,mmff_ka1,mmff_ka2
      real*8 mmff_ka3,mmff_ka4,mmff_ka5
      real*8 mmff_ka6,mmff_ka7,mmff_ka8
      real*8 mmff_ang0,mmff_ang01,mmff_ang02
      real*8 mmff_ang03,mmff_ang04,mmff_ang05
      real*8 mmff_ang06,mmff_ang07,mmff_ang08
      common /merck3/ mmff_ka(0:100,100,0:100),
     &                mmff_ka1(0:100,100,0:100),
     &                mmff_ka2(0:100,100,0:100),
     &                mmff_ka3(0:100,100,0:100),
     &                mmff_ka4(0:100,100,0:100),
     &                mmff_ka5(0:100,100,0:100),
     &                mmff_ka6(0:100,100,0:100),
     &                mmff_ka7(0:100,100,0:100),
     &                mmff_ka8(0:100,100,0:100),
     &                mmff_ang0(0:100,100,0:100),
     &                mmff_ang01(0:100,100,0:100),
     &                mmff_ang02(0:100,100,0:100),
     &                mmff_ang03(0:100,100,0:100),
     &                mmff_ang04(0:100,100,0:100),
     &                mmff_ang05(0:100,100,0:100),
     &                mmff_ang06(0:100,100,0:100),
     &                mmff_ang07(0:100,100,0:100),
     &                mmff_ang08(0:100,100,0:100)
c
c
c     Stretch-Bend Type 0
c     stbn_abc     stretch-bend parameters for A-B-C atom classes
c     stbn_cba     stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 1  (A-B is Bond Type 1)
c     stbn_abc1    stretch-bend parameters for A-B-C atom classes
c     stbn_cba1    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 2  (B-C is Bond Type 1) 
c     stbn_abc2    stretch-bend parameters for A-B-C atom classes
c     stbn_cba2    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type = 3  (A-B and B-C are Bond Type 1) 
c     stbn_abc3    stretch-bend parameters for A-B-C atom classes
c     stbn_cba3    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 4  (both Bond Types 0, 4-membered ring)
c     stbn_abc4    stretch-bend parameters for A-B-C atom classes
c     stbn_cba4    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 5  (both Bond Types 0, 3-membered ring)
c     stbn_abc5    stretch-bend parameters for A-B-C atom classes
c     stbn_cba5    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 6  (A-B is Bond Type 1, 3-membered ring)
c     stbn_abc6    stretch-bend parameters for A-B-C atom classes
c     stbn_cba6    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 7  (B-C is Bond Type 1, 3-membered ring)
c     stbn_abc7    stretch-bend parameters for A-B-C atom classes
c     stbn_cba7    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 8  (both Bond Types 1, 3-membered ring)
c     stbn_abc8    stretch-bend parameters for A-B-C atom classes
c     stbn_cba8    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 9  (A-B is Bond Type 1, 4-membered ring)
c     stbn_abc9    stretch-bend parameters for A-B-C atom classes
c     stbn_cba9    stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 10  (B-C is Bond Type 1, 4-membered ring)
c     stbn_abc10   stretch-bend parameters for A-B-C atom classes
c     stbn_cba10   stretch-bend parameters for C-B-A atom classes
c     Stretch-Bend Type 11  (both Bond Types 1, 4-membered ring)
c     stbn_abc11   stretch-bend parameters for A-B-C atom classes
c     stbn_cba11   stretch-bend parameters for C-B-A atom classes
c     defstbn_abc  default stretch-bend parameters for A-B-C classes
c     defstbn_cba  default stretch-bend parameters for C-B-A classes
c
c
      real*8 stbn_abc,stbn_cba
      real*8 stbn_abc1,stbn_cba1
      real*8 stbn_abc2,stbn_cba2
      real*8 stbn_abc3,stbn_cba3
      real*8 stbn_abc4,stbn_cba4
      real*8 stbn_abc5,stbn_cba5
      real*8 stbn_abc6,stbn_cba6
      real*8 stbn_abc7,stbn_cba7
      real*8 stbn_abc8,stbn_cba8
      real*8 stbn_abc9,stbn_cba9
      real*8 stbn_abc10,stbn_cba10
      real*8 stbn_abc11,stbn_cba11
      real*8 defstbn_abc,defstbn_cba
      common /merck4/ stbn_abc(100,100,100),stbn_cba(100,100,100),
     &                stbn_abc1(100,100,100),stbn_cba1(100,100,100),
     &                stbn_abc2(100,100,100),stbn_cba2(100,100,100),
     &                stbn_abc3(100,100,100),stbn_cba3(100,100,100),
     &                stbn_abc4(100,100,100),stbn_cba4(100,100,100),
     &                stbn_abc5(100,100,100),stbn_cba5(100,100,100),
     &                stbn_abc6(100,100,100),stbn_cba6(100,100,100),
     &                stbn_abc7(100,100,100),stbn_cba7(100,100,100),
     &                stbn_abc8(100,100,100),stbn_cba8(100,100,100),
     &                stbn_abc9(100,100,100),stbn_cba9(100,100,100),
     &                stbn_abc10(100,100,100),stbn_cba10(100,100,100),
     &                stbn_abc11(100,100,100),stbn_cba11(100,100,100),
     &                defstbn_abc(0:4,0:4,0:4),defstbn_cba(0:4,0:4,0:4)
c
c
c     t1_1     torsional parameters for 1-fold, MMFF Torsion Type 1
c     t1_2     torsional parameters for 1-fold, MMFF Torsion Type 2
c     t2_1     torsional parameters for 2-fold, MMFF Torsion Type 1
c     t2_2     torsional parameters for 2-fold, MMFF Torsion Type 2
c     t3_1     torsional parameters for 3-fold, MMFF Torsion Type 1
c     t3_2     torsional parameters for 3-fold, MMFF Torsion Type 2
c     kt_1     string of classes for torsions, MMFF Torsion Type 1
c     kt_2     string of classes for torsions, MMFF Torsion Type 2
c
c
      real*8 t1_1,t2_1,t3_1
      real*8 t1_2,t2_2,t3_2
      character*16 kt_1,kt_2
      common /merck5/ t1_1(2,0:2000),t2_1(2,0:2000),t3_1(2,0:2000),
     &                t1_2(2,0:2000),t2_2(2,0:2000),t3_2(2,0:2000),
     &                kt_1(0:2000),kt_2(0:2000)
c
c
c     g        scale factors for calculation of MMFF eps
c     alph     atomic polarizabilities for calculation of MMFF eps
c     nn       effective number of valence electrons for MMFF eps
c     da       donor/acceptor atom classes
c
c
      real*8 g,alph,nn
      character*1 da
      common /merck6/ g(maxclass),alph(maxclass),nn(maxclass),
     &                da(maxclass)
c
c
c     bci      bond charge increments for building atom charges
c     bci_1    bond charge increments for Bond Type 1   
c     pbci     partial BCI for building missing BCI's
c     fcadj    formal charge adjustment factor
c
c
      real*8 bci,bci_1
      real*8 pbci,fcadj
      common /merck7/ bci(100,100),bci_1(100,100),pbci(maxclass),
     &                fcadj(maxclass)
