/* abstree.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#define OpCode(x) (((x)->type)&0x0f)

/* Operator Types */
#define OpAnd            0x01
#define OpOr             0x02
#define OpNot            0x03
#define OpEqual          0x04
#define OpNotEq          0x05
#define OpLess           0x06
#define OpMore           0x07
#define OpLessEq         0x08
#define OpMoreEq         0x09
#define OpConst          0x0a
#define OpWithin         0x0b
#define OpMember         0xac

#define OpLftProp        0x10
#define OpLftVal         0x20
#define OpRgtProp        0x40
#define OpRgtVal         0x80

/* Property fields */
#define PropIdent        1
#define PropXCord        2
#define PropYCord        3
#define PropZCord        4
#define PropTemp         5
#define PropRad          6
#define PropResId        7
#define PropName         8
#define PropChain        9
#define PropResName      10
#define PropSelect       11
#define PropElemNo       12
#define PropModel        13

#define PredAbsOrd(x)    ((x)-20)
#define PredAbsChr(x)    ((x)+20)

#define PredAlpha        20
#define PredAmino        21
#define PredAT           22
#define PredBonded       23
#define PredCG           24
#define PredCystine      25
#define PredDNA          26
#define PredHelix        27
#define PredHetero       28
#define PredHydrogen     29
#define PredIon          30
#define PredLigand       31
#define PredMainChain    32
#define PredNucleic      33
#define PredProtein      34
#define PredPurine       35
#define PredPyrimidine   36
#define PredRNA          37
#define PredSelected     38 /* Unused! */
#define PredSheet        39
#define PredSidechain    40
#define PredSolvent      41
#define PredTurn         42
#define PredWater	 43

#define PredAcidic       44
#define PredAcyclic      45
#define PredAliphatic    46
#define PredAromatic     47
#define PredBasic        48
#define PredBuried       49
#define PredCharged      50
#define PredCyclic       51
#define PredHydrophobic  52
#define PredLarge        53
#define PredMedium       54
#define PredNeutral      55
#define PredPolar        56
#define PredSmall        57
#define PredSurface      58



#define SetSize     10
typedef struct _AtomSet {
	struct _AtomSet __far *next;
	Atom __far *data[SetSize];
        int count;
        } AtomSet;
        
typedef union {
	AtomSet __far *set;
	struct _Expr *ptr;
        Long limit;
	int val;
	} Branch;

typedef struct _Expr {
	int type;
        Branch rgt;
        Branch lft;
	} Expr;

/* CPK Colour Indices
 *  0 Light Grey    1 Sky Blue      2 Red           3 Yellow
 *  4 White         5 Pink          6 Golden Rod    7 Blue
 *  8 Orange        9 Dark Grey    10 Brown        11 Purple
 * 12 Deep Pink    13 Green        14 Fire Brick   15 Mid Green
 */

#define MAXELEMNO  104
typedef struct {
           char symbol[2];
           int covalrad;
           int vdwrad;
           int cpkcol;
           char *name;
        } ElemStruct;

/* Structures with Implicit Hydrogens */
#define VDWCarbon    468
#define VDWNitrogen  375
#define VDWOxygen    350
#define VDWSulphur   462

#ifdef ABSTREE
ElemStruct Element[MAXELEMNO] =  {
    { { ' ', ' ' }, 180, 360, 12, ""             },  /*   0 */
    { { 'H', ' ' },  80, 275,  4, "HYDROGEN"     },  /*   1 */
    { { 'H', 'e' }, 400, 550,  5, "HELIUM"       },  /*   2 */
    { { 'L', 'i' }, 170, 305, 14, "LITHIUM"      },  /*   3 */
    { { 'B', 'e' },  88, 157, 12, "BERYLLIUM"    },  /*   4 */
    { { 'B', ' ' }, 208, 387, 13, "BORON"        },  /*   5 */
    { { 'C', ' ' }, 180, 387,  0, "CARBON"       },  /*   6 */
    { { 'N', ' ' }, 170, 350,  1, "NITROGEN"     },  /*   7 */
    { { 'O', ' ' }, 170, 337,  2, "OXYGEN"       },  /*   8 */
    { { 'F', ' ' }, 160, 325,  6, "FLUORINE"     },  /*   9 */
    { { 'N', 'e' }, 280, 505, 12, "NEON"         },  /*  10 */
    { { 'N', 'a' }, 243, 550,  7, "SODIUM"       },  /*  11 */
    { { 'M', 'g' }, 275, 375, 15, "MAGNESIUM"    },  /*  12 */
    { { 'A', 'l' }, 338, 375,  9, "ALUMINIUM"    },  /*  13 */
    { { 'S', 'i' }, 300, 550,  6, "SILICON"      },  /*  14 */
    { { 'P', ' ' }, 259, 470,  8, "PHOSPHORUS"   },  /*  15 */  /* 262? */
    { { 'S', ' ' }, 255, 452,  3, "SULPHUR"      },  /*  16 */
    { { 'C', 'l' }, 250, 437, 13, "CHLORINE"     },  /*  17 */
    { { 'A', 'r' }, 392, 692, 12, "ARGON"        },  /*  18 */
    { { 'K', ' ' }, 332, 597, 12, "POTASSIUM"    },  /*  19 */
    { { 'C', 'a' }, 248, 487,  9, "CALCIUM"      },  /*  20 */
    { { 'S', 'c' }, 360, 330, 12, "SCANDIUM"     },  /*  21 */
    { { 'T', 'i' }, 368, 487,  9, "TITANIUM"     },  /*  22 */
    { { 'V', ' ' }, 332, 265, 12, "VANADIUM"     },  /*  23 */
    { { 'C', 'r' }, 338, 282,  9, "CHROMIUM"     },  /*  24 */
    { { 'M', 'n' }, 338, 297,  9, "MANGANESE"    },  /*  25 */
    { { 'F', 'e' }, 335, 487,  8, "IRON"         },  /*  26 */
    { { 'C', 'o' }, 332, 282, 12, "COBALT"       },  /*  27 */
    { { 'N', 'i' }, 405, 310, 10, "NICKEL"       },  /*  28 */  /* >375! */
    { { 'C', 'u' }, 380, 287, 10, "COPPER"       },  /*  29 */
    { { 'Z', 'n' }, 362, 287, 10, "ZINC"         },  /*  30 */
    { { 'G', 'a' }, 305, 387, 12, "GALLIUM"      },  /*  31 */
    { { 'G', 'e' }, 292, 999, 12, "GERMANIUM"    },  /*  32 */  /* 1225? */
    { { 'A', 's' }, 302, 207, 12, "ARSENIC"      },  /*  33 */
    { { 'S', 'e' }, 305, 225, 12, "SELENIUM"     },  /*  34 */
    { { 'B', 'r' }, 302, 437, 10, "BROMINE"      },  /*  35 */
    { { 'K', 'r' }, 400, 475, 12, "KRYPTON"      },  /*  36 */
    { { 'R', 'b' }, 368, 662, 12, "RUBIDIUM"     },  /*  37 */
    { { 'S', 'r' }, 280, 505, 12, "STRONTIUM"    },  /*  38 */
    { { 'Y', ' ' }, 445, 402, 12, "YTTRIUM"      },  /*  39 */
    { { 'Z', 'r' }, 390, 355, 12, "ZIRCONIUM"    },  /*  40 */
    { { 'N', 'b' }, 370, 332, 12, "NIOBIUM"      },  /*  41 */
    { { 'M', 'o' }, 368, 437, 12, "MOLYBDENUM"   },  /*  42 */
    { { 'T', 'c' }, 338, 450, 12, "TECHNETIUM"   },  /*  43 */
    { { 'R', 'u' }, 350, 300, 12, "RUTHENIUM"    },  /*  44 */
    { { 'R', 'h' }, 362, 305, 12, "RHODIUM"      },  /*  45 */
    { { 'P', 'd' }, 375, 360, 12, "PALLADIUM"    },  /*  46 */
    { { 'A', 'g' }, 398, 387,  9, "SILVER"       },  /*  47 */
    { { 'C', 'd' }, 422, 437, 12, "CADMIUM"      },  /*  48 */
    { { 'I', 'n' }, 408, 362, 12, "INDIUM"       },  /*  49 */
    { { 'S', 'n' }, 365, 417, 12, "TIN",         },  /*  50 */
    { { 'S', 'b' }, 365, 280, 12, "ANTIMONY"     },  /*  51 */
    { { 'T', 'e' }, 368, 315, 12, "TELLURIUM"    },  /*  52 */
    { { 'I', ' ' }, 350, 437, 11, "IODINE"       },  /*  53 */
    { { 'X', 'e' }, 425, 525, 12, "XENON"        },  /*  54 */
    { { 'C', 's' }, 418, 752, 12, "CAESIUM"      },  /*  55 */
    { { 'B', 'a' }, 335, 602,  8, "BARIUM"       },  /*  56 */
    { { 'L', 'a' }, 468, 457, 12, "LANTHANUM"    },  /*  57 */
    { { 'C', 'e' }, 458, 465, 12, "CERIUM"       },  /*  58 */
    { { 'P', 'r' }, 455, 405, 12, "PRASEODYMIUM" },  /*  59 */
    { { 'N', 'd' }, 452, 447, 12, "NEODYMIUM"    },  /*  60 */
    { { 'P', 'm' }, 450, 440, 12, "PROMETHIUM"   },  /*  61 */
    { { 'S', 'm' }, 450, 435, 12, "SAMARIUM"     },  /*  62 */
    { { 'E', 'u' }, 498, 490, 12, "EUROPIUM"     },  /*  63 */
    { { 'G', 'd' }, 448, 422, 12, "GADOLINIUM"   },  /*  64 */
    { { 'T', 'b' }, 440, 415, 12, "TERBIUM"      },  /*  65 */
    { { 'D', 'y' }, 438, 407, 12, "DYSPROSIUM"   },  /*  66 */
    { { 'H', 'o' }, 435, 402, 12, "HOLMIUM"      },  /*  67 */
    { { 'E', 'r' }, 432, 397, 12, "ERBIUM"       },  /*  68 */
    { { 'T', 'm' }, 430, 392, 12, "THULIUM"      },  /*  69 */
    { { 'Y', 'b' }, 485, 385, 12, "YTTERBIUM"    },  /*  70 */
    { { 'L', 'u' }, 430, 382, 12, "LUTETIUM"     },  /*  71 */
    { { 'H', 'f' }, 392, 350, 12, "HAFNIUM"      },  /*  72 */
    { { 'T', 'a' }, 358, 305, 12, "TANTALUM"     },  /*  73 */
    { { 'W', ' ' }, 342, 315, 12, "TUNGSTEN"     },  /*  74 */
    { { 'R', 'e' }, 338, 325, 12, "RHENIUM"      },  /*  75 */
    { { 'O', 's' }, 342, 395, 12, "OSMIUM"       },  /*  76 */
    { { 'I', 'r' }, 330, 305, 12, "IRIDIUM"      },  /*  77 */
    { { 'P', 't' }, 375, 387, 12, "PLATINUM"     },  /*  78 */
    { { 'A', 'u' }, 375, 362,  6, "GOLD"         },  /*  79 */
    { { 'H', 'g' }, 425, 495, 12, "MERCURY"      },  /*  80 */
    { { 'T', 'l' }, 388, 427, 12, "THALLIUM"     },  /*  81 */
    { { 'P', 'b' }, 385, 540, 12, "LEAD"         },  /*  82 */
    { { 'B', 'i' }, 385, 432, 12, "BISMUTH"      },  /*  83 */
    { { 'P', 'o' }, 420, 302, 12, "POLONIUM"     },  /*  84 */
    { { 'A', 't' }, 302, 280, 12, "ASTATINE"     },  /*  85 */
    { { 'R', 'n' }, 475, 575, 12, "RADON"        },  /*  86 */
    { { 'F', 'r' }, 450, 810, 12, "FRANCIUM"     },  /*  87 */
    { { 'R', 'a' }, 358, 642, 12, "RADIUM"       },  /*  88 */
    { { 'A', 'c' }, 295, 530, 12, "ACTINIUM"     },  /*  89 */
    { { 'T', 'h' }, 255, 460, 12, "THORIUM"      },  /*  90 */
    { { 'P', 'a' }, 222, 400, 12, "PROTACTINIUM" },  /*  91 */
    { { 'U', ' ' }, 242, 437, 12, "URANIUM"      },  /*  92 */
    { { 'N', 'p' }, 238, 427, 12, "NEPTUNIUM"    },  /*  93 */
    { { 'P', 'u' }, 232, 417, 12, "PLUTONIUM"    },  /*  94 */
    { { 'A', 'm' }, 230, 415, 12, "AMERICIUM"    },  /*  95 */
    { { 'C', 'm' }, 228, 412, 12, "CURIUM"       },  /*  96 */
    { { 'B', 'k' }, 225, 410, 12, "BERKELIUM"    },  /*  97 */
    { { 'C', 'f' }, 222, 407, 12, "CALIFORNIUM"  },  /*  98 */
    { { 'E', 's' }, 220, 405, 12, "EINSTEINIUM"  },  /*  99 */
    { { 'F', 'm' }, 218, 402, 12, "FERMIUM"      },  /* 100 */
    { { 'M', 'd' }, 215, 400, 12, "MENDELEVIUM"  },  /* 101 */
    { { 'N', 'o' }, 212, 397, 12, "NOBELIUM"     },  /* 102 */
    { { 'L', 'r' }, 210, 395, 12, "LAWRENCIUM"   }   /* 103 */ /* Lw? */
        };


Expr *QueryExpr;
Chain __far *QChain;
Group __far *QGroup;
Atom __far *QAtom;

#else
extern ElemStruct Element[MAXELEMNO];

extern Expr *QueryExpr;
extern Chain __far *QChain;
extern Group __far *QGroup;
extern Atom __far *QAtom;

#ifdef FUNCPROTO
Expr *AllocateNode();
void DeAllocateExpr( Expr* );
int EvaluateExpr( Expr* );
int DefineSetExpr( char*, Expr* );
Expr *LookUpSetExpr( char* );
AtomSet __far *BuildAtomSet( Expr* );
void DeleteAtomSet( AtomSet __far* );
Expr *LookUpElement( char* );

int ElemVDWRadius( int );
int ParsePrimitiveExpr( char** );
int GetElemNumber( Group __far*, Atom __far* );
void FormatLabel( Chain __far*, Group __far*, Atom __far*, char*, char* );
void InitialiseAbstree();
void ResetSymbolTable();

double CalcTorsion( Atom __far*, Atom __far*, Atom __far*, Atom __far* );
double CalcAngle( Atom __far*, Atom __far*, Atom __far* );
double CalcDistance( Atom __far*, Atom __far* );

#else /* non-ANSI C compiler */
Expr *AllocateNode();
void DeAllocateExpr();
int EvaluateExpr();
int DefineSetExpr();
Expr *LookUpSetExpr();
AtomSet __far *BuildAtomSet();
void DeleteAtomSet();
Expr *LookUpElement();

int ElemVDWRadius();
int ParsePrimitiveExpr();
int GetElemNumber();
void FormatLabel();
void InitialiseAbstree();
void ResetSymbolTable();

double CalcTorsion();
double CalcAngle();
double CalcDistance();
#endif
#endif

