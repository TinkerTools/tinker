/* molecule.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#define MAXMASK 40
#define MAXELEM 256
#define MINELEM 29
#define MAXRES  100
#define MINRES  53

/*=================*/
/*  Database Flags */
/*=================*/

#define SelectFlag      0x01
#define DrawBondFlag    0x0e
#define AllAtomFlag     0x1c

/* Atom Flags */
#define SphereFlag      0x02     /* Sphere representation */
#define HydrogenFlag    0x08     /* Hydrogen atom         */
#define NormAtomFlag    0x10
#define NonBondFlag     0x20
#define BreakFlag       0x40     /* Break in backbone     */

/* Bond Flags */
#define WireFlag        0x02     /* Depth-cued wireframe         */
#define DashFlag        0x04     /* Dashed Depth-cued wireframe  */
#define CylinderFlag    0x08     /* Line/Cylinder representation */

#define HydrBondFlag    0x00     /* Hydrogen bond [place keeper] */
#define NormBondFlag    0x10
#define DoubBondFlag    0x20
#define TripBondFlag    0x40
#define AromBondFlag    0x80

/*=====================*/
/*  Molecule Database  */
/*=====================*/

typedef struct _Atom {
        struct _Atom __far *anext;        /* Linked list of atoms  */
        struct _Atom __far *bucket;       /* Sphere Y-Bucket       */
        struct _Atom __far *next;         /* Active Object List    */
        Long   xorg, yorg, zorg;          /* World Co-ordinates    */
        short  x, y, z;                   /* Image Co-ordinates    */
        short  radius;                    /* World Radius          */
        short  col;                       /* Atom Colour           */
        Long   serno;                     /* Atom Serial Number    */
        void   *label;                    /* Atom Label Structure  */
        Byte   elemno;                    /* Atomic Number         */
        Byte   refno;                     /* ElemDesc index number */
        Byte   flag;                      /* Database flags        */
        char   altl;                      /* Alternate Location    */
        short  irad;                      /* Image Radius          */
        short  mbox;                      /* Shadow Casting NOnce  */
    } Atom;


typedef struct _Bond {
        struct _Bond __far *bnext;       /* Linked list of bonds  */
        Atom __far *srcatom;             /* Source Atom Ptr       */
        Atom __far *dstatom;             /* Destination Atom Ptr  */
        short radius;                    /* World Radius          */
        short irad;                      /* Image Radius          */
        short col;                       /* Bond Colour           */
        Byte  flag;                      /* Database flags        */
    } Bond;

typedef struct _Group {
        struct _Group __far *gnext;       /* Linked list of groups */
        Atom __far *alist;                /* Linked list of atoms  */
        short serno;                      /* Group serial number   */
        char  insert;                     /* PDB insertion code    */
        Byte  refno;                      /* Residue index number  */
        Byte  struc;                      /* Secondary Structure   */
        Byte  flag;                       /* Database flags        */
    } Group;
 
#ifdef APPLEMAC
/* Avoid Name Clash! */
#define Chain ChainSeg
#endif

typedef struct _ChainSeg {
        struct _ChainSeg __far *cnext;    /* Linked list of chains     */
        Group __far *glist;               /* Linked list of groups     */
        Bond __far *blist;                /* Linked list of back bonds */
        char ident;                       /* Chain identifier          */
        Byte model;                       /* NMR Model / Symmetry      */
    } Chain;

typedef struct _AtomRef {
        Chain __far *chn;
        Group __far *grp;
        Atom  __far *atm;
    } AtomRef;

typedef struct _Molecule {
        Chain __far *clist;               /* Linked list of chains    */
        Bond __far *blist;                /* Linked list of bonds     */
    } Molecule;

/*========================*/
/* Other Consts & Structs */
/*========================*/

#define SourceNone   0
#define SourcePDB    1
#define SourceCalc   2
 
#define SerNoFlag 0x01
#define ResNoFlag 0x02

typedef struct {
        short radius;
        char  mask[19];
        Byte  flags;
        Byte  r;
        Byte  g;
        Byte  b;
        } MaskDesc;

typedef struct _IntCoord {
        struct _IntCoord __far *inext;
        short na,nb,nc;
        short refno;
        Real dihed;
        Real angle;
        Real dist;
    } IntCoord;

typedef struct _InfoStruct {
        char filename[256];
        char moleculename[80];
        char classification[42];
        char identcode[6];

        char spacegroup[11];
        Real cellalpha, cellbeta, cellgamma;
        Real cella, cellb, cellc;

        Long bondcount;
        int chaincount;
        int ssbondcount;
        int hbondcount;

        int structsource;
        int laddercount;
        int helixcount;
        int turncount;
    } InfoStruct;


#ifdef MOLECULE
/* Avoid SGI Compiler Warnings! */
char Residue[MAXRES][4] = {

    /*===============*/
    /*  Amino Acids  */
    /*===============*/

/* Ordered by Cumulative Frequency in Brookhaven *
 * Protein Databank, December 1991               */

          "ALA", /* 8.4% */     "GLY", /* 8.3% */
          "LEU", /* 8.0% */     "SER", /* 7.5% */
          "VAL", /* 7.1% */     "THR", /* 6.4% */
          "LYS", /* 5.8% */     "ASP", /* 5.5% */
          "ILE", /* 5.2% */     "ASN", /* 4.9% */
          "GLU", /* 4.9% */     "PRO", /* 4.4% */
          "ARG", /* 3.8% */     "PHE", /* 3.7% */
          "GLN", /* 3.5% */     "TYR", /* 3.5% */
          "HIS", /* 2.3% */     "CYS", /* 2.0% */
          "MET", /* 1.8% */     "TRP", /* 1.4% */

          "ASX", "GLX", "PCA", "HYP",

    /*===================*/
    /*  DNA Nucleotides  */
    /*===================*/
          "  A", "  C", "  G", "  T",

    /*===================*/
    /*  RNA Nucleotides  */
    /*===================*/
          "  U", " +U", "  I", "1MA", 
          "5MC", "OMC", "1MG", "2MG", 
          "M2G", "7MG", "OMG", " YG", 
          "H2U", "5MU", "PSU",

    /*=================*/
    /*  Miscellaneous  */ 
    /*=================*/
          "UNK", "ACE", "FOR", "HOH",
          "DOD", "SO4", "PO4", "NAD",
          "COA", "NAP"  };

/* Avoid SGI Compiler Warnings! */
char ElemDesc[MAXELEM][4] = {
    { ' ', 'N', ' ', ' ' },  /* 0*/
    { ' ', 'C', 'A', ' ' },  /* 1*/
    { ' ', 'C', ' ', ' ' },  /* 2*/
    { ' ', 'O', ' ', ' ' },  /* 3*/   /* 0-3   Amino Acid Backbone    */
    { ' ', 'C', '\'', ' ' }, /* 4*/
    { ' ', 'O', 'T', ' ' },  /* 5*/
    { ' ', 'S', ' ', ' ' },  /* 6*/
    { ' ', 'P', ' ', ' ' },  /* 7*/   /* 4-7   Shapely Amino Backbone */
    { ' ', 'O', '1', 'P' },  /* 8*/
    { ' ', 'O', '2', 'P' },  /* 9*/
    { ' ', 'O', '5', '*' },  /*10*/
    { ' ', 'C', '5', '*' },  /*11*/
    { ' ', 'C', '4', '*' },  /*12*/
    { ' ', 'O', '4', '*' },  /*13*/
    { ' ', 'C', '3', '*' },  /*14*/
    { ' ', 'O', '3', '*' },  /*15*/
    { ' ', 'C', '2', '*' },  /*16*/
    { ' ', 'O', '2', '*' },  /*17*/
    { ' ', 'C', '1', '*' },  /*18*/   /* 7-18  Nucleic Acid Backbone  */
    { ' ', 'C', 'A', '2' },  /*19*/   /* 19    Shapely Special        */
    { ' ', 'S', 'G', ' ' },  /*20*/   /* 20    Cysteine Sulphur       */
    { ' ', 'N', '1', ' ' },  /*21*/
    { ' ', 'N', '2', ' ' },  /*22*/
    { ' ', 'N', '3', ' ' },  /*23*/
    { ' ', 'N', '4', ' ' },  /*24*/
    { ' ', 'N', '6', ' ' },  /*25*/
    { ' ', 'O', '2', ' ' },  /*26*/
    { ' ', 'O', '4', ' ' },  /*27*/
    { ' ', 'O', '6', ' ' }   /*28*/   /* 21-28 Nucleic Acid H-Bonding */
    };

InfoStruct Info;
int MainGroupCount;
Long MainAtomCount; 

Long MinX, MinY, MinZ;
Long MaxX, MaxY, MaxZ;

int MMinMaxFlag;
int MinMainRes,  MaxMainRes;

Molecule __far *CurMolecule;
Chain __far *CurChain;
Group __far *CurGroup;
Atom __far *CurAtom;

IntCoord __far *IntList;
Molecule __far *Database;
MaskDesc UserMask[MAXMASK];
Long MaxHBondDist;
Long MinBondDist,  MaxBondDist;
int ElemNo,ResNo;
int HasHydrogen;
int MaskCount;
int NMRModel;

#else
extern char Residue[MAXRES][4];
extern char ElemDesc[MAXELEM][4];
extern InfoStruct Info;

extern int MainGroupCount;
extern Long MainAtomCount;

extern Long MinX, MinY, MinZ;
extern Long MaxX, MaxY, MaxZ;

extern int MMinMaxFlag;
extern int MinMainRes,  MaxMainRes;

extern Molecule __far *CurMolecule;
extern Chain __far *CurChain;
extern Group __far *CurGroup;
extern Atom __far *CurAtom;

extern IntCoord __far *IntList;
extern Molecule __far *Database;
extern MaskDesc UserMask[MAXMASK];
extern Long MaxHBondDist;
extern Long MinBondDist,  MaxBondDist;
extern int ElemNo,ResNo;
extern int HasHydrogen;
extern int MaskCount;
extern int NMRModel;

#ifdef FUNCPROTO
void CreateChain( int );
void CreateGroup( int );
void ProcessGroup( int );
void CreateMolGroup();
int FindResNo( char* );

Atom __far *CreateAtom();
Atom __far *FindGroupAtom( Group __far*, int );
void ProcessAtom( Atom __far* );

int NewAtomType( char* );
int SimpleAtomType( char* );
int ComplexAtomType( char* );

Bond __far *ProcessBond( Atom __far*, Atom __far*, int );
void CreateBond( Long, Long, int );
void CreateBondOrder( Long, Long );
void CreateMoleculeBonds( int, int );
void CalcHydrogenBonds();

void InitInternalCoords();
IntCoord __far* AllocInternalCoord();
int ConvertInternal2Cartesian();
void FreeInternalCoords();

void DetermineStructure( int );
void RenumberMolecule( int );

void InitialiseDatabase();
void DescribeMolecule();
void DestroyDatabase();
void PurgeDatabase();

#else /* non-ANSI C compiler */

void CreateChain();
void CreateGroup();
void CreateMolGroup();
int FindResNo();

Atom __far *CreateAtom();
Atom __far *FindGroupAtom();
void ProcessAtom();

int NewAtomType();
int SimpleAtomType();
int ComplexAtomType();

Bond __far *ProcessBond();
void CreateBond();
void CreateBondOrder();
void CreateMoleculeBonds();
void CalcHydrogenBonds();

void InitInternalCoords();
IntCoord __far* AllocInternalCoord();
int ConvertInternal2Cartesian();
void FreeInternalCoords();

void DetermineStructure();
void RenumberMolecule();

void InitialiseDatabase();
void DescribeMolecule();
void DestroyDatabase();
void PurgeDatabase();

#endif
#endif
