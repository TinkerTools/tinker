/* molecule.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define MOLECULE
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"

#define BondPool    32
#define AtomPool    32

#define NoLadder     0x00
#define ParaLadder   0x01
#define AntiLadder   0x02

#define Cos70Deg     0.34202014332567

#define MaxHBondDist   ((Long)300*300)
#define MaxBondDist    ((Long)475*475)
#define MinBondDist    ((Long)100*100)
#define AbsMaxBondDist 600

#ifdef APPLEMAC
#define AllocSize   256
typedef struct _AllocRef {
        struct _AllocRef *next;
        void *data[AllocSize];
        int count;
        } AllocRef;
static AllocRef *AllocList;  
#endif

typedef struct {
          char name[4];
          int code;
      } SynonymTable;

#define RESSYNMAX 16
static SynonymTable ResSynonym[RESSYNMAX] = {
    { "ADE", 24 },  /*   A : Adenosine   */
    { "CPR", 11 },  /* PRO : Cis-proline */
    { "CSH", 17 },  /* CYS : Cystine     */
    { "CSM", 17 },  /* CYS : Cystine     */
    { "CYH", 17 },  /* CYS : Cystine     */
    { "CYT", 25 },  /*   C : Cytosine    */
    { "D2O", 47 },  /* DOD : Heavy Water */
    { "GUA", 26 },  /*   G : Guanosine   */
    { "H2O", 46 },  /* HOH : Solvent     */
    { "SOL", 46 },  /* HOH : Solvent     */
    { "SUL", 48 },  /* SO4 : Sulphate    */
    { "THY", 27 },  /*   T : Thymidine   */
    { "TIP", 46 },  /* HOH : Water       */
    { "TRY", 20 },  /* TRP : Tryptophan  */
    { "URI", 28 },  /*   U : Uridine     */
    { "WAT", 46 }   /* HOH : Water       */
        };

static Molecule __far *FreeMolecule;
static Chain __far *FreeChain;
static Group __far *FreeGroup;
static Atom __far *FreeAtom;
static Bond __far *FreeBond;

static IntCoord __far *IntPrev;
static int MemSize;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)

/* Forward Reference */
void DestroyDatabase();


#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif


static void FatalDataError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Database Error: %s!",ptr);
    RasMolFatalExit(buffer);
}



void DescribeMolecule()
{
    char buffer[40];

    if( CommandActive )
        WriteChar('\n');
    CommandActive=False;

    if( *Info.moleculename )
    {   WriteString("Molecule Name ....... ");
        WriteString(Info.moleculename);
        WriteChar('\n');
    }

    if( *Info.classification )
    {   WriteString("Classification ...... ");
        WriteString(Info.classification);
        WriteChar('\n');
    }

/*  Remove Chain/Group Output, JWP, June 2001

    if( Info.chaincount > 1 )
    {   sprintf(buffer,"Number of Chains .... %d\n",Info.chaincount);
        WriteString(buffer);
    }

    sprintf(buffer,"Number of Groups .... %d",MainGroupCount);
    WriteString(buffer);
    WriteChar('\n');

End of Remove Chain/Group Output  */

    sprintf(buffer,"Number of Atoms ..... %ld",(long)MainAtomCount);
    WriteString(buffer);
    WriteChar('\n');

    if( Info.bondcount )
    {   sprintf(buffer,"Number of Bonds ..... %ld\n",(long)Info.bondcount);
        WriteString(buffer);
    }
}


#ifdef APPLEMAC
/* Avoid System Memory Leaks! */
static void RegisterAlloc( data )
    void *data;
{
    register AllocRef *ptr;
    
    if( !AllocList || (AllocList->count==AllocSize) )
    {   ptr = (AllocRef *)_fmalloc( sizeof(AllocRef) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        
        ptr->next = AllocList;
        ptr->data[0] = data;
        ptr->count = 1;
        AllocList = ptr;
    } else AllocList->data[AllocList->count++] = data;
}
#else
#define RegisterAlloc(x)
#endif


/*==================================*/
/* Group & Chain Handling Functions */
/*==================================*/


void CreateChain( ident )
    int ident;
{
    register Chain __far *prev;

    if( !CurMolecule )
    {   if( !(CurMolecule = FreeMolecule) )
        {   MemSize += sizeof(Molecule);
            CurMolecule = (Molecule __far *)_fmalloc(sizeof(Molecule));
            if( !CurMolecule ) FatalDataError("Memory allocation failed");
            RegisterAlloc( CurMolecule );
        } else FreeMolecule = (void __far*)0;

        CurChain = (void __far*)0;
        CurMolecule->blist = (void __far*)0;
        CurMolecule->clist = (void __far*)0;
        Database = CurMolecule;
    }

    /* Handle chain breaks! */
    if( !(prev=CurChain) )
        if( (prev=CurMolecule->clist) )
            while( prev->cnext )
                prev = prev->cnext;

    if( !(CurChain = FreeChain) )
    {   MemSize += sizeof(Chain);
        CurChain = (Chain __far *)_fmalloc(sizeof(Chain));
        if( !CurChain ) FatalDataError("Memory allocation failed");
        RegisterAlloc( CurChain );
    } else FreeChain = FreeChain->cnext;

    if( prev )
    {   prev->cnext = CurChain;
    } else CurMolecule->clist = CurChain;
    CurChain->cnext = (void __far*)0;
     
    CurChain->ident = ident;
    CurChain->model = NMRModel;
    CurChain->glist = (void __far*)0;
    CurChain->blist = (void __far*)0;
    CurGroup = (void __far*)0;
    Info.chaincount++;
}


void CreateGroup( pool )
    int pool;
{
    register Group __far *ptr;
    register int i;

    if( !(ptr = FreeGroup) )
    {   MemSize += pool*sizeof(Group);
        ptr = (Group __far *)_fmalloc( pool*sizeof(Group) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        RegisterAlloc( ptr );
        for( i=1; i<pool; i++ )
        {   ptr->gnext = FreeGroup;
            FreeGroup = ptr++;
        } 
    } else FreeGroup = ptr->gnext;
    
    if( CurGroup )
    {   ptr->gnext = CurGroup->gnext;
        CurGroup->gnext = ptr;
    } else 
    {   ptr->gnext = CurChain->glist;
        CurChain->glist = ptr;
    }
    CurGroup = ptr;

    CurAtom = (void __far*)0;
    ptr->alist = (void __far*)0;
    ptr->insert = ' ';
    ptr->struc = 0;
    ptr->flag = 0;
}


int FindResNo( ptr )
    char *ptr;
{
    register int hi,lo;
    register int refno;
    register int flag;
    register int mid;

    for( refno=0; refno<ResNo; refno++ )
        if( !strncmp(Residue[refno],ptr,3) )
            return( refno );

    lo = 0;
    hi = RESSYNMAX;
    while( lo < hi )
    {   mid = (hi+lo)>>1;
        flag = strncmp(ResSynonym[mid].name,ptr,3);
        if( !flag ) return( ResSynonym[mid].code );

        /* Binary Search */
        if( flag<0 )
        {   lo = mid+1;
        } else hi = mid;
    }

    if( ResNo++ == MAXRES )
        FatalDataError("Too many new residues");
    Residue[refno][0] = *ptr++;
    Residue[refno][1] = *ptr++;
    Residue[refno][2] = *ptr;
    return( refno );
}


void ProcessGroup( heta )
    int heta;
{
    register int serno;

    serno = CurGroup->serno;

    MainGroupCount++;
    if( MMinMaxFlag )
    {   if( serno > MaxMainRes )
        {   MaxMainRes = serno;
        } else if( serno < MinMainRes )
            MinMainRes = serno;
    } else MinMainRes = MaxMainRes = serno;
}


void CreateMolGroup()
{
    strcpy(Info.filename,DataFileName);

    CreateChain( ' ' );
    CreateGroup( 1 );

    CurGroup->refno = FindResNo( "MOL" );
    CurGroup->serno = 1;
        
    MinMainRes = MaxMainRes = 1;
    MainGroupCount = 1;
}


/*=========================*/
/* Atom Handling Functions */
/*=========================*/


Atom __far *CreateAtom()
{
    register Atom __far *ptr;
    register int i;

    if( !(ptr = FreeAtom) )
    {   MemSize += AtomPool*sizeof(Atom);
        ptr = (Atom __far *)_fmalloc( AtomPool*sizeof(Atom) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        RegisterAlloc( ptr );
        for( i=1; i<AtomPool; i++ )
        {   ptr->anext = FreeAtom;
            FreeAtom = ptr++;
        } 
    } else FreeAtom = ptr->anext;

    if( CurAtom )
    {   ptr->anext = CurAtom->anext;
        CurAtom->anext = ptr;
    } else 
    {   ptr->anext = CurGroup->alist;
        CurGroup->alist = ptr;
    }
    CurAtom = ptr;

    SelectCount++;
    ptr->flag = SelectFlag | NonBondFlag;
    ptr->label = (void*)0;
    ptr->radius = 375;
    ptr->altl = ' ';
    ptr->mbox = 0;
    ptr->col = 0;

    return( ptr );
}


void ProcessAtom( ptr )
    Atom __far *ptr;
{
    ptr->elemno = GetElemNumber(CurGroup,ptr);
    if( ptr->elemno == 1 )
    {   ptr->flag |= HydrogenFlag;
        HasHydrogen = True;
    }
    if( !(ptr->flag&(HydrogenFlag)) )
        ptr->flag |= NormAtomFlag;

#ifdef INVERT
    ptr->yorg = -ptr->yorg;
#endif

    if( MMinMaxFlag )
    {   if( ptr->xorg < MinX ) 
        {   MinX = ptr->xorg;
        } else if( ptr->xorg > MaxX ) 
            MaxX = ptr->xorg;

        if( ptr->yorg < MinY ) 
        {   MinY = ptr->yorg;
        } else if( ptr->yorg > MaxY ) 
            MaxY = ptr->yorg;

        if( ptr->zorg < MinZ ) 
        {   MinZ = ptr->zorg;
        } else if( ptr->zorg > MaxZ ) 
            MaxZ = ptr->zorg;
    } else 
    {   MinX = MaxX = ptr->xorg;
        MinY = MaxY = ptr->yorg;
        MinZ = MaxZ = ptr->zorg;
    }
            
    MMinMaxFlag = True;
    MainAtomCount++;
}


Atom __far *FindGroupAtom( group, n )
    Group __far *group;  Byte n;
{
    register Atom __far *ptr;

    for( ptr=group->alist; ptr; ptr=ptr->anext )
        if( ptr->refno == n ) return( ptr );
    return( (Atom __far*)0 );
}


int NewAtomType( ptr )
    char *ptr;
{
    register int refno;
    register int i;

    for( refno=0; refno<ElemNo; refno++ )
        if( !strncmp(ElemDesc[refno],ptr,4) )
            return(refno);

    if( ElemNo++ == MAXELEM )
        FatalDataError("Too many new atom types");

    for( i=0; i<4; i++ )
        ElemDesc[refno][i] = ptr[i];
    return( refno );
}


int SimpleAtomType( type )
    char *type;
{
    char name[4];

    name[2] = name[3] = ' ';
    if( type[1] && (type[1]!=' ') )
    {   name[0] = ToUpper(type[0]);
        name[1] = ToUpper(type[1]);
    } else
    {   name[1] = ToUpper(type[0]);
        name[0] = ' ';
    }

    if( name[0]=='H' )
    {   name[1] = name[0];
        name[0] = ' ';
    } else if( name[0]=='C' && name[1]!='L' )
    {   name[1] = name[0];
        name[0] = ' ';
    } else if( name[0]=='N' )
    {   name[1] = name[0];
        name[0] = ' ';
    } else if( name[0]=='O' )
    {   name[1] = name[0];
        name[0] = ' ';
    } else if( name[0]=='P' )
    {   name[1] = name[0];
        name[0] = ' ';
    } else if( name[0]=='S' )
    {   name[1] = name[0];
        name[0] = ' ';
    }

    return( NewAtomType(name) );
}


/*===============================*/
/* Z-Matrix Conversion Functions */
/*===============================*/


#ifdef FUNCPROTO
static IntCoord __far* GetInternalCoord( int );
#endif


void InitInternalCoords()
{
    IntList = (IntCoord __far*)0;
    IntPrev = (IntCoord __far*)0;
}


IntCoord __far* AllocInternalCoord()
{
    register IntCoord __far *ptr;

    ptr = (IntCoord __far*)_fmalloc(sizeof(IntCoord));
    if( !ptr ) FatalDataError("Memory allocation failed");
    ptr->inext = (IntCoord __far*)0;

    if( IntPrev )
    {   IntPrev->inext = ptr;
    } else IntList = ptr;
    IntPrev = ptr;
    return( ptr );
}


static IntCoord __far* GetInternalCoord( index )
    int index;
{
    register IntCoord __far *ptr;

    ptr = IntList;
    while( (index>1) && ptr->inext )
    {   ptr = ptr->inext;
        index--;
    }
    return( ptr );
}


void FreeInternalCoords()
{
    register IntCoord __far *ptr;

    while( (ptr = IntList) )
    {    IntList = ptr->inext;
         _ffree( ptr );
    }
}


int ConvertInternal2Cartesian()
{
    register IntCoord __far *ptr;
    register IntCoord __far *na;
    register IntCoord __far *nb;
    register IntCoord __far *nc;
    register double cosph,sinph,costh,sinth,coskh,sinkh;
    register double cosa,sina,cosd,sind;
    register double dist,angle,dihed;

    register double xpd,ypd,zpd,xqd,yqd,zqd;
    register double xa,ya,za,xb,yb,zb;
    register double rbc,xyb,yza,temp;
    register double xpa,ypa,zqa;
    register double xd,yd,zd;
    register int flag;

    /* Atom #1 */
    ptr = IntList;
    ptr->dist  = 0.0;
    ptr->angle = 0.0;
    ptr->dihed = 0.0;

    if( !(ptr=ptr->inext) )
        return( True );

    /* Atom #2 */
    ptr->angle = 0.0;
    ptr->dihed = 0.0;

    if( !(ptr=ptr->inext) )
        return( True );

    /* Atom #3 */
    dist = ptr->dist;
    angle = Deg2Rad*ptr->angle;
    cosa = cos(angle);
    sina = sin(angle);
    if( ptr->na == 1 )
    {   na = IntList;
        ptr->dist = na->dist + cosa*dist;
    } else /* ptr->na == 2 */
    {   na = IntList->inext;
        ptr->dist = na->dist - cosa*dist;
    }
    ptr->angle = sina*dist;
    ptr->dihed = 0.0;

    while( (ptr=ptr->inext) )
    {   dist = ptr->dist;
        angle = Deg2Rad*ptr->angle;
        dihed = Deg2Rad*ptr->dihed;

        /* Optimise this access?? */
        na = GetInternalCoord(ptr->na);
        nb = GetInternalCoord(ptr->nb);
        nc = GetInternalCoord(ptr->nc);

        xb = nb->dist  - na->dist;
        yb = nb->angle - na->angle;
        zb = nb->dihed - na->dihed;

        rbc = xb*xb + yb*yb + zb*zb;
        if( rbc < 0.0001 )
            return( False );
        rbc= 1.0/sqrt(rbc);

        cosa = cos(angle);
        sina = sin(angle);


        if( fabs(cosa) >= 0.999999 )
        {   /* Colinear */
            temp = dist*rbc*cosa;
            ptr->dist  = na->dist  + temp*xb;
            ptr->angle = na->angle + temp*yb;
            ptr->dihed = na->dihed + temp*zb;
        } else
        {   xa = nc->dist  - na->dist;
            ya = nc->angle - na->angle;
            za = nc->dihed - na->dihed;

            sind = -sin(dihed);
            cosd = cos(dihed);

            xd = dist*cosa;
            yd = dist*sina*cosd;
            zd = dist*sina*sind;

            xyb = sqrt(xb*xb + yb*yb);
            if( xyb < 0.1 )
            {   /* Rotate about y-axis! */
                temp = za; za = -xa; xa = temp;
                temp = zb; zb = -xb; xb = temp;
                xyb = sqrt(xb*xb + yb*yb);
                flag = True;
            } else flag = False;

            costh = xb/xyb;
            sinth = yb/xyb;
            xpa = costh*xa + sinth*ya;
            ypa = costh*ya - sinth*xa;

            sinph = zb*rbc;
            cosph = sqrt(1.0 - sinph*sinph);
            zqa = cosph*za  - sinph*xpa;

            yza = sqrt(ypa*ypa + zqa*zqa);

            if( yza > 1.0E-10 )
            {   coskh = ypa/yza;
                sinkh = zqa/yza;

                ypd = coskh*yd - sinkh*zd;
                zpd = coskh*zd + sinkh*yd;
            } else
            {   /* coskh = 1.0; */
                /* sinkh = 0.0; */
                ypd = yd;
                zpd = zd;
            }

            xpd = cosph*xd  - sinph*zpd;
            zqd = cosph*zpd + sinph*xd;
            xqd = costh*xpd - sinth*ypd;
            yqd = costh*ypd + sinth*xpd;

            if( flag )
            {   /* Rotate about y-axis! */
                ptr->dist  = na->dist  - zqd;
                ptr->angle = na->angle + yqd;
                ptr->dihed = na->dihed + xqd;
            } else
            {   ptr->dist  = na->dist  + xqd;
                ptr->angle = na->angle + yqd;
                ptr->dihed = na->dihed + zqd;
            }
        }
    }
    return( True );
}


/*=========================*/
/* Bond Handling Functions */
/*=========================*/


#ifdef FUNCPROTO
Bond __far *ProcessBond( Atom __far*, Atom __far*, int );
#endif


Bond __far *ProcessBond( src, dst, flag )
    Atom __far *src, __far *dst;
    int flag;
{
    register Bond __far *ptr;
    register int i;

    if( !(ptr = FreeBond) )
    {   MemSize += BondPool*sizeof(Bond);
        ptr = (Bond __far *)_fmalloc( BondPool*sizeof(Bond) );
        if( !ptr ) FatalDataError("Memory allocation failed");
        RegisterAlloc( ptr );
        for( i=1; i<BondPool; i++ )
        {   ptr->bnext = FreeBond;
            FreeBond = ptr++;
        } 
    } else FreeBond = ptr->bnext;

    ptr->flag = flag | SelectFlag;
    ptr->srcatom = src;
    ptr->dstatom = dst;
    ptr->radius = 0;
    ptr->col = 0;

    return( ptr );
}


void CreateBond( src, dst, flag )
    Long src, dst; int flag;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Atom __far *sptr;
    register Atom __far *dptr;
    register Bond __far *bptr;
    register int done;


    if( src == dst )
        return;

    sptr = (void __far*)0;
    dptr = (void __far*)0;

    done = False;
    for( chain=Database->clist; chain && !done; chain=chain->cnext )
        for( group=chain->glist; group && !done; group=group->gnext )
            for( aptr=group->alist; aptr; aptr=aptr->anext )
            {   if( aptr->serno == src )
                {   sptr = aptr;
                    if( dptr )
                    {   done = True;
                        break;
                    }
                } else if( aptr->serno == dst )
                {   dptr = aptr;
                    if( sptr )
                    {   done = True;
                        break;
                    }
                }
            }


    /* Both found! */
    if( done ) 
    {   if( flag != HydrBondFlag )
        {   /* Reset Non-bonded flags! */
            sptr->flag &= ~NonBondFlag;
            dptr->flag &= ~NonBondFlag;

            bptr = ProcessBond( sptr, dptr, flag );
            bptr->bnext = CurMolecule->blist;
            CurMolecule->blist = bptr;
            Info.bondcount++;

        }
    }
}


static void TestBonded( sptr, dptr, flag )
    Atom __far *sptr, __far *dptr; 
    int flag;
{
    register Bond __far *bptr;
    register Long dx, dy, dz;
    register Long max, dist;

    if( flag )
    {    /* Sum of covalent radii with 0.56A tolerance */
         dist = Element[sptr->elemno].covalrad + 
                Element[dptr->elemno].covalrad + 140;
         max = dist*dist;  
    } else 
    {    /* Fast Bio-Macromolecule Bonding Calculation */
         if( (sptr->flag|dptr->flag) & HydrogenFlag )
         {      max = MaxHBondDist;
         } else max = MaxBondDist;
    }

    dx = sptr->xorg-dptr->xorg;   if( (dist=dx*dx)>max ) return;
    dy = sptr->yorg-dptr->yorg;   if( (dist+=dy*dy)>max ) return;
    dz = sptr->zorg-dptr->zorg;   if( (dist+=dz*dz)>max ) return;

    if( dist > MinBondDist )
    {   /* Reset Non-bonded flags! */
        sptr->flag &= ~NonBondFlag;
        dptr->flag &= ~NonBondFlag;

        bptr = ProcessBond(sptr,dptr,NormBondFlag);
        bptr->bnext = CurMolecule->blist;
        CurMolecule->blist = bptr;
        Info.bondcount++;
    }
}


static void ReclaimBonds( ptr )
    Bond __far *ptr;
{
    register Bond __far *temp;

    if( (temp = ptr) )
    {   while( temp->bnext )
            temp=temp->bnext;
        temp->bnext = FreeBond;
        FreeBond = ptr;
    }
}


static Bond __far *ExtractBonds( ptr )
    Bond __far *ptr;
{
    register Bond __far *result;
    register Bond __far *temp;

    result = (Bond __far*)0;

    while( (temp = ptr) )
    {   ptr = temp->bnext;
        if( temp->flag & NormBondFlag )
        {   temp->bnext = FreeBond;
            FreeBond = temp;
        }
    }
    return( result );
}


static void InsertBonds( list, orig )
    Bond __far **list;  Bond __far *orig;
{
    register Atom __far *src;
    register Atom __far *dst;
    register Bond __far *temp;
    register Bond __far *ptr;

    while( (ptr=orig) )
    {   orig = ptr->bnext;
        src = ptr->srcatom;
        dst = ptr->dstatom;
        for( temp=*list; temp; temp=temp->bnext )
            if( ((temp->srcatom==src)&&(temp->dstatom==dst)) ||
                ((temp->srcatom==dst)&&(temp->dstatom==src)) )
                break;

        if( temp )
        {   temp->flag = ptr->flag;
            ptr->bnext = FreeBond;
            FreeBond = ptr;
        } else
        {   ptr->bnext = *list;
            *list = ptr;
        }
    }
}


void CreateMoleculeBonds( info, flag )
    int info, flag;
{
    register int i, x, y, z;
    register Long tx, ty, tz;
    register Long mx, my, mz; 
    register Long dx, dy, dz;
    register int lx, ly, lz, ux, uy, uz;
    register Atom __far *aptr, __far *dptr;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *list;
    char buffer[40];


    if( !Database ) 
        return;

    dx = (MaxX-MinX)+1;
    dy = (MaxY-MinY)+1;
    dz = (MaxZ-MinZ)+1;

    list = ExtractBonds( CurMolecule->blist );
    CurMolecule->blist = (Bond __far*)0;
    Info.bondcount = 0;

    ResetVoxelData();

    for( chain=Database->clist; chain; chain=chain->cnext )
    {   ResetVoxelData();
        for( group=chain->glist; group; group=group->gnext )
            for( aptr=group->alist; aptr; aptr=aptr->anext )
            {   /* Initially non-bonded! */
                aptr->flag |= NonBondFlag;

                mx = aptr->xorg-MinX;
                my = aptr->yorg-MinY;
                mz = aptr->zorg-MinZ;

                tx = mx-AbsMaxBondDist;  
                ty = my-AbsMaxBondDist;  
                tz = mz-AbsMaxBondDist;  

                lx = (tx>0)? (int)((VOXORDER*tx)/dx) : 0;
                ly = (ty>0)? (int)((VOXORDER*ty)/dy) : 0;
                lz = (tz>0)? (int)((VOXORDER*tz)/dz) : 0;

                tx = mx+AbsMaxBondDist;  
                ty = my+AbsMaxBondDist;  
                tz = mz+AbsMaxBondDist;

                ux = (tx<dx)? (int)((VOXORDER*tx)/dx) : VOXORDER-1;
                uy = (ty<dy)? (int)((VOXORDER*ty)/dy) : VOXORDER-1;
                uz = (tz<dz)? (int)((VOXORDER*tz)/dz) : VOXORDER-1;

                for( x=lx; x<=ux; x++ )
                {   i = VOXORDER2*x + VOXORDER*ly;
                    for( y=ly; y<=uy; y++ )
                    {   for( z=lz; z<=uz; z++ )
                            if( (dptr = (Atom __far*)HashTable[i+z]) )
                                do { TestBonded(aptr,dptr,flag);
                                } while( (dptr = dptr->next) );
                        i += VOXORDER;
                    }
                }
                
                x = (int)((VOXORDER*mx)/dx);
        	y = (int)((VOXORDER*my)/dy);
                z = (int)((VOXORDER*mz)/dz);

                i = VOXORDER2*x + VOXORDER*y + z;
                aptr->next = (Atom __far*)HashTable[i];
                HashTable[i] = (void __far*)aptr;
            }
        VoxelsClean = False;
    }

    InsertBonds(&CurMolecule->blist,list);

    if( info )
    {   if( CommandActive )
            WriteChar('\n');
        CommandActive=False;
        sprintf(buffer,"Number of Bonds ..... %ld\n\n",(long)Info.bondcount);
        WriteString(buffer);
    }
}


void RenumberMolecule( start )
    int start;
{
    register Chain __far *chain;
    register Group __far *group;
    register int hinit, minit;
    register int resno;

    if( !Database )
        return;

    hinit = minit = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
    {   resno = start;
        for( group=chain->glist; group; group=group->gnext )
        {   if( minit )
            {   if( resno > MaxMainRes )
                {   MaxMainRes = resno;
                } else if( resno < MinMainRes )
                    MinMainRes = resno;
            } else MinMainRes = MaxMainRes = resno;
            minit = True;
            group->serno = resno++;
        }
    }
}


/*===============================*/
/* Molecule Database Maintenance */
/*===============================*/


static void ReclaimAtoms( ptr )
    Atom __far *ptr;
{
    register Atom __far *temp;

    if( (temp = ptr) )
    {   while( temp->anext )
	    temp=temp->anext;
	temp->anext = FreeAtom;
	FreeAtom = ptr;
    }
}


static void ResetDatabase()
{
    Database = CurMolecule = (void __far*)0;
    MainGroupCount = 0;
    Info.chaincount = 0;
    Info.bondcount = 0;
    MainAtomCount = 0;  
    SelectCount = 0;

    Info.structsource = SourceNone;

    CurGroup = (void __far*)0;
    CurChain = (void __far*)0;
    CurAtom = (void __far*)0;

    MinX = MinY = MinZ = 0;
    MaxX = MaxY = MaxZ = 0;

    MinMainRes = MaxMainRes = 0;

    *Info.moleculename = 0;
    *Info.classification = 0;
    *Info.spacegroup = 0;
    *Info.identcode = 0;
    *Info.filename = 0;

    VoxelsClean = False;
    MMinMaxFlag = False;
    HasHydrogen = False;
    ElemNo = MINELEM;
    ResNo = MINRES;
    MaskCount = 0;
}


void DestroyDatabase()
{
    register void __far *temp;
    register Group __far *gptr;

    if( Database )
    {   ReclaimBonds( Database->blist );

	while( Database->clist )
	{   ReclaimBonds(Database->clist->blist);
	    if( (gptr = Database->clist->glist) )
	    {   ReclaimAtoms(gptr->alist);
		while( gptr->gnext )
		{   gptr = gptr->gnext;
		    ReclaimAtoms(gptr->alist);
		}
		gptr->gnext = FreeGroup;
		FreeGroup = Database->clist->glist;
	    }
	    temp = Database->clist->cnext;
	    Database->clist->cnext = FreeChain;
	    FreeChain = Database->clist;
	    Database->clist = temp;
	}

	FreeMolecule = Database;
	Database = (void __far*)0;
    }
    ResetDatabase();
}


void PurgeDatabase()
{
#ifdef APPLEMAC
    register AllocRef *ptr;
    register AllocRef *tmp;
    register int i;
    
    /* Avoid Memory Leaks */
    for( ptr=AllocList; ptr; ptr=tmp )
    {   for( i=0; i<ptr->count; i++ )
	    _ffree( ptr->data[i] );
	tmp = ptr->next;
	_ffree( ptr );
    }
#endif
}


void InitialiseDatabase()
{
    FreeMolecule = (void __far*)0;
    FreeChain = (void __far*)0;
    FreeGroup = (void __far*)0;
    FreeAtom = (void __far*)0;
    FreeBond = (void __far*)0;

#ifdef APPLEMAC
    AllocList = (void*)0;
#endif

    ResetDatabase();
}
