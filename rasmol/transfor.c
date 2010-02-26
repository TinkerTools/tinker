/* transfor.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#include <Errors.h>
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#endif
#include <stdio.h>
#include <math.h>

#define TRANSFORM
#include "molecule.h"
#include "abstree.h"
#include "transfor.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "graphics.h"

typedef struct {
                short col;
                short shade;
                unsigned char r;
                unsigned char g;
                unsigned char b;
              } ShadeRef;

#define CPKMAX  16
static ShadeRef CPKShade[] = {
     { 0, 0, 200, 200, 200 },       /*  0 Light Grey   */
     { 0, 0, 143, 143, 255 },       /*  1 Sky Blue     */
     { 0, 0, 240,   0,   0 },       /*  2 Red          */
     { 0, 0, 255, 200,  50 },       /*  3 Yellow       */
     { 0, 0, 255, 255, 255 },       /*  4 White        */
     { 0, 0, 255, 192, 203 },       /*  5 Pink         */
     { 0, 0, 218, 165,  32 },       /*  6 Golden Rod   */
     { 0, 0,   0,   0, 255 },       /*  7 Blue         */
     { 0, 0, 255, 165,   0 },       /*  8 Orange       */
     { 0, 0, 128, 128, 144 },       /*  9 Dark Grey    */
     { 0, 0, 165,  42,  42 },       /* 10 Brown        */
     { 0, 0, 160,  32, 240 },       /* 11 Purple       */
     { 0, 0, 255,  20, 147 },       /* 12 Deep Pink    */
     { 0, 0,   0, 255,   0 },       /* 13 Green        */
     { 0, 0, 178,  34,  34 },       /* 14 Fire Brick   */
     { 0, 0,  34, 139,  34 } };     /* 15 Forest Green */

#define SHAPELYMAX  15
static ShadeRef Shapely[] = {
     { 0, 0, 255,   0,   0 },
     { 0, 0, 236,   0,  19 },
     { 0, 0, 217,   0,  38 },
     { 0, 0, 199,   0,  56 },
     { 0, 0, 181,   0,  74 },
     { 0, 0, 163,   0,  92 },
     { 0, 0, 145,   0, 110 },
     { 0, 0, 127,   0, 128 },
     { 0, 0, 110,   0, 145 },
     { 0, 0,  92,   0, 163 },
     { 0, 0,  74,   0, 181 },
     { 0, 0,  56,   0, 199 },
     { 0, 0,  38,   0, 217 },
     { 0, 0,  19,   0, 236 },
     { 0, 0,   0,   0, 255 } };

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(ptr=group->alist;ptr;ptr=ptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext) 
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)

#define MatchChar(a,b)   (((a)=='#')||((a)==(b)))

static ShadeRef ScaleRef[LastShade];
static int MaskColour[MAXMASK];
static int MaskShade[MAXMASK];
static int ScaleCount;

static Real LastRX,LastRY,LastRZ;
static Real Zoom;


void DetermineClipping()
{
    register int temp;
    register int max;

    max = 0;
    if( DrawAtoms && (MaxAtomRadius>max) )  max = MaxAtomRadius;
    if( DrawBonds && (MaxBondRadius>max) )  max = MaxBondRadius;
       
    temp = ImageRadius + max;
    if( (YOffset>=temp) && (XOffset>=temp) && (YOffset+temp<YRange) )
    {   if( UseStereo )
        {   UseScreenClip = (XOffset+temp) >= (XRange>>1);
        } else UseScreenClip = (XOffset+temp) >= XRange;
    } else UseScreenClip = True;
}


void SetRadiusValue( rad )
    int rad;
{
    register int irad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    irad = (int)(Scale*rad);
    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag & SelectFlag )
        {   if( irad>MaxAtomRadius )
                MaxAtomRadius = irad;
            ptr->flag |= SphereFlag;
            ptr->radius = rad;
            ptr->irad = irad;
            change = True;
        } else if( ptr->flag & SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void SetVanWaalRadius()
{
    register int rad,change;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    change = False;

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   rad = ElemVDWRadius(ptr->elemno);
            ptr->irad = (int)(Scale*rad);
            ptr->radius = rad;
            change = True;

            ptr->flag |=SphereFlag;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        } else if( ptr->flag&SphereFlag )
        {   DrawAtoms = True;
            if( ptr->irad>MaxAtomRadius )
                MaxAtomRadius = ptr->irad;
        }

    if( change )
    {   DrawAtoms = True;
        DetermineClipping();
        VoxelsClean = False;
        BucketFlag = False;
    }
}


void DisableSpacefill()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database || !DrawAtoms )
        return;

    MaxAtomRadius = 0;
    DrawAtoms = False;
    
    ForEachAtom
        if( !(ptr->flag&SelectFlag) )
        {   if( ptr->flag&SphereFlag )
            {   if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
                DrawAtoms = True;
            }
        } else if( ptr->flag&SphereFlag )
            ptr->flag &= ~SphereFlag;

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


void EnableWireframe( mask, rad )
    int mask, rad;
{
    register Bond __far *bptr;
    register int flag, irad;

    if( !Database )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;
    irad = (int)(Scale*rad);

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   DrawBonds = True;
            bptr->flag &= ~DrawBondFlag;
            bptr->flag |= mask;
            if( mask == CylinderFlag )
            {   if( irad>MaxBondRadius )
                    MaxBondRadius = irad;
                bptr->radius = rad;
                bptr->irad = irad;
            }
        } else if( bptr->flag&DrawBondFlag )
        {    DrawBonds = True;
             if( bptr->flag&CylinderFlag )
                 if( bptr->irad>MaxBondRadius )
                     MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


void DisableWireframe()
{
    register Bond __far *bptr;
    register int flag;

    if( !Database || !DrawBonds )
        return;

    DrawBonds = False;
    MaxBondRadius = 0;

    ForEachBond
    {   flag = ZoneBoth? bptr->dstatom->flag & bptr->srcatom->flag
                       : bptr->dstatom->flag | bptr->srcatom->flag;

        if( flag&SelectFlag )
        {   bptr->flag &= ~DrawBondFlag;
        } else if( bptr->flag&DrawBondFlag )
        {   DrawBonds = True;
            if( bptr->flag&CylinderFlag )
                if( bptr->irad>MaxBondRadius )
                    MaxBondRadius = bptr->irad;
        }
    }
    DetermineClipping();
}


/*===========================*/
/* Atom Selection Functions! */
/*===========================*/

static void DisplaySelectCount()
{
    char buffer[40];

    if( FileDepth == -1 )
    {   if( CommandActive )
           WriteChar('\n');
        CommandActive=False;

        if( SelectCount==0 )
        {   WriteString("No atoms selected!\n");
        } else if( SelectCount>1 )
        {   sprintf(buffer,"%ld atoms selected!\n",(long)SelectCount);
            WriteString(buffer);
        } else WriteString("1 atom selected!\n");
    }

    if( DisplayMode )
        ReDrawFlag |= RFRefresh;
    AdviseUpdate(AdvSelectCount);
}


void SelectZone( mask )
    int mask;
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( !Database )
        return;

    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;
        } else ptr->flag &= ~SelectFlag;
    DisplaySelectCount();

    if( ZoneBoth )
    {   ForEachBond
           if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
    } else
        ForEachBond
           if( (bptr->srcatom->flag|bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;

}


void RestrictZone( mask )
    int mask;
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int flag;

    if( !Database )
        return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;
    DrawLabels = False;
    
    SelectCount = 0;
    ForEachAtom
        if( ptr->flag & mask )
        {   ptr->flag |= SelectFlag;
            SelectCount++;

            if( ptr->flag & SphereFlag )
            {   DrawAtoms = True;
                if( ptr->irad>MaxAtomRadius )
                    MaxAtomRadius = ptr->irad;
            }

            if( ptr->label )
                DrawLabels = True;
        } else 
        {   ptr->flag &= ~(SelectFlag|SphereFlag);
            if( ptr->label )
            {   DeleteLabel( (Label*)ptr->label );
                ptr->label = (void*)0;
            }
        }
    DisplaySelectCount();
    
    ForEachBond
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( flag & SelectFlag )
        {   bptr->flag |= SelectFlag;
            if( bptr->flag&DrawBondFlag )
            {   DrawBonds = True;
                if( bptr->flag & CylinderFlag )
                    if( bptr->irad>MaxBondRadius )
                        MaxBondRadius = bptr->irad;
            } 
        } else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    ForEachBack
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( !(flag&SelectFlag) )
            bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


void SelectZoneExpr( expr )
    Expr *expr;
{
    register Bond __far *bptr;

    if( !Database )
        return;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
        for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
            for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                if( EvaluateExpr(expr) )
                {   QAtom->flag |= SelectFlag;
                    SelectCount++;
                } else QAtom->flag &= ~SelectFlag;
    DisplaySelectCount();

    if( ZoneBoth )
    {   ForEachBond
           if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
    } else
        ForEachBond
           if( (bptr->srcatom->flag|bptr->dstatom->flag) & SelectFlag )
           {   bptr->flag |= SelectFlag;
           } else bptr->flag &= ~SelectFlag;
}


void RestrictZoneExpr( expr )
    Expr *expr;
{
    register Bond __far *bptr;
    register Chain __far *chain;
    register int flag;

    if( !Database )
        return;

    DrawAtoms = False;   MaxAtomRadius = 0;
    DrawBonds = False;   MaxBondRadius = 0;
    DrawLabels = False;

    SelectCount = 0;
    for( QChain=Database->clist; QChain; QChain=QChain->cnext )
        for( QGroup=QChain->glist; QGroup; QGroup=QGroup->gnext )
            for( QAtom=QGroup->alist; QAtom; QAtom=QAtom->anext )
                if( EvaluateExpr(expr) )
                {   QAtom->flag |= SelectFlag;
                    SelectCount++;

                    if( QAtom->flag & SphereFlag )
                    {   DrawAtoms = True;
                        if( QAtom->irad>MaxAtomRadius )
                            MaxAtomRadius = QAtom->irad;
                    }
                    if( QAtom->label )
                        DrawLabels = True;

                }  else 
                {   QAtom->flag &= ~(SelectFlag|SphereFlag);
                    if( QAtom->label )
                    {   DeleteLabel( (Label*)QAtom->label );
                        QAtom->label = (void*)0;
                    }
                }
    DisplaySelectCount();

    ForEachBond
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( flag & SelectFlag )
        {   bptr->flag |= SelectFlag;
            if( bptr->flag & CylinderFlag )
            {   DrawBonds = True;
                if( bptr->irad>MaxBondRadius )
                    MaxBondRadius = bptr->irad;
            } else if( bptr->flag&WireFlag )
                DrawBonds = True;
        } else bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    ForEachBack
    {   /* Ignore ZoneBoth setting! */
        flag = bptr->dstatom->flag & bptr->srcatom->flag;
        if( !(flag&SelectFlag) )
            bptr->flag &= ~(SelectFlag|DrawBondFlag);
    }

    DetermineClipping();
    VoxelsClean = False;
    BucketFlag = False;
}


int DefineShade( r, g, b )
    unsigned char r, g, b;
{
    register int d,dr,dg,db;
    register int dist,best;
    register int i;

    /* Already defined! */
    for( i=0; i<LastShade; i++ )
        if( Shade[i].refcount )
            if( (Shade[i].r==r)&&(Shade[i].g==g)&&(Shade[i].b==b) )
                return(i);

    /* Allocate request */
    for( i=0; i<LastShade; i++ )
         if( !Shade[i].refcount )
         {   Shade[i].r = r;
             Shade[i].g = g;
             Shade[i].b = b;
             Shade[i].refcount = 0;
             return(i);
         }

    if( CommandActive )
        WriteChar('\n');
    WriteString("Warning: Unable to allocate shade!\n");
    CommandActive = False;

    /* To avoid lint warning! */
    best = dist = 0;

    /* Nearest match */
    for( i=0; i<LastShade; i++ )
    {   dr = Shade[i].r - r;
        dg = Shade[i].g - g;
        db = Shade[i].b - b;
        d = dr*dr + dg*dg + db*db;
        if( !i || (d<dist) )
        {   dist = d;
            best = i;
        }
    }
    return( best );
}


static void SetLutEntry( i, r, g, b )
    int i, r, g, b;
{
    ULut[i] = True;
    RLut[i] = r;
    GLut[i] = g;
    BLut[i] = b;

#ifdef EIGHTBIT
    Lut[i] = i;
#else
    Lut[i] = ( (Card)((r<<8)|g)<<8 ) | b;
#endif
}


static Real Power( x, y )
    Real x; int y;
{
    register Real result;

    result = x;
    while( y>1 )
    {   if( y&1 ) { result *= x; y--; }
        else { result *= result; y>>=1; }
    }
    return( result );
}


void DefineColourMap()
{
    register Real diffuse,fade;
    register Real temp,inten;
    register int col,r,g,b;
    register int i,j,k;

    for( i=0; i<LutSize; i++ )
        ULut[i] = False;

    if( !DisplayMode )
    {   SetLutEntry(BackCol,BackR,BackG,BackB);
        SetLutEntry(LabelCol,LabR,LabG,LabB);
        SetLutEntry(BoxCol,BoxR,BoxG,BoxB);
    } else SetLutEntry(BackCol,80,80,80);


    diffuse = 1.0 - Ambient;
    if( DisplayMode )
    {   for( i=0; i<ColourDepth; i++ )
        {   temp = (Real)i/ColourMask;
            inten = diffuse*temp + Ambient;

            /* Unselected [40,40,255] */
            /* Selected   [255,160,0]  */
            r = (int)(255*inten);
            g = (int)(160*inten);
            b = (int)(40*inten);

            /* Avoid Borland Compiler Warning! */
            /* Shade2Colour(0) == FirstCol     */
            SetLutEntry( FirstCol+i, b, b, r );
            SetLutEntry( Shade2Colour(1)+i, r, g, 0 );
        }
    } else
        for( i=0; i<ColourDepth; i++ )
        {   temp = (Real)i/ColourMask;
            inten = diffuse*temp + Ambient;
            fade = 1.0-inten;

            if( FakeSpecular )
            {   temp = Power(temp,SpecPower);
                k = (int)(255*temp);
                temp = 1.0 - temp;
                inten *= temp;
                fade *= temp;
            }

            for( j=0; j<LastShade; j++ )
                if( Shade[j].refcount )
                {   col = Shade2Colour(j);
                    if( UseBackFade )
                    {   temp = 1.0-inten;
                        r = (int)(Shade[j].r*inten + fade*BackR); 
                        g = (int)(Shade[j].g*inten + fade*BackG);
                        b = (int)(Shade[j].b*inten + fade*BackB);
                    } else
                    {   r = (int)(Shade[j].r*inten); 
                        g = (int)(Shade[j].g*inten);
                        b = (int)(Shade[j].b*inten);
                    }

                    if( FakeSpecular )
                    {   r += k;
                        g += k;
                        b += k;
                    }
                    SetLutEntry( col+i, r, g, b );
                }
        }

    if( Interactive )
        AllocateColourMap();
}


void ResetColourMap()
{
    register int i;

#ifdef EIGHTBIT
    for( i=0; i<256; i++ )
        ULut[i] = False;
#endif

    SpecPower = 8;
    FakeSpecular = False;
    Ambient = DefaultAmbient;
    UseBackFade = False;

    BackR = BackG = BackB = 0;
    BoxR = BoxG = BoxB = 255;
    LabR = LabG = LabB = 255;

    for( i=0; i<LastShade; i++ )
        Shade[i].refcount = 0;
    ScaleCount = 0;
}


void ColourBondNone()
{
    register Bond __far *bptr;

    if( Database )
        ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
            {   Shade[Colour2Shade(bptr->col)].refcount--;
                bptr->col = 0;
            }
}


void ColourBondAttrib( r, g, b )
    int r, g, b;
{
    register Bond __far *bptr;
    register int shade,col;

    if( Database )
    {   ForEachBond
            if( (bptr->flag&SelectFlag) && bptr->col )
                Shade[Colour2Shade(bptr->col)].refcount--;

        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        ForEachBond
            if( bptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                bptr->col = col;
            }
    }
}


void ColourMonitNone()
{
    register Monitor *ptr;
    register int flag;

    if( Database )
        for( ptr=MonitList; ptr; ptr=ptr->next )
            if( ptr->col )
            {   flag = ZoneBoth? ptr->src->flag & ptr->dst->flag
                               : ptr->src->flag | ptr->dst->flag;
                if( flag & SelectFlag )
                {   Shade[Colour2Shade(ptr->col)].refcount--;
                    ptr->col = 0;
                }
            }
}


void ColourMonitAttrib( r, g, b )
    int r, g, b;
{
    register Monitor *ptr;
    register int shade,col;
    register int flag;

    if( !Database )
        return;

    ColourMonitNone();
    shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
    col = Shade2Colour(shade);

    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   flag = ZoneBoth? ptr->src->flag & ptr->dst->flag 
                       : ptr->src->flag | ptr->dst->flag;
        if( flag & SelectFlag )
        {   Shade[shade].refcount++;
            ptr->col = col;
        }
    }
}


static void ResetColourAttrib()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    ForEachAtom
        if( (ptr->flag&SelectFlag) && ptr->col )
            Shade[Colour2Shade(ptr->col)].refcount--;
}


void MonoColourAttrib( r, g, b )
    int r, g, b;
{
    register int shade,col;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;

    if( Database )
    {   ResetColourAttrib();
        shade = DefineShade((Byte)r,(Byte)g,(Byte)b);
        col = Shade2Colour(shade);

        ForEachAtom
            if( ptr->flag&SelectFlag )
            {   Shade[shade].refcount++;
                ptr->col = col;
            }
    }
}


void CPKColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;

    if( !Database ) return;
    for( i=0; i<CPKMAX; i++ )
        CPKShade[i].col = 0;
    ResetColourAttrib();

    ForEachAtom
        if( ptr->flag&SelectFlag )
        {   ref = CPKShade + Element[ptr->elemno].cpkcol;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


void ShapelyColourAttrib()
{
    register ShadeRef *ref;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register int i;
    register int atoms;
    register int offset;

    if( !Database ) return;
    for( i=0; i<SHAPELYMAX; i++ )
        Shapely[i].col = 0;
    ResetColourAttrib();

    atoms=0;
    ForEachAtom
        if( ptr->flag&SelectFlag )  atoms=atoms+1;

    i=0;
    ForEachAtom
        if( ptr->flag&SelectFlag )

        {   offset = (SHAPELYMAX*i)/atoms;
            ref = Shapely+offset;
            i=i+1;

            if( !ref->col )
            {   ref->shade = DefineShade( ref->r, ref->g, ref->b );
                ref->col = Shade2Colour(ref->shade);
            }
            Shade[ref->shade].refcount++;
            ptr->col = ref->col;
        }
}


int IsCPKColour( ptr )
    Atom __far *ptr;
{
    register ShadeRef *cpk;
    register ShadeDesc *col;

    cpk = CPKShade + Element[ptr->elemno].cpkcol;
    col = Shade + Colour2Shade(ptr->col);
    return( (col->r==cpk->r) && 
            (col->g==cpk->g) && 
            (col->b==cpk->b) );
}


int IsVDWRadius( ptr )
    Atom __far *ptr;
{
    register int rad;

    if( ptr->flag & SphereFlag )
    {   rad = ElemVDWRadius( ptr->elemno );
        return( ptr->radius == rad );
    } else return( False );
}


void DefaultRepresentation()
{
    if( Database )
    {   ReDrawFlag |= RFRefresh | RFColour;
        if( Info.bondcount < 1 )
        {   SetVanWaalRadius();

        } else if( MainAtomCount<256 )
               {   EnableWireframe(CylinderFlag,40);
               } else EnableWireframe(CylinderFlag,80);

        CPKColourAttrib();
    }
}


void InitialTransform()
{
    register Card dist,max;
    register double fdist,fmax;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *ptr;
    register Card ax, ay, az;
    register Long dx, dy, dz;


    dx = MaxX-MinX;   OrigCX = (dx>>1)+MinX;
    dy = MaxY-MinY;   OrigCY = (dy>>1)+MinY;
    dz = MaxZ-MinZ;   OrigCZ = (dz>>1)+MinZ;

    MaxX -= OrigCX;   MinX -= OrigCX;
    MaxY -= OrigCY;   MinY -= OrigCY;
    MaxZ -= OrigCZ;   MinZ -= OrigCZ;

    SideLen = MaxFun(dx,dy);
    if( dz>SideLen ) SideLen = dz;
    SideLen += 1500;  Offset = SideLen>>1;
    XOffset = WRange;  YOffset = HRange;
    ZOffset = 10000;

    ForEachAtom
    {   ptr->xorg -= OrigCX;
        ptr->yorg -= OrigCY;
        ptr->zorg -= OrigCZ;
    }

    if( Offset > 37836 )
    {   fmax = 0.0;
        ForEachAtom
        {   ax = (Card)AbsFun(ptr->xorg);
            ay = (Card)AbsFun(ptr->yorg);
            az = (Card)AbsFun(ptr->zorg);
            fdist = (double)ax*ax + 
                    (double)ay*ay + 
                    (double)az*az;
            if( fdist > fmax )
                fmax = fdist;
        }
    } else
    {   max = 1;
        ForEachAtom
        {   ax = (Card)AbsFun(ptr->xorg);
            ay = (Card)AbsFun(ptr->yorg);
            az = (Card)AbsFun(ptr->zorg);
            dist = ax*ax + ay*ay + az*az;
            if( dist > max )
                max = dist;
        }
        fmax = (double)max;
    }


    WorldRadius = (Card)sqrt(fmax);
    WorldSize = WorldRadius<<1;
    DScale = 1.0/(WorldSize+1500);

    /* Code should match ReSizeScreen() */
    /* MaxZoom*DScale*Range*750 == 252  */
    MaxZoom = 0.336*(WorldSize+1500)/Range;
    if( MaxZoom < 1.0 )
    {   DScale *= MaxZoom;
        MaxZoom = 1.0;
    }
    ZoomRange = Range;
    MaxZoom -= 1.0;
}


void ReviseInvMatrix()
{
    /* The inverse of a rotation matrix
     * is its transpose, and the inverse
     * of Scale is 1.0/Scale [IScale]!
     */
    InvX[0] = IScale*RotX[0];
    InvX[1] = IScale*RotY[0];
    InvX[2] = IScale*RotZ[0];

    InvY[0] = IScale*RotX[1];
    InvY[1] = IScale*RotY[1];
    InvY[2] = IScale*RotZ[1];

    InvZ[0] = IScale*RotX[2];
    InvZ[1] = IScale*RotY[2];
    InvZ[2] = IScale*RotZ[2];
    ShadowTransform();
}


void PrepareTransform()
{
    register Real theta, temp;
    register Real cost, sint;
    register Real x, y, z;
    register Real ncost;

    if( (ReDrawFlag&RFRotateX) && (DialValue[0]!=LastRX) )
    {   theta = PI*(DialValue[0]-LastRX);
        cost = cos(theta);  sint = sin(theta);
        LastRX = DialValue[0];

        y=RotY[0]; z=RotZ[0];
        RotY[0]=cost*y+sint*z; 
        RotZ[0]=cost*z-sint*y;

        y=RotY[1]; z=RotZ[1];
        RotY[1]=cost*y+sint*z;
        RotZ[1]=cost*z-sint*y;

        y=RotY[2]; z=RotZ[2];
        RotY[2]=cost*y+sint*z;
        RotZ[2]=cost*z-sint*y;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotZ[2]);
            temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateY) && (DialValue[1]!=LastRY) )
    {   theta = PI*(DialValue[1]-LastRY);
        cost = cos(theta);  sint = sin(theta);
        LastRY = DialValue[1];

        x=RotX[0]; z=RotZ[0];
        RotX[0]=cost*x+sint*z;
        RotZ[0]=cost*z-sint*x;

        x=RotX[1]; z=RotZ[1];
        RotX[1]=cost*x+sint*z;
        RotZ[1]=cost*z-sint*x;

        x=RotX[2]; z=RotZ[2];
        RotX[2]=cost*x+sint*z;
        RotZ[2]=cost*z-sint*x;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] + sint*RotZ[0]);
            temp += CenY*(ncost*RotX[1] + sint*RotZ[1]);
            temp += CenZ*(ncost*RotX[2] + sint*RotZ[2]);
            temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;
        }
    }

    if( (ReDrawFlag&RFRotateZ) && (DialValue[2]!=LastRZ) )
    {   theta = PI*(DialValue[2]-LastRZ);
        cost = cos(theta);  sint = sin(theta);
        LastRZ = DialValue[2];

        x=RotX[0]; y=RotY[0];
        RotX[0]=cost*x-sint*y;
        RotY[0]=cost*y+sint*x;

        x=RotX[1]; y=RotY[1];
        RotX[1]=cost*x-sint*y;
        RotY[1]=cost*y+sint*x;

        x=RotX[2]; y=RotY[2];
        RotX[2]=cost*x-sint*y;
        RotY[2]=cost*y+sint*x;

        if( CenX || CenY || CenZ )
        {   ncost = 1.0-cost;
            temp =  CenX*(ncost*RotX[0] - sint*RotY[0]);
            temp += CenY*(ncost*RotX[1] - sint*RotY[1]);
            temp += CenZ*(ncost*RotX[2] - sint*RotY[2]);
            temp = DialValue[4] - (Scale*temp)/XRange;

            if( temp < -1.0 )
            {   DialValue[4] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[4] = 1.0;
            } else DialValue[4] = temp;

            ncost = 1.0-cost;
            temp =  CenX*(ncost*RotY[0] + sint*RotX[0]);
            temp += CenY*(ncost*RotY[1] + sint*RotX[1]);
            temp += CenZ*(ncost*RotY[2] + sint*RotX[2]);
            temp = DialValue[5] - (Scale*temp)/YRange;

            if( temp < -1.0 )
            {   DialValue[5] = -1.0;
            } else if( temp > 1.0 )
            {   DialValue[5] = 1.0;
            } else DialValue[5] = temp;
        }
    }
}


void ApplyTransform()
{
    register int temp;
    register Real x, y, z;
    register int oldx,oldy;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *ptr;


    if( ReDrawFlag & RFMagnify )
    {   if( DialValue[3] <= 0.0 )
        {   Zoom = DialValue[3]+1.0;
            if( Zoom<0.1 ) Zoom=0.1;
        } else Zoom = (DialValue[3]*MaxZoom) + 1.0;

        Scale = Zoom*DScale*Range;
        ImageSize = (int)(Scale*WorldSize);
        if( ImageSize < 2 )
        {   ImageRadius = 1;
            ImageSize = 2;
        } else 
            ImageRadius = ImageSize>>1;
        IScale = 1.0/Scale;

        MaxAtomRadius = 0;
        MaxBondRadius = 0;
    }

    if( ReDrawFlag & RFRotate )
    {   PrepareTransform();
        if( UseShadow )
            ShadowTransform();
    }

    if( ReDrawFlag & (RFRotate|RFMagnify) )
    {   MatX[0] = Scale*RotX[0]; 
        MatX[1] = Scale*RotX[1];
        MatX[2] = Scale*RotX[2];

        MatY[0] = Scale*RotY[0];
        MatY[1] = Scale*RotY[1];
        MatY[2] = Scale*RotY[2];

        MatZ[0] = Scale*RotZ[0];
        MatZ[1] = Scale*RotZ[1];
        MatZ[2] = Scale*RotZ[2];

        if( UseShadow )
        {   InvX[0] = IScale*RotX[0]; 
            InvX[1] = IScale*RotY[0];
            InvX[2] = IScale*RotZ[0];

            InvY[0] = IScale*RotX[1];
            InvY[1] = IScale*RotY[1];
            InvY[2] = IScale*RotZ[1];

            InvZ[0] = IScale*RotX[2];
            InvZ[1] = IScale*RotY[2];
            InvZ[2] = IScale*RotZ[2];
        }
    }

    oldx = XOffset;
    oldy = YOffset;
    XOffset = WRange + (int)(DialValue[4]*XRange);
    YOffset = HRange + (int)(DialValue[5]*YRange);
    if( UseStereo ) XOffset /= 2;

    /* Zoom dependent Translation! */
    /* XOffset = WRange + (int)(DialValue[4]*ImageSize); */
    /* YOffset = HRange + (int)(DialValue[5]*ImageSize); */


    switch( ReDrawFlag )
    {   case(RFTransX):
                if( XOffset != oldx ) 
                {   temp = XOffset - oldx;
                    ForEachAtom ptr->x += temp;
                }
                break;

        case(RFTransY):
                if( YOffset != oldy ) 
                {   temp = YOffset - oldy;
                    ForEachAtom ptr->y += temp;
                }
                break;

        case(RFRotateX):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
            }
            break;

        case(RFRotateY):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
            }
            break;

        case(RFRotateZ):
            ForEachAtom
            {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
            }
            break;

        default:
            /* This condition scales atomic radii! */
            if( DrawAtoms && (ReDrawFlag&RFMagnify) )
            {   ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
                    if( ptr->flag&SphereFlag )
                    {   ptr->irad = (int)(Scale*ptr->radius);
                        if( ptr->irad>MaxAtomRadius )
                            MaxAtomRadius = ptr->irad;
                    }
                }
            } else
                ForEachAtom 
                {   x = ptr->xorg; y = ptr->yorg; z = ptr->zorg;
                    ptr->x = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
                    ptr->y = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                    ptr->z = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
                }

            if( ReDrawFlag & RFMagnify )
            {   if( DrawBonds )
                    ForEachBond
                        if( bptr->flag&CylinderFlag )
                        {   bptr->irad = (int)(Scale*bptr->radius);
                            if( bptr->irad>MaxBondRadius )
                            MaxBondRadius = bptr->irad;
                        }

                ForEachBack
                    if( bptr->flag&CylinderFlag )
                        bptr->irad = (int)(Scale*bptr->radius);
            }
    }

    DetermineClipping();
    if( UseScreenClip || ReDrawFlag!=RFRotateY )
        BucketFlag = False;
}


void ResetTransform()
{
    RotX[0] = 1.0;  RotX[1] = 0.0;  RotX[2] = 0.0;
    RotY[0] = 0.0;  RotY[1] = 1.0;  RotY[2] = 0.0;
    RotZ[0] = 0.0;  RotZ[1] = 0.0;  RotZ[2] = 1.0;
    LastRX = LastRY = LastRZ = 0.0;
    CenX = CenY = CenZ = 0;
}


void InitialiseTransform()
{
    ResetColourMap();
    ResetTransform();

    ZoneBoth = True;
    Hydrogens = True;
}
