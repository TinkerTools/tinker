/* script.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#define SCRIPT
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
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
#ifndef sun386
#include <stdlib.h>
#endif
 
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "script.h"
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "repres.h"
#include "graphics.h"
#include "pixutils.h"

#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif

#define Round(x)       ((int)(x))
#define DatWirFlag  (Long)0x10000
#define DatDasFlag  (Long)0x20000
#define DatCylFlag  (Long)0x40000

typedef struct {
        Long datum;
        Long count;
    } FreqEntry;

#define FREQSIZE  8
static FreqEntry Freq[FREQSIZE];


static Atom __far *MagePrev;
static char *MageCol;
static FILE *OutFile;
static int SelectAll;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)


static void FatalScriptError( ptr )
    char *ptr;
{
    if( CommandActive ) WriteChar('\n');
    WriteString("Script Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
    CommandActive = False;
}


#ifdef FUNCPROTO
static void IncFreqTable( Long );
static Long GetBondDatum( Bond __far* );
#endif


static void ResetFreqTable()
{
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        Freq[i].count = 0;
}


static void IncFreqTable( datum )
    Long datum;
{
    register Long count;
    register int i;

    for( i=0; i<FREQSIZE; i++ )
        if( !Freq[i].count )
        {   Freq[i].datum = datum;
            Freq[i].count = 1;
            return;
        } else if( Freq[i].datum == datum )
        {   count = Freq[i].count+1;
            while( i && (Freq[i-1].count<=count) )
            {   Freq[i] = Freq[i-1];  
                i--;
            }
            Freq[i].datum = datum;
            Freq[i].count = count;
            return;
        }

    /* Replace Singletons! */
    if( Freq[FREQSIZE-1].count == 1 )
        Freq[FREQSIZE-1].datum = datum;
}


static Long GetBondDatum( bptr )
    Bond __far *bptr;
{
    if( bptr->flag & CylinderFlag )
    {   return( DatCylFlag | bptr->radius );
    } else if( bptr->flag & WireFlag )
    {   return( DatWirFlag );
    } else if( bptr->flag & DashFlag )
    {   return( DatDasFlag );
    } else return( (Long)0 );
}


#ifdef FUNCPROTO
static void WriteScriptDatum( char*, Long );
static void WriteScriptSelectBond( Atom __far*, Atom __far* );
#endif


static void WriteScriptAll()
{
    if( !SelectAll )
    {   fputs("select all\n",OutFile);
        SelectAll = True;
    }
}


static void WriteScriptColour( ptr, col )
    char *ptr;  int col;
{
    register ShadeDesc *shade;
    
    if( col )
    {   shade = Shade + Colour2Shade(col);
        fprintf(OutFile,"colour %s [%d,%d,%d]\n",ptr,
                shade->r,shade->g,shade->b);
    } else fprintf(OutFile,"colour %s none\n",ptr);
}


static void WriteScriptBetween( lo, hi )
    int lo, hi;
{
    if( lo != hi )
    {   fprintf(OutFile,"select (atomno>=%d) and (atomno<=%d)\n",lo,hi);
    } else fprintf(OutFile,"select atomno=%d\n",lo);
    SelectAll = False;
}


static void WriteScriptSelectBond( src, dst )
    Atom __far *src, __far *dst;
{
    fprintf(OutFile,"select (atomno=%d) or (atomno==%d)\n",
                    src->serno, dst->serno);
    SelectAll = False;
}


static void WriteScriptAtoms()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Long first,last;
    register int same,init;
    register int cpk,vdw;
    register int col,rad;

    fputs("\n# Atoms\n",OutFile);

    same = True;
    init = False;
    ForEachAtom
        if( !init )
        {   first = last = aptr->serno;
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            init = True;
        } else if( cpk && IsCPKColour(aptr) )
        {   last = aptr->serno;
            if( aptr->col != col )
                col = 0;
        } else if( aptr->col == col )
        {   last = aptr->serno;
            cpk = False;
        } else if( aptr->col != col )
        {   WriteScriptBetween( first, last );
            if( !col )
            {   fputs("colour atoms cpk\n",OutFile);
            } else WriteScriptColour("atoms",col);
                
            first = last = aptr->serno;
            cpk = IsCPKColour( aptr );
            col = aptr->col;
            same = False;
        } else last = aptr->serno; 
        
    if( init )
    {   if( !same )
        {   WriteScriptBetween(first,last);
        } else WriteScriptAll();

        if( !col )
        {   fputs("colour atoms cpk\n",OutFile);
        } else WriteScriptColour("atoms",col);
    }

    if( DrawAtoms )
    {   same = True;
        init = False;
        ForEachAtom
            if( !init )
            {   rad = aptr->flag&SphereFlag? aptr->radius : 0;
                first = last = aptr->serno;
                vdw = IsVDWRadius( aptr );
                init = True;
            } else if( rad == ((aptr->flag&SphereFlag)? aptr->radius : 0) )
            {   if( vdw ) vdw = IsVDWRadius( aptr );
                last = aptr->serno;
            } else if( vdw && IsVDWRadius(aptr) )
            {   last = aptr->serno;
                rad = -1;
            } else 
            {   WriteScriptBetween(first,last);
                if( rad == -1 )
                {   fputs("spacefill on\n",OutFile);
                } else if( rad )
                {   fprintf(OutFile,"spacefill %d\n",rad);
                } else fputs("spacefill off\n",OutFile); 

                rad = aptr->flag&SphereFlag? aptr->radius : 0;
                first = last = aptr->serno;
                vdw = IsVDWRadius( aptr );
                same = False;
            }

        if( !same )
        {   WriteScriptBetween(first,last);
        } else WriteScriptAll();

        if( rad == -1 )
        {   fputs("spacefill on\n",OutFile);
        } else if( rad )
        {   fprintf(OutFile,"spacefill %d\n",rad);
        } else fputs("spacefill off\n",OutFile); 

        if( UseShadow )
        {   fputs("set shadow on\n",OutFile);
        } else fputs("set shadow off\n",OutFile);

    } else
    {   WriteScriptAll();
        fputs("spacefill off\n",OutFile);
    }
        
}


static void WriteScriptDatum( ptr, datum )
    char *ptr;  Long datum;
{
    if( datum & DatCylFlag )
    {   fprintf(OutFile,"%s %d\n",ptr,(int)(datum-DatCylFlag));
    } else if( datum & DatWirFlag )
    {   fprintf(OutFile,"%s on\n",ptr);
    } else if( datum & DatDasFlag )
    {   fprintf(OutFile,"%s dash\n",ptr);
    } else fprintf(OutFile,"%s off\n",ptr);
}


static void WriteScriptBonds()
{
    register Bond __far *bptr;
    register Long defdat;
    register Long datum;
    register int col;

    fputs("\n# Bonds\n",OutFile);

    ResetFreqTable();
    ForEachBond
        IncFreqTable(GetBondDatum(bptr));

    WriteScriptAll();
    defdat = Freq[0].datum;
    WriteScriptDatum("wireframe",defdat);

    if( Freq[1].count )
    {   ForEachBond
        {   datum = GetBondDatum(bptr);
            if( datum != defdat )
            {    WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
                 WriteScriptDatum("wireframe",datum);
            }
        }
    } else if( !defdat )
        return;

    ResetFreqTable();
    ForEachBond
        IncFreqTable(bptr->col);

    col = (int)Freq[0].datum;
    if( col )
    {   WriteScriptAll();
        WriteScriptColour("bonds",col);
    }

    if( Freq[1].count )
        ForEachBond
            if( bptr->col != col )
            {   WriteScriptSelectBond(bptr->srcatom,bptr->dstatom);
                WriteScriptColour("bonds",bptr->col);
            }
}


static void WriteScriptLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Long first,last;
    register Label *label;

    fputs("\n# Labels\n",OutFile);
    WriteScriptAll();
    fputs("labels off\n",OutFile);
    if( !DrawLabels ) return;

    if( UseLabelCol )
    {   fprintf(OutFile,"colour labels [%d,%d,%d]\n",LabR,LabG,LabB);
    } else fputs("colour labels none\n",OutFile);
    fprintf(OutFile,"set fontsize %d\n",FontSize);

    label = (Label*)0;
    ForEachAtom
        if( aptr->label != label )
        {   if( label )
            {   WriteScriptBetween(first,last);
                fprintf(OutFile,"label \"%s\"\n",label->label);
            }
            label = (Label*)aptr->label;
            first = last = aptr->serno;
        } else last = aptr->serno;

    if( label )
    {   WriteScriptBetween(first,last);
        fprintf(OutFile,"label \"%s\"",label->label);
    }
}


static void WriteScriptMonitors()
{
    register Monitor *ptr;
    register int col;

    fputs("\n# Monitors\n",OutFile);
    if( !MonitList )
    {   fputs("monitors off\n",OutFile);
        return;
    }

    fprintf(OutFile,"set monitors %s\n",DrawMonitDistance?"on":"off");

    ResetFreqTable();
    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   fprintf(OutFile,"monitor %d %d\n",ptr->src->serno,ptr->dst->serno);
        IncFreqTable(ptr->col);
    }

    col = (int)Freq[0].datum;
    if( col )
    {   WriteScriptAll();
        WriteScriptColour("monitors",col);
    }

    if( Freq[1].count )
        for( ptr=MonitList; ptr; ptr=ptr->next )
            if( ptr->col != col )
            {   WriteScriptSelectBond(ptr->src,ptr->dst);
                WriteScriptColour("monitor",ptr->col);
            }
}


int WriteScriptFile( name )
    char *name;
{
    register int theta,phi,psi;
    register char *ptr;
    register int temp;

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalScriptError(name);
        return(False);
    }

    fprintf(OutFile,"#!rasmol -script\n# File: %s\n",name);
    fputs("# Creator: TINKER RasMol\n\n",OutFile);
    fputs("zap\n",OutFile);
    fprintf(OutFile,"background [%d,%d,%d]\n",BackR,BackG,BackB);

    if( !Database )
    {   /* No Molecule! */
        fclose(OutFile);
#ifdef APPLEMAC
        SetFileInfo(name,'RSML','RSML',133);
#endif
        return(True);
    }

    /* Molecule File Name */
    switch( DataFileFormat )
    {   default:
        case(FormatXYZ):      ptr = "xyz";      break;
    }
    fprintf(OutFile,"load %s \"%s\"\n",ptr,Info.filename);

    /* Colour Details */
    fprintf(OutFile,"set ambient %d\n", (int)(100*Ambient) );
    fputs("set specular ",OutFile);
    if( FakeSpecular )
    {   fprintf(OutFile,"on\nset specpower %d\n",SpecPower);
    } else fputs("off\n",OutFile);
    putc('\n',OutFile);

    /* Transformation */
    fputs("reset\n",OutFile);
    if( UseSlabPlane )
    {   temp = (int)(50.0*DialValue[7]);
        if( temp )
        {   fprintf(OutFile,"slab %d\n",temp+50);
        } else fputs("slab on\n",OutFile);

        fputs("set slabmode ",OutFile);
        switch( SlabMode )
        {   default:            
            case(SlabClose):    ptr = "solid";    break;
            case(SlabReject):   ptr = "reject";   break;
            case(SlabHalf):     ptr = "half";     break;
            case(SlabHollow):   ptr = "hollow";   break;
            case(SlabSection):  ptr = "section";
        }
        fputs(ptr,OutFile);
        putc('\n',OutFile);
    } else fputs("slab off\n",OutFile);

    phi = Round(Rad2Deg*asin(RotX[2]));
    if( phi == 90 )
    {   theta = -Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else if( phi == -90 )
    {   theta = Round(Rad2Deg*atan2(RotY[0],RotY[1]));
        psi = 0;
    } else /* General Case! */
    {   theta = Round(Rad2Deg*atan2(RotY[2],RotZ[2]));
        psi =  Round(-Rad2Deg*atan2(RotX[1],RotX[0]));
    }

    if( psi )   fprintf(OutFile,"rotate z %d\n",InvertY(-psi));
    if( phi )   fprintf(OutFile,"rotate y %d\n",phi);
    if( theta ) fprintf(OutFile,"rotate x %d\n",InvertY(-theta));

    temp = (int)(100.0*DialValue[4]);
    if( temp ) fprintf(OutFile,"translate x %d\n",temp);
    temp = (int)(100.0*DialValue[5]);
    if( temp ) fprintf(OutFile,"translate y %d\n",InvertY(-temp));

    if( DialValue[3] != 0.0 )
    {   if( DialValue[3]<0.0 )
        {   temp = (int)(100*DialValue[3]);
        } else temp = (int)(100*MaxZoom*DialValue[3]);
        fprintf(OutFile,"zoom %d\n",temp+100);
    }
    putc('\n',OutFile);

    /* Rendering */
    if( DrawAxes || DrawBoundBox || DrawUnitCell )
        fprintf(OutFile,"colour axes [%d,%d,%d]\n",BoxR,BoxG,BoxB);

    fprintf(OutFile,"set axes %s\n", DrawAxes? "on":"off" );
    fprintf(OutFile,"set boundingbox %s\n", DrawBoundBox? "on":"off" );
    fprintf(OutFile,"set unitcell %s\n", DrawUnitCell? "on":"off" );

    fputs("set bondmode and\ndots off\n\n",OutFile); 
    fputs("\n# Avoid Colour Problems!\nselect all\n",OutFile);
    fputs("colour bonds none\ncolour backbone none\n",OutFile);
    fputs("colour hbonds none\ncolour ssbonds none\n",OutFile);
    SelectAll = True;

    WriteScriptAtoms();
    if( UseSlabPlane && (SlabMode==SlabSection) )
    {   /* Section Mode Slabbing! */
        fclose(OutFile);
#ifdef APPLEMAC
        SetFileInfo(name,'RSML','RSML',133);
#endif
        return(True);
    }

    WriteScriptBonds();
    WriteScriptLabels();
    WriteScriptMonitors();
    fputc('\n',OutFile);
    
    WriteScriptAll();

    fclose(OutFile);
#ifdef APPLEMAC
    SetFileInfo(name,'RSML','RSML',133);
#endif
    return( True );
}
