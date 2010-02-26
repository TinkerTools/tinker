/* repres.c
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

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define REPRES
#include "molecule.h"
#include "graphics.h"
#include "repres.h"
#include "render.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "pixutils.h"

#define RootSix          2.44948974278

#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif

/* These define light source position */
#define LightDot(x,y,z)  ((x)+(y)+(z)+(z))
#define LightLength      RootSix
 
static Atom __far *Exclude;
static Monitor *FreeMonit;
static Label *FreeLabel;

#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)
 

static void FatalRepresError(ptr)
    char *ptr;
{
    char buffer[80];
 
    sprintf(buffer,"Renderer Error: Unable to allocate %s!",ptr);
    RasMolFatalExit(buffer);
}
 

/*============================*/
/*  Label Handling Functions  */
/*============================*/


static void ResetLabels()
{
    register Label *ptr;
 
    while( LabelList )
    {   ptr = LabelList;
        LabelList = ptr->next;
        ptr->next = FreeLabel;
        free(ptr->label);
        FreeLabel = ptr;
    }
}
 

void DeleteLabel( label )
    Label *label;
{
    register Label **ptr;
 
    if( label->refcount == 1 )
    {   ptr = &LabelList;
        while( *ptr != label )
        ptr = &(*ptr)->next;
 
        *ptr = label->next;
        label->next = FreeLabel;
        free(label->label);
        FreeLabel = label;
    } else label->refcount--;
}
 

int DeleteLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int result;
 
    if( !Database )
        return( True );
 
    result = True;
 
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   if( aptr->label )
            {   DeleteLabel( (Label*)aptr->label );
                aptr->label = (void*)0;
            }
            result = False;
        }
    DrawLabels = LabelList? True : False;
    return( result );
}


Label *CreateLabel( text, len )
    char *text;  int len;
{
    register Label *ptr;
 
    /* Test for existing label */
    for( ptr=LabelList; ptr; ptr=ptr->next )
        if( !strcmp(ptr->label,text) )
            return( ptr );
 
    if( FreeLabel )
    {   ptr = FreeLabel;  FreeLabel = ptr->next;
    } else if( !(ptr=(Label*)malloc(sizeof(Label))) )
        FatalRepresError("label");
 
    ptr->label = (char*)malloc(len+1);
    if( !ptr->label ) FatalRepresError("label");
    strcpy(ptr->label,text);
 
    ptr->next = LabelList;
    ptr->refcount = 0;
    LabelList = ptr;
    return( ptr );
}
 
 
void DefineLabels( label )
    char *label;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *ptr;
    register char *cptr;
    register int len;
 
 
    if( !Database ) return;
    if( DeleteLabels() )
        return;
 
    len = 0;
    for( cptr=label; *cptr; cptr++ )
        len++;
 
    /* Strip trailing spaces */
    while( len && cptr[-1]==' ' )
    {   cptr--;  len--;
        *cptr = '\0';
    }
 
    if( !len )
        return;
 
    ptr = CreateLabel(label,len);
    DrawLabels = True;
 
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   aptr->label = ptr;
            ptr->refcount++;
        }
}
 
 
void DefaultLabels( enable )
    int enable;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label1;
    register Label *label2;
 
    if( !Database )
        return;
 
    label1 = (Label*)0;
    label2 = (Label*)0;
 
/*  Start of Original Version
    ForEachAtom
        if( (aptr->flag&SelectFlag) && (aptr->elemno!=6)
                                        && (aptr->elemno!=1) )
        {   if( enable )
            {   if( !label1 )
                    label1 = CreateLabel("%e",2);
                aptr->label = label1;
                label1->refcount++;
            } else if( aptr->label )
            {   DeleteLabel( (Label*)aptr->label );
                aptr->label = (Label*)0;
            }
            ReDrawFlag |= RFRefresh;
        }
End of Original Version  */
 
/*  Start of New Version; Added by JWP, June 2001  */
    ForEachAtom
        if( aptr->flag&SelectFlag )
        {   if( enable )
            {   if( !label1 )
                    label1 = CreateLabel("%e%i",4);
                aptr->label = label1;
                label1->refcount++;
            } else if( aptr->label )
            {   DeleteLabel( (Label*)aptr->label );
                aptr->label = (Label*)0;
            }
            ReDrawFlag |= RFRefresh;
        }
 
    DrawLabels = LabelList? True : False;
}


void DisplayLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label;
    register int col,z;
 
    auto char buffer[256];
 
 
    if( !Database )
        return;
 
    if( !UseSlabPlane )
    {   z = ImageRadius + ZOffset;
    } else z = SlabValue - 1;
 
    ForEachAtom
        if( aptr->label )
        {   /* Peform Label Slabbing! */
            if( !ZValid(aptr->z) )
                continue;
 
            label = (Label*)aptr->label;
            FormatLabel(chain,group,aptr,label->label,buffer);
 
            if( !UseLabelCol )
            {   /* Depth-cue atom labels */
                /* col = aptr->col + (ColorDepth*                  */
                /*       (aptr->z+ImageRadius-ZOffset))/ImageSize; */
                col = aptr->col + (ColourMask>>1);
            } else col = LabelCol;
 
            /* (aptr->z+2) + ((aptr->flag & SphereFlag)?aptr->irad:0); */
            DisplayString(aptr->x+4,aptr->y,z,buffer,col);
        }
}


/*==============================*/
/*  Monitor Handling Functions  */
/*==============================*/


#ifdef FUNCPROTO
/* Function Prototype */
void AddMonitors( Atom __far*, Atom __far* );
#endif


void DeleteMonitors()
{
    register Monitor *ptr;
 
    while( MonitList )
    {   ptr = MonitList;
        if( ptr->col )
            Shade[Colour2Shade(ptr->col)].refcount--;
 
        MonitList = ptr->next;
        ptr->next = FreeMonit;
        FreeMonit = ptr;
    }
}
 

void AddMonitors( src, dst )
    Atom __far *src, __far *dst;
{
    register Monitor **prev;
    register Monitor *ptr;
    register Long dx,dy,dz;
    register Long dist;
 
    /* Delete an already existing monitor! */
    for( prev=&MonitList; (ptr=*prev); prev=&ptr->next )
         if( ((ptr->src==src) && (ptr->dst==dst)) ||
             ((ptr->src==dst) && (ptr->dst==src)) )
         {   if( ptr->col )
                 Shade[Colour2Shade(ptr->col)].refcount--;
 
             *prev = ptr->next;
             ptr->next = FreeMonit;
             FreeMonit = ptr;
             return;
         }
 
 
    /* Create a new monitor! */
    if( FreeMonit )
    {   ptr = FreeMonit;  FreeMonit = ptr->next;
    } else if( !(ptr=(Monitor*)malloc(sizeof(Monitor))) )
        FatalRepresError("monitor");
 
    dx = src->xorg - dst->xorg;
    dy = src->yorg - dst->yorg;
    dz = src->zorg - dst->zorg;
 
    /* ptr->dist = 100.0*CalcDistance(src,dst) */
    dist = isqrt( dx*dx + dy*dy + dz*dz );
    ptr->dist = (unsigned short)((dist<<1)/5);
 
    ptr->src = src;
    ptr->dst = dst;
    ptr->col = 0;
 
    ptr->next = MonitList;
    MonitList = ptr;
}


void CreateMonitor( src, dst )
    Long src, dst;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Atom __far *sptr;
    register Atom __far *dptr;
    register int done;
    char buffer[20];
 
    if( src == dst )
    {   if( !CommandActive )
            WriteChar('\n');
        WriteString("Error: Duplicate atom serial numbers!\n");
        CommandActive = False;
        return;
    }
 
    done = False;
    sptr = (Atom __far*)0;
    dptr = (Atom __far*)0;
 
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
 
    if( !done )
    {   if( !CommandActive )
            WriteChar('\n');
        WriteString("Error: Atom serial number");
        if( sptr )
        {   sprintf(buffer," %d",dst);
        } else if( dptr )
        {   sprintf(buffer," %d",src);
        } else sprintf(buffer,"s %d and %d",src,dst);
        WriteString(buffer); WriteString(" not found!\n");
        CommandActive = False;
 
    } else AddMonitors( sptr, dptr );
}
 
 
void DisplayMonitors()
{
    register Atom __far *s;
    register Atom __far *d;
    register Monitor *ptr;
    register int x,y,z;
    register int sc,dc;
    register int col;
 
    register char *cptr;
    register int dist;
    char buffer[10];
 
    if( !Database )
        return;
 
    if( !UseSlabPlane )
    {   z = ImageRadius + ZOffset;
    } else z = SlabValue-1;
    buffer[9] = '\0';
    buffer[6] = '.';
 
    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   s = ptr->src;
        d = ptr->dst;
 
        if( !ptr->col )
        {   sc = s->col;
            dc = d->col;
        } else sc = dc = ptr->col;
 
        ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
 
        if( DrawMonitDistance )
            if( ZValid( (s->z+d->z)/2 ) )
            {   x = (s->x+d->x)/2;
                y = (s->y+d->y)/2;
 
                if( !UseLabelCol )
                {   /* Use Source atom colour! */
                    col = sc + (ColourMask>>1);
                } else col = LabelCol;
 
                dist = ptr->dist;
                buffer[8] = (dist%10)+'0';  dist /= 10;
                buffer[7] = (dist%10)+'0';
                cptr = &buffer[5];
 
                if( dist > 9 )
                {   do {
                       dist /= 10;
                       *cptr-- = (dist%10)+'0';
                    } while( dist > 9 );
                    cptr++;
                } else *cptr = '0';
 
                DisplayString(x+4,y,z,cptr,col);
            }
    }
}
 

void ResetRepres()
{
    DeleteMonitors();

    DrawLabels = False;
    ResetLabels();

    DrawMonitDistance = True;
}


void InitialiseRepres()
{
    MonitList = (Monitor __far*)0;
    LabelList = (void*)0;

    FreeMonit = (Monitor __far*)0;
    FreeLabel = (void*)0;

    ResetRepres();
}
