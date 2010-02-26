/* render.c
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

#define RENDER
#include "molecule.h"
#include "graphics.h"
#include "render.h"
#include "repres.h"
#include "abstree.h"
#include "transfor.h"
#include "command.h"
#include "pixutils.h"

/* Avoid PowerPC Errors! */
#ifdef INFINITY
#undef INFINITY
#endif

#define PoolSize       16
#define ApproxZero     1.0E-3
#define RootSix        2.44948974278
#define INFINITY       200000
#define FUDGEFACTOR    1000

#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif

/* These define light source position */
#define LightDot(x,y,z)  ((x)+(y)+(z)+(z))
#define LightLength      RootSix
#define LightXComp       1
#define LightYComp       1
#define LightZComp       2

#if !defined(IBMPC) && !defined(APPLEMAC)
Card ColConstTable[MAXRAD];
#endif

typedef struct { Real h,l; } Interval;

static Atom __far * __far *YBucket;
static Atom __far * __far *IBuffer;
static int BuckY,ItemX;
static int FBufX,FBufY;
static int DBClear;

static Atom __far *SBuffer;
static Atom __far *Exclude;
static Real ShadowI, ShadowJ, ShadowK;
static int  ShadowX, ShadowY, ShadowZ;
static int deltax, deltay, deltaz;
static int xcord, ycord, zcord;
static int xflag, yflag, zflag;
static int xhash, yhash, zhash;
static int RayCount;

static Item __far *FreeItem;
static Real VoxRatio;
static int VoxelCount,InVoxCount;
static int ProbeCount;
static int VoxelsDone;

/* Identified Atom Info */
static AtomRef PickHist[4];
static Long IdentDist;
static int IdentFound;
static int IdentDepth;
static int PickCount;
static int PickMode;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(bptr=chain->blist;bptr;bptr=bptr->bnext)


static void FatalRenderError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Renderer Error: Unable to allocate %s!",ptr);
    RasMolFatalExit(buffer);
}


int isqrt( val )
    Card val;
{
#ifndef sun386
    register int i,result;
    register Card temp;
    register Card rem;

    i = 16;
    while( !(val&((Card)3<<30)) && i )
    {   val <<= 2;
        i--;
    }

    if( i )
    {   rem = (val>>30)-1;
        val <<= 2;
        result = 1;
        i--;

        while( i )
        {   rem = (rem<<2) | (val>>30);
            result <<= 1;
            val <<= 2;

            temp = result<<1;
            if( rem > temp )
            {   rem -= temp|1;
                result |= 1;
            }
            i--;
        }
        return( result );
    } else return( 0 );
#else
    return( (int)sqrt((double)val) );
#endif
}


/*=============================*/
/*  ClearBuffers Subroutines!  */
/*=============================*/
 

#ifdef IBMPC
/* Windows NT vs Microsoft C vs Borland Turbo C */
static void ClearMemory( register char __huge*, register long );

static void ClearMemory( ptr, count )
    register char __huge *ptr;
    register long count;
{
#ifndef _WIN32
#ifdef __TURBOC__
    register long left;
    register int off;

    off = OFFSETOF(ptr);
    if( off )
    {   left = 65536 - off;
        if( count < left )
        {   _fmemset(ptr,0,(size_t)count);
            return;
        } else
        {   _fmemset(ptr,0,(size_t)left);
            count -= left;
            ptr += left;
        }
    }

    while( count > 65535 )
    {   _fmemset(ptr,0,(size_t)65535);
        count -= 65536;
        ptr += 65535;
        *ptr++ = '\0';
    }

    if( count )
        _fmemset(ptr,0,(size_t)count);
#else /* Microsoft C/C++ */
    while( count > 65534 )
    {   _fmemset(ptr,0,(size_t)65534);
        count -= 65534;
	ptr += 65534;
    }
    if( count )
        _fmemset(ptr,0,(size_t)count);
#endif
#else  /* Windows NT  */
    memset(ptr,0,(size_t)count);
#endif /* Windows NT */
}


void ClearBuffers()
{
    register char __huge *ptr;

    if( !FBClear )
    {   FBClear = True;
	ptr = (Pixel __huge*)GlobalLock(FBufHandle);
        ClearMemory(ptr,(Long)XRange*YRange);
	GlobalUnlock(FBufHandle);
    }

    if( !DBClear )
    {   DBClear = True;
	ptr = (char __huge*)GlobalLock(DBufHandle);
        ClearMemory(ptr,(Long)XRange*YRange*sizeof(short));
	GlobalUnlock(DBufHandle);
    }
}
#else
#if defined(__sgi)        /* memset */
void ClearBuffers()
{
#ifndef EIGHTBIT
    register Long *ptr;
    register Long *end;
    register Long fill;

    if( !FBClear )
    {   FBClear = True;
	fill = Lut[BackCol];
	ptr = (Long*)FBuffer;
	end = (Long*)(FBuffer+(Long)XRange*YRange);
	do { *ptr++=fill; *ptr++=fill;
	     *ptr++=fill; *ptr++=fill;
	} while( ptr<end );
    }
#else

    if( !FBClear )
    {   FBClear = True;
        memset(FBuffer,Lut[BackCol],(Long)XRange*YRange);
    }
#endif

    if( !DBClear )
    {   DBClear = True;
        memset(DBuffer,0,(Long)XRange*YRange*sizeof(short));
    }
}
#else /* !memset */
void ClearBuffers()
{
    register Long *ptr;
    register Long *end;
    register Long fill;

    if( !FBClear )
    {   FBClear = True;
	fill = Lut[BackCol];
#ifdef EIGHTBIT
	fill |= fill<<8;
	fill |= fill<<16;
#endif
	ptr = (Long*)FBuffer;
	end = (Long*)(FBuffer+(Long)XRange*YRange);
	do { *ptr++=fill; *ptr++=fill;
	     *ptr++=fill; *ptr++=fill;
	} while( ptr<end );
    }

    if( !DBClear )
    {   DBClear = True;
	ptr = (Long*)DBuffer;
	end = (Long*)(DBuffer+(Long)XRange*YRange);
	do { *ptr++=0; *ptr++=0;
	     *ptr++=0; *ptr++=0;
	} while( ptr<end );
    }
}
#endif /* !memset    */
#endif /* UNIX */


void ReAllocBuffers()
{
    register Atom __far * __far *iptr;
    register int index,len;
    register Long temp;

    temp = (Long)XRange*YRange*sizeof(short)+32;
#ifdef IBMPC
    if( DBufHandle ) GlobalFree(DBufHandle);
    DBufHandle = GlobalAlloc(GMEM_MOVEABLE,temp);
    if( !DBufHandle ) FatalRenderError("depth buffer");
#else
    if( DBuffer ) _ffree( DBuffer );
    DBuffer = (short*)_fmalloc( temp );
    if( !DBuffer ) FatalRenderError("depth buffer");
#endif
    DBClear=False;

    if( YBucket && (BuckY<YRange) )
    {   _ffree(YBucket); 
	YBucket=(void __far*)0; 
    }

    if( !YBucket )
    {   len = YRange*sizeof(Atom __far*);
	YBucket = (Atom __far* __far*)_fmalloc( len );
	if( !YBucket ) FatalRenderError("Y buckets");
	BuckY = YRange;
    }

    if( IBuffer && (ItemX<XRange) )
    {   _ffree(IBuffer); 
	IBuffer=(void __far*)0; 
    }

    if( !IBuffer )
    {   len = (XRange+4)*sizeof(Atom __far*);
	IBuffer = (Atom __far* __far*)_fmalloc(len);
	if( !IBuffer ) FatalRenderError("item buffer");
	len = XRange>>2;  iptr = IBuffer;
	for( index=0; index<=len; index++ )
	{   *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
	    *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
	}
	ItemX = XRange;
    }
}


void ReSizeScreen()
{
    register Real orig;

    if( Range != ZoomRange )
    {   orig = MaxZoom;
        /* Code should match InitialTransform() */
        /* MaxZoom*DScale*Range*750 == 252      */
	MaxZoom = 0.336*(WorldSize+1500)/Range;
	ZoomRange = Range;  MaxZoom -= 1.0;

	/* Handle Change in MaxZoom */
	if( DialValue[3]>0.0 )
	{   DialValue[3] *= orig/MaxZoom;
	    if( DialValue[3]>1.0 )
		DialValue[3] = 1.0;
	}
    }

#ifdef IBMPC
    if( !FBufHandle || (FBufX!=XRange) || (FBufY!=YRange) )
#else /* UNIX */
    if( !FBuffer || (FBufX!=XRange) || (FBufY!=YRange) )
#endif
    {   if( !CreateImage() )
	    FatalRenderError("frame buffer");

	BucketFlag = False;
	FBufX=XRange;  FBufY=YRange;  FBClear = False;
	ReAllocBuffers();
	ClearBuffers();
    }
}


static void PrepareYBucket()
{
    register Atom __far * __far *temp;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int scan;
    register int rad;

    temp = YBucket;
    for( scan=0; scan<BuckY; scan++ )
	*temp++ = (void __far*)0;

    if( UseClipping )
    {   ForEachAtom
	    if( aptr->flag&SphereFlag )
	    {   rad = aptr->irad;
		if( (aptr->x-rad>=XRange) || 
		    (aptr->x+rad<0) || (aptr->y+rad<0) )
		    continue;
		if( (scan=aptr->y-rad) > BuckY ) continue;

		if( scan>0 )
		{   aptr->bucket = YBucket[scan];
		    YBucket[scan] = aptr;
		} else
		{   aptr->bucket = *YBucket;
		    *YBucket = aptr;
		}
	    }
    } else
	ForEachAtom
	    if( aptr->flag&SphereFlag )
	    {   scan = aptr->y-aptr->irad;
		aptr->bucket = YBucket[scan];
		YBucket[scan] = aptr;
	    }
    BucketFlag = True;
}


#ifdef FUNCPROTO
/* Function Prototypes */
static void SqrInterval( Interval __far* );
static void VoxelInsert( Atom __far*, int );
static int AtomInter( Atom __far* );
#endif


static void SqrInterval( ival )
    Interval __far *ival;
{   register Real l,h;

    l = ival->l;
    h = ival->h;

    if( l>=0.0 )
    {   ival->l = l*l;
	ival->h = h*h;
    } else if( h<0.0 )
    {   ival->l = h*h;
	ival->h = l*l;
    } else
    {   ival->h = (-l>h)? l*l : h*h;
	ival->l = 0.0;
    }
}


static void VoxelInsert( ptr, ref )
    Atom __far *ptr;
    int ref;
{
    register Item __far *datum;
    register int i;

    if( !FreeItem )
    {   datum = (Item __far*)_fmalloc( PoolSize*sizeof(Item) );
	if( !datum ) FatalRenderError("voxel item");
	for( i=1; i<PoolSize; i++ )
	{   datum->list = FreeItem;
	    FreeItem = datum++;
	}
    } else
    {   datum = FreeItem;
	FreeItem = datum->list;
    }
    datum->data = ptr;
    InVoxCount++;

    if( !HashTable[ref] ) VoxelCount++;
    datum->list = (Item __far*)HashTable[ref];
    HashTable[ref] = (void __far*)datum;
}


void ResetVoxelData()
{
    register Item __far *datum;
    register int i;

    if( VoxelsDone )
    {   for( i=0; i<VOXSIZE; i++ )
	    if( HashTable[i] )
	    {   datum = (Item __far*)HashTable[i];
                while( datum->list ) datum = datum->list;
		datum->list = FreeItem;
		FreeItem = (Item __far*)HashTable[i];
		HashTable[i] = (void __far*)0;
	    }
	VoxelsDone = False;
    } else for( i=0; i<VOXSIZE; i++ )
	HashTable[i] = (void __far*)0;
    VoxelsClean = True;
}


void CreateVoxelData( flag )
    int flag;
{
    static Interval ix, iy, iz;
    register int lvx, lvy, lvz;
    register int hvx, hvy, hvz;
    register Long mx, my, mz;
    register int px, py, pz;
    register int i, rad;
    register Real radius2;

    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;


    ResetVoxelData();
    ProbeCount = InVoxCount = VoxelCount = 0;
    VoxRatio = (Real)SideLen/VOXORDER;
    IVoxRatio = 1.0/VoxRatio;
    VoxelsDone = True;

    ForEachAtom
    if( aptr->flag & flag )
    {   mx = aptr->xorg+Offset;
	my = aptr->yorg+Offset;
	mz = aptr->zorg+Offset;
	if( flag != SphereFlag )
	{   if( !ProbeRadius )
	    {   rad = ElemVDWRadius(aptr->elemno);
	    } else rad = ProbeRadius;
	} else rad = aptr->radius;  
	radius2 = (Long)rad*rad;

	lvx = (int)((mx-rad)*IVoxRatio);  hvx = (int)((mx+rad)*IVoxRatio);
	lvy = (int)((my-rad)*IVoxRatio);  hvy = (int)((my+rad)*IVoxRatio);
	lvz = (int)((mz-rad)*IVoxRatio);  hvz = (int)((mz+rad)*IVoxRatio);


	for( px=lvx; px<=hvx; px++ )
	{   ix.l=px*VoxRatio-mx;
	    ix.h=ix.l+VoxRatio;  
	    SqrInterval(&ix);
	    i = VOXORDER2*px + VOXORDER*lvy;
       
	    for( py=lvy; py<=hvy; py++ )
	    {   iy.l=py*VoxRatio-my;
		iy.h=iy.l+VoxRatio;
		SqrInterval(&iy);
		
		for( pz=lvz; pz<=hvz; pz++ )
		{   iz.l=pz*VoxRatio-mz; 
		    iz.h=iz.l+VoxRatio;
		    SqrInterval(&iz);

#ifdef ORIG
                    /* Test for surface voxels */
		    if( ((ix.l+iy.l+iz.l)<radius2) && 
			((ix.h+iy.h+iz.h)>radius2) )
			VoxelInsert( aptr, i+pz );
#else
                    /* Test for contained voxels */
                    if( ix.l+iy.l+iz.l < radius2 )
			VoxelInsert( aptr, i+pz );
#endif
		} /*pz*/
		i += VOXORDER;
	    } /*py*/
	} /*px*/
    }
}


void ShadowTransform()
{
    ShadowI = (LightXComp*LightDot(RotX[0],-RotY[0],RotZ[0]))/LightLength;
    ShadowJ = (LightYComp*LightDot(RotX[0],-RotY[0],RotZ[0]))/LightLength;
    ShadowK = (LightZComp*LightDot(RotX[0],-RotY[0],RotZ[0]))/LightLength;

    if( ShadowI>ApproxZero )
    {   deltax =  (int)(FUDGEFACTOR/ShadowI); xhash =  VOXORDER2; xflag =  1;
    } else if( ShadowI<-ApproxZero )
    {   deltax = -(int)(FUDGEFACTOR/ShadowI); xhash = -VOXORDER2; xflag = -1;
    } else xflag = 0;

    if( ShadowJ>ApproxZero )
    {   deltay =  (int)(FUDGEFACTOR/ShadowJ); yhash =  VOXORDER; yflag =  1;
    } else if( ShadowJ<-ApproxZero )
    {   deltay = -(int)(FUDGEFACTOR/ShadowJ); yhash = -VOXORDER; yflag = -1;
    } else yflag = 0;

    if( ShadowK>ApproxZero )
    {   deltaz =  (int)(FUDGEFACTOR/ShadowK); zhash = zflag =  1;
    } else if( ShadowK<-ApproxZero )
    {   deltaz = -(int)(FUDGEFACTOR/ShadowK); zhash = zflag = -1;
    } else zflag = 0;
}


static int AtomInter( ptr )
    Atom __far *ptr;
{
    register Long modv,radius2;
    register int vx, vy, vz;
    register Real tca;

    if( ptr->mbox == RayCount )
	return( False );
    ptr->mbox = RayCount;

    vx = (int)ptr->xorg-ShadowX;
    vy = (int)ptr->yorg-ShadowY;
    vz = (int)ptr->zorg-ShadowZ;

    tca = vx*ShadowI + vy*ShadowJ + vz*ShadowK;
    if( tca<0.0 ) return( False );
    
    radius2 = ptr->radius+10;  radius2 = radius2*radius2;
    modv = (Long)vx*vx + (Long)vy*vy + (Long)vz*vz - radius2;
    return( modv<tca*tca );
}


static int ShadowRay()
{
    register Item __far * __far *ident;
    register Item __far *ptr;
    register Real ex, ey, ez;
    register Long dx, dy, dz;
    register int ref;

   
    RayCount++;
    if( SBuffer )
    {   if( (SBuffer!=Exclude) && AtomInter(SBuffer) )
	    return( True );
	SBuffer = (void __far*)0;
    }

    ex = IVoxRatio*(ShadowX+Offset);  xcord = (int)ex;
    ey = IVoxRatio*(ShadowY+Offset);  ycord = (int)ey;
    ez = IVoxRatio*(ShadowZ+Offset);  zcord = (int)ez;

    ref = VOXORDER2*xcord+VOXORDER*ycord+zcord;
    ident = (Item __far* __far*)(HashTable+ref);

    if( xflag==1 ) 
    {   dx = (Long)(((xcord+1)-ex)*deltax);
    } else if( xflag == -1 )
    {   dx = (Long)((ex-xcord)*deltax); 
    } else dx = INFINITY;

    if( yflag==1 ) 
    {   dy = (Long)(((ycord+1)-ey)*deltay);
    } else if( yflag == -1 )
    {   dy = (Long)((ey-ycord)*deltay); 
    } else dy = INFINITY;

    if( zflag==1 ) 
    {   dz = (Long)(((zcord+1)-ez)*deltaz);
    } else if( zflag == -1 )
    {   dz = (Long)((ez-zcord)*deltaz); 
    } else dz = INFINITY;

    
    while( True )
    {   for( ptr = *ident; ptr; ptr = ptr->list )
	    if( (ptr->data!=Exclude) && AtomInter(ptr->data) )
	    {   SBuffer = ptr->data;
		return( True );
	    }

	if( (dx<=dy) && (dx<=dz) )
	{   xcord += xflag;
	    if( (xcord<0) || (xcord>=VOXORDER) ) return( False );
	    ident += xhash; dx += deltax;
	} else if( dy<=dz  ) /*(dy<=dx)*/
	{   ycord += yflag;
	    if( (ycord<0) || (ycord>=VOXORDER) ) return( False );
	    ident += yhash; dy += deltay;
	} else /* (dz<=dx) && (dz<=dy) */
	{   zcord += zflag;
	    if( (zcord<0) || (zcord>=VOXORDER) ) return( False );
	    ident += zhash; dz += deltaz;
	}
    }
}


#define UpdateScanAcross \
	if( depth>*dptr )   \
	{   *dptr = depth;  \
	    iptr[dx] = ptr; \
	} dptr++; dx++;


/* ScanLine for Shadows! */
static void ScanLine()
{
    static Atom __far *list;
    register Atom __far *ptr;
    register Atom __far * __far *iptr;
    register Atom __far * __far *prev;
    register short __huge *dbase;
    register short __huge *dptr;
    register Pixel __huge *fptr;
    register Byte __far *tptr;

    register int pos,depth,inten;
    register int lastx,wide,scan;
    register int dx,dy,dz;

    fptr = FBuffer;
    dbase = DBuffer;
    list = (void __far*)0;  

    wide = XRange>>2;  iptr = IBuffer;
    for( pos=0; pos<=wide; pos++ )
    {   *iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
	*iptr++ = (void __far*)0;  *iptr++ = (void __far*)0;
    }


    for( scan=0; scan<YRange; scan++ )
    {   for( ptr = YBucket[scan]; ptr; ptr = ptr->bucket )
	{    ptr->next = list; list = ptr; }

	prev = &list;
	for( ptr=list; ptr; ptr=ptr->next )
	{   dy = scan - ptr->y;
	    wide = LookUp[ptr->irad][AbsFun(dy)];
	    lastx = (XRange-1)-ptr->x;
	    if( wide<lastx ) lastx=wide;
	    dx = - MinFun(wide,ptr->x);

	    iptr = IBuffer+ptr->x;
	    tptr = LookUp[wide];

	    dptr = dbase+ptr->x+dx;
	    while( dx<=lastx )
	    {   depth = tptr[AbsFun(dx)]+ptr->z;
		UpdateScanAcross;
	    }

	    /* Remove completed atoms */
	    if( dy == ptr->irad )
	    {   *prev = ptr->next;
	    } else prev = &ptr->next;
	} /*ptr*/


	/* Process visible scanline */
	prev = (Atom __far* __far*)IBuffer;
	SBuffer = (void __far*)0;
	dptr = dbase; 

	for( pos=0; pos<XRange; pos++ )
	{   if( *prev )
	    {   ptr = *prev;
                dz = *dptr-ptr->z;
                inten = LightDot(pos-ptr->x,InvertY(scan-ptr->y),dz);
		if( inten>0 )
		{   inten = (int)( (inten*ColConst[ptr->irad])>>ColBits);
		    dz = *dptr-ZOffset;
		    dx = pos-XOffset;
		    dy =   scan-YOffset;

		    ShadowX = (int)(dx*InvX[0]+dy*InvX[1]+dz*InvX[2]);
		    ShadowY = (int)(dx*InvY[0]+dy*InvY[1]+dz*InvY[2]);
		    ShadowZ = (int)(dx*InvZ[0]+dy*InvZ[1]+dz*InvZ[2]);

		    Exclude = ptr;
		    if( ShadowRay() )
		    {   *fptr = Lut[ptr->col+(inten>>2)];
		    } else *fptr = Lut[ptr->col+inten];
		} else *fptr = Lut[ptr->col];
		*prev = (void __far*)0;
	    }
	    dptr++; fptr++; prev++;
	}
	dbase = dptr;
    } /*scan*/
}


static void DisplaySpaceFill()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;

    if( UseShadow )
    {   if( !BucketFlag )
	    PrepareYBucket();
	ScanLine();
    } else if( UseClipping )
    {   ForEachAtom
	    if( aptr->flag&SphereFlag )
		ClipSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
    } else 
	ForEachAtom
	    if( aptr->flag&SphereFlag )
		DrawSphere(aptr->x,aptr->y,aptr->z,aptr->irad,aptr->col);
}


static void DisplayWireframe()
{
    register Bond __far *bptr;
    register Atom __far *s;
    register Atom __far *d;
    register int sc,dc;

    if( UseClipping )
    {   ForEachBond
           if( bptr->flag & DrawBondFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
	       if( !bptr->col ) 
	       {   sc = s->col;  dc = d->col;
	       } else sc = dc = bptr->col;

	       if( bptr->flag&WireFlag )
	       {   ClipTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
	       } else if( bptr->flag&CylinderFlag )
	       {   if( bptr->irad>0 )
	           {  ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
                                   sc,dc,bptr->irad);
	           } else ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,
				   sc+ColourMask,dc+ColourMask);
               } else /* bptr->flag & DashFlag */
                   ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
	   }
    } else
	ForEachBond
           if( bptr->flag & DrawBondFlag )
           {   s = bptr->srcatom; d = bptr->dstatom;
               if( !bptr->col )
               {   sc = s->col;  dc = d->col;
               } else sc = dc = bptr->col;

               if( bptr->flag&WireFlag )
               {      DrawTwinVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
	       } else if( bptr->flag&CylinderFlag )
	       {   if( bptr->irad>0 )
	           {  DrawCylinder(s->x,s->y,s->z,d->x,d->y,d->z,
                                   sc,dc,bptr->irad);
	           } else DrawTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,
				   sc+ColourMask,dc+ColourMask);
	       } else ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
           }
}


static void DisplayBoxes()
{
    register Real lena, lenb, lenc;
    register Real tmpx, tmpy, tmpz;
    register Real cosa, cosb, cosg;
    register Real temp, sing;

    register int xorg, yorg, zorg;
    register int dxx,dxy,dxz;
    register int dyx,dyy,dyz;
    register int dzx,dzy,dzz;
    register int x, y, z;


    if( DrawAxes  || DrawBoundBox )
    {   dxx = (int)(MaxX*MatX[0]);
	dxy = (int)(MaxX*MatY[0]);
	dxz = (int)(MaxX*MatZ[0]);

	dyx = (int)(MaxY*MatX[1]);
	dyy = (int)(MaxY*MatY[1]);
	dyz = (int)(MaxY*MatZ[1]);

	dzx = (int)(MaxZ*MatX[2]);
	dzy = (int)(MaxZ*MatY[2]);
	dzz = (int)(MaxZ*MatZ[2]);

	if( DrawAxes )
	{   /* Line (MinX,0,0) to (MaxX,0,0) */
            x = XOffset+dxx;  y = YOffset+dxy;  z = ZOffset+dxz;
            if( ZValid(z) ) DisplayString(x+2,y,z,"X",BoxCol);
	    ClipTwinLine(XOffset-dxx,YOffset-dxy,ZOffset-dxz,
                         x,y,z,BoxCol,BoxCol);

	    /* Line (0,MinY,0) to (0,MaxY,0) */
            x = XOffset+dyx;  y = YOffset+dyy;  z = ZOffset+dyz;
            if( ZValid(z) ) DisplayString(x+2,y,z,"Y",BoxCol);
	    ClipTwinLine(XOffset-dyx,YOffset-dyy,ZOffset-dyz, 
			 x,y,z,BoxCol,BoxCol);


	    /* Line (0,0,MinZ) to (0,0,MaxZ) */
            x = XOffset-dzx;  y = YOffset-dzy;  z = ZOffset-dzz;
            if( ZValid(z) ) DisplayString(x+2,y,z,"Z",BoxCol);
	    ClipTwinLine(XOffset+dzx,YOffset+dzy,ZOffset+dzz, 
			 x,y,z,BoxCol,BoxCol);

	}

	if( DrawBoundBox )
	{   /* Line (MinX,MinY,MinZ) to (MaxX,MinY,MinZ) */
	    x=XOffset-dyx-dzx;  y=YOffset-dyy-dzy;  z=ZOffset-dyz-dzz;
	    ClipTwinLine(x-dxx,y-dxy,z-dxz,x+dxx,y+dxy,z+dxz,BoxCol,BoxCol);

	    /* Line (MaxX,MinY,MinZ) to (MaxX,MaxY,MinZ) */
	    x=XOffset+dxx-dzx;  y=YOffset+dxy-dzy;  z=ZOffset+dxz-dzz;
	    ClipTwinLine(x-dyx,y-dyy,z-dyz,x+dyx,y+dyy,z+dyz,BoxCol,BoxCol);

	    /* Line (MaxX,MaxY,MinZ) to (MinX,MaxY,MinZ) */
	    x=XOffset+dyx-dzx;  y=YOffset+dyy-dzy;  z=ZOffset+dyz-dzz;
	    ClipTwinLine(x+dxx,y+dxy,z+dxz,x-dxx,y-dxy,z-dxz,BoxCol,BoxCol);

	    /* Line (MinX,MaxY,MinZ) to (MinX,MinY,MinZ) */
	    x=XOffset-dxx-dzx;  y=YOffset-dxy-dzy;  z=ZOffset-dxz-dzz;
	    ClipTwinLine(x+dyx,y+dyy,z+dyz,x-dyx,y-dyy,z-dyz,BoxCol,BoxCol);


	    /* Line (MinX,MinY,MinZ) to (MinX,MinY,MaxZ) */
	    x=XOffset-dxx-dyx;  y=YOffset-dxy-dyy;  z=ZOffset-dxz-dyz;
	    ClipTwinLine(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);

	    /* Line (MaxX,MinY,MinZ) to (MaxX,MinY,MaxZ) */
	    x=XOffset+dxx-dyx;  y=YOffset+dxy-dyy;  z=ZOffset+dxz-dyz;
	    ClipTwinLine(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);

	    /* Line (MaxX,MaxY,MinZ) to (MaxX,MaxY,MaxZ) */
	    x=XOffset+dxx+dyx;  y=YOffset+dxy+dyy;  z=ZOffset+dxz+dyz;
	    ClipTwinLine(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);

	    /* Line (MinX,MaxY,MinZ) to (MinX,MaxY,MaxZ) */
	    x=XOffset-dxx+dyx;  y=YOffset-dxy+dyy;  z=ZOffset-dxz+dyz;
	    ClipTwinLine(x-dzx,y-dzy,z-dzz,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);


	    /* Line (MinX,MinY,MaxZ) to (MaxX,MinY,MaxZ) */
	    x=XOffset-dyx+dzx;  y=YOffset-dyy+dzy;  z=ZOffset-dyz+dzz;
	    ClipTwinLine(x-dxx,y-dxy,z-dxz,x+dxx,y+dxy,z+dxz,BoxCol,BoxCol);

	    /* Line (MaxX,MinY,MaxZ) to (MaxX,MaxY,MaxZ) */
	    x=XOffset+dxx+dzx;  y=YOffset+dxy+dzy;  z=ZOffset+dxz+dzz;
	    ClipTwinLine(x-dyx,y-dyy,z-dyz,x+dyx,y+dyy,z+dyz,BoxCol,BoxCol);

	    /* Line (MaxX,MaxY,MaxZ) to (MinX,MaxY,MaxZ) */
	    x=XOffset+dyx+dzx;  y=YOffset+dyy+dzy;  z=ZOffset+dyz+dzz;
	    ClipTwinLine(x+dxx,y+dxy,z+dxz,x-dxx,y-dxy,z-dxz,BoxCol,BoxCol);

	    /* Line (MinX,MaxY,MaxZ) to (MinX,MinY,MaxZ) */
	    x=XOffset-dxx+dzx;  y=YOffset-dxy+dzy;  z=ZOffset-dxz+dzz;
	    ClipTwinLine(x+dyx,y+dyy,z+dyz,x-dyx,y-dyy,z-dyz,BoxCol,BoxCol);
	}
    }

    if( DrawUnitCell && *Info.spacegroup )
    {   /* Calculate Unit Cell! */
	lena = 250.0*Info.cella;
	lenb = 250.0*Info.cellb;
	lenc = 250.0*Info.cellc;

	cosa = cos(Info.cellalpha);
	cosb = cos(Info.cellbeta);
	cosg = cos(Info.cellgamma);  
        sing = sin(Info.cellgamma);

	temp = cosa*cosa + cosb*cosb + cosg*cosg - 2.0*cosa*cosb*cosg;
	tmpx = cosb; 
	tmpy = (cosa - cosb*cosg)/sing;
	tmpz = -sqrt(1.0-temp)/sing;

	dxx = (int)(lena*MatX[0]);
	dxy = (int)(lena*MatY[0]);
	dxz = (int)(lena*MatZ[0]);

	dyx = (int)(lenb*(cosg*MatX[0] + sing*MatX[1]));
	dyy = (int)(lenb*(cosg*MatY[0] + sing*MatY[1]));
	dyz = (int)(lenb*(cosg*MatZ[0] + sing*MatZ[1]));

	dzx = (int)(lenc*(tmpx*MatX[0] + tmpy*MatX[1] + tmpz*MatX[2]));
	dzy = (int)(lenc*(tmpx*MatY[0] + tmpy*MatY[1] + tmpz*MatY[2]));
	dzz = (int)(lenc*(tmpx*MatZ[0] + tmpy*MatZ[1] + tmpz*MatZ[2]));

	xorg = XOffset - (int)(OrigCX*MatX[0]+OrigCY*MatX[1]+OrigCZ*MatX[2]);
	yorg = YOffset - (int)(OrigCX*MatY[0]+OrigCY*MatY[1]+OrigCZ*MatY[2]);
	zorg = ZOffset + (int)(OrigCX*MatZ[0]+OrigCY*MatZ[1]+OrigCZ*MatZ[2]);


	/* Draw Unit Cell! */
	x = xorg;
	y = yorg;
	z = zorg;
	ClipTwinLine(x,y,z,x+dxx,y+dxy,z+dxz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x+dyx,y+dyy,z+dyz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);

	x = xorg + dxx + dyx;
	y = yorg + dxy + dyy;
	z = zorg + dxz + dyz;
	ClipTwinLine(x,y,z,x-dxx,y-dxy,z-dxz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x-dyx,y-dyy,z-dyz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x+dzx,y+dzy,z+dzz,BoxCol,BoxCol);

	x = xorg + dxx + dzx;
	y = yorg + dxy + dzy;
	z = zorg + dxz + dyz;
	ClipTwinLine(x,y,z,x-dxx,y-dxy,z-dxz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x+dyx,y+dyy,z+dyz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x-dzx,y-dzy,z-dzz,BoxCol,BoxCol);

	x = xorg + dyx + dzx;
	y = yorg + dyy + dzy;
	z = zorg + dyz + dzz;
	ClipTwinLine(x,y,z,x+dxx,y+dxy,z+dxz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x-dyx,y-dyy,z-dyz,BoxCol,BoxCol);
	ClipTwinLine(x,y,z,x-dzx,y-dzy,z-dzz,BoxCol,BoxCol);
    }
}


static void DisplaySelected()
{
    register Atom __far *s, __far *d;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register int irad,sc,dc;
    register int col;

    irad = (int)(Scale*20);

    if( irad>0 )
    {   ForEachBond
	{   s = bptr->srcatom;  
	    col = (s->flag&SelectFlag)? 1 : 0;
	    sc = Shade2Colour(col);

	    d = bptr->dstatom;  
	    col = (d->flag&SelectFlag)? 1 : 0;
	    dc = Shade2Colour(col);
	    ClipCylinder(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc,irad);
	}
    } else ForEachBond
	{   s = bptr->srcatom;  
	    col = (s->flag&SelectFlag)? 1 : 0;
	    sc = Shade2Colour(col);

	    d = bptr->dstatom;  
	    col = (d->flag&SelectFlag)? 1 : 0;
	    dc = Shade2Colour(col);
	    ClipTwinLine(s->x,s->y,s->z,d->x,d->y,d->z,
			 sc+ColourMask,dc+ColourMask);
	}


    irad = (int)(Scale*50);
    ForEachAtom
	if( aptr->flag&NonBondFlag )
	{   col = Shade2Colour( (aptr->flag&SelectFlag)? 1 : 0 );
	    ClipSphere(aptr->x,aptr->y,aptr->z,irad,col);
	}
}


static void RenderFrame()
{
    if( !DisplayMode )
    {   if( DrawAtoms ) 
	    DisplaySpaceFill();

	if( !UseSlabPlane || (SlabMode != SlabSection) )
	{   if( DrawBonds ) DisplayWireframe();
	    if( DrawLabels ) DisplayLabels();
            if( MonitList ) DisplayMonitors();
	}
    } else DisplaySelected();
    DisplayBoxes();
}


void DrawFrame()
{
    register double temp;
    register int wide;

    if( !Database ) 
	return;

    ClearBuffers();
    if( !DisplayMode )
    {   if( UseShadow && DrawAtoms )
	    if( !VoxelsClean )
		CreateVoxelData( SphereFlag );
    }

    if( UseSlabPlane )
    {   SlabValue = (int)(DialValue[7]*ImageRadius)+ZOffset;
	SlabInten = (int)(ColourMask*LightZComp/LightLength);
	SliceValue = SlabValue+16;
	UseClipping = True;
    } else UseClipping = UseScreenClip;

#ifdef IBMPC
    /* Lock Buffers into Memory */
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
    DBuffer = (short __huge*)GlobalLock(DBufHandle);
#endif

    /* Common View Elements */
    View.yskip = XRange;
    View.ymax = YRange;

    if( UseStereo )
    {   temp = StereoAngle/180.0;
        wide = XRange>>1;

        /* Create 'Left' View structure */
        View.fbuf = FBuffer;
        View.dbuf = DBuffer;
        View.xmax = wide;

        DialValue[1] -= temp;
        ReDrawFlag |= RFRotateY;
        ApplyTransform();
        RenderFrame();

        /* Create 'Right' View structure */
        View.fbuf = FBuffer+wide;
        View.dbuf = DBuffer+wide;
        View.xmax = wide;

        DialValue[1] += temp;
        ReDrawFlag |= RFRotateY;
        ApplyTransform();
        RenderFrame();

    } else /* Mono */
    {   /* Create 'Mono' View structure */
        View.fbuf = FBuffer;
        View.dbuf = DBuffer;
        View.xmax = XRange;
        RenderFrame();
    }

#ifdef IBMPC
    /* Unlock Buffers */
    GlobalUnlock(FBufHandle);
    GlobalUnlock(DBufHandle);
#endif
    DBClear = False;
    FBClear = False;
}



#ifdef FUNCPROTO
/* Function Prototype */
static void TestAtomProximity( Atom __far *, int, int );
#endif


static void TestAtomProximity( ptr, xpos, ypos )
    Atom __far *ptr;
    int xpos, ypos;
{
    register Long dist;
    register int dx,dy;

    if( UseSlabPlane && (ptr->z>SlabValue) )
	return;

    dx = ptr->x - xpos;
    dy = ptr->y - ypos;

    dist = (Long)dx*dx + (Long)dy*dy;

    if( IdentFound )
    {   if( dist==IdentDist )
	{   if( ptr->z<IdentDepth )
		return;
	} else if( dist>IdentDist ) 
	    return;
    }

    IdentDepth = ptr->z;
    IdentFound = True;
    IdentDist = dist;
    QAtom = ptr;
}

static void IdentifyAtom( xpos, ypos )
    int xpos, ypos;
{
    register int rad, wide, dpth;
    register int new, dx, dy, dz;
    register Chain __far *chain;
    register Group __far *group;
    register Atom  __far *aptr;
    register Bond __far *bptr;

    /* Reset Search */
    QChain = (void __far*)0;
    QGroup = (void __far*)0;
    QAtom = (void __far*)0;
    IdentFound = False;

    if( !DisplayMode )
    {   if( !UseSlabPlane || (SlabMode != SlabSection) )
	{   if( DrawBonds )
		ForEachBond
		    if( bptr->flag&DrawBondFlag )
		    {   TestAtomProximity(bptr->srcatom,xpos,ypos);
			TestAtomProximity(bptr->dstatom,xpos,ypos);
		    }

	    ForEachBack
		if( bptr->flag&DrawBondFlag )
		{   TestAtomProximity(bptr->srcatom,xpos,ypos);
		    TestAtomProximity(bptr->dstatom,xpos,ypos);
		}
	}

	ForEachAtom
	{   /* Identify bond! */
	    if( aptr == QAtom )
	    {   QChain = chain;
		QGroup = group;
	    }

	    if( aptr->flag & SphereFlag )
	    {   dy = AbsFun(aptr->y-ypos);
		if( dy>aptr->irad ) continue;
		rad = LookUp[aptr->irad][dy];
		dx = AbsFun(aptr->x-xpos);
		if( dx>rad ) continue;

		new = False;
		dpth = aptr->z+LookUp[rad][dx];
		if( UseSlabPlane && (aptr->z+rad>=SlabValue) )
		{   dz = SlabValue-aptr->z;
		    if( SlabMode && (dz >= -rad) )
		    {   wide = LookUp[aptr->irad][AbsFun(dz)];
			if( (dy<=wide) && (dx<=(int)(LookUp[wide][dy])) )
			{   if( SlabMode == SlabFinal )
			    {   dpth = SliceValue;
				new = True;
			    } else if( SlabMode == SlabHollow )
			    {   dpth = aptr->z-LookUp[rad][dx];
				new = !IdentFound || (dpth>IdentDepth);
			    } else if( SlabMode != SlabHalf )
			    {   /* SlabClose, SlabSection */
				dpth = dx*dx+dy*dy+dz*dz+SliceValue;
				if( IdentFound )
				{   new = (IdentDepth<SliceValue) 
					  || (dpth<IdentDepth);
				} else new=True;
			    }
			} else if( (dz>0) && (SlabMode!=SlabSection) )
			    new = !IdentFound || (dpth>IdentDepth);
		    }
		} else if( !UseSlabPlane || (SlabMode != SlabSection) )
		    new = !IdentFound || IdentDist || (dpth>IdentDepth);

		if( new )
		{   IdentFound = True;
		    IdentDepth = dpth;
		    IdentDist = 0;

		    QChain = chain;
		    QGroup = group;
		    QAtom = aptr;
		}
	    } 
	}
    } else /* Display Mode */
    {   ForEachAtom
	{   TestAtomProximity(aptr,xpos,ypos);
	    /* Identify bond! */
	    if( aptr == QAtom )
	    {   QChain = chain;
		QGroup = group;
	    }
	}
    }


    if( !IdentFound || (IdentDist>=50) )
    {   /* Reset Pick Atom! */
        QChain = (void __far*)0;
	QGroup = (void __far*)0;
	QAtom = (void __far*)0;
    }
}


void SetPickMode( mode )
    int mode;
{
    PickMode = mode;
    PickCount = 0;
}


static void DescribeAtom( ptr, flag )
    AtomRef *ptr;  int flag;
{
    register char *str;
    register int i,ch;
    char buffer[40];

    str = Residue[ptr->grp->refno];
    for( i=0; i<3; i++ )
        if( str[i]!=' ' ) 
             WriteChar(str[i]);

    sprintf(buffer,"%d",ptr->grp->serno);
    WriteString(buffer);

    ch = ptr->chn->ident;
    if( ch != ' ' )
    {   if( isdigit(ch) )
            WriteChar(':');
        WriteChar(ch);
    }

    WriteChar('.');
    str = ElemDesc[ptr->atm->refno];
    for( i=0; i<3; i++ )
        if( str[i]!=' ' ) 
             WriteChar(str[i]);

    if( flag )
    {   sprintf(buffer," (%d)",ptr->atm->serno);
        WriteString(buffer);
    }
}


void PickAtom( shift, xpos, ypos )
    int shift, xpos, ypos;
{
    register AtomRef *ptr;
    register Label *label;
    register float temp;
    register char *str;
    register int len;

    char buffer[40];
    AtomRef ref;


    if( PickMode == PickNone )
        return;

    IdentifyAtom(xpos,ypos);

    if( QAtom )
    {   if( PickMode == PickIdent )
        {   if( CommandActive )
	        WriteChar('\n');
            CommandActive = False;

            WriteString("Atom: ");
            str = ElemDesc[QAtom->refno];
            if( str[0]!=' ' )   WriteChar(str[0]);
            WriteChar(str[1]);  WriteChar(str[2]);
            if( str[3]!=' ' )   WriteChar(str[3]);

	    sprintf(buffer," %d  ",QAtom->serno);
            WriteString(buffer);

/*  Remove Group/Chain Output, JWP, June 2001

            str = Residue[QGroup->refno];
            WriteString("Group: ");

            if( str[0]!=' ' )  WriteChar(str[0]);
            WriteChar(str[1]); WriteChar(str[2]);

            sprintf(buffer," %d",QGroup->serno);
            WriteString(buffer);

            if( QChain->ident!=' ' )
            {   WriteString("  Chain: ");
                WriteChar(QChain->ident);
            }

End of Remove Group/Chain Output  */

            WriteChar('\n');

        } else if( PickMode == PickLabel )
        {   if( !QAtom->label )
            {   if( MainGroupCount > 1 )
                {   strcpy(buffer,"%n%r");
                    str = buffer+4;
                    if( Info.chaincount > 1 )
                    {   if( isdigit(QChain->ident) )
                            *str++ = ':';
                        *str++ = '%';
                        *str++ = 'c';
                    }
                    strcpy(str,".%a");

                    len = (str-buffer) + 3;
                    label = CreateLabel(buffer,len);
                } else label = CreateLabel("%e%i",4);

                QAtom->label = label;
                label->refcount++;
            } else
            {   DeleteLabel( (Label*)QAtom->label );
                QAtom->label = (void*)0;
            }

            DrawLabels = LabelList? True : False;
            ReDrawFlag |= RFRefresh;

        } else if( PickMode == PickCentr )
        {   CenX = QAtom->xorg;
            CenY = QAtom->yorg;
            CenZ = QAtom->zorg;

            ref.chn = QChain;
            ref.grp = QGroup;
            ref.atm = QAtom;

            if( CommandActive )
	        WriteChar('\n');
            CommandActive = False;

            WriteString("Rotating about ");
            DescribeAtom(&ref,True);
            WriteChar('\n');

        } else if( PickMode == PickMonit )
        {   /* State Machine Implementation */

            if( PickCount == 0 )
            {   PickHist[0].atm = QAtom;
                PickCount = 1;
            } else if( PickCount == 1 )
            {   if( !shift )
                {   if( PickHist[0].atm != QAtom )
                    {   AddMonitors(PickHist[0].atm,QAtom);
                        ReDrawFlag |= RFRefresh;
                    }
                    PickCount = 2;
                } else PickHist[0].atm = QAtom;
            } else /* PickCount == 2 */
                if( !shift )
                {   PickHist[0].atm = QAtom;
                    PickCount = 1;
                } else if( PickHist[0].atm != QAtom )   
                {   AddMonitors(PickHist[0].atm,QAtom);
                    ReDrawFlag |= RFRefresh;
                }

        } else /* Distance, Angle or Torsion! */
        {   if( PickCount )
            {   if( shift )
                {   PickCount--;
                } else if( PickCount == PickMode )
                    PickCount = 0;
            }

            ptr = PickHist+PickCount;
            ptr->chn = QChain;
            ptr->grp = QGroup;
            ptr->atm = QAtom;
            PickCount++;

            if( CommandActive )
	        WriteChar('\n');
            CommandActive = False;

            WriteString("Atom #");
            WriteChar(PickCount+'0');
            WriteString(": ");
            DescribeAtom(ptr,True);
            WriteChar('\n');

            if( PickCount == PickMode )
            {   if( PickMode == PickDist )
                {   temp = (float)CalcDistance(PickHist[0].atm,
                                               PickHist[1].atm);

                    WriteString("Distance ");
                    DescribeAtom(PickHist,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+1,False);
                    sprintf(buffer,": %.3f\n\n",temp);
                    WriteString(buffer);

                } else if( PickMode == PickAngle )
                {   temp = (float)CalcAngle(PickHist[0].atm,
                                            PickHist[1].atm,
                                            PickHist[2].atm);

                    WriteString("Angle ");
                    DescribeAtom(PickHist,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+1,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+2,False);
                    sprintf(buffer,": %.1f\n\n",temp);
                    WriteString(buffer);

                } else /* PickMode == PickTorsn */
                {   temp = (float)CalcTorsion(PickHist[0].atm,
                                              PickHist[1].atm,
                                              PickHist[2].atm,
                                              PickHist[3].atm);

                    WriteString("Torsion ");
                    DescribeAtom(PickHist,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+1,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+2,False);
                    WriteChar('-');
                    DescribeAtom(PickHist+3,False);
                    sprintf(buffer,": %.1f\n\n",temp);
                    WriteString(buffer);
                }
            }
        }
    }
}


void SetStereoMode( enable )
    int enable;
{
    ReDrawFlag |= RFRefresh | RFTransX;
    StereoView = ViewLeft;
    UseStereo = enable;
    DetermineClipping();
}


void ResetRenderer()
{
    DrawAtoms = False;  MaxAtomRadius = 0;
    DrawBonds = False;  MaxBondRadius = 0;

    SlabMode = SlabClose;
    UseSlabPlane = False;
    UseLabelCol = False;
    UseShadow = False;

    UseDepthCue = False;

    DisplayMode = 0;

    DrawBoundBox = False;
    DrawUnitCell = False;
    DrawAxes = False;

    SetStereoMode(False);
    StereoAngle = 6.0;
}


static void InitialiseTables()
{
    register Byte __far *ptr;
    register unsigned int root,root2;
    register unsigned int i,rad,arg;

    ptr = Array;
    LookUp[0] = ptr;  *ptr++ = 0;
    LookUp[1] = ptr;  *ptr++ = 1;  *ptr++ = 0;
    
    for( rad=2; rad<MAXRAD; rad++ )
    {   LookUp[rad] = ptr;

        /* i == 0 */
        *ptr++ = (Byte)rad;  

        root = rad-1;
        root2 = root*root;

        arg = rad*rad;
	for( i=1; i<rad; i++ )
        {   /* arg = rad*rad - i*i */
            arg -= (i<<1)-1;

            /* root = isqrt(arg)   */
            while( arg < root2 )
            {   root2 -= (root<<1)-1;
                root--;
            }
            /* Thanks to James Crook */
            *ptr++ = ((arg-root2)<i)? root : root+1;
        }

        /* i == rad */
        *ptr++ = 0;    
    }
}


void InitialiseRenderer()
{
    register int rad,maxval;

    FBuffer = (void __huge*)0;  
    DBuffer = (void __huge*)0;
    IBuffer = (void __far*)0;   
    YBucket = (void __far*)0;

#if defined(IBMPC) || defined(APPLEMAC)
    FBufHandle = NULL;
    DBufHandle = NULL;
#endif

#if defined(IBMPC) || defined(APPLEMAC)
    /* Allocate tables on FAR heaps */ 
    Array = (Byte __far*)_fmalloc(MAXTABLE*sizeof(Byte));
    LookUp = (Byte __far* __far*)_fmalloc(MAXRAD*sizeof(Byte __far*));
    HashTable = (void __far* __far*)_fmalloc(VOXSIZE*sizeof(void __far*));
    ColConst = (Card __far*)_fmalloc(MAXRAD*sizeof(Card));
    
    if( !Array || !LookUp || !HashTable || !ColConst )
	FatalRenderError("tables");
#else
    ColConst = ColConstTable;
#endif

    InitialiseTables();

    /* Initialise ColConst! */
    for( rad=0; rad<MAXRAD; rad++ )
    {   maxval = (int)(LightLength*rad)+4;
	ColConst[rad] = ((Card)ColourDepth<<ColBits)/maxval;
    }

    FreeItem = (Item __far*)0;
    PickMode = PickIdent;

    VoxelsClean = False;
    VoxelsDone = False;
    BucketFlag = False;
    RayCount = 0;

    ResetRenderer();
    ReSizeScreen();
}
