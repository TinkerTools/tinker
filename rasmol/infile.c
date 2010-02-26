/* infile.c
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

#define INFILE
#include "infile.h"
#include "molecule.h"
#include "abstree.h"
#include "command.h"
#include "transfor.h"

#ifndef APPLEMAC
#ifndef IBMPC
#include <sys/types.h>
#include <sys/time.h>
#endif
#include <time.h>
#endif
 

#ifdef MMIO
#include "mmio.h"
#endif

#define GroupPool 8


static Atom __far *ConnectAtom;
static char Record[202];
static FILE *DataFile;

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


static void FatalInFileError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"InFile Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


/*================================*/
/* File/String Handling Functions */
/*================================*/

static int FetchRecord()
{
    register char *ptr;
    register int ch;

    if( feof(DataFile) )
    {   *Record = '\0';
        return( False );
    }

    ptr = Record;
    do {
        ch = getc(DataFile);
        if( ch == '\n' )
        {   *ptr = 0;
            return( True );
        } else if( ch == '\r' )
        {   ch = getc(DataFile);
            if( ch != '\n' )
                ungetc(ch,DataFile);
            *ptr = 0;
            return( True );
        } else if( ch == EOF )
        {   *ptr = 0;
            return( ptr != Record+1 );
        } else *ptr++ = ch;
    } while( ptr < Record+200 );

    /* skip to the end of the line! */
    do { ch = getc(DataFile);
    } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );

    if( ch == '\r' )
    {   ch = getc(DataFile);
        if( ch != '\n' )
            ungetc(ch,DataFile);
    }
    *ptr = 0;
    return( True );
}


static void ExtractString( len, src, dst )
    int len;  char *src, *dst;
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( *src )
	{   ch = *src++;
            *dst++ = ch;
            if( ch != ' ' ) 
		ptr = dst;
	} else break;
    }
    *ptr = 0;
}


static Long ReadValue( pos, len )
    int pos, len;
{
    register Long result;
    register char *ptr;
    register char ch;
    register int neg;

    result = 0;
    neg = False;
    ptr = Record+pos;
    while( len-- )
    {   ch = *ptr++;
	if( (ch>='0') && (ch<='9') )
	{   result = (10*result)+(ch-'0');
	} else if( ch=='-' )
	    neg = True;
    }
    return( neg? -result : result );
}


/*==============================*/
/* Molecule File Format Parsing */
/*==============================*/ 

int LoadXYZMolecule( fp )
    FILE *fp;
{
    auto char type[12];
    auto double xpos, ypos, zpos;
    auto double charge, u, v, w;
    auto long atoms;

    auto long atmnum;
 
    register Atom __far *ptr;
    register char *src,*dst;
    register int count;
    register Long i;
 
 
    DataFile = fp;
    /* Number of Atoms */
    FetchRecord();
    sscanf(Record,"%ld",&atoms);
 
    /* Molecule Description */
    src = Record;
    while( *src == ' ' )
        src++;
    while( *src >= '0' && *src <= '9' )
        src++;
    while( *src == ' ' )
        src++;
 
    dst = Info.moleculename;
    for( i=0; i<78; i++ )
        if( *src ) *dst++ = *src++;
    *dst = '\0';
 
    if( atoms )
    {   CreateMolGroup();
        for( i=0; i<atoms; i++ )
        {   FetchRecord();
            ptr = CreateAtom();
            ptr->serno = i+1;
 
            xpos = ypos = zpos = 0.0;
            count = sscanf(Record,"%ld %s %lf %lf %lf %lf %lf %lf %lf",
                  &atmnum, type, &xpos, &ypos, &zpos, &charge, &u, &v, &w );
 
            ptr->refno = SimpleAtomType(type);
            ptr->xorg =  (Long)(250.0*xpos);
            ptr->yorg =  (Long)(250.0*ypos);
            ptr->zorg = -(Long)(250.0*zpos);
 
            ProcessAtom( ptr );

            CreateGroup( GroupPool );
            CurGroup->refno = i+2;
            CurGroup->serno = i+2;
            ProcessGroup( False );
        }
    }
    return( True );
}


int SaveXYZMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}
