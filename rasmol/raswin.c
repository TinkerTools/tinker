/* raswin.c
 * RasWin Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#include <windows.h>
#include <shellapi.h>
#include <commdlg.h>
#include <direct.h>
#include <dde.h>

#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>

#define RASMOL
#include "rasmol.h"
#include "raswin.idm"
#include "molecule.h"
#include "abstree.h"
#include "graphics.h"
#include "pixutils.h"
#include "transfor.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "outfile.h"

/* Microsoft C vs Borland Turbo C */
#ifndef __TURBOC__
#define stricmp _stricmp
#define getcwd  _getcwd
#endif

/* Determine Mouse Sensitivity! */
#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))

#define DefaultTimeOut 10000
#define MaxAdviseNum   32
#define MaxConvNum     8

#define ColdLink       0x01
#define WarmLink       0x02
#define HotLink        0x03
#define AckLink        0x04

typedef struct {
	    HWND  server;
	    HWND  client;
	    Byte  closed;
	} DDEConv;
	
typedef struct {
	    HANDLE data;
	    HWND  server;
	    HWND  client;
	    ATOM  atom;
	    Byte  mode;
	    Byte  item;
	    Byte  wait;
       } DDEAdvise;

static int DDETimeOut;
static DDEAdvise AdviseData[MaxAdviseNum];
static DDEConv ConvData[MaxConvNum];
static int RasWinReady;
static int AdviseCount;
static int ConvCount;

static char *ItemName[ItemCount] = {
	"Pick",      /* AdvPickAtom    */
	"PickNo",    /* AdvPickNumber  */
	"Count",     /* AdvSelectCount */
	"Name",      /* AdvName        */
	"Ident",     /* AdvIdent       */
	"Class",     /* AdvClass       */
	"Image",     /* AdvImage       */
	"PickXYZ"    /* AdvPickCoord   */
    };

#define CmndSize   (CmndRows*CmndCols)
#define ScrlMax    80
#define CmndRows   160
#define CmndCols   80

static int CmndStart;
static int ScrlStart;
static int TermCursor;
static int CharWide,CharHigh;
static int TermXPos,TermYPos;
static int TermRows,TermCols;
static char __far *TermScreen;
static HFONT TermFont;
static HWND CmndWin;

static int PointX,PointY;
static int LabelOptFlag;
static int InitX,InitY;
static int HeldButton;
static int FileFormat;

static char snamebuf[128];
static char fnamebuf[128];
static char ifilters[512];
static char ofilters[512];
static OPENFILENAME ofn1;
static OPENFILENAME ofn2;
static HANDLE hInstance;
static char Text[256];

long FAR PASCAL MainCallB(HWND, UINT, WPARAM, LPARAM);
long FAR PASCAL CmndCallB(HWND, UINT, WPARAM, LPARAM);
long FAR PASCAL DDECallB( HWND, UINT, WPARAM, LPARAM); 
BOOL FAR PASCAL AboutCallB(HWND, unsigned, WORD, LONG);
BOOL FAR PASCAL InfoCallB(HWND, unsigned, WORD, LONG);


static void CloseDDELinks()
{
    register long alarm;
    register int i;
    MSG message;
    
    for( i=0; i<MaxConvNum; i++ )
	if( ConvData[i].server )
	{   ConvData[i].closed = True;
	    PostMessage( ConvData[i].client, WM_DDE_TERMINATE,
			 (WPARAM)ConvData[i].server, 0L );
	} 
   
    alarm = GetTickCount() + DDETimeOut;
    while( PeekMessage(&message,NULL,WM_DDE_FIRST,WM_DDE_LAST,PM_REMOVE) )
    {   DispatchMessage( &message );
	if( message.message == WM_DDE_TERMINATE )
	    if( !ConvCount ) break;
	    
	/* Terminate Time Out */
	if( (long)GetTickCount() > alarm )
	{   for( i=0; i<MaxConvNum; i++ )
		if( ConvData[i].server )
		    DestroyWindow( ConvData[i].server );
	    break;
	}
    }
}


void RasMolExit()
{
    DeleteObject(TermFont);
    CloseDDELinks();
    CloseDisplay();
    exit(0);
}


void RasMolFatalExit( msg )
    char *msg;
{
    MessageBox(NULL,msg,"RasMol Fatal Error!",
	MB_OK | MB_ICONEXCLAMATION | MB_APPLMODAL );
    
    /* PostQuitMessage(0); */
    DeleteObject(TermFont);
    CloseDDELinks();
    CloseDisplay();
    exit(1);    
}


static void LoadInitFile()
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;

    fname = "RASMOL.INI";
    initrc = fopen(fname,"r");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
	while( *src )
	    *dst++ = *src++;
	*dst++ = '\\';

	src = fname; fname = fnamebuf;
	while( *dst++ = *src++ );
	initrc = fopen(fname,"r");
    }

    if( initrc )
	LoadScriptFile(initrc,fname);
}


static void SetTermScroll( pos )
    int pos;
{
    SetScrollPos(CmndWin,SB_VERT,pos,True);
    InvalidateRect(CmndWin,NULL,True);
    ScrlStart = ScrlMax - pos;
}


void WriteChar( ch )
    char ch;
{
    register int i;
    RECT rect;

    /* Scroll to bottom! */
    if( ScrlStart )
	SetTermScroll( ScrlMax );
	
    switch( ch )
    {    case(0x07):  MessageBeep(0);
		      break;
		      
	 case(0x08):  if( TermXPos>0 )
		      {   TermXPos--;
			  if( TermCursor )
			      SetCaretPos(TermXPos*CharWide,
					  TermYPos*CharHigh);
		      }
		      break;
		      
	 case(0x0D):  if( TermXPos )
		      {   if( TermCursor )
			      SetCaretPos(0,TermYPos*CharHigh);
			  TermXPos=0;
		      }
		      break;
		      
	 case(0x0A):  if( TermYPos==TermRows-1 )
		      {   CmndStart++;
			  if( CmndStart == CmndRows )
			      CmndStart = 0;
			      
			  i = TermYPos + CmndStart;
			  if( i >= CmndRows ) i -= CmndRows;
			  _fmemset(TermScreen+i*CmndCols,' ',CmndCols);
			  InvalidateRect(CmndWin,NULL,FALSE);
		      } else TermYPos++;
		      TermXPos = 0;

		      if( TermCursor )
			  SetCaretPos(0,TermYPos*CharHigh);
		      UpdateWindow(CmndWin);
		      break;
		      
	 
	 default:     i = TermYPos + CmndStart;
		      if( i >= CmndRows ) i -= CmndRows;
		      TermScreen[i*CmndCols+TermXPos]=ch;
		      if( TermXPos < TermCols )
		      {   rect.top = TermYPos*CharHigh; 
			  rect.left = TermXPos*CharWide;
			  rect.bottom = rect.top+CharHigh;
			  rect.right = rect.left+CharWide;
			  InvalidateRect(CmndWin,&rect,FALSE);
		      }
		      
		      if( TermXPos==CmndCols-1 )
		      {   if( TermYPos==TermRows-1 )
			  {   CmndStart++;
			      if( CmndStart == CmndRows )
				  CmndStart = 0;
				  
			      i = TermYPos + CmndStart;
			      if( i >= CmndRows ) i -= CmndRows;
			      _fmemset(TermScreen+i*CmndCols,' ',CmndCols);
			      InvalidateRect(CmndWin,NULL,FALSE);
			  } else TermYPos++;
			  TermXPos=0;
		      } else TermXPos++;

		      if( TermCursor )
			  SetCaretPos(TermXPos*CharWide,
				      TermYPos*CharHigh);
		      break;
		      
    }
}


void WriteString( ptr )
    char *ptr;
{
    while( *ptr )
	WriteChar(*ptr++);
}


static int InitTerminal( instance )
    HANDLE instance;
{
    TEXTMETRIC Text;
    LOGFONT LogFont;
    long style;
    RECT rect;
    HDC hDC;

    
    LogFont.lfHeight     = 12;
    LogFont.lfWidth      = 8;
    LogFont.lfEscapement = 0;
    LogFont.lfWeight     = 0;
    LogFont.lfItalic     = 0;
    LogFont.lfUnderline  = 0;
    LogFont.lfStrikeOut  = 0;
    
    LogFont.lfCharSet        = OEM_CHARSET;
    LogFont.lfOutPrecision   = OUT_DEFAULT_PRECIS;
    LogFont.lfClipPrecision  = CLIP_DEFAULT_PRECIS;
    LogFont.lfQuality        = DEFAULT_QUALITY;
    LogFont.lfPitchAndFamily = FIXED_PITCH | FF_MODERN;
    LogFont.lfFaceName[0]    = '\0';
    TermFont = CreateFontIndirect(&LogFont);

    /* TermFont = GetStockObject(ANSI_FIXED_FONT); */

    /* Default Window Size */
    TermCols = 80;  TermRows = 24;
    ScrlStart = CmndStart = 0;
    TermXPos = TermYPos = 0;
    
    TermScreen = (char __far*)_fmalloc(CmndSize*sizeof(char));
    if( !TermScreen ) return( False );
    _fmemset(TermScreen,' ',CmndSize);
    TermCursor = False;


    hDC = GetDC(NULL);
    SelectObject(hDC,TermFont);
    GetTextMetrics(hDC,&Text);  
    ReleaseDC(NULL,hDC);
    
    CharWide = Text.tmAveCharWidth;
    CharHigh = Text.tmHeight + Text.tmExternalLeading;

    rect.top  = 0;   rect.bottom = TermRows*CharHigh;
    rect.left = 0;   rect.right  = TermCols*CharWide;
    
    style = WS_OVERLAPPED | WS_CAPTION | WS_THICKFRAME | 
	    WS_MAXIMIZEBOX | WS_MINIMIZEBOX | WS_VSCROLL;
    
    AdjustWindowRect(&rect,style,False);
#ifdef _WIN32
    CmndWin = CreateWindowEx(WS_EX_CLIENTEDGE,
                           "RasCliClass", "RasMol Command Line",
			   style, CW_USEDEFAULT, CW_USEDEFAULT,
			   rect.right-rect.left, rect.bottom-rect.top,
			   NULL, NULL, instance, NULL );
#else
    CmndWin = CreateWindow("RasCliClass", "RasMol Command Line",
			   style, CW_USEDEFAULT, CW_USEDEFAULT,
			   rect.right-rect.left, rect.bottom-rect.top,
			   NULL, NULL, instance, NULL );
#endif
			   
    if( !CmndWin ) return( False );
   
    SetScrollRange(CmndWin,SB_VERT,0,ScrlMax,FALSE); 
    SetScrollPos(CmndWin,SB_VERT,ScrlMax,FALSE);
    ShowWindow(CmndWin,SW_SHOWMINNOACTIVE);
    return( True );
}


static void PaintScreen()
{
    int SRow,ERow,SCol,ECol;
    register char __far *ptr;
    register int row,len;
    register int x,y;
    
    PAINTSTRUCT ps;
    HFONT font;
    RECT rect;
    HDC hDC;
    
    hDC = BeginPaint(CmndWin,&ps);
    font = SelectObject(hDC,TermFont);
    SetBkColor(hDC,GetSysColor(COLOR_WINDOW));
    SetTextColor(hDC,RGB(0,0,0));
    SetBkMode(hDC,OPAQUE);
    
    SRow = ps.rcPaint.top/CharHigh;
    if( SRow >= TermRows )
    {   SRow = TermRows-1;
    } else if( SRow < 0 )
	SRow = 0;

    ERow = ps.rcPaint.bottom/CharHigh;
    if( ERow >= TermRows )
    {   ERow = TermRows-1;
    } else if( ERow < 0 ) 
	ERow = 0;
    
    SCol = ps.rcPaint.left/CharWide;
    if( SCol >= TermCols )
    {   SCol = TermCols-1;
    } else if( SCol < 0 ) 
	SCol = 0;
    
    ECol = ps.rcPaint.right/CharWide;
    if( ECol >= TermCols )
    {   ECol = TermCols-1;
    } else if( ECol < 0 ) 
	ECol = 0;

    len = ECol-SCol+1;
    x = SCol*CharWide;
    y = SRow*CharHigh;
    
    rect.right = x+len*CharWide;
    rect.left = x;   

    SRow += CmndStart - ScrlStart;
    if( SRow >= CmndRows )
    {   SRow -= CmndRows;
    } else if( SRow < 0 )
	SRow += CmndRows;
	
    ERow += CmndStart - ScrlStart;
    if( ERow >= CmndRows )
    {   ERow -= CmndRows;
    } else if( ERow < 0 )
	ERow += CmndRows;
	
    row = SRow;
    ptr = TermScreen + CmndCols*row + SCol;
    while( True )
    {   rect.top = y;    
	rect.bottom = y+CharHigh;
	ExtTextOut(hDC,x,y,ETO_OPAQUE|ETO_CLIPPED,
		   &rect,ptr,len,NULL);
		   
	if( row != ERow )
	{   ptr += CmndCols;
	    row++;
	    if( row == CmndRows )
	    {   ptr = TermScreen + SCol;
		row = 0;
	    }
	} else break;
	y += CharHigh;
    }
    
    SelectObject(hDC,font);
    EndPaint(CmndWin,&ps);
    
    if( TermCursor )
    {   row = TermYPos + ScrlStart;
	if( row < TermRows )
	{   SetCaretPos(TermXPos*CharWide,row*CharHigh);
	    ShowCaret(CmndWin);
	} else HideCaret(CmndWin);
    }    
}


BOOL FAR PASCAL AboutCallB(hWin,uMsg,wArg,lArg)
    HWND hWin; unsigned uMsg; WORD wArg; LONG lArg;
{
#ifdef _WIN32
    auto SYSTEM_INFO sysinfo;
    register char *ptr;
#else
    register DWORD flags;
#endif
    register int len; 
   
    switch(uMsg)
    {   case(WM_INITDIALOG):  
#ifdef _WIN32
				GetSystemInfo(&sysinfo);
				if( sysinfo.dwNumberOfProcessors > 1 )
				{   len = sprintf(Text,"%d processor ",
					       sysinfo.dwNumberOfProcessors);
				} else len = 0;

				switch(sysinfo.dwProcessorType)
				{   case(386): ptr = "Intel 386";     break;
				    case(486): ptr = "Intel 486";     break;
				    case(586): ptr = "Intel Pentium"; break;
				    case(860): ptr = "Intel i860";    break;

				    case(2000):  ptr = "MIPS R2000"; break;
				    case(3000):  ptr = "MIPS R3000"; break;
				    case(4000):  ptr = "MIPS R4000"; break;

				    case(21064): ptr = "DEC Alpha";  break;
				    default:     ptr = "unrecognised";
				}
				sprintf(Text+len,"%s machine",ptr);
#else
				flags = GetWinFlags();
				if( flags & WF_CPU286 )
				{      len = sprintf(Text,"286 with");
				} else if( flags & WF_CPU386 )
				{      len = sprintf(Text,"386 with");
				} else len = sprintf(Text,"486 with");
			      
				if( !(flags&WF_80x87) )
				{   sprintf(Text+len,"out maths coprocessor");
				} else sprintf(Text+len," maths coprocessor"); 
#endif
				SetDlgItemText(hWin,IDD_HARDWARE,Text);         
				return(TRUE);
    
	case(WM_COMMAND):     if( wArg == IDOK )
			      {   EndDialog(hWin,TRUE);
				  return(TRUE);
			      }
			      break;
    }
    return(FALSE);
}


static void DisplayMoleculeInfo( hWin )
    HWND hWin;
{
    register int line;
    register int len;

    line = IDD_INFOTEXT;
    
    if( *Info.moleculename )
    {   sprintf(Text," Molecule Name ...... %s",Info.moleculename);
	SetDlgItemText(hWin,line++,Text);
    }
    
    if( *Info.classification )
    {   sprintf(Text," Classification ......... %s",Info.classification);
	SetDlgItemText(hWin,line++,Text);
    }
    
    if( *Info.identcode )
    {   sprintf(Text," Brookhaven code .. %s",Info.identcode);
	SetDlgItemText(hWin,line++,Text);
    }
    
    if( Info.chaincount>1 )
    {   sprintf(Text," Number of chains .. %d",Info.chaincount);
	SetDlgItemText(hWin,line++,Text);
    }
    
    len = sprintf(Text," Number of groups .. %d",MainGroupCount);
    SetDlgItemText(hWin,line++,Text);

    len = sprintf(Text," Number of atoms ... %ld",MainAtomCount);
    SetDlgItemText(hWin,line++,Text);

    sprintf(Text," Number of bonds ... %ld",Info.bondcount);
    SetDlgItemText(hWin,line++,Text);
}

	
BOOL FAR PASCAL InfoCallB(hWin,uMsg,wArg,lArg)
    HWND hWin; unsigned uMsg; WORD wArg; LONG lArg;
{
    switch(uMsg)
    {   case(WM_INITDIALOG):  DisplayMoleculeInfo(hWin);
			      return(TRUE);
			      
	case(WM_COMMAND):     if( wArg == IDOK )
			      {   EndDialog(hWin,TRUE);
				  return(TRUE);
			      }
			      break;
    }
    return(FALSE);
}


static char *GetItemName( item )
    int item;
{
    switch( item )
    {   case(-1):  return("Topics");
	case(-2):  return("SysItems");
	case(-3):  return("Formats");
	case(-4):  return("Status");
	case(-5):  return("Items");
    }

    if( item<=ItemCount )
	return( ItemName[item-1] );
    return( "" );
}


static HANDLE RenderClipboard( format )
    WPARAM format;
{
    register BITMAPINFO __far *bitmap;
    register char __huge *src;
    register char __huge *dst;
    register HANDLE result;
    register long size,len;
    register int i; 
   
    if( format==CF_PALETTE )
    {   if( ColourMap )
	{   return( CreatePalette(Palette) );
	} else return( NULL );
    }    
    
    if( !PixMap || (format!=CF_DIB) )
	return( NULL );

    len = (long)XRange*YRange*sizeof(Pixel);
    size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    if( !(result=GlobalAlloc(GHND,size+len)) ) return( NULL );
    
    bitmap = (BITMAPINFO __far *)GlobalLock(result);
    bitmap->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bitmap->bmiHeader.biWidth = XRange;
    bitmap->bmiHeader.biHeight = YRange;
    bitmap->bmiHeader.biPlanes = 1;
    bitmap->bmiHeader.biBitCount = 8;
    bitmap->bmiHeader.biCompression = BI_RGB;
    bitmap->bmiHeader.biSizeImage = len;
    bitmap->bmiHeader.biXPelsPerMeter = 0;
    bitmap->bmiHeader.biYPelsPerMeter = 0;
    bitmap->bmiHeader.biClrImportant = 0;
    bitmap->bmiHeader.biClrUsed = 0;
    
    for( i=0; i<256; i++ )
	if( ULut[i] )
	{   bitmap->bmiColors[Lut[i]].rgbBlue  = BLut[i];
	    bitmap->bmiColors[Lut[i]].rgbGreen = GLut[i];
	    bitmap->bmiColors[Lut[i]].rgbRed   = RLut[i];
	}
    
   
    src = (Pixel __huge*)GlobalLock(FBufHandle);
    dst = ((Pixel __huge*)bitmap)+size;
    
    /* Transfer the frame buffer */
    while( len-- ) *dst++ = *src++;
    
    GlobalUnlock(FBufHandle);
    GlobalUnlock(result);
    return( result );    
}


static void SendItemData( hSrc, hDst, mode, item, advise )
    HWND hSrc, hDst;  int mode, item, advise;
{
    DDEDATA FAR *data;
    HANDLE FAR *hImage;
    HANDLE hData;
    ATOM atom;

    register char __far *dest;
    register char *src, *dst;
    register char *name;
    register Long len;
    register int i;

    name = GetItemName(item);

    if( mode==WarmLink )
    {   atom = GlobalAddAtom(name);
	if( !PostMessage(hDst,WM_DDE_DATA,(WPARAM)hSrc,MAKELONG(0,atom)) )
	    GlobalDeleteAtom(atom);
	return;
    }

    dst = Text;
    if( item>0 )
	item--;
	
    switch( item )
    {   case(-1): /* Topics */
		  src="System\tRemoteControl"; 
		  while( *dst++ = *src++ ); break;

	case(-2): /* SysItems */
		  src = "Topics\tSysItems\tFormats\tStatus\tItems";
		  while( *dst++ = *src++ ); break;

	case(-3): /* Formats */
		  src = "DIB\tTEXT\tPalette\tLink";
		  while( *dst++ = *src++ ); break;

	case(-4): /* Status */
		  src = RasWinReady? "Ready" : "Busy";
		  while( *dst++ = *src++ ); break;

	case(-5): /* Items */
		  for( i=0; i<ItemCount; i++ )
		  {   if( i ) *dst++ = '\t';
		      src = GetItemName(i); 
		      while( *src )
			  *dst++ = *src++;
		  }
		  *dst = '\0';
		  break;

	case(AdvPickAtom):
		  if( QAtom )
		  {   src = Residue[QGroup->refno];
		      if( src[0]!=' ' ) *dst++ = src[0];
		      *dst++  = src[1]; *dst++ = src[2];

		      sprintf(dst,"%d",QGroup->serno);
		      for( dst=Text; *dst; dst++ );
		      if( QChain->ident!=' ' )
		      {   *dst++ = ':';
			  *dst++ = QChain->ident;
		      }
		      *dst++ = '.';
		      
		      src = ElemDesc[QAtom->refno];
		      if( src[0]!=' ' ) *dst++ = src[0];
		      *dst++  = src[1]; *dst++ = src[2];
		      if( src[3]!=' ' ) *dst++ = src[3];
		  } 
		  *dst = '\0';
		  break;

	case(AdvPickNumber):
		  if( QAtom )
		  { sprintf(dst,"%d",QAtom->serno);
		  } else *dst = '\0';
		  break;

	case(AdvSelectCount):
		  sprintf(dst,"%ld",SelectCount);
		  break;
		  
	case(AdvName):
		  src = Info.moleculename;
		  while( *dst++ = *src++ );
		  break;

	case(AdvPickCoord):
		  if( QAtom )
		  { sprintf( dst, "%ld\t%ld\t%ld",
			     QAtom->xorg, QAtom->yorg, QAtom->zorg);
		  } else *dst = '\0';
		  break;

	default:  *dst = '\0';
		  break;
    }

    len = sizeof(DDEDATA);
    if( item == AdvImage )
    {   len += sizeof(HANDLE);
    } else for( dst=Text; *dst; dst++ )
	len++;

    if( hData = GlobalAlloc(GHND|GMEM_DDESHARE,len) )
    {   if( data = (DDEDATA FAR*)GlobalLock(hData) )
	{   data->fResponse = (mode!=AckLink);
	    data->fAckReq = (mode==ColdLink);
	    data->fRelease = True;

	    if( item == AdvImage )
	    {   data->cfFormat = CF_DIB;
		hImage = (HANDLE __far*)&data->Value[0];
		*hImage = RenderClipboard(CF_DIB);

	    } else 
	    {   data->cfFormat = CF_TEXT;
		dest = (char __far*)&data->Value[0];
		for( src=Text; *src; *dest++ = *src++ );
		/* Correctly terminate the data string */
		/* *dest++ = '\r'; *dest++ = '\n';     */
		*dest = '\0';
	    }
	    
	    GlobalUnlock(hData);
	    atom = GlobalAddAtom(name);
	    if( PostMessage(hDst,WM_DDE_DATA,(WPARAM)hSrc,MAKELONG(hData,atom)) )
	    {   if( mode==AckLink )
		{   SetTimer( hSrc, (UINT)hDst, DDETimeOut, NULL );
		    AdviseData[advise].data = hData;
		    AdviseData[advise].atom = atom;
		    AdviseData[advise].wait = True;
		}
		return;
	    }
	    GlobalDeleteAtom(atom);
	}
	GlobalFree( hData );
    }
    return;
}


static int GetItemNumber( atom )
    ATOM atom;
{
    register int i;

    GlobalGetAtomName(atom,Text,240);

    for( i=1; i<6; i++ )
	if( !stricmp(Text,GetItemName(-i)) )
	    return( -i );

    for( i=0; i<ItemCount; i++ )
	if( !stricmp(Text,ItemName[i]) )
	    return( i+1 );
    return( 0 );
}


void AdviseUpdate( item )
    int item;
{
    register DDEAdvise *ptr;
    register int i;

    if( AdviseCount )
    {   if( item >= 0 ) item++;
	for( i=0; i<MaxAdviseNum; i++ )
	{   ptr = AdviseData + i;
	    if( ptr->server && (ptr->item==(Byte)item) )
		SendItemData(ptr->server,ptr->client,ptr->mode,item,i);
	}
    }
}


void RefreshScreen()
{
    ReDrawFlag &= ~(RFTransZ|RFPoint);

    if( ReDrawFlag )
    {   if( RasWinReady )
	{   RasWinReady = False;
	    AdviseUpdate( -4 );
	}
	
	if( ReDrawFlag & RFReSize )
	    ReSizeScreen();

	if( ReDrawFlag & RFColour )
	{   ClearImage();
	    DefineColourMap();
	}

	if( Database )
	{   BeginWait();
	    if( ReDrawFlag & RFApply ) 
		ApplyTransform();
	    DrawFrame();
	    TransferImage();
	    EndWait();
	} else
	{   ClearBuffers();
	    TransferImage();
	}
	ReDrawFlag = 0;
    }
}


long FAR PASCAL DDECallB(hWin,uMsg,wArg,lArg)
    HWND hWin; UINT uMsg; WPARAM wArg; LPARAM lArg;
{
    DDEADVISE FAR *options;
    HWND hDest;

    register DDEConv *ptr;
    register char __huge *cmnd;
    register int item, stat;
    register int format,i;
    register int flag;
    
    stat = IPC_Ok;
    flag = False;

    switch( uMsg )
    {   case( WM_TIMER ):    
		    /* Simulate DDE NAck */
		    lArg = 0L;  
	case( WM_DDE_ACK ):
		    KillTimer( hWin, wArg );
		    if( !(LOWORD(lArg)&0x8000) )
		    {   item = GetItemNumber( HIWORD(lArg) );
			for( i=0; i<MaxAdviseNum; i++ )
			    if( (AdviseData[i].server==hWin) &&
				(AdviseData[i].item==(Byte)item) &&
				 AdviseData[i].wait )
			    {   GlobalDeleteAtom(AdviseData[i].atom);
				GlobalFree(AdviseData[i].data);
				AdviseData[i].wait = False;
				break;
			    }
		    }

		    if( HIWORD(lArg) )
			GlobalDeleteAtom( HIWORD(lArg) );
		    return( 0L );

	case( WM_DDE_REQUEST ):
		    if( item=GetItemNumber(HIWORD(lArg)) )
		    {   format = (item==AdvImage+1)? CF_DIB : CF_TEXT; 
			if( format == (int)(LOWORD(lArg)) )
			{   SendItemData( hWin, wArg, ColdLink, item, 0 );
			    GlobalDeleteAtom( HIWORD(lArg) );
			    return( 0L );
			}
		    } 
		    break;

	case( WM_DDE_UNADVISE ):
		    if( HIWORD(lArg) )
		    {   item = GetItemNumber( HIWORD(lArg) );
			if( !item ) break;
		    } else item = 0;

		    for( i=0; i<MaxAdviseNum; i++ )
			if( (AdviseData[i].server==hWin) &&
			    ( !item || AdviseData[i].item==(Byte)item ) )
			{   if( AdviseData[i].wait )
			    {   GlobalDeleteAtom(AdviseData[i].atom);
				GlobalFree(AdviseData[i].data);
			    }
			    AdviseData[i].server = NULL;
			    AdviseCount--;
			    flag = True;
			}
		    break;

	case( WM_DDE_ADVISE ):
		    item = GetItemNumber( HIWORD(lArg) );
		    if( !item || (AdviseCount==MaxAdviseNum ) )
			break;

		    /* Check for established link! */
		    for( i=0; i<MaxAdviseNum; i++ )
			if( (AdviseData[i].server==hWin) &&
			    (AdviseData[i].item==(Byte)item) ) break;
		    if( i<MaxAdviseNum ) break;

		    options = (DDEADVISE FAR*)GlobalLock((HGLOBAL)LOWORD(lArg));
		    if( !options ) break;

		    format = (item==AdvImage+1)? CF_DIB : CF_TEXT;
		    if( options->cfFormat == format ) 
		    {   for( i=0; i<MaxConvNum; i++ )
			   if( ConvData[i].server==hWin )
			   {   hDest = ConvData[i].client;
			       break;
			   }

		       for( i=0; i<MaxAdviseNum; i++ )
			   if( !AdviseData[i].server )
			       break;

		       AdviseData[i].server = hWin;
		       AdviseData[i].client = hDest;
		       AdviseData[i].atom = HIWORD(lArg);
		       AdviseData[i].wait = False;
		       AdviseData[i].item = item;
		       AdviseCount++;

		       if( options->fDeferUpd )
		       {      AdviseData[i].mode = WarmLink;
		       } else if( options->fAckReq )
		       {      AdviseData[i].mode = AckLink;
		       } else AdviseData[i].mode = HotLink;

		       PostMessage( (HWND)wArg, WM_DDE_ACK, (WPARAM)hWin,
				 MAKELONG(0x8000,HIWORD(lArg)) ); 
		    /* SendItemData(hWin,hDest,AdviseData[i].mode,item,i); */
		    }
		    GlobalUnlock((HGLOBAL)LOWORD(lArg));
		    return( 0L );

	case( WM_DDE_EXECUTE ):  
		    if( cmnd=(char __huge*)GlobalLock((HANDLE)HIWORD(lArg)) )
		    {   flag = stat = ExecuteIPCCommand( cmnd );
			GlobalUnlock((HANDLE)HIWORD(lArg));
		    }
		    break;

		    
	case( WM_DDE_TERMINATE ):
		    /* Destroy all Hot/Warm Links */
		    for( i=0; i<MaxAdviseNum; i++ )
			if( AdviseData[i].server == hWin )
			{   AdviseData[i].server = NULL;
			    if( AdviseData[i].wait )
			    {   GlobalDeleteAtom(AdviseData[i].atom);
				GlobalFree(AdviseData[i].data);
			    }
			    AdviseCount--;
			}

		    /* Remove the Conversation */
		    for( i=0; i<MaxConvNum; i++ )
			if( ConvData[i].server == hWin )
			{   ptr = ConvData+i;
			    if( !ptr->closed )
				PostMessage( ptr->client, WM_DDE_TERMINATE,
					     (WPARAM)ptr->server, 0L );
			    DestroyWindow( ptr->server );
			    ptr->server = NULL;
			    ConvCount--;
			    break;
			}
		    return( 0L );
		    
	default:  return( DefWindowProc(hWin,uMsg,wArg,lArg) );
    }

    /* Return a DDE acknowledgement */
    PostMessage( (HWND)wArg, WM_DDE_ACK, (WPARAM)hWin,
		 MAKELONG( (flag?0x8000:0), HIWORD(lArg)) ); 

    if( (stat==IPC_Quit) || (stat==IPC_Exit) )
        RasMolExit();

    if( ReDrawFlag ) 
	RefreshScreen();
    if( !CommandActive ) 
	ResetCommandLine(0);
    return( 0L );
}


static void ResizeTerminal( x, y )
    int x, y;
{
    register int rows, cols;
    register int sr, er;

    HBRUSH hBr;
    RECT rect;
    HDC hDC;
    
    if( x > CharWide )
    {   cols = x/CharWide;
    } else cols = 1;
    
    if( y > CharHigh )
    {   rows = y/CharHigh;
    } else rows = 1;

    /* Scroll to bottom! */
    if( ScrlStart )
	SetTermScroll( ScrlMax );

    if( rows < TermRows )
    {   if( TermYPos >= rows )
	{   CmndStart += (TermYPos - rows) + 1;
	    if( CmndStart >= CmndRows )
		CmndStart -= CmndRows;
	    TermYPos = rows - 1;
	
	    hDC = GetDC(CmndWin);
	    GetClientRect(CmndWin,&rect);
	    hBr = CreateSolidBrush(GetSysColor(COLOR_WINDOW));
	    FillRect(hDC,&rect,hBr);
	    ReleaseDC(CmndWin,hDC);
	} 

    } else if( rows > TermRows )
    {   sr = TermRows + CmndStart;
	if( sr >= CmndRows )
	    sr -= CmndRows;
	    
	er = CmndStart + rows;
	if( er >= CmndRows )
	    er -= CmndRows;
	    
	do {
	    _fmemset(TermScreen+sr*CmndCols,' ',CmndCols);
	    sr++; if( sr == CmndRows ) sr = 0;
	} while( sr != er );
    }
    
    InvalidateRect(CmndWin,NULL,False);
    if( cols > CmndCols )
    {   TermCols = CmndCols;
    } else TermCols = cols;
    TermRows = rows;
}


long FAR PASCAL CmndCallB(hWin,uMsg,wArg,lArg)
    HWND hWin; UINT uMsg; WPARAM wArg; LPARAM lArg;
{
    register int row;
    
    switch(uMsg)
    {   case(WM_CLOSE):       DestroyWindow(CanvWin);
			      DestroyWindow(CmndWin);
			      CommandActive = True;
			      ReDrawFlag = False;
			      break;

	case(WM_DESTROY):     /* Destroy RasWin */
			      PostQuitMessage(0);
			      break;

	 case(WM_SYSCHAR):    if( lArg & (1L<<29) )  /* ALT-key pressed? */
				  return(DefWindowProc(hWin,uMsg,wArg,lArg));

	 case(WM_CHAR):       if( ProcessCharacter(LOBYTE(wArg)) )
				  if( ExecuteCommand() )
				      RasMolExit();
			      break;

	 case(WM_PAINT):      PaintScreen();
			      return(0L);
			      
	case(WM_SYSKEYDOWN):
	case(WM_KEYDOWN):     switch(LOBYTE(wArg))
			      {   case(0x23): ProcessCharacter(0x05); break;
				  case(0x24): ProcessCharacter(0x01); break;
				  case(0x25): ProcessCharacter(0x02); break;
				  case(0x26): ProcessCharacter(0x10); break;
				  case(0x27): ProcessCharacter(0x06); break;
				  case(0x28): ProcessCharacter(0x0e); break;
				  case(0x2e): ProcessCharacter(0x04); break;
				  
				  default:
				  return(DefWindowProc(hWin,uMsg,wArg,lArg));
			      }
			      break;
	
	 case(WM_SETFOCUS):   if( !TermCursor )
			      {   CreateCaret(hWin,NULL,CharWide,CharHigh);
				  TermCursor = True;
			      }
			      
			      row = TermYPos + ScrlStart;
			      if( row < TermRows )
			      {   SetCaretPos(TermXPos*CharWide,row*CharHigh);
				  ShowCaret(hWin);
			      } else HideCaret(hWin);
			      return(0L);
			      
	 case(WM_SIZE):       if( wArg != SIZE_MINIMIZED )
				  ResizeTerminal(LOWORD(lArg),HIWORD(lArg));
			      return(0L);
			      
	 case(WM_KILLFOCUS):  if( TermCursor )
			      {   TermCursor=False;
				  HideCaret(hWin);
				  DestroyCaret();
			      }
			      return(0L);

	 case(WM_VSCROLL):    switch( wArg )
			      {  case(SB_TOP):    SetTermScroll(0);  break;
				 case(SB_BOTTOM): SetTermScroll(ScrlMax);  
						  break;
				 
				 case(SB_LINEUP):   
				     if( ScrlStart < ScrlMax )
					 SetTermScroll((ScrlMax-ScrlStart)-1);
				     break;
				     
				 case(SB_LINEDOWN):
				     if( ScrlStart > 0 )
					 SetTermScroll((ScrlMax-ScrlStart)+1);
				     break;
				     
				 case(SB_PAGEUP):
				     if( ScrlStart < (ScrlMax-10) )
				     {   SetTermScroll((ScrlMax-ScrlStart)-10);
				     } else SetTermScroll(0);
				     break;
				     
				 case(SB_PAGEDOWN):
				     if( ScrlStart > 10 )
				     {   SetTermScroll((ScrlMax-ScrlStart)+10);
				     } else SetTermScroll(ScrlMax);
				     break;
				     
				 case(SB_THUMBTRACK):
				 case(SB_THUMBPOSITION):
				     SetTermScroll(LOWORD(lArg));
				     break;
			      }
			      break;
							    
	 default:  return( DefWindowProc(hWin,uMsg,wArg,lArg) );
    }

    if( ReDrawFlag )
	RefreshScreen();
    if( !CommandActive )
	ResetCommandLine(0);
    return(0L);
}


static void LoadInputFile( format )
    int format;
{
    register char *ext;
    register int num;

    switch( format )
    {   case(FormatXYZ):      ext = "XYZ";  num = 5;  break;
    }

    ofn1.nFilterIndex = num;
    ofn1.lpstrDefExt = ext;
    *fnamebuf = '\0';

    if( GetOpenFileName(&ofn1) )
    {   switch( ofn1.nFilterIndex )
	{   case(2): FetchFile(FormatXYZ,False,fnamebuf);     break;
	}
        DefaultRepresentation();
    }
}


static void SaveOutputFile( format )
    int format;
{
    register char *ext;
    register int num;

    switch( format )
    {   case(IDM_GIF):   ext="GIF";  num=1;   break;
	case(IDM_EPSF):  ext="PS";   num=2;   break;
    }

    ofn2.nFilterIndex = num;
    ofn2.lpstrDefExt = ext;
    *fnamebuf = '\0';
    
/*  Default Filename   
 *  dst = fnamebuf;
 *  for( src="RASWIN."; *src; src++ ) *dst++ = *src;
 *  for( src=ext; *src; src++ ) *dst++ = *src;
 *  *dst++ = '\0';
 */
    
    if( GetSaveFileName(&ofn2) )    
	switch( ofn2.nFilterIndex )
	{   case(1):  WriteGIFFile(fnamebuf);             break;
	    case(2):  WriteEPSFFile(fnamebuf,True,True);  break;
	    case(3):  WriteEPSFFile(fnamebuf,False,True); break;
            case(4):  WriteVectPSFile(fnamebuf);          break;
            case(5):  WritePICTFile(fnamebuf);            break;
            case(6):  WriteIRISFile(fnamebuf);            break;
	}
}


static void HandlePrintSetUp()
{
    PRINTDLG pd;

    memset(&pd,0,sizeof(PRINTDLG));
    pd.lStructSize = sizeof(PRINTDLG);
    pd.hwndOwner = CanvWin;
    pd.Flags = PD_PRINTSETUP;

    PrintDlg(&pd);

    if( pd.hDevNames ) GlobalFree(pd.hDevNames);
    if( pd.hDevMode )  GlobalFree(pd.hDevMode);
}


static BOOL HandleMenu( option )
    WPARAM option;
{
    register char *src, *dst;
    register FARPROC lpProc;
    register int mask;
   
    switch(option)
    {   /* File Menu */
	case(IDM_OPEN):   if( !Database )
			      LoadInputFile(FormatXYZ);
			  break;
			  
	case(IDM_INFO):   lpProc = MakeProcInstance(InfoCallB,hInstance);
			  DialogBox(hInstance,"InfoBox",CanvWin,lpProc);
			  FreeProcInstance(lpProc);
			  break;
			  
	case(IDM_CLOSE):  ZapDatabase();
			  break;

	case(IDM_PRINT):  if( !PrintImage() )
			  {   if( CommandActive )
				  WriteChar('\n');
			      WriteString("Warning: No suitable printer!\n");
			      CommandActive = False;
			  }
			  break;
	
	case(IDM_SETUP):  HandlePrintSetUp();
                          break;

	case(IDM_EXIT):   PostMessage(CanvWin,WM_CLOSE,0,0L);
			  break;
			  

        /* Edit Menu */
        case(IDM_SELECT): mask = NormAtomFlag;
                          if( Hydrogens )  mask |= HydrogenFlag;
                          SelectZone(mask);
                          break;

	case(IDM_COPY):   if( !ClipboardImage() )
			  {   if( CommandActive )
				  WriteChar('\n');
			      WriteString("Unable to copy to clipboard!\n");
			      CommandActive = False;
			  }
			  break;
	
	/* Help Menu */
	case(IDM_ABOUT):  lpProc = MakeProcInstance(AboutCallB,hInstance);
			  DialogBox(hInstance,"AboutBox",CanvWin,lpProc);
			  FreeProcInstance(lpProc);
			  break;

	case(IDM_HELP):   if( getcwd(fnamebuf,100) )
			  {   dst = fnamebuf;
			      while( *dst ) dst++;
			      if( *(dst-1) != '\\' ) 
				  *dst++ = '\\';
				  
			      src = "RASWIN.HLP";    
			      while( *dst++ = *src++ );
			      WinHelp(CanvWin,fnamebuf,HELP_INDEX,0L);
			  }
			  break;
       

	/* Display Menu */
	case(IDM_WIREFRAME):  DisableSpacefill();
			      EnableWireframe(WireFlag,0);
			      ReDrawFlag |= RFRefresh;
			      break;

	case(IDM_STICKS):     DisableSpacefill();
			      if( MainAtomCount<256 )
			      {   EnableWireframe(CylinderFlag,40);
			      } else EnableWireframe(CylinderFlag,80);
			      ReDrawFlag |= RFRefresh;
			      break;

	case(IDM_SPHERES):    SetVanWaalRadius();
			      DisableWireframe();
			      ReDrawFlag |= RFRefresh;
			      break;

	case(IDM_BALLSTICK):  SetRadiusValue(120);
			      EnableWireframe(CylinderFlag,40);
			      ReDrawFlag |= RFRefresh;
			      break;

	/* Colours Menu */
	case(IDM_MONO):     MonoColourAttrib(255,255,255);
			    ReDrawFlag |= RFColour;  break;

	case(IDM_CPK):      CPKColourAttrib();
			    ReDrawFlag |= RFColour;  break;

	case(IDM_SHAPELY):  ShapelyColourAttrib();
			    ReDrawFlag |= RFColour;  break;
       
	/* Options Menu */
	case(IDM_SLAB):      ReDrawFlag |= RFRefresh;
			     UseSlabPlane = !UseSlabPlane;
			     if( UseSlabPlane )
				 UseShadow = False;
			     break;

	case(IDM_HYDROGEN):  mask = NormAtomFlag;
			     Hydrogens = !Hydrogens;
			     ReDrawFlag |= RFRefresh;

			     if( Hydrogens )
			     {      SelectZone(mask|HydrogenFlag);
			     } else RestrictZone(mask);
			     break;
	
	case(IDM_SPECULAR):  FakeSpecular = !FakeSpecular;
			     ReDrawFlag |= RFColour;
			     break;
	
	case(IDM_SHADOW):    ReDrawFlag |= RFRefresh;
			     UseShadow = !UseShadow;
			     if( UseShadow )
			     {   ReviseInvMatrix();
				 VoxelsClean = False;
				 UseSlabPlane = False;
				 ReAllocBuffers();
			     }
			     break;

	case(IDM_STEREO):    /* Stereo */
                             if( UseStereo )
                             {   SetStereoMode(False);
                             } else SetStereoMode(True);
                             ReDrawFlag |= RFRefresh;
			     break;

        case(IDM_LABELS):    /* Labels */
                             LabelOptFlag = !LabelOptFlag;
                             DefaultLabels(LabelOptFlag);
                             ReDrawFlag |= RFRefresh;
                             break;

	/* Save Menu */
	case(IDM_GIF):

	case(IDM_EPSF):    SaveOutputFile( option ); 
			   break;
	
	default:  return(FALSE);
    }
    return(TRUE);
}    


static void InitiateServer( hWinCli, lParam )
    HWND hWinCli;  LONG lParam;
{
    HWND hWinServ;
    ATOM aTopicIn, aTopicOut;
    ATOM aApplIn, aApplOut;
    
    char TopicName[16];
    char ApplName[16];
    register int i;

    if( ConvCount == MaxConvNum )
	return;
	    
    if( aApplIn = LOWORD(lParam) )
    {   GlobalGetAtomName(aApplIn,ApplName,14);
	if( stricmp(ApplName,"RasWin") )
	    return;
    } else return;
    
    if( aTopicIn = HIWORD(lParam) )
    {   GlobalGetAtomName(aTopicIn,TopicName,14);
	/* Test for Valid Topic */
	/* if( _stricmp(Topic,"System") &&
	 *     _stricmp(Topic,"RemoteControl") )
	 * return;
	 */
    } else *TopicName = '\0';
    
   
    hWinServ = CreateWindow("RasDDEClass","RasWinDDE",
			    WS_CHILD, 0, 0, 0, 0,
			    CanvWin, NULL, hInstance, NULL );
    if( !hWinServ ) return;
	 
    for( i=0; i<MaxConvNum; i++ )
	if( !ConvData[i].server )
	    break;
	    
    ConvData[i].server = hWinServ;
    ConvData[i].client = hWinCli;
    ConvData[i].closed = False;
    ConvCount++;       
	  
	 
    /* Main DDE Server */       
    aTopicOut = (ATOM)NULL;
    
    aApplOut = GlobalAddAtom("RasWin");
    SendMessage( hWinCli, WM_DDE_ACK, (WPARAM)hWinServ,
		 MAKELONG(aApplOut,aTopicOut) ); 
}


static void ClampDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;

    if( temp > 1.0 )
    {   DialValue[dial] = 1.0;
    } else if( temp < -1.0 )
    {   DialValue[dial] = -1.0;
    } else DialValue[dial] = temp;
}


static void WrapDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;
    while( temp < -1.0 )  temp += 2.0;
    while( temp > 1.0 )   temp -= 2.0;
    DialValue[dial] = temp;
}


static void MouseMove( status, dx, dy )
    int status, dx, dy;
{
    if( MouseMode == MMRasMol )
    {   if( status & MK_SHIFT )
	{   if( status & MK_LBUTTON )
	    {   if( dy ) /* Zoom Vertical */
		{   ClampDial( 3, (Real)dy/HRange );
		    ReDrawFlag |= RFZoom;
		}
	    } else if( status & (MK_MBUTTON|MK_RBUTTON) )
		if( dx ) /* Z Rotation Horizontal */
		{   WrapDial( 2, (Real)-dx/WRange );
		    ReDrawFlag |= RFRotateZ;
		}
	} else if( status & MK_CONTROL )
	{   if( status & MK_LBUTTON )
	    {   if( dy ) /* Slab Vertical */
		{   ClampDial( 7, (Real)dy/YRange );
		    ReDrawFlag |= RFSlab;
		}
	    }

	} else /* Unmodified! */
	    if( status & MK_LBUTTON )
	    {   if( dx ) /* Rotate Y Horizontal */
		{   WrapDial( 1, (Real)dx/WRange );
		    ReDrawFlag |= RFRotateY;
		}

		if( dy ) /* Rotate X Vertical */
		{   WrapDial( 0, (Real)-dy/HRange );
		    ReDrawFlag |= RFRotateX;
		}
		UpdateScrollBars();
	    } else if( status & (MK_MBUTTON|MK_RBUTTON) )
	    {   if( dx ) /* Translate X Horizontal */
		{   ClampDial( 4, (Real)dx/XRange );
		    ReDrawFlag |= RFTransX;
		}

		if( dy ) /* Translate Y Vertical */
		{   ClampDial( 5, (Real)-dy/YRange );
		    ReDrawFlag |= RFTransY;
		}
	    }
    } else if( MouseMode == MMQuanta )
    {   if( status & MK_SHIFT )
	{   if( status & MK_LBUTTON )
	    {   if( dy ) /* Slab Vertical */
		{   ClampDial( 7, (Real)dy/YRange );
		    ReDrawFlag |= RFSlab;
		}
	    } else if( status & (MK_MBUTTON|MK_RBUTTON) )
	    {   if( dx ) /* Translate X Horizontal */
		{   ClampDial( 4, (Real)dx/XRange );
		    ReDrawFlag |= RFTransX;
		}

		if( dy ) /* Translate Y Vertical */
		{   ClampDial( 5, (Real)-dy/YRange );
		    ReDrawFlag |= RFTransY;
		}
	    } else /* No Mouse Buttons */
		if( dy ) /* Zoom Vertical */
		{   ClampDial( 3, (Real)dy/HRange );
		    ReDrawFlag |= RFZoom;
		}
	} else if( status & (MK_LBUTTON|MK_MBUTTON) )
	{   if( dx ) /* Rotate Y Horizontal */
	    {   WrapDial( 1, (Real)dx/WRange );
		ReDrawFlag |= RFRotateY;
	    }

	    if( dy ) /* Rotate X Vertical */
	    {   WrapDial( 0, (Real)-dy/HRange );
		ReDrawFlag |= RFRotateX;
	    }
	    UpdateScrollBars();
	} else if( status & MK_RBUTTON )
	    if( dx ) /* Z Rotation Horizontal */
	    {   WrapDial( 2, (Real)-dx/WRange );
		ReDrawFlag |= RFRotateZ;
	    }
	
    } else /* MMInsight */
	switch( status & (MK_LBUTTON|MK_MBUTTON|MK_RBUTTON) )
	{   case( MK_LBUTTON ):
		    if( dx ) /* Rotate Y Horizontal */
		    {   WrapDial( 1, (Real)dx/WRange );
			ReDrawFlag |= RFRotateY;
		    }

		    if( dy ) /* Rotate X Vertical */
		    {   WrapDial( 0, (Real)dy/HRange );
			ReDrawFlag |= RFRotateX;
		    }
		    break;

	    case( MK_MBUTTON ):
		    if( dx ) /* Translate X Horizontal */
		    {   ClampDial( 4, (Real)dx/XRange );
			ReDrawFlag |= RFTransX;
		    }

		    if( dy ) /* Translate Y Vertical */
		    {   ClampDial( 5, (Real)dy/YRange );
			ReDrawFlag |= RFTransY;
		    }
		    break;

	    case( MK_LBUTTON|MK_MBUTTON ):
		    ClampDial( 3, (Real)dx/WRange - (Real)dy/HRange );
		    ReDrawFlag |= RFZoom;
		    break;

	    case( MK_LBUTTON|MK_RBUTTON ):
		    WrapDial( 2, (Real)dx/WRange - (Real)dy/HRange );
		    ReDrawFlag |= RFRotateZ;
		    break;

	    case( MK_LBUTTON|MK_MBUTTON|MK_RBUTTON ):
		    ClampDial( 7, (Real)dx/XRange - (Real)dy/YRange );
		    ReDrawFlag |= RFSlab;
		    break;
	}
}


long FAR PASCAL MainCallB(hWin,uMsg,wArg,lArg)
    HWND hWin; UINT uMsg; WPARAM wArg; LPARAM lArg;
{
    register int pos,status;
    register int dx, dy;

    register COLORREF BackColRef;    
    register HPALETTE hCMap;
    register HANDLE hand;
    register HMENU hMenu;
    register HDC hMemDC;
    register HDC hDC;

    PAINTSTRUCT ps;
    RECT rc;
    

    CanvWin = hWin;
    
    switch(uMsg)
    {   case(WM_DROPFILES):   /* Disable Drag & Drop */
                              if( IsPaused ) break;

                              ZapDatabase();
			      *fnamebuf = '\0';
			      DragQueryFile((HDROP)wArg,0,fnamebuf,127);
			      FetchFile(FormatXYZ,False,fnamebuf);
                              DefaultRepresentation();
			      DragFinish((HDROP)wArg);
			      break;

	case(WM_DESTROY):     /* Destroy RasWin */
			      DragAcceptFiles(CanvWin,FALSE);
			      PostQuitMessage(0);
			      break;

	case(WM_CLOSE):       DestroyWindow(CanvWin);
			      DestroyWindow(CmndWin);
			      break;

	case(WM_ACTIVATE):    if( !wArg ) break;
	case(WM_QUERYNEWPALETTE):
			      if( ColourMap )
			      {   hDC = GetDC(hWin);
				  hCMap = SelectPalette(hDC,ColourMap,False);
				  status = RealizePalette(hDC);
				  if( hCMap ) SelectPalette(hDC,hCMap,False);
				  ReleaseDC(hWin,hDC);
				  
				  if( status )
				  {   InvalidateRect(hWin,NULL,True);
				      return True;
				  }
			      }
			      return(0L);
			      
	case(WM_PALETTECHANGED):
			      if( ColourMap && ((HWND)wArg != hWin) )
			      {   hDC = GetDC(hWin);
				  hCMap = SelectPalette(hDC,ColourMap,False);
				  if( RealizePalette(hDC) )
				      InvalidateRect(hWin,NULL,True);
				  if( hCMap ) SelectPalette(hDC,hCMap,False);
				  ReleaseDC(hWin,hDC);
			      }
			      return(0L);
			      
			     
	case(WM_INITMENUPOPUP):  /* Initialise Checks */
			      if( lArg == 4 )
			      {   /* Options Menu */
				  hMenu = (HMENU)wArg;

				  status = UseSlabPlane ? MF_CHECKED 
							: MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_SLAB,status);

				  status = Hydrogens ? MF_CHECKED 
						     : MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_HYDROGEN,status);

				  status = FakeSpecular ? MF_CHECKED 
							: MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_SPECULAR,status);

				  status = UseShadow ? MF_CHECKED 
						     : MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_SHADOW,status);

				  status = UseStereo ? MF_CHECKED
						     : MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_STEREO,status);

				  status = LabelOptFlag ? MF_CHECKED 
						        : MF_UNCHECKED;
				  CheckMenuItem(hMenu,IDM_LABELS,status);
			      }
			      return( 0L );                      

	case(WM_SIZE):        if( wArg != SIZE_MINIMIZED )
			      {   GetClientRect(hWin,&rc);
				  YRange = rc.bottom;
				  XRange = rc.right;
			      
				  /* Ensure Long Aligned */
				  if( dx = XRange%4 )
				      XRange += 4-dx;
			      
				  Range = MinFun(XRange,YRange);
				  ReDrawFlag |= RFReSize;
				  HRange = YRange>>1;
				  WRange = XRange>>1;
				  ClearImage();
			      }
			      break;

	case(WM_HSCROLL):     /* Horizontal Scroll */
			      pos = GetScrollPos(hWin,SB_HORZ);
			      switch( wArg )
			      {   case(SB_LINEDOWN):  pos += 5;   break;
				  case(SB_PAGEDOWN):  pos += 10;  break;
				  case(SB_PAGEUP):    pos -= 10;  break;
				  case(SB_LINEUP):    pos -= 5;   break;
				  default:            return(0L);

				  case(SB_THUMBTRACK):
				  case(SB_THUMBPOSITION):
					     pos = LOWORD(lArg);
					     break;
			      }
			      
			      if( pos>100 ) 
			      {   pos -= 100;
			      } else if( pos<0 ) 
				  pos += 100; 
			     
			      SetScrollPos(hWin,SB_HORZ,pos,TRUE);
			      DialValue[1] = (pos/50.0)-1.0;
			      ReDrawFlag |= RFRotateY;
			      break;                      

	case(WM_VSCROLL):     /* Vertical Scroll */
			      pos = GetScrollPos(hWin,SB_VERT);
			      switch( wArg )
			      {   case(SB_LINEDOWN):  pos += 5;   break;
				  case(SB_PAGEDOWN):  pos += 10;  break;
				  case(SB_PAGEUP):    pos -= 10;  break;
				  case(SB_LINEUP):    pos -= 5;   break;
				  default:            return(0L);

				  case(SB_THUMBTRACK):
				  case(SB_THUMBPOSITION):
					     pos = LOWORD(lArg);
					     break;
			      }
			      
			      if( pos>100 ) 
			      {   pos -= 100;
			      } else if( pos<0 ) 
				  pos += 100; 
			     
			      SetScrollPos(hWin,SB_VERT,pos,TRUE);
			      DialValue[0] = 1.0-(pos/50.0);
			      ReDrawFlag |= RFRotateX;
			      break;                      

	case(WM_LBUTTONDOWN): 
	case(WM_MBUTTONDOWN):
	case(WM_RBUTTONDOWN): InitX = PointX = LOWORD(lArg);
			      InitY = PointY = HIWORD(lArg);
			      HeldButton = True;
			      SetCapture(hWin);
			      break;
	
	case(WM_LBUTTONUP):
	case(WM_MBUTTONUP):
	case(WM_RBUTTONUP): /* Mouse Buttons */
			      HeldButton = False;
			      ReleaseCapture();

			      if( Database )
			      {   PointX = dx = LOWORD(lArg);
				  PointY = dy = HIWORD(lArg);
				  if( IsClose(dx,InitX) && IsClose(dy,InitY) )
				  {   if( wArg & (MK_SHIFT|MK_CONTROL) )
                                      {      PickAtom(True,dx,YRange-dy);
                                      } else PickAtom(False,dx,YRange-dy);
				      AdviseUpdate(AdvPickNumber);
				      AdviseUpdate(AdvPickAtom);
				  }
			      }
			      break;


	case(WM_MOUSEMOVE):   /* Mouse Movement */
			      if( !HeldButton )
			      {   if( (MouseMode==MMQuanta) &&
				      (wArg & MK_SHIFT) )
				  {   InitX = PointX = LOWORD(lArg);
				      InitY = PointY = HIWORD(lArg);
				      HeldButton = True;
				      SetCapture(hWin);
				  }
				  break;
			      }

			      if( IsClose((int)LOWORD(lArg),InitX) &&
				  IsClose((int)HIWORD(lArg),InitY) )
				  break;

			      if( wArg & (MK_LBUTTON|MK_MBUTTON|
					  MK_RBUTTON|MK_SHIFT) )
			      {   dx = (int)LOWORD(lArg)-PointX;
				  dy = (int)HIWORD(lArg)-PointY;
				  MouseMove( wArg, dx, dy );
				  PointX = LOWORD(lArg);
				  PointY = HIWORD(lArg);
			      } else  /* No Buttons! */
			      {   HeldButton = False;
				  ReleaseCapture();
			      }
			      break;

				      
	case(WM_SETFOCUS):    /* Obtain Window Focus */ 
	case(WM_KILLFOCUS):   /* Release Window Focus */
			      SendMessage(CmndWin,uMsg,wArg,lArg);     
			      return( 0L );
	
	case(WM_PAINT):       hDC = BeginPaint(hWin,&ps);
			      SetBkMode(hDC,TRANSPARENT);
			      if( PixMap )
			      {   hCMap = SelectPalette(hDC,ColourMap,False);
				  RealizePalette(hDC);
#ifdef _WIN32
				  SetWindowOrgEx(hDC,0,0,NULL);
#else
				  SetWindowOrg(hDC,0,0);
#endif
				  hMemDC = CreateCompatibleDC(hDC);
				  SelectObject(hMemDC,PixMap);
				  BitBlt(hDC,0,0,XRange,YRange,
					 hMemDC,0,0,SRCCOPY);
					 
				  SelectPalette(hDC,hCMap,False);      
				  DeleteDC(hMemDC);
			      } else /* Erase Update Region */
			      {    if( ColourMap )
				   {   hCMap=SelectPalette(hDC,ColourMap,0);
				       RealizePalette(hDC);
				   }
				   BackColRef = RGB(BackR,BackG,BackB);
				   hand = CreateSolidBrush(BackColRef);
				   GetUpdateRect(hWin,&rc,False);
				   FillRect( hDC, &rc, hand );
				   if( ColourMap && hCMap )
				       SelectPalette(hDC,hCMap,False);
				   DeleteObject(hand);
			      }
			      EndPaint(hWin,&ps);
			      if( !RasWinReady )
			      {   RasWinReady = True;
				  AdviseUpdate(-4);
			      }
			      return( 0L );
	
	case(WM_SYSCHAR):     if( lArg & (1L<<29) )  /* ALT-key pressed? */
				  return(DefWindowProc(hWin,uMsg,wArg,lArg));
	case(WM_CHAR):        if( ProcessCharacter(LOBYTE(wArg)) )
				  if( ExecuteCommand() )
				      RasMolExit();
			      break;

	case(WM_SYSKEYDOWN):
	case(WM_KEYDOWN):     switch(LOBYTE(wArg))
			      {   case(0x23): ProcessCharacter(0x05); break;
				  case(0x24): ProcessCharacter(0x01); break;
				  case(0x25): ProcessCharacter(0x02); break;
				  case(0x26): ProcessCharacter(0x10); break;
				  case(0x27): ProcessCharacter(0x06); break;
				  case(0x28): ProcessCharacter(0x0e); break;
				  case(0x2e): ProcessCharacter(0x04); break;
				  
				  default:
				  return(DefWindowProc(hWin,uMsg,wArg,lArg));
			      }
			      break;

	case(WM_RENDERALLFORMATS):
			      OpenClipboard(hWin);
			      SendMessage(hWin,WM_RENDERFORMAT,CF_DIB,0L);
			      SendMessage(hWin,WM_RENDERFORMAT,CF_PALETTE,0L);
			      CloseClipboard();
			      return( 0L );
			      
	case(WM_RENDERFORMAT):
			      if( hand = RenderClipboard(wArg) )
				  SetClipboardData(wArg,hand);
			      return( 0L );
			      

	case(WM_DDE_INITIATE): /* DDE Server Connection */
			      InitiateServer((HWND)wArg,lArg);
			      return( 0L );
					      
	case(WM_COMMAND):     if( !IsPaused && HandleMenu(wArg) )
				  break;
			      
	default:              return( DefWindowProc(hWin,uMsg,wArg,lArg) );

    }

	
    if( ReDrawFlag )
	RefreshScreen();
    if( !CommandActive )
	ResetCommandLine(0);
	
    return(0L);
}


static int InitialiseApplication()
{
    WNDCLASS wc;
    
    wc.hIcon = LoadIcon(hInstance,"RasWinIcon");
    wc.hInstance = hInstance;
    wc.cbWndExtra = 0;
    wc.cbClsExtra= 0;

    /* Canvas Window Class */
    wc.style = 0;
    wc.lpfnWndProc = MainCallB;
    wc.hbrBackground = CreateSolidBrush(RGB(0,0,0));
    wc.hCursor = LoadCursor(hInstance,"RasWinCursor");
    wc.lpszClassName = "RasWinClass";
    wc.lpszMenuName = NULL;

    if( !RegisterClass(&wc) )
	return( False );

    /* Terminal Window Class */
    wc.style = CS_NOCLOSE;
    wc.lpfnWndProc = CmndCallB;
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wc.hCursor = LoadCursor(NULL,IDC_ARROW);
    wc.lpszClassName = "RasCliClass";
    wc.lpszMenuName = NULL;

    if( !RegisterClass(&wc) )
	return( False );

    /* DDE Server Window Class */
    wc.lpfnWndProc = DDECallB;
    wc.lpszClassName = "RasDDEClass";
    wc.hbrBackground = NULL;
    wc.hCursor = NULL;
    wc.hIcon = NULL;

    return( RegisterClass(&wc) );
}


static char *RegisterFormat( buffer, desc, ext )
    char *buffer, *desc, *ext;
{
    while( *buffer++ = *desc++ );
    while( *buffer++ = *ext++ );
    return( buffer );
}


static void InitDefaultValues()
{
    register int i;

    Interactive = True;

    ConvCount = 0;
    AdviseCount = 0;
    for( i=0; i<MaxConvNum; i++ )
        ConvData[i].server = NULL;
    for( i=0; i<MaxAdviseNum; i++ )
        AdviseData[i].server = NULL;
    LabelOptFlag = False;
    RasWinReady = True;
    HeldButton = False;

    fnamebuf[0] = '\0';
    snamebuf[0] = '\0';

    FileFormat = FormatXYZ;
    CalcBondsFlag = True;
}


#define FORMATOPTMAX   1
static struct {
        char *ident;
        int format, len;
    } FormatOpt[FORMATOPTMAX] = {
            { "xyz",        FormatXYZ,       3 }
                                };
   

static int ProcessOptions( ptr )
    char __far *ptr;
{
    register char *dst;
    register int i;

    while( *ptr )
    {   if( (*ptr==' ') || (*ptr=='=') )
	{   ptr++;
	} else if( (*ptr=='/') || (*ptr=='-') )
	{   ptr++;
            for( i=0; i<FORMATOPTMAX; i++ )
	        if( !_fstrnicmp(ptr,FormatOpt[i].ident,FormatOpt[i].len) )
                    break;

            if( i < FORMATOPTMAX )
	    {   FileFormat = FormatOpt[i].format;
                ptr += FormatOpt[i].len;
	    } else if( !_fstrnicmp(ptr,"sybyl",5) )
	    {   FileFormat = FormatMol2;     ptr += 4;

	    } else if( !_fstrnicmp(ptr,"script",6) )
	    {   ptr += 6;
		while( *ptr && (*ptr==' ') )
		    ptr++;

		if( *ptr )
		{   dst = snamebuf;
		    while( *ptr && (*ptr!=' ') )
			*dst++ = *ptr++;
		    *dst = '\0';
		} else return( False );
	    } else return( False );

	} else if( !*fnamebuf )
	{   dst = fnamebuf;
	    while( *ptr && (*ptr!=' ') )
		*dst++ = *ptr++;
	    *dst = '\0';
	} else return( False );
    }
    return( True );
}


int PASCAL WinMain(hCurrent,hPrevious,lpCmdLine,nCmdShow)
    HANDLE hCurrent,hPrevious;
    LPSTR lpCmdLine;
    int nCmdShow;
{
    register char *dst;
    register FILE *fp;
    MSG event;

    hInstance = hCurrent;
    if( !hPrevious && !InitialiseApplication() )
	return(False);

    dst = ifilters;
    dst = RegisterFormat(dst,"TINKER XYZ Format","*.XYZ","*.xyz");
    *dst = '\0';

    /* Load File Common Dialog Box */
    ofn1.lStructSize=sizeof(OPENFILENAME);
    ofn1.Flags = OFN_NOCHANGEDIR | OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
    ofn1.lpstrFilter = ifilters;
    ofn1.lpstrTitle = "Select Molecular Coordinate File";
    ofn1.lpstrFile = fnamebuf;
    ofn1.nMaxFile = 128;

    ofn1.lpstrCustomFilter = NULL;
    ofn1.lpstrInitialDir = NULL;
    ofn1.lpstrFileTitle = NULL;
    ofn1.hwndOwner = NULL;
    
    dst = ofilters;
    dst = RegisterFormat(dst,"GIF Image","*.GIF");
    dst = RegisterFormat(dst,"Colour PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Mono PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Vector PostScript","*.PS;*.EPS");
    dst = RegisterFormat(dst,"Apple Macintosh PICT","*.PIC");
    dst = RegisterFormat(dst,"Silicon Graphics RGB","*.RGB");
    *dst = '\0';
    
    /* Save File Common Dialog Box */
    ofn2.lStructSize=sizeof(OPENFILENAME);
    ofn2.Flags = OFN_NOCHANGEDIR | OFN_HIDEREADONLY | OFN_NOREADONLYRETURN
	       | OFN_OVERWRITEPROMPT | OFN_PATHMUSTEXIST;
    ofn2.lpstrFilter = ofilters;
    ofn2.lpstrTitle = "Select Graphics Ouptut File";
    ofn2.lpstrFile = fnamebuf;
    ofn2.nMaxFile = 128;

    ofn2.lpstrCustomFilter = NULL;
    ofn2.lpstrInitialDir = NULL;
    ofn2.lpstrFileTitle = NULL;
    ofn2.hwndOwner = NULL;

    if( !GetProfileString("extensions","xyz","",fnamebuf,128) )
	WriteProfileString("extensions","xyz","raswin.exe ^.pdb");
    DDETimeOut = GetPrivateProfileInt("RasWin","DDETimeOut",
				      DefaultTimeOut,"RASWIN.INI"); 

    InitDefaultValues();
    ProcessOptions(lpCmdLine);

    /* Avoid Windows NT problems! */
    ReDrawFlag = RFInitial;
    CommandActive = True;

    if( !InitTerminal(hInstance) ||
	!OpenDisplay(hInstance,nCmdShow) )
       return(False);

    WriteString("TINKER RasMol Molecule Viewer\n");
    WriteString("Version 4.0, October 2002\n");
    WriteString("Based on RasMol 2.6 by Roger Sayle\n\n");
	    
    InitialiseCommand(); 
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();
    InitialiseRepres();
   
    if( *fnamebuf && FetchFile(FileFormat,True,fnamebuf) )
        DefaultRepresentation();

    ResetCommandLine(1);
    LoadInitFile();

    if( *snamebuf )
    {   if( !(fp=fopen(snamebuf,"r")) )
	{   if( CommandActive )
		WriteChar('\n');
	    WriteString("Error: File '");
	    WriteString(snamebuf);
	    WriteString("' not found!\n");
	    CommandActive = False;
	} else LoadScriptFile(fp,snamebuf);
    }

    RefreshScreen();
    if( !CommandActive )
	ResetCommandLine(0);

    DragAcceptFiles(CanvWin,TRUE);
    while( GetMessage(&event,NULL,0,0) )
    {   TranslateMessage(&event);
	DispatchMessage(&event);
    }
    DeleteObject(TermFont);
    CloseDDELinks();
    CloseDisplay();

    return( event.wParam );
}
