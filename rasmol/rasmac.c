/* rasmac.c
 * RasMol Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#include <string.h>
#include <stdio.h>

#include <Errors.h>
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif

#include <StandardFile.h>
#include <AppleEvents.h>
#include <GestaltEqu.h>
#include <ToolUtils.h>
#include <Resources.h>
#include <OSEvents.h>
#include <DiskInit.h>
#include <OSUtils.h>
#include <SegLoad.h>
#include <Events.h>
#include <Memory.h>
#include <Finder.h>
#include <Files.h>
#include <Types.h>
#include <Scrap.h>
#include <Desk.h>
#include <EPPC.h>

#include <Quickdraw.h>
#include <Controls.h>
#include <Palettes.h>
#include <TextEdit.h>
#include <Windows.h>
#include <Dialogs.h>
#include <Menus.h>
#include <Fonts.h>

#define RASMOL
#include "rasmol.h"
#include "molecule.h"
#include "abstree.h"
#include "graphics.h"
#include "pixutils.h"
#include "transfor.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "outfile.h"


/* Determine Mouse Sensitivity! */
#define IsClose(u,v)  (((u)>=(v)-1) && ((u)<=(v)+1))

#define CmndSize   (CmndRows*CmndCols)
#define ScrlMax    80
#define CmndRows   160
#define CmndCols   80

/* Terminal Emulation Variables */
static char TermScreen[CmndSize];

static int CharSkip;
static int CmndStart;
static int ScrlStart;
static int TermCursor;
static int CharWide,CharHigh;
static int TermXPos,TermYPos;
static int TermRows,TermCols;
static Rect CmndRect;

static ControlHandle CmndScroll;
static short RasMolResFile;

static int DialogFormat;
static char Filename[256];
static Point MousePrev;
static int PointX,PointY;
static int InitX,InitY;
static int LabelOptFlag;
static int HeldButton;


#ifdef __CONDITIONALMACROS__
#define ArrowCursor   SetCursor(&qd.arrow)
#else
#define ArrowCursor   SetCursor(&arrow)
#endif


/* Forwardly Compatable Definitions */
#ifndef gestaltNativeCPUtype
#define gestaltNativeCPUtype 'cput'
#define gestaltCPU68000    0x000
#define gestaltCPU68010    0x001
#define gestaltCPU68020    0x002
#define gestaltCPU68030    0x003
#define gestaltCPU68040    0x004
#define gestaltCPU601      0x101
#define gestaltCPU603      0x103
#define gestaltCPU604      0x104
#endif

#ifndef __CONDITIONALMACROS__
/* Undocumented System Functions */
pascal OSErr SetDialogDefaultItem( DialogPtr theDialog,
        short newItem ) = {0x303C,0x0304,0xAA68};
#endif


#ifdef __CONDITIONALMACROS__
/* Global RoutineDescriptors */
AEEventHandlerUPP HandleAEIgnorePtr;
AEEventHandlerUPP HandleAEOpenDocPtr;
AEEventHandlerUPP HandleAEQuitAppPtr;
ControlActionUPP CanvScrollProcPtr;
ControlActionUPP CmndScrollProcPtr;
DlgHookYDUPP OpenDlgHookPtr;
DlgHookYDUPP SaveDlgHookPtr;
#endif


/* Function Prototypes */
static void PaintScreen();


void RasMolExit()
{
    /* Restore System Event Mask */
    SetEventMask( everyEvent & ~keyUp );
    
    /* Free System Memory Resources */
    if( FBuffer ) _ffree( FBuffer );
    if( DBuffer ) _ffree( DBuffer );
    PurgeDatabase();
    CloseDisplay();
    ExitToShell();
}


void RasMolFatalExit( msg )
    char *msg;
{
    register char *ptr;
    register int len;
    
    Str255 buffer;
    DialogPtr dlg;
    Handle hand;
    short item;
    Rect rect;
    
    dlg = GetNewDialog(171,(Ptr)0,(WindowPtr)-1);
    SetDialogDefaultItem(dlg,1);
    
    /* Error Message Text! */
    ptr = (char*)&buffer[1];
    while( *msg ) 
        *ptr++ = *msg++;
    
    len = ptr - (char*)&buffer[1];
    buffer[0] = (unsigned char)len;
    GetDItem(dlg,4,&item,&hand,&rect);
    SetIText(hand,buffer);
    ShowWindow(dlg);

    do {
#ifdef __CONDITIONALMACROS__
        ModalDialog((ModalFilterUPP)0,&item);
#else
        ModalDialog((ProcPtr)0,&item);
#endif
    } while( item != 1 );
    DisposDialog(dlg);

    /* Restore System Event Mask */
    SetEventMask( everyEvent & ~keyUp );
    
    /* Free System Memory Resources */
    if( FBuffer ) _ffree( FBuffer );
    if( DBuffer ) _ffree( DBuffer );
    PurgeDatabase();
    CloseDisplay();
    ExitToShell();
}


/* Externally Visible Functions! */
static void CopyResource( type, file, id )
    ResType type;  short file;  short id;
{
    register Handle hand;
    Str255 name;

    UseResFile(RasMolResFile);
    hand = GetResource(type,id);
    if( hand )
    {   GetResInfo(hand,&id,&type,name);
        DetachResource(hand);

        UseResFile(file);
        AddResource(hand,type,kCustomIconResource,name);
        if( !ResError() ) WriteResource(hand);
        ReleaseResource(hand);
        UseResFile(RasMolResFile);
    }
}

void SetFileInfo( ptr, appl, type, icon )
    char *ptr;  OSType appl, type;  short icon;
{
    register char *dst;
    register short file;
    register int len;
    
    Str255 buffer;
    FSSpec fss;
    FInfo info;
    
    /* Create Pascal Filename */
    dst = (char*)&buffer[1];
    while( *ptr ) *dst++ = *ptr++;
    len = dst - (char*)&buffer[1];
    buffer[0] = (unsigned char)len;
    
    FSMakeFSSpec(0,0,buffer,&fss);
    
    if( icon )
    {   /* smSystemScript not always defined! */
        FSpCreateResFile(&fss,appl,type,0);
        file = FSpOpenResFile(&fss,fsRdWrPerm);
        CopyResource('icl4',file,icon);
        CopyResource('icl8',file,icon);
        CopyResource('ICN#',file,icon);
        CopyResource('ics#',file,icon);
        CopyResource('ics4',file,icon);
        CopyResource('ics8',file,icon);
        CloseResFile(file);
    }
    
    if( !FSpGetFInfo(&fss,&info) )
    {   info.fdType = type;
        info.fdCreator = appl;
        if( icon )
            info.fdFlags |= kHasCustomIcon;
        FSpSetFInfo(&fss,&info);
    }
    FlushVol(NULL,fss.vRefNum);
}


static void SetTermScroll( pos )
    int pos;
{
    GrafPtr savePort;
    
    GetPort(&savePort);
    SetPort(CmndWin);
    SetCtlValue(CmndScroll,pos);
    ScrlStart = ScrlMax - pos;
    InvalRect(&CmndWin->portRect);
    SetPort(savePort);
}


#define ClipCaret  (CaretX<CmndWin->portRect.right-15)

static long LastCaretTime;
static int CaretX,CaretY;
static int CaretFlag;


static void ShowCaret()
{
    GrafPtr savePort;
    
    if( !CaretFlag && ClipCaret )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        PenMode(patCopy);
        MoveTo(CaretX,CaretY);
        Line(0,CharSkip);
        SetPort(savePort);
    }
    LastCaretTime = TickCount();
    CaretFlag = True;
}

static void HideCaret()
{
    GrafPtr savePort;
    
    if( CaretFlag && ClipCaret )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        PenMode(patBic);
        MoveTo(CaretX,CaretY);
        Line(0,CharSkip);
        SetPort(savePort);
    }
    CaretFlag = False;
}

static void SetCaretPos( x, y )
    int x, y;
{
    if( CaretFlag )
    {   HideCaret();
        CaretX = x;
        CaretY = y;
        ShowCaret();
    } else
    {   CaretX = x;
        CaretY = y;
    }
}

static void HandleCaret()
{
    register long ticks;
    
    ticks = TickCount();
    if(  ticks > LastCaretTime + GetCaretTime() )
    {   if( CaretFlag )
       {   LastCaretTime = ticks;
           HideCaret();
       } else ShowCaret();
    }
}


void WriteChar( ch )
    int ch;
{
    register int i;
    
    GrafPtr savePort;
    
    /* Scroll to bottom! */
    if( ScrlStart )
        SetTermScroll( ScrlMax );

    switch( ch )
    {   case(0x07):  SysBeep(15);
                     break;
                     
        case(0x08):  if( TermXPos>0 )
                     {   TermXPos--;
                         if( TermCursor )
                             SetCaretPos(TermXPos*CharWide,
                                         TermYPos*CharSkip);
                     }
                     break;
                     
#ifdef ORIG
        case(0x0D):  if( TermXPos )
                     {   if( TermCursor )
                             SetCaretPos(0,TermYPos*CharSkip);
                         TermXPos = 0;
                     }
                     break;
#else
        case(0x0D):
#endif
                     
        case(0x0A):  if( TermYPos==TermRows-1 )
                     {   CmndStart++;
                         if( CmndStart == CmndRows )
                             CmndStart = 0;
                             
                         i = TermYPos + CmndStart;
                         if( i >= CmndRows ) i -= CmndRows;
                         memset(TermScreen+i*CmndCols,32,CmndCols);
                         
                         GetPort(&savePort);
                         SetPort(CmndWin);
                         /* InvalRect(&CmndWin->portRect); */
                         PaintScreen();
                         SetPort(savePort);
                     } else TermYPos++;
                     TermXPos = 0;
                     
                     if( TermCursor )
                         SetCaretPos(0,TermYPos*CharSkip);
                     break;
                     
        default:     i = TermYPos + CmndStart;
                     if( i >= CmndRows ) i -= CmndRows;
                     TermScreen[i*CmndCols+TermXPos] = ch;
                     if( TermXPos < TermCols )
                     {   GetPort(&savePort);
                         SetPort(CmndWin);
                         MoveTo(TermXPos*CharWide+1,TermYPos*CharSkip+CharHigh);
                         DrawChar(ch);
                         SetPort(savePort);
                     }
                     
                     if( TermXPos==CmndCols-1 )
                     {   if( TermYPos==TermRows-1 )
                         {  CmndStart++;
                            if( CmndStart == CmndRows )
                                CmndStart = 0;
                             
                            i = TermYPos + CmndStart;
                            if( i >= CmndRows ) i-= CmndRows;
                            memset(TermScreen+i*CmndCols,32,CmndCols);
                            
                            GetPort(&savePort);
                            SetPort(CmndWin);
                            /* InvalRect(&CmndWin->portRect); */
                            PaintScreen();
                            SetPort(savePort);
                         } else TermYPos++;
                         TermXPos = 0;
                     } else TermXPos++;
                     
                     if( TermCursor )
                         SetCaretPos(TermXPos*CharWide,
                                     TermYPos*CharSkip);
                     break;
    }
}


void WriteString( ptr )
    char *ptr;
{
    while( *ptr )
        WriteChar( *ptr++ );
}


#ifndef topLeft
#define topLeft(r) (*((Point*)(&(r))))
#endif


static void InitTerminal()
{
    register WindowPeek wpeek;
    register WStateData *wsdp;
    register Region *rgn;
    register int mb,tb;
    register int x, y;
    FontInfo finfo;
    Point pnt;
    Rect rect;
    
    TermCursor = False;
    CaretFlag = False;
    
    /* Default Window Size */
    TermCols = 80;  TermRows = 24;
    ScrlStart = CmndStart = 0;
    TermXPos = TermYPos = 0;

    CmndWin = GetNewWindow(151,0,CanvWin);
    SetPort(CmndWin);
    
    /* Font Style */
    TextMode(srcCopy);
    TextFont(monaco);
    TextSize(9);
    TextFace(0);
    
    GetFontInfo(&finfo);
    CharSkip = finfo.ascent + finfo.descent + finfo.leading;
    CharHigh = finfo.ascent;
    CharWide = finfo.widMax;
    
    /* Set Initial Terminal Size */
    x = TermCols*CharWide;
    y = TermRows*CharSkip;
    SizeWindow(CmndWin,x+16,y,true);
    ShowWindow(CmndWin);

    /* Create Scroll Bar */
    rect.left = x+1;   rect.right = x+17;
    rect.top = -1;     rect.bottom = y-14;
    CmndScroll = NewControl(CmndWin,&rect,"\p",false,ScrlMax,
                            0,ScrlMax,scrollBarProc,0L);
    
    wpeek = (WindowPeek)CmndWin;
    pnt = topLeft(CmndWin->portRect);
    rgn = *(wpeek->strucRgn);
    LocalToGlobal(&pnt);

    /* Title & Menu Bar Heights */
    tb = pnt.v - rgn->rgnBBox.top;
    mb = GetMBarHeight()+tb+4;

#ifdef __CONDITIONALMACROS__    
    y = (qd.screenBits.bounds.bottom-(mb+1))/CharSkip;
#else
    y = (screenBits.bounds.bottom-(mb+1))/CharSkip;
#endif
    if( y>CmndRows ) y = CmndRows;

    CmndRect.bottom = y*CharSkip+1;
    CmndRect.right = CmndCols*CharWide+17;
    CmndRect.left = CharWide+49;
    CmndRect.top = 62;   

    /* Set Zoom-In Size */
    wsdp = (WStateData*) *(wpeek->dataHandle);
    wsdp->stdState.bottom = CmndRect.bottom+mb;
    wsdp->stdState.right = CmndRect.right+4;
    wsdp->stdState.left = 4;
    wsdp->stdState.top = mb;
    
    memset(TermScreen,32,CmndSize);
}


static void DrawCmndGrowIcon()
{
    register RgnHandle saveRgn;
    GrafPtr savePort;
    Rect rect;
    
    GetPort(&savePort);
    SetPort(CmndWin);
    saveRgn = NewRgn();
    GetClip(saveRgn);
    
    rect = CmndWin->portRect;
    rect.left = rect.right-15;
    ClipRect(&rect);
    DrawGrowIcon(CmndWin);
    
    SetClip(saveRgn);
    DisposeRgn(saveRgn);
    SetPort(savePort);
}


static void PaintScreen()
{
    register char *ptr;
    register WindowPtr win;
    register int SRow,ERow;
    register int SCol,ECol;
    register int row,len;
    register int x,y;
    Rect rect;
    
    /* Toolbox Components */
    DrawControls(CmndWin);
    DrawCmndGrowIcon();
    
    /* Determine Invalid Rect? */
    SRow = 0;   ERow = TermRows-1;
    SCol = 0;   ECol = TermCols-1;
    
    len = ECol-SCol+1;
    x = SCol*CharWide + 1;
    y = SRow*CharSkip + CharHigh;
    
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
    ptr = TermScreen + SRow*CmndCols + SCol;
    while( True )
    {   MoveTo(x,y);
        DrawText(ptr,0,len);
        if( row != ERow )
        {   ptr += CmndCols;
            row++;
            if( row == CmndRows )
            {   ptr = TermScreen + SCol;
                row = 0;
            }
        } else break;
        y += CharSkip;
    }
    
    /* Erase Screen Edges! */
    rect = CmndWin->portRect;
    rect.right -= 15;
    
    x = TermCols*CharWide+1;
    y = TermRows*CharSkip;
    if( x<rect.right )
    {   rect.left = x;
        EraseRect(&rect);
    }
    
    if( y<rect.bottom )
    {   rect.left = 0;
        rect.top = y;
        EraseRect(&rect);
    }
    
    win = FrontWindow();
    if( (win==CanvWin) || (win==CmndWin) )
    {   row = TermYPos + ScrlStart;
        if( row < TermRows )
        {   SetCaretPos(TermXPos*CharWide,row*CharSkip);
            TermCursor = True;
            ShowCaret();
        } else 
        {   TermCursor = False;
            HideCaret();
        }
    } else
    {   TermCursor = False;
        HideCaret();
    }
    
    /* Validate CmndWin */
    ValidRect(&CmndWin->portRect);
}


void AdviseUpdate( item )
    int item;
{
    if( item == AdvName )
        EnableMenus( !DisableMenu );
}


void RefreshScreen()
{
    ReDrawFlag &= ~(RFTransZ|RFPoint);
    
    if( ReDrawFlag )
    {   if( ReDrawFlag & RFReSize )
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
            ClearImage();
        }
        ReDrawFlag = 0;
    }
}


static void ConvertFilename( fss )
    FSSpec *fss;
{
    register char *src;
    register char *dst;
    register int i;
    char buffer[256];
    
    Str255 dirname;
    DirInfo dinfo;
    
    src = buffer;
    dinfo.ioDrParID = fss->parID;
    dinfo.ioNamePtr = dirname;
    do {
        dinfo.ioVRefNum = fss->vRefNum;
        dinfo.ioFDirIndex = -1;
        dinfo.ioDrDirID = dinfo.ioDrParID;
        PBGetCatInfo((CInfoPBPtr)&dinfo,0);
        
        *src++ = ':';
        for( i=dirname[0]; i; i-- )
            *src++ = dirname[i];
    } while( dinfo.ioDrDirID != 2 );
    
    /* Reverse the file path! */
    dst = Filename;
    while( src != buffer )
        *dst++ = *(--src);
    for( i=1; i<=fss->name[0]; i++ )
        *dst++ = fss->name[i];
    *dst = '\0';
}


static void HandleAboutDialog()
{
    register char *fpu;
    register char *src;
    register char *dst;
    register int len;
    
    Str255 temp;
    Str255 buffer;
    DialogPtr dlg;
    Handle hand;
    long reply;
    short item;
    Rect rect;
    
    dlg = GetNewDialog(170,(Ptr)0,(WindowPtr)-1);
    SetDialogDefaultItem(dlg,1);
    
    /* System Information! */
    dst = (char*)&buffer[1];
    Gestalt(gestaltPhysicalRAMSize,&reply);
    reply = (long)reply>>10;
   
    if( reply >= (long)1024 )
    {   len = sprintf(dst,"%ldMb ",reply>>10);
    } else len = sprintf(dst,"%ldKb ",reply);
    dst += len;
    
    Gestalt(gestaltMachineType,&reply);
    GetIndString(temp,kMachineNameStrID,(short)reply);
    for( len=1; len<=temp[0]; len++ ) 
        *dst++ = temp[len];        
    
    if( Gestalt(gestaltAUXVersion,&reply) )
    {   for( src=" (System "; *src; src++ )
            *dst++ = *src;
            
        Gestalt(gestaltSystemVersion,&reply);
        *dst++ = (char)((reply>>8)&0x0f) + '0';  *dst++ = '.';
        *dst++ = (char)((reply>>4)&0x0f) + '0';  *dst++ = '.';
        *dst++ = (char)(reply&0x0f) + '0';       *dst++ = ')';
    } else /* A/UX */
        for( src=" (A/UX)"; *src; src++ )
            *dst++ = *src;
            
    len = dst - (char*)&buffer[1];
    buffer[0] = (unsigned char)len;
    GetDItem(dlg,6,&item,&hand,&rect);
    SetIText(hand,buffer);
    
    /* Machine Information! */
    dst = (char*)&buffer[1];
    if( Gestalt(gestaltNativeCPUtype,&reply) )
    {   /* NativeCPUtype not available! */
        Gestalt(gestaltProcessorType,&reply);
        switch( reply )
        {   case(gestalt68000):  src = "MC68000"; break;
            case(gestalt68010):  src = "MC68010"; break;
            case(gestalt68020):  src = "MC68020"; break;
            case(gestalt68030):  src = "MC68030"; break;
            case(gestalt68040):  src = "MC68040"; break;
            default:   src = "Unknown processor";
        }
    } else switch( reply )
        {   case(gestaltCPU68000):  src = "MC68000"; break;
            case(gestaltCPU68010):  src = "MC68010"; break;
            case(gestaltCPU68020):  src = "MC68020"; break;
            case(gestaltCPU68030):  src = "MC68030"; break;
            case(gestaltCPU68040):  src = "MC68040"; break;
            case(gestaltCPU601):    src = "PPC601";  break;
            case(gestaltCPU603):    src = "PPC603";  break;
            case(gestaltCPU604):    src = "PPC604";  break;
            default:   src = "Unknown processor";
        }

    if( *src != 'P' )
    {   Gestalt(gestaltFPUType,&reply);
        switch( reply )
        {   default:               fpu=" unknown";  break;
            case(gestalt68881):    fpu=" 68881";    break;
            case(gestalt68882):    fpu=" 68882";    break;
            case(gestalt68040FPU): fpu=" internal"; break;
            case(gestaltNoFPU):    
                        fpu="out";
                        if( !strcmp(src,"MC68040") )
                            src = "MC68LC040";       
                        break;
        }
    } else fpu = " internal";
    
    while( *src ) *dst++ = *src++;
    for( src=" with"; *src; src++ )
        *dst++ = *src;
    while( *fpu ) *dst++ = *fpu++;
    for( src=" maths coprocessor"; *src; src++ )
        *dst++ = *src;
        
    len = dst - (char*)&buffer[1];
    buffer[0] = (unsigned char)len;
    GetDItem(dlg,7,&item,&hand,&rect);
    SetIText(hand,buffer);

    /* Display Dialog Box! */    
    ShowWindow(dlg);
    do {
#ifdef __CONDITIONALMACROS__
        ModalDialog((ModalFilterUPP)0,&item);
#else
        ModalDialog((ProcPtr)0,&item);
#endif
    } while( item != 1 );
    DisposDialog(dlg);
}


static void PasteCommandText()
{
    register Handle hand;
    register long i,len;
    register char *ptr;
    register char ch;
    long offset;
    
    hand = NewHandle(0);
    len = GetScrap(hand,'TEXT',&offset);
    if( len > 0 )
    {   HLock(hand);
        ptr = (char*)*hand;
        for( i=0; i<len; i++ )
        {   ch = *ptr++;
            if( ch >= ' ' )
            {   ProcessCharacter(ch);
            } else if( (ch==0x0d) || (ch==0x0a) )
            {   ProcessCharacter(ch);
                if( !ExecuteCommand() )
                {   ResetCommandLine(0);
                } else RasMolExit();
            } else if( ch>=' ' )
                ProcessCharacter(ch);
        }
        HUnlock(hand);
    }
}


pascal short OpenDlgHook( item, dialog, data )
    short item;  DialogPtr dialog;  void *data;
{
    Handle hand;
    short type;
    Rect rect;

    if( ((DialogPeek)dialog)->window.refCon != sfMainDialogRefCon )
        return( item );
 
    if( item == sfHookFirstCall )
    {   GetDItem(dialog,10,&type,&hand,&rect);
        SetCtlValue((ControlHandle)hand,DialogFormat);
        item = sfHookNullEvent;
    } else if( item == 10 )
    {   GetDItem(dialog,10,&type,&hand,&rect);
        DialogFormat = GetCtlValue((ControlHandle)hand);
        item = sfHookNullEvent;
    }
    return( item );
}

pascal short SaveDlgHook( item, dialog, data )
    short item;  DialogPtr dialog;  void *data;
{
    Handle hand;
    short type;
    Rect rect;

    if( ((DialogPeek)dialog)->window.refCon != sfMainDialogRefCon )
        return( item );
    
    if( item == sfHookFirstCall )
    {   GetDItem(dialog,13,&type,&hand,&rect);
        SetCtlValue((ControlHandle)hand,DialogFormat);
        item = sfHookNullEvent;
    } else if( item == 13 )
    {   GetDItem(dialog,13,&type,&hand,&rect);
        DialogFormat = GetCtlValue((ControlHandle)hand);
        item = sfHookNullEvent;
    }
    return( item );
}


static void HandleFileOpen()
{
    register int format;
    StandardFileReply reply;
    SFTypeList types;
    Point pnt;
    

    /* File Types */
    types[0]='TEXT';
    types[1]='RSML';
    types[2]='mMOL';

    DialogFormat = 1;
    pnt.v = pnt.h = -1;
    CustomGetFile( NULL, 3, types, &reply, 172, pnt,
#ifdef __CONDITIONALMACROS__
                   OpenDlgHookPtr,
#else
                   (DlgHookYDProcPtr)OpenDlgHook,
#endif 
                   NULL, 0, NULL, NULL );
    
    if( reply.sfGood )
    {   ConvertFilename( &reply.sfFile );
        switch( DialogFormat )
        {   case(1):  format = FormatXYZ;      break;
        }
        
        FetchFile(format,True,Filename);
        DefaultRepresentation();
    }
}

static void HandleFileSave()
{
    StandardFileReply reply;
    Point pnt;
    
    DialogFormat = 1;
    pnt.v = pnt.h = -1;
    CustomPutFile("\pSave Coordinate File:","\p",&reply,173,pnt,
#ifdef __CONDITIONALMACROS__
                   SaveDlgHookPtr,
#else
                   (DlgHookYDProcPtr)SaveDlgHook, 
#endif
                   NULL, 0, NULL, NULL );
    
    if( reply.sfGood )
    {   ConvertFilename( &reply.sfFile );
        switch( DialogFormat )
        {   default:
            case(1): SaveXYZMolecule(Filename);      break;
        }
    }
}

static void HandleExportMenu( item )
    int item;
{
    StandardFileReply reply;
    register unsigned char *ptr;
    register int resid;
    Point pnt;
    
    switch( item )
    {   case(1):  resid=0;   ptr="\pSave GIF Image:";          break;
        case(2):  resid=174; ptr="\pSave PostScript File:";    break;
        case(3):  resid=175; ptr="\pSave Portable PixMap:";    break;
        case(4):  resid=176; ptr="\pSave SUN Rasterfile:";     break;
        case(5):  resid=0;   ptr="\pSave Microsoft BMP File:"; break;
        case(6):  resid=0;   ptr="\pSave PICT Image:";         break;
        default:  return;
    }
    
    if( resid )
    {   DialogFormat = 1;
        pnt.v = pnt.h = -1;
        CustomPutFile( ptr, "\p", &reply, resid, pnt,
#ifdef __CONDITIONALMACROS__
                       SaveDlgHookPtr,
#else
                      (DlgHookYDProcPtr)SaveDlgHook,
#endif 
                      NULL, 0, NULL, NULL );
    } else StandardPutFile( ptr, "\p", &reply );
    if( !reply.sfGood ) return;
    
    ConvertFilename( &reply.sfFile );
    switch( item )
    {   case(1):  WriteGIFFile(Filename);  break;

        case(2):  if( DialogFormat == 1 )
                  {   WriteEPSFFile(Filename,True,True);
                  } else if( DialogFormat == 2 )
                  {   WriteEPSFFile(Filename,False,True); 
                  } else /* DialogFormat == 3 */
                      WriteVectPSFile(Filename);
                  break;
                  
        case(3):  WritePICTFile(Filename); break;
    }
}


static void HandleAppleMenu( item )
    int item;
{
    GrafPtr port;
    Str255 name;
    
    GetPort(&port);
    GetItem( GetMHandle(140), item, name );
    OpenDeskAcc(name);
    SetPort(port);
}


static void HandleMenu( hand )
    long hand;
{
    register int i,mask;
    register int menu;
    register int item;
    
    menu = HiWord(hand);
    if( menu )
    {   item = LoWord(hand);
        switch( menu )
        {   case(140):  /* Apple Menu */
                        if( item == 1 )
                        {   HandleAboutDialog();
                        } else if( item>2 )
                            HandleAppleMenu(item);
                        break;
                      
            case(141):  /* File Menu */
                        switch( item )
                        {   case(1):  /* Open */
                                      if( !Database )
                                          HandleFileOpen();
                                      break;
                            case(2):  /* Save As */
                                      if( Database )
                                          HandleFileSave();
                                      break;
                            case(3):  /* Close */
                                      ZapDatabase();
                                      break;
                                      
                            case(5):  /* Page Setup */
                                      PrOpen();
                                      if( !PrintHand )
                                      {   PrintHand = (THPrint)
                                              NewHandle(sizeof(TPrint));
                                          PrintDefault(PrintHand);
                                      }
                                      PrStlDialog(PrintHand);
                                      PrClose();
                                      break;
                                      
                            case(6):  /* Print */
                                      PrintImage();
                                      break;
                                      
                            case(8):  /* Quit */
                                      RasMolExit();
                        }
                        break;
                      
            case(142):  /* Edit Menu */
                        switch( item )
                        {   case(1):  /* Undo */
                                      for( i=0; i<8; i++ )
                                          DialValue[i] = 0.0;
                                      ReDrawFlag |= RFDials;
                                      ResetTransform();
                                      UpdateScrollBars();
                                      break;
                                      
                            case(3):  /* Copy */
                                      ClipboardImage();
                                      break;
                            
                            case(4):  /* Paste */
                                      PasteCommandText();
                                      break;
                                      
                            case(7):  /* Select All */
                                      mask = NormAtomFlag;
                                      if( Hydrogens )  mask |= HydrogenFlag;
                                      SelectZone(mask);
                                      break;
                        }
                        break;
                        
            case(143):  /* Display Menu */
                        switch( item )
                        {   case(1):  /* Wireframe */
                                      DisableSpacefill();
                                      EnableWireframe(WireFlag,0);
                                      ReDrawFlag |= RFRefresh;
                                      break;
                                      
                            case(2):  /* Sticks */
                                      DisableSpacefill();
                                      if( MainAtomCount<256 )
                                      {   EnableWireframe(CylinderFlag,40);
                                      } else EnableWireframe(CylinderFlag,80);
                                      ReDrawFlag |= RFRefresh;
                                      break;
                                      
                            case(3):  /* Spheres */
                                      SetVanWaalRadius();
                                      DisableWireframe();
                                      ReDrawFlag |= RFRefresh;
                                      break;
                                      
                            case(4):  /* Ball & Stick */
                                      SetRadiusValue(120);
                                      EnableWireframe(CylinderFlag,40);
                                      ReDrawFlag |= RFRefresh;
                                      break;
                        }
                        break;              
                            
            case(144):  /* Colours Menu */
                        switch( item )
                        {   case(1):  /* Monochrome */
                                      MonoColourAttrib(255,255,255);
                                      ReDrawFlag |= RFColour;  break;

                            case(2):  /* CPK */
                                      CPKColourAttrib();
                                      ReDrawFlag |= RFColour;  break;

                            case(3):  /* Shapely */
                                      ShapelyColourAttrib();
                                      ReDrawFlag |= RFColour;  break;
                        }
                        break;
                                                    
            case(145):  /* Option Menu */
                        switch( item )
                        {   case(1):  /* Slabbing */
                                      ReDrawFlag |= RFRefresh;
                                      UseSlabPlane = !UseSlabPlane;
                                      if( UseSlabPlane )
                                          UseShadow = False;
                                      break;
                                      
                            case(2):  /* Hydrogens */
                                      mask = NormAtomFlag;
                                      Hydrogens = !Hydrogens;
                                      ReDrawFlag |= RFRefresh;
                                      
                                      if( Hydrogens )
                                      {   SelectZone(mask|HydrogenFlag);
                                      } else RestrictZone(mask);
                                      break;
                                      
                            case(3):  /* Specular */
                                      FakeSpecular = !FakeSpecular;
                                      ReDrawFlag |= RFColour;
                                      break;
                                      
                            case(4):  /* Shadows */
                                      ReDrawFlag |= RFRefresh;
                                      UseShadow = !UseShadow;
                                      if( UseShadow )
                                      {   ReviseInvMatrix();
                                          VoxelsClean = False;
                                          UseSlabPlane = False;
                                          ReAllocBuffers();
                                      }
                                      break;

                            case(5):  /* Stereo */
                                      if( UseStereo )
                                      {   SetStereoMode(False);
                                      } else SetStereoMode(True);
                                      ReDrawFlag |= RFRefresh;
                                      break;

                            case(6):  /* Labels */
                                      LabelOptFlag = !LabelOptFlag;
                                      DefaultLabels(LabelOptFlag);
                                      ReDrawFlag |= RFRefresh;
                                      break;
                        }
                        break;
                        
            case(146):  /* Export Menu */
                        if( Database )
                            HandleExportMenu( item );
                        break;
                        
            case(147):  /* Windows Menu */
                        if( item==1 )
                        {   SelectWindow(CanvWin);
                            ShowWindow(CanvWin);
                        } else /* item==2 */
                        {   SelectWindow(CmndWin);
                            ShowWindow(CmndWin);
                        }
                        break;
        }
        HiliteMenu(0);
    }
}


static void AdjustMenus()
{
    register MenuHandle menu;
    register WindowPtr win;
    
    /* Refresh Menus */
    EnableMenus(!DisableMenu);
    
    /* Options Menu */
    menu = GetMHandle(145);
    CheckItem(menu,1,UseSlabPlane);
    CheckItem(menu,2,Hydrogens);
    CheckItem(menu,3,FakeSpecular);
    CheckItem(menu,4,UseShadow);
    CheckItem(menu,5,UseStereo);
    CheckItem(menu,6,LabelOptFlag);

    /* Windows Menu */
    win = FrontWindow();
    menu = GetMHandle(147);
    CheckItem(menu,1,(win==CanvWin));
    CheckItem(menu,2,(win==CmndWin));
}

static void ReSizeCanvWin()
{
    register int x,y;

    x = CanvWin->portRect.right - CanvWin->portRect.left;
    y = CanvWin->portRect.bottom - CanvWin->portRect.top;

    /* ClipRect(&CanvWin->portRect); */
    InvalRect(&CanvWin->portRect);
    
    HidePen();
    HideControl(HScroll);
    HideControl(VScroll);
    MoveControl(HScroll,-1,y-15);  SizeControl(HScroll,x-14,16);    
    MoveControl(VScroll,x-15,-1);  SizeControl(VScroll,16,y-14);
    ShowPen();
    
    ShowControl(HScroll);
    ShowControl(VScroll);
    DrawGrowIcon(CanvWin);

    XRange = x-15;   WRange = XRange>>1;
    YRange = y-15;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);
    ReDrawFlag |= RFReSize;
    ClearImage();
}


static void GrowCanvWin( pos )
    Point pos;
{
    register long size;
    register int x,y,dx;
    GrafPtr savePort;
    Rect rect;
    
#ifdef __CONDITIONALMACROS__
    rect.bottom = qd.screenBits.bounds.bottom;
    rect.right = qd.screenBits.bounds.right;
#else
    rect.bottom = screenBits.bounds.bottom;
    rect.right = screenBits.bounds.right;
#endif

    rect.left = 128;
    rect.top = 128;
    
    size = GrowWindow(CanvWin,pos,&rect);
    if( !size ) return; /* No Change! */
    
    GetPort(&savePort);
    SetPort(CanvWin);
    
    x = LoWord(size)-15;
    y = HiWord(size);
    
    /* Ensure Long Aligned */
    dx = x%4;
    if( dx ) x += 4-dx;

    SizeWindow(CanvWin,x+15,y,false);
    ReSizeCanvWin();
    SetPort(savePort);
}

static void ReSizeCmndWin()
{
    register int rows,cols;
    register int sr,er;
    register int x, y;
    
    x = CmndWin->portRect.right - CmndWin->portRect.left;
    y = CmndWin->portRect.bottom - CmndWin->portRect.top;
    
    InvalRect(&CmndWin->portRect);
    MoveControl(CmndScroll,x-15,-1);
    SizeControl(CmndScroll,16,y-14);    
    
    cols = (x-16)/CharWide;
    rows = y/CharSkip;
    
    /* Scroll to bottom */
    if( ScrlStart )
        SetTermScroll( ScrlMax );

    if( rows < TermRows )
    {   if( TermYPos >= rows )
        {   CmndStart += (TermYPos-rows) + 1;
            if( CmndStart > CmndRows )
                CmndStart -= CmndRows;
            TermYPos = rows-1;
        }
        
    } else if( rows > TermRows )
    {   sr = TermRows + CmndStart;
        if( sr >= CmndRows )
            sr -= CmndRows;
            
        er = CmndStart + rows;
        if( er >= CmndRows )
            er -= CmndRows;
            
        do {
            memset(TermScreen+sr*CmndCols,32,CmndCols);
            sr++; if( sr == CmndRows ) sr = 0;
        } while( sr != er );
    }
    
    if( cols > CmndCols )
    {   TermCols = CmndCols;
    } else TermCols = cols;
    TermRows = rows;
}


static void GrowCmndWin( pos )
    Point pos;
{
    register long size;
    GrafPtr savePort;
    
    size = GrowWindow(CmndWin,pos,&CmndRect);
    if( !size ) return;  /* No Change! */
    
    GetPort(&savePort);
    SetPort(CmndWin);
    SizeWindow(CmndWin,LoWord(size),HiWord(size),false);
    ReSizeCmndWin();
    SetPort(savePort);
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


#define ShiftModifier   0x2200
#define CntrlModifier   0x9000
#define LButtModifier   0x0080
#define MButtModifier   0x0100  /* [Option]  */
#define RButtModifier   0x4800  /* [Command] */

static void MouseMove( status, dx, dy )
    int status, dx, dy;
{
    if( MouseMode == MMRasMol )
    {   if( status & ShiftModifier )
        {   if( status & (MButtModifier|RButtModifier) )
            {   if( dx )  /* Z Rotation Horizontal */
                {   WrapDial( 2, (Real)dx/WRange );
                    ReDrawFlag |= RFRotateZ;
                }
            } else
                if( dy )  /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
        
        } else if( status & CntrlModifier )
        {   if( dy )   /* Slab Vertical */
            {   ClampDial( 7, (Real)dy/YRange );
                ReDrawFlag |= RFSlab;
            }
        
        } else /* Unmodified */
            if( status & (MButtModifier|RButtModifier) )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }
                if( dy ) /* Translate Y Vertical   */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            } else
            {   if( dx ) /* Rotate Y Horizontal */
                {   WrapDial( 1, (Real)dx/WRange );
                    ReDrawFlag |= RFRotateY;
                }
                if( dy ) /* Rotate X Vertical */
                {   WrapDial( 0, (Real)dy/HRange );
                    ReDrawFlag |= RFRotateX;
                }
                UpdateScrollBars();
            }

    } else if( MouseMode == MMQuanta )
    {   if( status & ShiftModifier )
        {   if( status & LButtModifier )
            {   if( dy ) /* Slab Vertical */
                {   ClampDial( 7, (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                }
            } else if( status & MButtModifier )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }
                if( dy ) /* Translate Y Vertical   */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            } else if( !(status & RButtModifier) )
                if( dy )  /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
        } else if( status & MButtModifier )
        {   if( dx ) /* Rotate Y Horizontal */
            {   WrapDial( 1, (Real)dx/WRange );
                ReDrawFlag |= RFRotateY;
            }
            if( dy ) /* Rotate X Vertical */
            {   WrapDial( 0, (Real)dy/HRange );
                ReDrawFlag |= RFRotateX;
            }
            UpdateScrollBars();
        } else if( status & RButtModifier )
            if( dx )  /* Z Rotation Horizontal */
            {   WrapDial( 2, (Real)dx/WRange );
                ReDrawFlag |= RFRotateZ;
            }
        
    } else /* MMInsight */
        if( status & LButtModifier )
        {   if( status & MButtModifier )
            {   if( status & RButtModifier )
                {   ClampDial( 7, (Real)dx/XRange - 
                                  (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                } else
                {   ClampDial( 3, (Real)dx/WRange - 
                                  (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
            } else if( status & RButtModifier )
            {   WrapDial( 2, (Real)dx/WRange - 
                             (Real)dy/HRange );
                ReDrawFlag |= RFRotateZ;
            } else
            {   if( dx ) /* Rotate Y Horizontal */
                {   WrapDial( 1, (Real)dx/WRange );
                    ReDrawFlag |= RFRotateY;
                }
                if( dy ) /* Rotate X Vertical */
                {   WrapDial( 0, (Real)dy/HRange );
                    ReDrawFlag |= RFRotateX;
                }
                UpdateScrollBars();
            }
        } else if( status & MButtModifier )
        {   if( dx ) /* Translate X Horizontal */
            {   ClampDial( 4, (Real)dx/XRange );
                ReDrawFlag |= RFTransX;
            }
            if( dy ) /* Translate Y Vertical   */
            {   ClampDial( 5, (Real)dy/YRange );
                ReDrawFlag |= RFTransY;
            }
        } 
}


pascal void CanvScrollProc( cntrl, code )
    ControlHandle cntrl;  short code;
{   
    register int pos;
    
    pos = GetCtlValue(cntrl); 
    switch( code )
    {   case(inUpButton):    pos -= 5;   break;
        case(inDownButton):  pos += 5;   break;
        case(inPageUp):      pos -= 10;  break;
        case(inPageDown):    pos += 10;  break;
        default:             return;
    }
    
    if( pos>100 )
    {   pos -= 100;
    } else if( pos<0 )
        pos += 100;
        
    SetCtlValue(cntrl,pos);
    if( cntrl == HScroll )
    {   DialValue[1] = (pos/50.0)-1.0;
        ReDrawFlag |= RFRotateY;
    } else /* cntrl == VScroll */
    {   DialValue[0] = (pos/50.0)-1.0;
        ReDrawFlag |= RFRotateX; 
    }
    RefreshScreen();
}


static void ClickCanvWin( ptr )
    EventRecord *ptr;
{
    register int code,pos;

    ControlHandle hand;
    GrafPtr savePort;
    
    if( CanvWin == FrontWindow() )
    {   GetPort(&savePort);
        SetPort(CanvWin);
        
        GlobalToLocal(&ptr->where);
        code = FindControl(ptr->where,CanvWin,&hand);
        if( !code )
        {   InitX = PointX = ptr->where.h;
            InitY = PointY = ptr->where.v;
            HeldButton = True;
            
        } else if( code == inThumb )   
        {   TrackControl(hand,ptr->where,0L);
            pos = GetCtlValue(hand);
            if( hand == HScroll )
            {   DialValue[1] = (pos/50.0)-1.0;
                ReDrawFlag |= RFRotateY;
            } else /* hand == VScroll */
            {   DialValue[0] = (pos/50.0)-1.0;
                ReDrawFlag |= RFRotateX; 
            }
            RefreshScreen();
        } else TrackControl(hand,ptr->where,
#ifdef __CONDITIONALMACROS__
                         CanvScrollProcPtr);
#else
                         (ProcPtr)CanvScrollProc );
#endif
        SetPort(savePort);
    } else SelectWindow( CanvWin );
}



pascal void CmndScrollProc( cntrl, code )
    ControlHandle cntrl;  short code;
{
    switch( code )
    {   case(inUpButton):    if( ScrlStart < ScrlMax )
                             {   SetTermScroll((ScrlMax-ScrlStart)-1);
                                 PaintScreen();
                             }
                             break;
                             
        case(inDownButton):  if( ScrlStart > 0 )
                             {   SetTermScroll((ScrlMax-ScrlStart)+1);
                                 PaintScreen();
                             }
                             break;
                                 
        case(inPageUp):      if( ScrlStart < (ScrlMax-10) )
                             {   SetTermScroll((ScrlMax-ScrlStart)-10);
                                 PaintScreen();
                             }
                             break;
                             
        case(inPageDown):    if( ScrlStart > 10 )
                             {   SetTermScroll((ScrlMax-ScrlStart)+10);
                                 PaintScreen();
                             }
                             break;
    }
}


static void ClickCmndWin( ptr )
    EventRecord *ptr;
{
    register int code;
    ControlHandle hand;
    GrafPtr savePort;
    
    if( CmndWin == FrontWindow() )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        
        GlobalToLocal(&ptr->where);
        code = FindControl(ptr->where,CmndWin,&hand);
        if( code == inThumb )   
        {   TrackControl(CmndScroll,ptr->where,0L);
            SetTermScroll(GetCtlValue(CmndScroll));
        } else if( code )
            TrackControl(CmndScroll,ptr->where,
#ifdef __CONDITIONALMACROS__
                         CmndScrollProcPtr);
#else
                         (ProcPtr)CmndScrollProc );
#endif
        SetPort(savePort);
    } else SelectWindow( CmndWin );
}


static void ZoomCanvWin( pos, code )
    Point pos;  int code;
{
    GrafPtr savePort;
    
    if( TrackBox(CanvWin,pos,code) )
    {   GetPort(&savePort);
        SetPort(CanvWin);
        /* EraseRect(&CanvWin->portRect); */
        ZoomWindow(CanvWin,code,true);
        ReSizeCanvWin();
        SetPort(savePort);
    }   
}


static void ZoomCmndWin( pos, code )
    Point pos;  int code;
{
    GrafPtr savePort;
                                    
    if( TrackBox(CmndWin,pos,code) )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        EraseRect(&CmndWin->portRect);
        ZoomWindow(CmndWin,code,true);
        ReSizeCmndWin();
        SetPort(savePort);
     }
}


static void HandleMouseDownEvent( ptr )
    EventRecord *ptr;
{
    register long hand;
    register int code;
    WindowPtr win;
    
    code = FindWindow(ptr->where,&win);
    switch( code )
    {   case(inMenuBar):    AdjustMenus();
                            hand = MenuSelect(ptr->where);
                            if( !IsPaused )
                                HandleMenu( hand );
                            break;
                            
        case(inSysWindow):  SystemClick(ptr,win);
                            break;

        case(inContent):    if( win == CanvWin )
                            {   ClickCanvWin( ptr );
                            } else if( win == CmndWin )
                            {   ClickCmndWin( ptr );
                            }
                            break;
        
        case(inDrag):       if( (win==CanvWin) || (win==CmndWin) )
#ifdef __CONDITIONALMACROS__
                                DragWindow(win,ptr->where,&qd.screenBits.bounds);
#else
                                DragWindow(win,ptr->where,&screenBits.bounds);
#endif
                            break;
                             
        case(inGrow):       if( win==CanvWin )
                            {   GrowCanvWin( ptr->where );
                            } else if( win==CmndWin )
                                GrowCmndWin( ptr->where );
                            break;
        
        case(inGoAway):     if( (win==CanvWin) || (win==CmndWin) )
                                if( TrackGoAway(win,ptr->where) )
                                    HideWindow(win);
                            break;
                            
        case(inZoomIn):     
        case(inZoomOut):    if( win==CanvWin )
                            {   ZoomCanvWin(ptr->where,code);
                            } else if( win==CmndWin )
                                ZoomCmndWin(ptr->where,code);
                            break;
    }
}



static void HandleMoveEvent( status )
    int status;
{
    register WindowPtr win;
    register int dx,dy;
    GrafPtr savePort;
    Point pos;
    
    win = FrontWindow();
    if( win==CanvWin )
    {   GetPort(&savePort);
        SetPort(CanvWin);
        GetMouse(&pos);
        
        /* Invert Button Bit */
        status ^= LButtModifier;
        
        if( (status & (LButtModifier|MButtModifier|RButtModifier))
            || ((MouseMode==MMQuanta) && (status&ShiftModifier)) )
        {   if( !HeldButton )
            {   InitX = PointX = pos.h;
                InitY = PointY = pos.v;
                HeldButton = True;
            }
        } else HeldButton = False;

        if( HeldButton && !IsClose(pos.h,InitX) 
                       && !IsClose(pos.v,InitY) )
        {   dx = pos.h - PointX;
            dy = pos.v - PointY;
            MouseMove(status,dx,dy);
            PointX = pos.h;
            PointY = pos.v;
        }   
        
        if( HeldButton || ((pos.h>0) && (pos.v>0) &&
              (pos.v<CanvWin->portRect.bottom-15) &&
              (pos.h<CanvWin->portRect.right-15)) )
        {   SetCursor(*CanvCursor);
        } else ArrowCursor;
        SetPort(savePort);
    } else if( win==CmndWin )
    {   GetPort(&savePort);
        SetPort(CmndWin);
        GetMouse(&pos);
        
        if( (pos.h>0) && (pos.v>0) &&
            (pos.v<CmndWin->portRect.bottom) &&
            (pos.h<CmndWin->portRect.right-15) )
        {   SetCursor(*CmndCursor);
        } else ArrowCursor;
        SetPort(savePort);
    }
}


static int MacKeyMap[32] = { 0x00, 0x01, 0x00, 0x0d, 0x05, 0x00, 0x00, 0x00,
                             0x08, 0x00, 0x00, 0x00, 0x00, 0x0d, 0x00, 0x00,
                             0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                             0x00, 0x00, 0x00, 0x00, 0x02, 0x06, 0x10, 0x0e };
                              
static void HandleEvents()
{
    register int key,row;
    register long hand;

    EventRecord event;
    GrafPtr savePort;
    WindowPtr win;

    SystemTask();
    if( GetNextEvent(everyEvent,&event) )
    {   switch( event.what )
        {   case(mouseDown):    HeldButton = False;
                                HandleMouseDownEvent(&event);
                                break;
                                
            case(mouseUp):      if( HeldButton && Database )
                                {   GetPort(&savePort);
                                    SetPort(CanvWin);
                                    GlobalToLocal(&event.where);
                                    SetPort(savePort);
                                    
                                    PointX = event.where.h;
                                    PointY = event.where.v;
                                    
                                    if( IsClose(PointX,InitX) &&
                                        IsClose(PointY,InitY) )
                                    {   if( event.modifiers & 
                                            (ShiftModifier|CntrlModifier) )
                                        {      PickAtom(True,PointX,PointY);
                                        } else PickAtom(False,PointX,PointY);
                                        /* AdviseUpdate(AdvPickNumber); */
                                        /* AdviseUpdate(AdvPickAtom);   */
                                    }
                                }
                                HeldButton = False;
                                break;
            
            case(autoKey):
            case(keyDown):      key = (char)(event.message & charCodeMask);
                                if( event.modifiers & cmdKey )
                                {   AdjustMenus();
                                    hand = MenuKey(key);
                                    if( !IsPaused )
                                        HandleMenu(hand);
                                } else if( key<32 )
                                {   if( !(event.modifiers & controlKey) )
                                        key = MacKeyMap[key];
                                    if( ProcessCharacter(key) )
                                        if( ExecuteCommand() )
                                            RasMolExit();
                                } else if( key==127 )
                                {   /* Forward Delete */
                                    ProcessCharacter(0x04);
                                } else ProcessCharacter(key);
                                break;
            
            case(keyUp):        break;
            
            case(updateEvt):    GetPort(&savePort);
                                win = (WindowPtr)event.message;
                                if( win == CanvWin )
                                {   SetPort(CanvWin);
                                    BeginUpdate(CanvWin);
                                    /* ActivatePalette(CanvWin); */
                                    DrawGrowIcon(CanvWin);
                                    
                                    /* Macintosh 128K ROMs and later.       */ 
                                    /* UpdtControl(CanvWin,CanvWin->visRgn); */
                                    DrawControls(CanvWin);
                                    
                                    /* CopyBits(PixMap);   */
                                    /* EraseRect(updRect);*/
                                    if( Database )
                                    {   TransferImage();
                                    } else ClearImage();
                                    
                                    EndUpdate(CanvWin);
                                } else if( win == CmndWin )
                                {   SetPort(CmndWin);
                                    BeginUpdate(CmndWin);
                                    PaintScreen();
                                    EndUpdate(CmndWin);
                                }
                                SetPort(savePort);
                                break;
            
            case(activateEvt):  HiliteMenu(0);
                                HeldButton = False;
                                
                                win = (WindowPtr)event.message;
                                if( win==CanvWin )
                                {   DrawGrowIcon(CanvWin);
                                    if( event.modifiers & activeFlag )
                                    {   ShowControl(HScroll);
                                        ShowControl(VScroll);
                                    } else /* deactiveEvt */
                                    {   HideControl(HScroll);
                                        HideControl(VScroll);
                                    }
                                } else /* win==CmndWin */
                                {   DrawCmndGrowIcon();
                                    if( event.modifiers & activeFlag )
                                    {   ShowControl(CmndScroll);
                                    } else /* deactiveEvt */
                                        HideControl(CmndScroll);
                                }
                                
                                /* Caret Handling! */
                                win = FrontWindow();
                                if( (win==CanvWin) || (win==CmndWin) )
                                {   row = TermYPos + ScrlStart;
                                    if( row < TermRows )
                                    {   SetCaretPos( TermXPos*CharWide,
                                                     row*CharSkip );
                                        TermCursor = True;
                                        ShowCaret();
                                    } else
                                    {   TermCursor = False;
                                        HideCaret();
                                    }      
                                } else
                                {   TermCursor = False;
                                    HideCaret();
                                }
                                
                                /* Mouse Move Handling! */
                                if( win == CanvWin )
                                {   GetPort(&savePort);
                                    SetPort(CanvWin);
                                    GetMouse(&MousePrev);
                                    SetPort(savePort);
                                }
                                break;
                                
            case(kHighLevelEvent):
                                AEProcessAppleEvent(&event);
                                break;
                                
        }
    } else if( TermCursor )
        HandleCaret();
        
    HandleMoveEvent( event.modifiers );
    
    if( ReDrawFlag )
        RefreshScreen();
    if( !CommandActive )
        ResetCommandLine(0);
}


static int HandleFileSpec( fss )
    FSSpec *fss;
{
    register int format;
    register FILE *fp;
    FInfo info;
    
    FSpGetFInfo(fss,&info);
    if( info.fdType == 'RSML' )
    {   fp = fopen(Filename,"r");
        if( fp )
        {   LoadScriptFile(fp,Filename);
            return( True );
        } else return( False );
    }
    
    if( Database )
        ZapDatabase();

    format = FormatXYZ;
        
    FetchFile(format,True,Filename);
    if( Database )
    {   DefaultRepresentation();
        return( True );
    }
    return( False );
}


/* Apple Event Handler Prototypes! */
pascal OSErr HandleAEOpenDoc( AppleEvent*, AppleEvent*, long );
pascal OSErr HandleAEQuitApp( AppleEvent*, AppleEvent*, long );
pascal OSErr HandleAEIgnore( AppleEvent*, AppleEvent*, long );


pascal OSErr HandleAEOpenDoc( event, reply, ref )
    AppleEvent *event; AppleEvent *reply; long ref;
{
    register OSErr stat;
    register long i;
    
    AEDescList list;
    AEKeyword keywd;
    DescType dtype;
    FSSpec fss;
    long count;
    Size size;

    /* Disable event while paused! */
    if( IsPaused ) return( noErr );

    stat = AEGetParamDesc(event,keyDirectObject,
                          typeAEList,&list);
    if( stat ) return( stat );
    
    stat = AEGetAttributePtr(event,keyMissedKeywordAttr,
                             typeWildCard, &dtype, 0, 0, &size );
    if( stat != errAEDescNotFound ) 
    {   AEDisposeDesc( &list );
        return( stat? stat : errAEEventNotHandled );
    }
    
    AECountItems( &list, &count );
    for( i=1; i<=count; i++ )
    {   stat = AEGetNthPtr(&list,i,typeFSS,&keywd,
                           &dtype,(Ptr)&fss,sizeof(fss),
                           &size);
        if( !stat )
        {   ConvertFilename(&fss);
            if( HandleFileSpec(&fss) )
            {   RefreshScreen();
                if( Database && ref )
                    PrintImage();
                break;
            }
        }
    }
    AEDisposeDesc( &list );         
    return noErr;
}

pascal OSErr HandleAEQuitApp( event, reply, ref )
    AppleEvent *event; AppleEvent *reply; long ref;
{
    RasMolExit();
    return noErr;
}

pascal OSErr HandleAEIgnore( event, reply, ref )
    AppleEvent *event; AppleEvent *reply; long ref;
{
    return noErr;
}


static void InitialiseApplication()
{
    register PScrapStuff ptr;
    
    /* Init Mac ToolBox */
#ifdef __CONDITIONALMACROS__
    InitGraf(&qd.thePort);
#else
    InitGraf(&thePort);
#endif

    InitCursor();
    InitFonts();
    InitWindows();
    InitMenus();
    TEInit();
    InitDialogs(0L);

    RasMolResFile = CurResFile();
    FlushEvents(everyEvent,0);
    MaxApplZone();
    
    /* Initialise Clipboard */
    ptr = InfoScrap();
    if( ptr && (ptr->scrapState<0) )
        ZeroScrap();
    
    /* Enable KeyUp Events  */
    SetEventMask(everyEvent);
    
#ifdef __CONDITIONALMACROS__
    /* Create Routine Descriptors */
    HandleAEIgnorePtr = NewAEEventHandlerProc(HandleAEIgnore);
    HandleAEOpenDocPtr = NewAEEventHandlerProc(HandleAEOpenDoc);
    HandleAEQuitAppPtr = NewAEEventHandlerProc(HandleAEQuitApp);
    CanvScrollProcPtr = NewControlActionProc(CanvScrollProc);
    CmndScrollProcPtr = NewControlActionProc(CmndScrollProc);
    OpenDlgHookPtr = NewDlgHookYDProc(OpenDlgHook);
    SaveDlgHookPtr = NewDlgHookYDProc(SaveDlgHook);

    /* Install Required Event Handlers */
    AEInstallEventHandler(kCoreEventClass,kAEOpenApplication,
                          HandleAEIgnorePtr, 0, false);
    AEInstallEventHandler(kCoreEventClass,kAEOpenDocuments,
                          HandleAEOpenDocPtr, 0, false);
    AEInstallEventHandler(kCoreEventClass,kAEPrintDocuments,
                          HandleAEOpenDocPtr, 1, false);
    AEInstallEventHandler(kCoreEventClass,kAEQuitApplication,
                          HandleAEQuitAppPtr, 0, false);

#else
    /* Install Required Event Handlers */
    AEInstallEventHandler(kCoreEventClass,kAEOpenApplication,
                          HandleAEIgnore, 0, false);
    AEInstallEventHandler(kCoreEventClass,kAEOpenDocuments,
                          HandleAEOpenDoc, 0, false);
    AEInstallEventHandler(kCoreEventClass,kAEPrintDocuments,
                          HandleAEOpenDoc, 1, false);
    AEInstallEventHandler(kCoreEventClass,kAEQuitApplication,
                          HandleAEQuitApp, 0, false);
#endif
}


static void InitDefaultValues()
{
    Interactive = True;

    ReDrawFlag = RFInitial;
    LabelOptFlag = False;
    CalcBondsFlag = True;
}


int main()
{
    GrafPtr savePort;
    
    InitialiseApplication();
    InitDefaultDefaultValues();

    OpenDisplay(DefaultWide,DefaultHigh);
    InitTerminal();
    
    WriteString("TINKER RasMol Molecule Viewer\n");
    WriteString("Version 4.0, October 2002\n");
    WriteString("Based on RasMol 2.6 by Roger Sayle\n\n");

#ifdef __powerc
    WriteString("[PowerPC Native]\n");
#endif

#ifdef EIGHTBIT
    WriteString("[8bit version]\n\n");
#else
    WriteString("[24bit version]\n\n");
#endif

    InitialiseCommand();
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();
    InitialiseRepres();

    /* LoadInitFile(); */
    
    ResetCommandLine(1);
    RefreshScreen();
    
    GetPort(&savePort);
    SetPort(CanvWin);
    GetMouse(&MousePrev);
    SetPort(savePort);
    
    while( True )
        HandleEvents();
}
