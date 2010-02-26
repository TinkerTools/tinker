/* x11win.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>

#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "bitmaps.h"
#include "command.h"

/* Menu Definitions */
#define mbEnable    0x01
#define mbOption    0x02
#define mbCheck     0x04
#define mbSepBar    0x08
#define mbAccel     0x10

typedef struct _MenuItem {
            char *text;
            int flags;
            int pos;
            int len;
        } MenuItem;

static MenuItem FilMenu[5] = {
    { "Open...",    0x11, 0,  7 },
    { "Save As...", 0x11, 0, 10 },
    { "Close",      0x11, 0,  5 },
    { "",           0x08, 0,  0 },
    { "Exit",       0x11, 0,  4 } };

static MenuItem DisMenu[4] = {
    { "Wireframe",    0x11, 0,  9 },
    { "Sticks",       0x11, 1,  6 },
    { "Spacefill",    0x11, 0,  9 },
    { "Ball & Stick", 0x11, 0, 12 } };

static MenuItem ColMenu[3] = {
    { "Monochrome",  0x11, 0, 10 },
    { "CPK",         0x11, 0,  3 },
    { "Shapely",     0x11, 0,  7 } };

static MenuItem OptMenu[6] = {
    { "Slab Mode",    0x13, 0,  9 },
    { "Hydrogens",    0x17, 1,  9 },
    { "Specular",     0x13, 1,  8 },
    { "Shadows",      0x13, 1,  7 },
    { "Stereo",       0x13, 1,  6 },
    { "Labels",       0x13, 0,  6 } };

static MenuItem ExpMenu[4] = {
    { "GIF...",        0x11, 0,  6 },
    { "PostScript...", 0x11, 0, 13 },
    { "IRIS RGB...",   0x11, 5, 11 },
    { "PICT...",       0x11, 1,  7 } };

static MenuItem HelMenu[2] = {
    { "About RasMol...",  0x10, 0, 15 },
    { "User Manual...",   0x10, 0, 14 } };

typedef struct _BarItem {
            MenuItem *menu;
            char *text;
            int count;
            int flags;
            int len;
        } BarItem;

#define MenuBarMax 6
static BarItem MenuBar[MenuBarMax] = { 
    { FilMenu, "File",    5, 0x01, 4 },
    { DisMenu, "Display", 4, 0x01, 7 },
    { ColMenu, "Colours", 3, 0x01, 7 },
    { OptMenu, "Options", 6, 0x01, 7 },
    { ExpMenu, "Export",  4, 0x01, 6 },
    { HelMenu, "Help",    2, 0x01, 4 } };

static int MenuFocus;
static int ItemFocus;
static int MenuItemSelect;
static int MenuBarSelect;
static int MenuBarCount;
static int PopUpWide;
static int PopUpHigh;
static int PopUpFlag;
static int ItemFlag;

#ifdef DIALBOX
#include <X11/extensions/XInput.h>

static char *DialLabel[] = { "ROTATE X", "ROTATE Y", "ROTATE Z", "  ZOOM  ",
                             "TRANS X ", "TRANS Y ", "TRANS Z ", "  SLAB  " };
static int *DialMap;
static int ESVDialMap[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
static int SGIDialMap[8] = { 3, 7, 2, 6, 1, 5, 0, 4 };

static Real DialRes[8];
static int DialPrev[8];
static int DialMode;

static int UseDialLEDs;
static XDevice *Dials;
static int DialEvent;
static int UseDials;
#endif

#ifdef MITSHM
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>

XShmSegmentInfo xshminfo;
int SharedMemOption;
int SharedMemFlag;
#endif

#define XScrlDial  1 /*1*/
#define YScrlDial  0 /*0*/
#define XScrlSkip  8
#define YScrlSkip  8

typedef union {
    Long longword;
    Byte bytes[4];
    } ByteTest;

/* Determine Mouse Sensitivity! */
#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))

static int MenuHigh;
static int FontHigh;

static Cursor cross;
static Cursor arrow;
static Cursor hglass;
static Pixmap Scrl;
static Pixmap tilepix;
static Pixmap uppix, dnpix;
static Pixmap lfpix, rgpix;
static XFontStruct *MenuFont;
static XSetWindowAttributes attr;
static Window XScrlWin, YScrlWin;
static Window PopUpWin;
static Window MainWin;
static Window CanvWin;
static Window MenuWin;
static Window RootWin;
static XWMHints hints;
static Colormap cmap;
static Colormap lmap;
static XImage *image;
static Display *dpy;
static Visual *vis;
static GC gcon;

#ifdef EIGHTBIT
static unsigned long Ident[256];
static int IdentCount;
#endif

static int InitX, InitY;
static int HeldButton;
static int HeldStep;

static Byte Intensity[LutSize];
static Pixel WhiteCol;
static Pixel BlackCol;
static int Monochrome;

#ifdef THIRTYTWOBIT
static int SwapBytes;
#endif

static int MaxWidth, MaxHeight;
static int MinWidth, MinHeight;
static int MainWide, MainHigh;
static int ScrlX,NewScrlX;
static int ScrlY,NewScrlY;
static int PixDepth;
static int LocalMap;

/* WM_PROTOCOLS */
static char TkInterp[10];
static Atom AppNameAtom;
static Atom DelWinXAtom;
static Atom ProtoXAtom;
static Atom InterpAtom;
static Atom CommAtom;

/* Routine in rasmol.c! */
extern int ProcessCommand();

/* Forward Declarations */
static int HandleMenuLoop();


static void FatalGraphicsError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Graphics Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


void AllocateColourMap()
{
#ifdef EIGHTBIT
    static XColor Col;
    register int i,j;

    if( Monochrome )
    {   for( i=0; i<LutSize; i++ )
            Intensity[i] = (Byte)((int)(20*RLut[i]+32*GLut[i]+12*BLut[i])>>6);
        return;
    }

    if( LocalMap )
    {   XSetWindowColormap(dpy,MainWin,cmap);
        XSetWindowColormap(dpy,CanvWin,cmap);
        XUninstallColormap(dpy,lmap);
        XFreeColormap(dpy,lmap);
        LocalMap = False;
    } else if( IdentCount )
        XFreeColors(dpy,cmap,Ident,IdentCount,(long)0);
    IdentCount = 0;


    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            Col.flags = DoRed | DoGreen | DoBlue;
            if( !XAllocColor(dpy,cmap,&Col) )
                break;
            Ident[IdentCount++] = Col.pixel;
            Lut[i] = Col.pixel;
        }

    if( i<LutSize )
    {   lmap = XCopyColormapAndFree(dpy,cmap);
        LocalMap = True;

        for( j=0; j<5; j++ )
        {   Col.red   = RLut[j]<<8 | RLut[j];
            Col.green = GLut[j]<<8 | GLut[j];
            Col.blue  = BLut[j]<<8 | BLut[j];
            XAllocColor(dpy,cmap,&Col);
            Lut[i] = Col.pixel;
        }

        for( j=i; j<LutSize; j++ )
            if( ULut[j] )
            {   Col.red   = RLut[j]<<8 | RLut[j];
                Col.green = GLut[j]<<8 | GLut[j];
                Col.blue  = BLut[j]<<8 | BLut[j];
                XAllocColor(dpy,lmap,&Col);
                Lut[j] = Col.pixel;
            }
        XSetWindowColormap(dpy,MainWin,lmap);
        XSetWindowColormap(dpy,CanvWin,lmap);
        XInstallColormap(dpy,lmap);
    }
#else
    static XColor Col;
    static ByteTest buf;
    register Byte temp;
    register int i;

    for( i=0; i<LutSize; i++ )
        if( ULut[i] )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            XAllocColor(dpy,cmap,&Col);
#ifdef THIRTYTWOBIT
            if( SwapBytes )
            {   buf.longword = (Long)Col.pixel;
                temp = buf.bytes[0];
                buf.bytes[0] = buf.bytes[3];
                buf.bytes[3] = temp;

                temp = buf.bytes[1];
                buf.bytes[1] = buf.bytes[2];
                buf.bytes[2] = temp;
                Lut[i] = buf.longword;
            } else Lut[i] = (Long)Col.pixel;
#else
            Lut[i] = Col.pixel;
#endif
       }
#endif
    XSetWindowBackground(dpy,CanvWin,(unsigned long)Lut[5]);
}


static void OpenCanvas( x, y )
    int x, y;
{
    register unsigned long mask;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.cursor = cross;                           mask |= CWCursor;
    attr.background_pixel = Lut[0];                mask |= CWBackPixel;

    CanvWin = XCreateWindow(dpy, MainWin, 14, MenuHigh+14, x, y, 0, 
                            CopyFromParent, InputOutput, vis, mask, &attr );
}


static void OpenFonts()
{
    static char *fontname[] = { "-*-helvetica-bold-o-normal-*-14-*",
                                     "-*-serf-bold-o-normal-*-14-*",
                                        "-*-*-bold-o-normal-*-14-*" };
    register int i;

    cross = XCreateFontCursor(dpy,XC_tcross);
    arrow = XCreateFontCursor(dpy,XC_top_left_arrow);

    for( i=0; i<3; i++ )
        if( (MenuFont=XLoadQueryFont(dpy,fontname[i])) ) 
            break;

    if( !MenuFont )
        FatalGraphicsError("Unable to find suitable font");
    FontHigh = MenuFont->max_bounds.descent +
               MenuFont->max_bounds.ascent + 1;
    MenuHigh = FontHigh+6;
}


static void OpenCursors()
{
    Pixmap source,mask;
    XColor black,white;

    white.red = 65535;     black.red = 0;      
    white.green = 65535;   black.green = 0;
    white.blue = 65535;    black.blue = 0;
     
    white.flags = DoRed | DoGreen | DoBlue;
    black.flags = DoRed | DoGreen | DoBlue;

    source = XCreateBitmapFromData(dpy,MainWin,(char*)HGlassData,16,16);
    mask   = XCreateBitmapFromData(dpy,MainWin,(char*)HGlassMask,16,16);
    hglass = XCreatePixmapCursor(dpy,source,mask,&black,&white,7,7);
}


static void OpenColourMap()
{
    static XColor Col;
    register int i;

#ifdef EIGHTBIT
    if( !Monochrome )
    {   Col.flags = DoRed | DoGreen | DoBlue;

        for( i=0; i<5; i++ )
        {   Col.red   = RLut[i]<<8 | RLut[i];
            Col.green = GLut[i]<<8 | GLut[i];
            Col.blue  = BLut[i]<<8 | BLut[i];
            if( !XAllocColor(dpy,cmap,&Col) )
            {   cmap = XCopyColormapAndFree(dpy,cmap);
                XAllocColor(dpy,cmap,&Col);
            } 
            Lut[i] = Col.pixel;
        }
        Lut[5] = Lut[0];
    } else /* Black & White */
    {   Lut[0] = BlackCol;
        Lut[1] = BlackCol;
        Lut[2] = WhiteCol;
        Lut[3] = BlackCol;
        Lut[4] = WhiteCol;

        Intensity[5] = 0;
        Lut[5] = 5;         
    }

    LocalMap = False;
    IdentCount = 0;
#else
    Col.flags = DoRed | DoGreen | DoBlue;

    for( i=0; i<5; i++ )
    {   Col.red   = RLut[i]<<8 | RLut[i];
        Col.green = GLut[i]<<8 | GLut[i];
        Col.blue  = BLut[i]<<8 | BLut[i];
        XAllocColor(dpy,cmap,&Col);
        Lut[i] = Col.pixel;
    }
    Lut[5] = Lut[0];
#endif
}


static int RegisterInterpName( name )
    char *name;
{
    static unsigned char *registry;
    static unsigned long len,left;
    static char buffer[32];
    static int format;
    static Atom type;

    register int result;
    register char *ptr;

    registry = NULL;
    result = XGetWindowProperty(dpy, RootWindow(dpy,0), InterpAtom,
                                0, 100000, False, XA_STRING, &type,
                                &format, &len, &left, &registry );

    if( (result!=Success) || (format!=8) || (type!=XA_STRING) )
    {   if( (type!=None) && registry ) XFree( (char*)registry );

        sprintf(buffer,"%x %s",(int)MainWin,name);
        XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                         8, PropModeReplace, (unsigned char*)buffer, 
                         strlen(buffer)+1 );
        return( True );
    }

    ptr = (char*)registry;
    while( *ptr )
    {   /* Skip Window ID */
        while( *ptr++ != ' ' )
            if( !*ptr ) break;

        /* Compare Interp Name */
        if( !strcmp(ptr,name) )
        {   XFree( (char*)registry );
            return(False);
        }

        while( *ptr++ );
    }

    XFree( (char*)registry );
    sprintf(buffer,"%x %s",(int)MainWin,name);
    XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING, 
                     8, PropModeAppend, (unsigned char*)buffer, 
                     strlen(buffer)+1 );
    return( True );
}


static void DeRegisterInterpName( name )
    char *name;
{
    static unsigned char *registry;
    static unsigned long len,left;
    static int format;
    static Atom type;

    register char *src, *dst;
    register int result;

    registry = NULL;
    result = XGetWindowProperty(dpy, RootWindow(dpy,0), InterpAtom,
                                0, 100000, False, XA_STRING, &type,
                                &format, &len, &left, &registry );
    if( type==None )
        return;

    if( (result!=Success) || (format!=8) || (type!=XA_STRING) )
    {   XDeleteProperty( dpy, RootWindow(dpy,0), InterpAtom );
        if( registry ) XFree( (char*)registry );
        return;
    }

    dst = (char*)registry;
    while( *dst )
    {   /* Skip Window ID */
        src = dst;
        while( *src++ != ' ' )
            if( !*src ) break;

        /* Compare Interp Name */
        if( strcmp(src,name) )
        {   while( *dst++ );
        } else break;
    }

    if( *dst )
    {   /* Skip Interp Name */
        while( *src++ );
        
        /* Shuffle Registry */
        while( *src )
            while( (*dst++ = *src++) );
        *dst = 0;

        XChangeProperty( dpy, RootWindow(dpy,0), InterpAtom, XA_STRING,
                         8, PropModeReplace, registry, dst-(char*)registry );
    }
    XFree( (char*)registry );
}


static void OpenIPCComms()
{
    auto char buffer[16];
    register int i;

    CommAtom = XInternAtom( dpy, "Comm", False );
    InterpAtom = XInternAtom( dpy, "InterpRegistry", False );
    AppNameAtom = XInternAtom(dpy, "TK_APPLICATION", False );
    DelWinXAtom = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
    /* XSetWMProtocols(dpy,MainWin,&DelWinXAtom,True); */
    if( (ProtoXAtom = XInternAtom(dpy,"WM_PROTOCOLS",False)) )
        XChangeProperty( dpy, MainWin, ProtoXAtom, XA_ATOM, 32, 
                        PropModeReplace, (Byte*)&DelWinXAtom, True );

    i = 0;
    XGrabServer( dpy );
    if( !RegisterInterpName("rasmol") )
    {   strcpy(TkInterp,"rasmol #0");
        for( i=1; i<10; i++ )
        {    TkInterp[8] = i+'0';
             if( RegisterInterpName(TkInterp) )
                 break;
        }

        if( i < 10 ) 
        {   /* Tk4.0 and later! */
            strcpy(buffer,"{rasmol #0}");  buffer[9] = i+'0';
            XChangeProperty( dpy, MainWin, AppNameAtom, XA_STRING, 
                             8, PropModeReplace, (Byte*)buffer, 12 );
        } else *TkInterp = 0;
    } else 
    {   XChangeProperty( dpy, MainWin, AppNameAtom, XA_STRING,
                         8, PropModeReplace, (Byte*)"rasmol", 7 );
        strcpy(TkInterp,"rasmol");
    }
    XUngrabServer( dpy );
}


static void DrawUpBox( wdw, x1, y1, x2, y2 )
    Drawable wdw;  int x1,y1,x2,y2;
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[3]);
    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x1,y1,x1,y2);
    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,lx,ly,lx,uy);

    XSetForeground(dpy,gcon,(unsigned long)Lut[1]);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x2,y2);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,ux,uy);
}


static void DrawDnBox( wdw, x1, y1, x2, y2 )
    Drawable wdw;  int x1,y1,x2,y2;
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[1]);
    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x1,y1,x1,y2);
    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,lx,ly,lx,uy);

    XSetForeground(dpy,gcon,(unsigned long)Lut[3]);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x2,y2);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,ux,uy);
}


static void DrawNoBox( wdw, x1, y1, x2, y2 )
    Drawable wdw;  int x1,y1,x2,y2;
{
    register int lx,ly,ux,uy;

    lx = x1+1;  ly = y1+1;
    ux = x2-1;  uy = y2-1;

    XSetForeground(dpy,gcon,(unsigned long)Lut[2]);

    XDrawLine(dpy,wdw,gcon,x1,y1,x2,y1);
    XDrawLine(dpy,wdw,gcon,x2,y1,x2,y2);
    XDrawLine(dpy,wdw,gcon,x2,y2,x1,y2);
    XDrawLine(dpy,wdw,gcon,x1,y2,x1,y1);

    XDrawLine(dpy,wdw,gcon,lx,ly,ux,ly);
    XDrawLine(dpy,wdw,gcon,ux,ly,ux,uy);
    XDrawLine(dpy,wdw,gcon,ux,uy,lx,uy);
    XDrawLine(dpy,wdw,gcon,lx,uy,lx,ly);
}


static void OpenMenuBar()
{
    register unsigned long mask;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonReleaseMask;
    MenuWin = XCreateWindow( dpy, MainWin, 2, 2, XRange+49, FontHigh+5, 0,
                             CopyFromParent, InputOnly, vis, mask, &attr );


    /* Create Unmapped PopUp Window! */
    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonReleaseMask | 
                      KeyPressMask;
    attr.background_pixel = Lut[2];     mask |= CWBackPixel;
    attr.border_pixel = Lut[2];         mask |= CWBorderPixel;
    attr.override_redirect = True;      mask |= CWOverrideRedirect;
    attr.save_under = True;             mask |= CWSaveUnder;
    attr.colormap = cmap;               mask |= CWColormap;

    PopUpWin = XCreateWindow(dpy, RootWin, 0, 0, 100, 100, 0, 
                             PixDepth, InputOutput, vis,
                             mask, &attr );
    MenuFocus = False;
    PopUpFlag = False;
}


static void OpenScrollBars()
{
    register unsigned long mask;

    Scrl = XCreatePixmap( dpy, MainWin, 16, 16, PixDepth );
    XSetForeground(dpy,gcon,(unsigned long)Lut[2]); 
    XFillRectangle(dpy,Scrl,gcon,0,0,15,15);
    XSetForeground(dpy,gcon,(unsigned long)Lut[0]); 
    XDrawRectangle(dpy,Scrl,gcon,0,0,15,15);
    DrawUpBox( Scrl, 1, 1, 14, 14 );

    tilepix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)ScrlTile, 8, 8,
                                           (unsigned long)Lut[0], 
                                           (unsigned long)Lut[2], PixDepth );

    mask = CWEventMask;
    attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                    | ButtonReleaseMask;
    attr.background_pixmap = tilepix;              mask |= CWBackPixmap;

    XScrlWin = XCreateWindow(dpy,MainWin,14,YRange+MenuHigh+24,XRange,16, 
                             0,CopyFromParent,InputOutput,vis,mask,&attr);
    lfpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)LfArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );
    rgpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)RgArrow, 16, 16, 
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );


    YScrlWin = XCreateWindow(dpy,MainWin,XRange+24,MenuHigh+14,16,YRange, 
                             0,CopyFromParent,InputOutput,vis,mask,&attr);
    uppix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)UpArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );
    dnpix = XCreatePixmapFromBitmapData( dpy, MainWin, (char*)DnArrow, 16, 16,
                                         (unsigned long)Lut[0], 
                                         (unsigned long)Lut[2], PixDepth );

    ScrlX = (XRange/2)-8;
    ScrlY = (YRange/2)-8;
}


static void DrawXScroll()
{
    XCopyArea(dpy,rgpix,XScrlWin,gcon,0,0,16,16,XRange-16,0);
    XCopyArea(dpy,Scrl ,XScrlWin,gcon,0,0,16,16,ScrlX,0);
    XCopyArea(dpy,lfpix,XScrlWin,gcon,0,0,16,16,0,0);
}


static void DrawYScroll()
{
    XCopyArea(dpy,dnpix,YScrlWin,gcon,0,0,16,16,0,YRange-16);
    XCopyArea(dpy,Scrl ,YScrlWin,gcon,0,0,16,16,0,ScrlY);
    XCopyArea(dpy,uppix,YScrlWin,gcon,0,0,16,16,0,0);
}


void UpdateScrollBars()
{
    register int temp;

    temp = (DialValue[YScrlDial]+1.0)*(YRange-48);  
    NewScrlY = (temp>>1)+16;

    if( NewScrlY != ScrlY )
    {   XClearArea(dpy,YScrlWin,0,ScrlY,16,16,False);
        XCopyArea(dpy,Scrl,YScrlWin,gcon,0,0,16,16,0,NewScrlY);
        ReDrawFlag |= (1<<YScrlDial);
        ScrlY = NewScrlY; 
    }

    temp = (DialValue[XScrlDial]+1.0)*(XRange-48);  
    NewScrlX = (temp>>1)+16;

    if( NewScrlX != ScrlX )
    {   XClearArea(dpy,XScrlWin,ScrlX,0,16,16,False);
        XCopyArea(dpy,Scrl,XScrlWin,gcon,0,0,16,16,NewScrlX,0);
        ReDrawFlag |= (1<<XScrlDial);
        ScrlX = NewScrlX;
    }
    XFlush(dpy);
}



#ifdef DIALBOX
static void SetDialLabel( num, ptr )
    int num; char *ptr;
{
    static XStringFeedbackControl ctrl;
    static KeySym text[8];
    register int length;

    length = 0;
    while( *ptr )
       text[length++] = *ptr++;

    ctrl.id = num;
    ctrl.num_keysyms = length;
    ctrl.class = ValuatorClass;
    ctrl.syms_to_display = text;
    XChangeFeedbackControl(dpy,Dials,DvString,
                           (XFeedbackControl*)&ctrl);
}


static void GetDialState()
{
    register XValuatorState *ptr;
    register XDeviceState *stat;
    register int i,j,max;

    stat = XQueryDeviceState(dpy,Dials);
    ptr = (XValuatorState*)stat->data;
    for( i=0; i<stat->num_classes; i++ )
    {   if( ptr->class == ValuatorClass )
        {   if( ptr->mode & 0x01 )
            {   DialMode = Absolute;
                max = MinFun(ptr->num_valuators,8);
                for( j=0; j<max; j++ )
                    DialPrev[j] = ptr->valuators[j];
            } else DialMode = Relative;
            break;
        } else ptr = (XValuatorState*)(((char*)ptr) + 
                                       ptr->length);
    }
    XFreeDeviceState(stat);
}


static void OpenDialsBox()
{
    register XValuatorInfo *valptr;
    register XFeedbackState *list;
    register XFeedbackState *feed;
    register XDeviceInfo *devlist;
    register XDeviceInfo *ptr;
    register Atom devtype;
    register int i,j,max;

    static XEventClass dclass;
    static int count;


    UseDials = False;
    /* Avoid X Server's without the extension */
    if( !XQueryExtension(dpy,"XInputExtension",
                         &count,&count,&count) )
        return;
    
    devlist = XListInputDevices(dpy,&count);
    devtype = XInternAtom(dpy,XI_KNOB_BOX,True );
    if( (devtype==None) || !devlist ) return;

    ptr = devlist;
    for( i=0; i<count; i++ )
        if( (ptr->use==IsXExtensionDevice) && (ptr->type==devtype) )
        {   valptr = (XValuatorInfo*)ptr->inputclassinfo;
            for( j=0; j<ptr->num_classes; j++ )
            {   if( valptr->class == ValuatorClass )
                    if( (Dials=XOpenDevice(dpy,ptr->id)) )
                    {   UseDials = True;
                        break;
                    }
                valptr = (XValuatorInfo*)(((char*)valptr) +
                                          valptr->length);
            }
            if( UseDials ) break;
        } else ptr++;
    /* XFreeDeviceList(devlist); */

    if( UseDials ) 
    {   /* Determine Dial Mapping! */
        if( !strcmp(ServerVendor(dpy),"Silicon Graphics") )
        {      DialMap = SGIDialMap;
        } else DialMap = ESVDialMap;

        DialMode = valptr->mode;
        max = MinFun(valptr->num_axes,8);
        for( i=0; i<max; i++ )
            DialRes[i] = (Real)valptr->axes[i].resolution;
        GetDialState();
    } else return;

    UseDialLEDs = 0;
    feed = list = XGetFeedbackControl( dpy, Dials, &count );
    for( i=0; i<count; i++ )
    {   if( feed->class == StringFeedbackClass ) UseDialLEDs++;
        feed = (XFeedbackState*)(((char*)feed) + feed->length);
    }
    XFreeFeedbackList( list );

    if( UseDialLEDs >= 8 )
    {   for( i=0; i<8; i++ )
            SetDialLabel(i,DialLabel[DialMap[i]]);
    } else UseDialLEDs = False;

    DeviceMotionNotify( Dials, DialEvent, dclass );
    XSelectExtensionEvent( dpy, MainWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, MenuWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, CanvWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, XScrlWin, &dclass, 1 );
    XSelectExtensionEvent( dpy, YScrlWin, &dclass, 1 );
}


static void HandleDialEvent( ptr )
    XDeviceMotionEvent *ptr;
{
    register double temp;
    register int count;
    register int value;
    register int index;
    register int num;

    /* Limit Number of Dials */
    count = 8 - ptr->first_axis;
    if( count > (int)ptr->axes_count )
        count = (int)ptr->axes_count;

    for( index=0; index<count; index++ )
    {   num = ptr->first_axis+index;
        if( DialMode == Absolute )
        {   value = ptr->axis_data[index] - DialPrev[num];
            DialPrev[num] = ptr->axis_data[index];
        } else value = ptr->axis_data[index];

        if( value )
        {   temp = (Real)value/DialRes[num];
            num = DialMap[num];
            temp += DialValue[num];
            ReDrawFlag |= (1<<num);

            if( num<3 )
            {   while( temp<-1.0 ) temp += 2.0;
                while( temp>1.0 )  temp -= 2.0;
            } else
            {   if( temp<-1.0 ) temp = -1.0;
                if( temp>1.0 )  temp = 1.0;
            }
            DialValue[num] = temp;

            if( num==YScrlDial )
            {   value = (temp+1.0)*(YRange-48);
                NewScrlY = (value>>1)+16;
            }

            if( num==XScrlDial )
            {   value = (temp+1.0)*(XRange-48);
                NewScrlX = (value>>1)+16;
            }
        }
    }
}
#endif


static void DrawMainWin()
{
    register int temp;

    DrawUpBox(MainWin,0,0,MainWide,MainHigh);
    DrawUpBox(MainWin,0,0,MainWide-2,FontHigh+7);

    temp = YRange+MenuHigh;
    DrawDnBox(MainWin,12,MenuHigh+12,XRange+16,temp+16);
    DrawDnBox(MainWin,XRange+22,MenuHigh+12,XRange+41,temp+16);
    DrawDnBox(MainWin,12,temp+22,XRange+16,temp+41);
}


/********************/
/* Menu Bar Display */
/********************/


static void DisplayMenuBarText( ptr, x, y )
    BarItem *ptr;  int x, y;
{
    register unsigned long col;
    register int under,wide;

    if( ptr->flags&mbEnable && !DisableMenu )
    {      col = Lut[0];
    } else col = Lut[1];
    XSetForeground( dpy, gcon, col );

    XDrawString( dpy, MainWin, gcon, x, y, ptr->text, ptr->len );

    under = y + MenuFont->descent;
    wide = XTextWidth( MenuFont, ptr->text, 1 );
    XDrawLine( dpy, MainWin, gcon, x, under, x+wide, under );
}


static void DrawMenuBar()
{
    register BarItem *ptr;
    register int wide;
    register int x,y;
    register int i;

    x = 6; y = MenuFont->ascent+4;
    XSetFont( dpy, gcon, MenuFont->fid );

    for( i=0; i<MenuBarMax; i++ )
    {   ptr = MenuBar+i;
        wide = XTextWidth( MenuFont, ptr->text, ptr->len );
        if( x+wide+24 > MainWide ) break;

        /* Right Justify "Help" */
        if( i == MenuBarMax-1 )
            x = MainWide - (wide+24);

        DisplayMenuBarText( ptr, x+8, y );

        if( MenuFocus && (i==MenuBarSelect) )
        {      DrawUpBox( MainWin, x, 2, x+wide+16, FontHigh+5 );
        } else DrawNoBox( MainWin, x, 2, x+wide+16, FontHigh+5 );
        x += wide+24;
    }
    MenuBarCount = i;
    /* XSync(dpy,False); */
    XFlush(dpy);
}


/***********************/
/* Pop-up Menu Display */
/***********************/


static void DisplayPopUpText( ptr, x, y )
    MenuItem *ptr; int x, y;
{
    register unsigned long col;
    register int pos, wide;
    register int i,under;
    register int index;

    col = (ptr->flags&mbEnable)? Lut[0] : Lut[1];
    XSetForeground( dpy, gcon, col );

    XDrawString( dpy, PopUpWin, gcon, x, y, ptr->text, ptr->len );

    if( ptr->flags & mbAccel )
    {   under = y + MenuFont->descent;

        pos = x;
        for( i=0; i<ptr->pos; i++ )
        {   index = ptr->text[i] - MenuFont->min_char_or_byte2;
            pos += MenuFont->per_char[index].width;
        }

        index = ptr->text[ptr->pos] - MenuFont->min_char_or_byte2;
        wide = pos+MenuFont->per_char[index].rbearing;
        pos += MenuFont->per_char[index].lbearing;

        XDrawLine( dpy, PopUpWin, gcon, pos, under, wide, under );
    }
}


static void DrawPopUpMenu()
{
    register MenuItem *ptr;
    register int count;
    register int x,y;
    register int i;

    DrawUpBox(PopUpWin,0,0,PopUpWide,PopUpHigh);

    ptr = MenuBar[MenuBarSelect].menu;
    count = MenuBar[MenuBarSelect].count;

    y = 2;  x = 2;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   DisplayPopUpText( ptr, x+8, y+MenuFont->ascent+2 );

            if( ItemFlag && (i==MenuItemSelect) )
            {      DrawUpBox(PopUpWin,2,y,PopUpWide-2,y+FontHigh+3);
            } else DrawNoBox(PopUpWin,2,y,PopUpWide-2,y+FontHigh+3);
            y += FontHigh+4;
        } else
        {   XSetForeground( dpy, gcon, (unsigned long)Lut[1] );
            XDrawLine(dpy,PopUpWin,gcon,2,y,PopUpWide-2,y);
            XSetForeground( dpy, gcon, (unsigned long)Lut[3] );
            XDrawLine(dpy,PopUpWin,gcon,2,y+1,PopUpWide-2,y+1);
            y += 2;
        }
        ptr++;
    }
    /* XSync(dpy,False); */
    XFlush(dpy);
}


static void DisplayPopUpMenu( i, x )
    int i, x;
{
    register int wide, count;
    register MenuItem *ptr;
    register int flag;

    static int xpos, ypos;
    static Window win;


    MenuBarSelect = i;
    DrawMenuBar();

    ptr = MenuBar[i].menu;
    count = MenuBar[i].count;

    flag = False;
    PopUpHigh = 4;
    PopUpWide = 4;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   if( ptr->flags & mbOption ) flag = True;
            wide = XTextWidth(MenuFont,ptr->text,ptr->len);
            if( wide+28 > PopUpWide ) PopUpWide = wide+28;
            PopUpHigh += FontHigh+4;
        } else PopUpHigh += 2;
        ptr++;
    }

    /* Determine pop-up menu position! */
    XTranslateCoordinates(dpy,MainWin,RootWin,x,FontHigh+6,
                          &xpos, &ypos, &win );

    if( ypos+PopUpHigh > MaxHeight )
        ypos -= (PopUpHigh+FontHigh+6);
    if( xpos+PopUpWide > MaxWidth )
        xpos = MaxWidth-PopUpWide;
    if( xpos < 0 ) xpos = 0;

    XUnmapWindow(dpy,PopUpWin);
    XMoveResizeWindow(dpy,PopUpWin,xpos,ypos,PopUpWide+1,PopUpHigh+1);
    XRaiseWindow(dpy,PopUpWin);
    XMapWindow(dpy,PopUpWin);
    PopUpFlag = True;
    DrawPopUpMenu();
}


/******************************/
/* Pop-Up Menu Event Handling */
/******************************/


static void HandleItemClick( x, y )
    int x, y;
{
    register MenuItem *ptr;
    register int count,i;

    static int xpos, ypos;
    static Window win;

    XTranslateCoordinates(dpy,MenuWin,PopUpWin,x,y,
                          &xpos,&ypos,&win);


    /* Ignore by not setting ItemFocus! */
    if( (xpos<0) || (xpos>PopUpWide) ) return;
    if( (ypos<0) || (ypos>PopUpHigh) ) return;
    ItemFocus = True;

    ptr = MenuBar[MenuBarSelect].menu;
    count = MenuBar[MenuBarSelect].count;

    y = 2;
    for( i=0; i<count; i++ )
    {   if( !(ptr->flags&mbSepBar) )
        {   if( (ypos>=y) && (ypos<=y+FontHigh+3) )
            {   if( ptr->flags & mbEnable )
                {   if( !ItemFlag || (MenuItemSelect!=i) )
                    {   /* Avoid Flickering */
                        MenuItemSelect = i;
                        ItemFlag = True;
                        DrawPopUpMenu();
                    }
                    return;
                } else break;
            }
            y += FontHigh+4;
        } else y += 2;
        ptr++;
    }

    if( ItemFlag )
    {   ItemFlag = False;
        DrawPopUpMenu();
    }
}


static void HandleItemMove( x, y )
    int x, y;
{
    register MenuItem *ptr;
    register int count,i;

    static int xpos, ypos;
    static Window win;

    XTranslateCoordinates(dpy,MenuWin,PopUpWin,x,y,
                          &xpos,&ypos,&win);

    if( (xpos>=0) && (xpos<=PopUpWide) )
    {   ptr = MenuBar[MenuBarSelect].menu;
        count = MenuBar[MenuBarSelect].count;

        y = 2;
        for( i=0; i<count; i++ )
        {   if( !(ptr->flags&mbSepBar) )
            {   if( (ypos>=y) && (ypos<=y+FontHigh+3) )
                {   if( !ItemFlag || (MenuItemSelect!=i) )
                    {   /* Avoid Flicker! */
                        MenuItemSelect = i;
                        ItemFlag = True;
                        DrawPopUpMenu();
                    }
                    ItemFocus = True;
                    return;
                }
                y += FontHigh+4;
            } else y += 2;
            ptr++;
        }
    }

    if( ItemFlag )
    {   /* Avoid Flicker! */
        ItemFlag = False;
        DrawPopUpMenu();
    }
}


static int HandleItemKey( key )
    int key;
{
    register MenuItem *ptr;
    register int count;
    register int item;
    register int ch;
    register int i;

    key = ToUpper( key );
    item = MenuItemSelect;
    ptr = &MenuBar[MenuBarSelect].menu[item];
    count = MenuBar[MenuBarSelect].count;
    for( i=0; i<count; i++ )
    {   if( (ptr->flags&(mbEnable|mbAccel)) && 
           !(ptr->flags&mbSepBar) )
        {   ch = ptr->text[ptr->pos];
            if( ToUpper(ch) == key )
                return( (MenuBarSelect<<8)+item+1 );
        }

        /* Advance to next item! */
        if( item == count-1 )
        {   ptr = MenuBar[MenuBarSelect].menu;
            item = 0;
        } else 
        {   item++;
            ptr++;
        }
    }
    return( 0 );
}


static void SelectFirstItem( menu )
    int menu;
{
    register MenuItem *ptr;
    register int count;
    register int i;

    count = MenuBar[menu].count;
    ptr = MenuBar[menu].menu;

    ItemFlag = False;
    for( i=0; i<count; i++ )
        if( (ptr->flags&mbEnable) &&
           !(ptr->flags&mbSepBar) )
        {   MenuItemSelect = i;
            ItemFlag = True;
            break;
        } else ptr++;
}


static void SelectPrevItem()
{
    register BarItem *ptr;
    register int flags;
    register int item;
    register int i;

    if( !ItemFlag )
        return;

    item = MenuItemSelect;
    ptr = MenuBar + MenuBarSelect;
    for( i=0; i<ptr->count; i++ )
    {   if( !item )
        {   item = ptr->count-1;
        } else item--;

        flags = ptr->menu[item].flags;
        if( (flags&mbEnable) && !(flags&mbSepBar) )
            break;
    }

    if( item != MenuItemSelect )
    {   MenuItemSelect = item;
        DrawPopUpMenu();
    }
}


static void SelectNextItem()
{
    register BarItem *ptr;
    register int flags;
    register int item;
    register int i;

    if( !ItemFlag )
        return;

    item = MenuItemSelect;
    ptr = MenuBar + MenuBarSelect;
    for( i=0; i<ptr->count; i++ )
    {   if( item == ptr->count-1 )
        {   item = 0;
        } else item++;

        flags = ptr->menu[item].flags;
        if( (flags&mbEnable) && !(flags&mbSepBar) )
            break;
    }

    if( item != MenuItemSelect )
    {   MenuItemSelect = item;
        DrawPopUpMenu();
    }
}


/***************************/
/* Menu Bar Event Handling */
/***************************/


static void SelectMenu( menu )
    int menu;
{
    register BarItem *ptr;
    register int wide;
    register int i,x;


    if( !PopUpFlag )
    {   MenuBarSelect = menu;
        DrawMenuBar();
        return;
    }

    if( menu != MenuBarMax-1 )
    {   x = 6;
        for( i=0; i<menu; i++ )
        {   ptr = MenuBar+i;
            wide = XTextWidth(MenuFont,ptr->text,ptr->len);
            x += wide+24;
        }
    } else 
    {   ptr = MenuBar+menu;
        wide = XTextWidth(MenuFont,ptr->text,ptr->len);
        x = MainWide - (wide+24);
    }

    SelectFirstItem( menu );
    DisplayPopUpMenu( menu, x );
    ItemFocus = False;
}


static int HandleMenuClick( pos )
    int pos;
{
    register BarItem *ptr;
    register int wide;
    register int x,i;

    x = 6;
    for( i=0; i<MenuBarCount; i++ )
    {   ptr = MenuBar+i;
        wide = XTextWidth( MenuFont, ptr->text, ptr->len );
        if( i == MenuBarMax-1 ) x = MainWide - (wide+24);

        if( (pos>=x) && (pos<=x+wide+16) )
        {   if( !PopUpFlag || (MenuBarSelect!=i) )
            {   ItemFlag = False;
                DisplayPopUpMenu(i,x);
            } else if( ItemFlag )
            {   ItemFlag = False;
                DrawPopUpMenu();
            }
            ItemFocus = True;
            return( True );
        } else x += wide+24;
    }
    return( False );
}


static int HandleMenuKey( key )
    char key;
{
    register int i;

    key = ToUpper(key);
    for( i=0; i<MenuBarCount; i++ )
        if( MenuBar[i].text[0] == key )
        {   if( !PopUpFlag || (MenuBarSelect!=i) )
            {   PopUpFlag = True;
                SelectMenu( i );
            }
            return( True );
        }
    return( False );
}


void EnableMenus( flag )
    int flag;
{
    DisableMenu = !flag;
    if( Interactive )
        DrawMenuBar();
}


static void ReSizeWindow( wide, high )
    int wide, high;
{
    register Real xpos;
    register Real ypos;
    register int dx;

    xpos = (XRange>48)? (Real)(ScrlX-16)/(XRange-48) : 0.0;
    ypos = (YRange>48)? (Real)(ScrlY-16)/(YRange-48) : 0.0;

    YRange = high-(MenuHigh+53);
    XRange = wide-53;

    if( (dx = XRange%4) )
        XRange += 4-dx;

    MainHigh = YRange+(MenuHigh+53);  HRange = YRange>>1;
    MainWide = XRange+53;             WRange = XRange>>1;
    Range = MinFun(XRange,YRange);

    XResizeWindow( dpy, CanvWin, XRange, YRange);
    XResizeWindow( dpy, MenuWin, XRange+49, FontHigh+5 );
    XMoveResizeWindow( dpy, XScrlWin, 14, YRange+MenuHigh+24, XRange, 16 );
    XMoveResizeWindow( dpy, YScrlWin, XRange+24, MenuHigh+14, 16, YRange );

    NewScrlX = ScrlX = (xpos*(XRange-48))+16;
    NewScrlY = ScrlY = (ypos*(YRange-48))+16;

    XClearWindow( dpy, MainWin );
    XClearWindow( dpy, CanvWin );

    DrawXScroll();
    DrawYScroll();
    DrawMainWin();
    DrawMenuBar();

    ReDrawFlag |= RFReSize;
    XSync(dpy,True);
}


int FatalXError( ptr )
    Display *ptr;
{
    dpy = (Display*)NULL;
    RasMolFatalExit("*** Fatal X11 I/O Error! ***");
    /* Avoid Compilation Warnings! */
    return( (int)ptr );
}


int OpenDisplay( x, y )
    int x, y;
{
    register unsigned long mask;
    register int i,num;
    register char *ptr;
 
#ifdef THIRTYTWOBIT
    static ByteTest test;
#endif
    static XVisualInfo visinfo;
    static XClassHint xclass;
    static XSizeHints size;
    static Pixmap icon;
    static int temp;


    image = (XImage*)NULL;

    MouseMode = MMRasMol;
    UseHourGlass = True;
    DisableMenu = False;
    Monochrome = False;
    HeldButton = -1;

    for( i=0; i<8; i++ )
         DialValue[i] = 0.0;

    RLut[0]=0;   GLut[0]=0;   BLut[0]=0;    ULut[0]=True;
    RLut[1]=100; GLut[1]=100; BLut[1]=100;  ULut[1]=True;
    RLut[2]=150; GLut[2]=150; BLut[2]=150;  ULut[2]=True;
    RLut[3]=200; GLut[3]=200; BLut[3]=200;  ULut[3]=True;
    RLut[4]=255; GLut[4]=255; BLut[4]=255;  ULut[4]=True;

    XRange = x;  WRange = XRange>>1;
    YRange = y;  HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    if( !Interactive ) return( False );
    if( (dpy=XOpenDisplay(NULL)) == NULL )
        return( 0 );

    num = DefaultScreen(dpy);
    RootWin = RootWindow(dpy,num);
    XSetIOErrorHandler( FatalXError );

#ifdef EIGHTBIT
    if( !(XMatchVisualInfo(dpy,num,8,PseudoColor,&visinfo) ||
          XMatchVisualInfo(dpy,num,8,GrayScale,&visinfo)) )
    {   /* Attempt to use Monochrome Mode! */
        if( !(XMatchVisualInfo(dpy,num,1,StaticColor,&visinfo) ||
              XMatchVisualInfo(dpy,num,1,StaticGray,&visinfo)) )
        {   XCloseDisplay(dpy);
            return( 0 );
        }
        Monochrome = True;
        PixDepth = 1;
    } else PixDepth = 8;
#else
#ifdef THIRTYTWOBIT
    if( XMatchVisualInfo(dpy,num,32,TrueColor,&visinfo) ||
        XMatchVisualInfo(dpy,num,32,DirectColor,&visinfo) )
    {   PixDepth = 32;
    } else if( XMatchVisualInfo(dpy,num,24,TrueColor,&visinfo) ||
               XMatchVisualInfo(dpy,num,24,DirectColor,&visinfo) )
    {   PixDepth = 24;
    } else /* No suitable display! */
    {   XCloseDisplay(dpy);
        return(0);
    }
#else /* SIXTEENBIT */
    if( !XMatchVisualInfo(dpy,num,16,TrueColor,&visinfo) &&
        !XMatchVisualInfo(dpy,num,16,DirectColor,&visinfo) )
    {   XCloseDisplay(dpy);
        return(0);
    } else PixDepth = 16;
#endif
#endif

    if( !Monochrome )
    {   vis = visinfo.visual;
        if( vis != DefaultVisual(dpy,num) )
        {   cmap = XCreateColormap(dpy,RootWin,vis,AllocNone);
        } else cmap = DefaultColormap(dpy,num);
    } else /* Black & White */
    {   /* PixDepth = DefaultDepth(dpy,num); */
        /* vis = DefaultVisual(dpy,num);     */
        cmap = DefaultColormap(dpy,num);

        BlackCol = BlackPixel(dpy,num) & 1;
        WhiteCol = WhitePixel(dpy,num) & 1;
    }

    OpenFonts();
    OpenColourMap();

    MaxHeight = DisplayHeight(dpy,num);  MinHeight = MenuHigh+101;
    MaxWidth = DisplayWidth(dpy,num);    MinWidth = 101;

    MainHigh = YRange+MenuHigh+53;
    MainWide = XRange+53;

    mask = CWEventMask;
    attr.event_mask = ExposureMask | KeyPressMask | StructureNotifyMask
                    | EnterWindowMask | LeaveWindowMask | PropertyChangeMask;
    attr.background_pixel = Lut[2];     mask |= CWBackPixel;
    attr.border_pixel = Lut[2];         mask |= CWBorderPixel;
    attr.colormap = cmap;               mask |= CWColormap;
    attr.cursor = arrow;                mask |= CWCursor;

    MainWin = XCreateWindow(dpy, RootWin, 0, 0, MainWide, MainHigh, 2,
			    PixDepth, InputOutput, vis, mask, &attr );

    gcon = XCreateGC(dpy,MainWin,0L,NULL);
    /* DefaultGC(dpy,num) */

    XSetGraphicsExposures(dpy,gcon,False);
    icon = XCreateBitmapFromData(dpy,MainWin,(char*)icon_bits,
                                 icon_width,icon_height );

    size.flags = PMinSize | PMaxSize;
    size.min_width = MinWidth;    size.max_width = MaxWidth;
    size.min_height = MinHeight;  size.max_height = MaxHeight;
    XSetStandardProperties(dpy, MainWin, "TINKER RasMol",
                           "RasMol", icon, NULL, 0, &size );

    xclass.res_name = "rasmol";
    xclass.res_class = "RasMol";
    XSetClassHint(dpy,MainWin,&xclass);

    hints.icon_pixmap = icon;       
    hints.flags = IconPixmapHint;
    XSetWMHints(dpy,MainWin,&hints);

    OpenCanvas( XRange, YRange );
    OpenScrollBars();
    OpenMenuBar();
    OpenCursors();
    OpenIPCComms();

#ifdef DIALBOX
    OpenDialsBox();
#endif

#ifdef MITSHM
    ptr = DisplayString(dpy);
    if( !ptr || (*ptr==':') || !strncmp(ptr,"localhost:",10) || 
        !strncmp(ptr,"unix:",5) || !strncmp(ptr,"local:",6) )
    {   SharedMemOption = XQueryExtension(dpy,"MIT-SHM",&temp,&temp,&temp);
        if( Monochrome && (PixDepth!=1) ) SharedMemOption = False;
    } else SharedMemOption = False;
    SharedMemFlag = False;
#endif

#ifdef THIRTYTWOBIT
    /* Determine Byte Ordering */
    test.longword = (Long)0x000000ff;
    if( ImageByteOrder(dpy) == MSBFirst )
    {      SwapBytes = test.bytes[0];
    } else SwapBytes = test.bytes[3];
#endif

    XMapSubwindows(dpy,MainWin);
    XMapWindow(dpy,MainWin);

    DrawXScroll();
    DrawYScroll();
    DrawMainWin();
    DrawMenuBar();

    XClearWindow( dpy, CanvWin );
    XSync(dpy,False);

    num = 1<<ConnectionNumber(dpy);
    return( num );
}


int CreateImage()
{
    register int format,depth;
    register Long size, temp;
    register Pixel *ptr;

    if( !Interactive )
    {   if( FBuffer ) free(FBuffer);
        size = (Long)XRange*YRange*sizeof(Pixel);
        FBuffer = (Pixel*)malloc( size+32 );
        return( (int)FBuffer );
    }

    if( Monochrome )
    {   format = XYPixmap;
        depth = 1;
    } else /* Colour */
    {   format = ZPixmap;
        depth = PixDepth;
    }


    if( image ) 
    {   /* Monochrome Mode Frame Buffer! */
        if( FBuffer && (FBuffer!=(Pixel*)image->data) )
            free(FBuffer);
#ifdef MITSHM
        if( SharedMemFlag )
        {   XShmDetach( dpy, &xshminfo );
            image->data = (char*)NULL;
            shmdt( xshminfo.shmaddr );
        }
#endif
        XDestroyImage( image );
        image = (XImage*)NULL;
    }


    if( Monochrome )
    {   /* Monochrome Mode Frame Buffer! */
        size = (Long)XRange*YRange*sizeof(Pixel);
        FBuffer = (Pixel*)malloc( size+32 );
        if( !FBuffer ) return( False );

        /* Bit per Pixel ScanLines! */
        temp = ((XRange+31)>>5)<<2;
        size = (Long)temp*YRange + 32;
    } else 
        size = (Long)XRange*YRange*sizeof(Pixel) + 32;

#ifdef MITSHM
    if( SharedMemOption )
    {   SharedMemFlag = False;
        image = XShmCreateImage( dpy, vis, depth, format,
                                 NULL, &xshminfo, XRange, YRange );

        if( image )
        {   temp = (Long)image->bytes_per_line * image->height;
            if( temp > size ) size = temp;
            xshminfo.shmid = shmget( IPC_PRIVATE, size, IPC_CREAT|0777 );
            if( xshminfo.shmid != -1 ) 
            {   xshminfo.shmaddr = (char*)shmat(xshminfo.shmid,0,0);
                if( xshminfo.shmaddr != (char*)-1 )
                {   image->data = xshminfo.shmaddr;
                    if( !Monochrome )
                        FBuffer = (Pixel*)image->data;
                    xshminfo.readOnly = True;

                    SharedMemFlag = XShmAttach( dpy, &xshminfo );
                    XSync(dpy,False);
                }
                /* Always Destroy Shared Memory Ident */
                shmctl( xshminfo.shmid, IPC_RMID, 0 );
            }

            if( SharedMemFlag )
            {   if( Monochrome )
                {   if( BlackCol )
                    {      memset((void*)image->data,255,size);
                    } else memset((void*)image->data,255,size);
                }
                return( True );
            } else 
            {   XDestroyImage( image );
                image = (XImage*)NULL;
            }
        }
    }
#endif

    /* Allocate Frame Buffer! */
    ptr = (Pixel*)malloc( size );
    if( !ptr ) return( False );

    if( !Monochrome ) FBuffer = ptr;
    image = XCreateImage( dpy, vis, depth, format, 0, (char*)ptr, 
                          XRange, YRange, ((PixDepth>8)?32: 8) , 0 );
    return( (int)image );
}


static void DitherImage()
{
    register Card bits;
    register Card *dst;
    register Pixel *src;
    register int xmax,ymax;
    register int count,x,y;
    register int error;

    register Card bmask,wmask;
    register Card bhigh,whigh;
    register int wlen;
    register int blen;
    register int len;

    src = (Pixel*)FBuffer;
    dst = (Card*)image->data;

    wlen = XRange>>5;
    blen = XRange&31;
    if( blen )
    {   wmask = WhiteCol << (blen-1);
        bmask = BlackCol << (blen-1);
        len = wlen+1;
    } else len = wlen;

    whigh = WhiteCol << 31;
    bhigh = BlackCol << 31;

    /* Allow Compiler Optimisation */
    xmax = XRange;  ymax = YRange;

    error = 0;
    for( y=0; y<ymax; y++ )
    {    for( x=0; x<wlen; x++ )
         {   for( count=0; count<32; count++ )
             {   error += Intensity[*src++];
                 bits <<= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= WhiteCol;
                 } else bits |= BlackCol;
             }
             *dst++ = bits;
         }

         if( blen )
         {   for( count=0; count<blen; count++ )
             {   error += Intensity[*src++];
                 bits <<= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= WhiteCol;
                 } else bits |= BlackCol;
             }
             *dst++ = bits;
           
             /* Asymmetric Loop Unrolling! */
             if( ++y == ymax ) break;
             src += xmax;
             dst += wlen;

             bits = 0;
             for( count=0; count<blen; count++ )
             {   error += Intensity[*(--src)];
                 bits >>= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= wmask;
                 } else bits |= bmask;
             }
             *(--dst) = bits;
         } else
         {   /* Asymmetric Loop Unrolling! */
             if( ++y == ymax ) break;
             src += xmax;
             dst += len;
         }

         for( x=wlen-1; x>=0; x-- )
         {   for( count=0; count<32; count++ )
             {   error += Intensity[*(--src)];
                 bits >>= 1;
                 if( error >= 128  )
                 {   error -= 255;
                        bits |= whigh;
                 } else bits |= bhigh;
             }
             *(--dst) = bits;
         }
         src += xmax;
         dst += len;
    }
}


void TransferImage()
{
    if( Monochrome )
        DitherImage();

#ifdef MITSHM
    if( SharedMemFlag )
    {   XShmPutImage(dpy,CanvWin,gcon,image,0,0,0,0,XRange,YRange,False);
        XSync(dpy,False);
    } else
    {   XPutImage( dpy, CanvWin, gcon,image,0,0,0,0,XRange,YRange);
        XFlush(dpy);
    }
#else
    XPutImage( dpy, CanvWin, gcon, image, 0, 0, 0, 0, XRange, YRange );
    XFlush(dpy);
#endif
}


void ClearImage()
{
    XClearWindow( dpy, CanvWin );
    XFlush(dpy);
}


int PrintImage()
{
    return( False );
}


int ClipboardImage()
{
    return( False );
}


static int HandleIPCError( disp, ptr )
    Display *disp;  XErrorEvent *ptr;
{
    return( 0 );
}


static void HandleIPCCommand()
{
    static unsigned long len,left;
    static unsigned char *command;
    static Window source;
    static int serial;
    static int format;
    static Atom type;
    char buffer[32];
    int (*handler)();

    register int rlen;
    register int result;
    register char *cmnd;
    register char *ptr;

    command = NULL;
    result = XGetWindowProperty( dpy, MainWin, CommAtom, 0, 1024, True, 
                                 XA_STRING, &type, &format, &len, &left,
                                 &command );
    if( (result!=Success) || (type!=XA_STRING) || (format!=8) )
    {   if( command ) XFree( (char*)command );
        return;
    }

    result = 0;
    ptr = (char*)command;
    if( !*ptr )
    {   /* Tcl/Tk4.0 and later */

        ptr++;
        while( ptr < (char*)command+len )
        {    if( (ptr[0]=='c') && (ptr[1]=='\0') )
             {   ptr += 2;
                 cmnd = (char*)NULL;
                 source = serial = 0;
                 while( (ptr<(char*)command+len) && (*ptr=='-') )
                 {   if( (ptr[1]=='r') && (ptr[2]==' ') )
                     {   sscanf(ptr+3,"%x %d\n",(int*)&source,&serial);
                     } else if( (ptr[1]=='s') && (ptr[2]==' ') )
                         cmnd = ptr+3;
                     while( *ptr ) ptr++;
                     ptr++;
                 }

                 if( !cmnd ) continue;
                 result = ExecuteIPCCommand(cmnd);
                 if( !source || !serial ) continue;

                 buffer[0]='\0';
                 buffer[1]='r';
                 buffer[2]='\0';
                 buffer[3]='-';
                 buffer[4]='r';
                 buffer[5]=' ';
                 buffer[6]= result? '1' : '0';
                 buffer[7]='\0';
                 sprintf(buffer+8,"-s %d",serial);
                 rlen = strlen(buffer+8)+9;

                 /* Return Tcl/Tk v4.0 result! */
                 handler = XSetErrorHandler( HandleIPCError );
                 XChangeProperty(dpy,source,CommAtom, XA_STRING, 8,
                                 PropModeAppend,(unsigned char*)buffer,rlen);
                 XSync(dpy,False);
                 XSetErrorHandler(handler);
             } else /* Unrecognised command! */
             {   while( *ptr ) ptr++;
                 ptr++;
             }
        }

    } else while( *ptr )
    {   /* Tcl/Tk3.0 and later */
        if( *ptr=='C' )
        {   sscanf(ptr+1,"%x %x\n",(int*)&source,&serial);
            while( *ptr && (*ptr!='|') ) ptr++;
            if( *ptr=='|' )
            {   result = ExecuteIPCCommand(ptr+1);
            } else result = 0;

            sprintf(buffer,"R %x 0 %d",serial,result);
            handler = XSetErrorHandler( HandleIPCError );
            XChangeProperty( dpy, source, CommAtom, XA_STRING, 8,
                             PropModeAppend, (unsigned char*)buffer, 
                             strlen(buffer)+1 );
            XSync(dpy,False);
            XSetErrorHandler(handler);
        } 

        /* Next Command! */
        while( *ptr++ );
    }
    XFree( (char*)command );

    if( (result==IPC_Quit) || (result==IPC_Exit) )
        RasMolExit();
}


static int CropRange( val, min, max )
    int val, min, max;
{
    if( val<min ) return( min );
    if( val>max ) return( max );
    return( val );
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


void SetMouseMode( mode )
    int mode;
{
    if( mode==MouseMode )
        return;

    if( (mode==MMQuanta) || (MouseMode==MMQuanta) )
    {   /* Enable/Disable Pointer Motion Events! */
        attr.event_mask = ExposureMask | ButtonPressMask | ButtonMotionMask 
                        | ButtonReleaseMask;
        if( mode==MMQuanta ) attr.event_mask |= PointerMotionMask;
        XChangeWindowAttributes( dpy, CanvWin, CWEventMask, &attr );
    }
    MouseMode = mode;
}


static void MouseMove( status, dx, dy )
    int status, dx, dy;
{
    register int index;

    if( MouseMode == MMRasMol )
    {   if( status & ShiftMask )
        {   if( status & Button1Mask ) 
            {   if( dy ) /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
            } else if( status & (Button2Mask|Button3Mask) )
                if( dx ) /* Z Rotation Horizontal */
                {   WrapDial( 2, (Real)dx/WRange );
                    ReDrawFlag |= RFRotateZ;
                }
        } else if( status & ControlMask )
        {   if( status & Button1Mask )
            {   if( dy ) /* Slab Vertical */
                {   ClampDial( 7, (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                }
            }

        } else /* Unmodified! */
            if( status & Button1Mask )
            {   if( dx ) /* Rotate Y Horizontal */
                {   WrapDial( 1, (Real)dx/WRange );
                    index = (DialValue[1]+1.0)*(XRange-48);
                    NewScrlX = (index>>1)+16;
                    ReDrawFlag |= RFRotateY;
                }

                if( dy ) /* Rotate X Vertical */
                {   WrapDial( 0, (Real)dy/HRange );
                    index = (DialValue[0]+1.0)*(YRange-48);
                    NewScrlY = (index>>1)+16;
                    ReDrawFlag |= RFRotateX;
                }
            } else if( status & (Button2Mask|Button3Mask) )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }

                if( dy ) /* Translate Y Vertical */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            }
    } else if( MouseMode==MMQuanta )
    {   if( status & ShiftMask )
        {   if( status & Button1Mask )
            {   if( dy ) /* Slab Vertical */
                {   ClampDial( 7, (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                }
            } else if( status & Button2Mask )
            {   if( dx ) /* Translate X Horizontal */
                {   ClampDial( 4, (Real)dx/XRange );
                    ReDrawFlag |= RFTransX;
                }

                if( dy ) /* Translate Y Vertical */
                {   ClampDial( 5, (Real)dy/YRange );
                    ReDrawFlag |= RFTransY;
                }
            } else if( !(status & Button3Mask) )
                if( dy ) /* Zoom Vertical */
                {   ClampDial( 3, (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                }
        } else if( status & Button2Mask )
        {   if( dx ) /* Rotate Y Horizontal */
            {   WrapDial( 1, (Real)dx/WRange );
                index = (DialValue[1]+1.0)*(XRange-48);
                NewScrlX = (index>>1)+16;
                ReDrawFlag |= RFRotateY;
            }

            if( dy ) /* Rotate X Vertical */
            {   WrapDial( 0, (Real)dy/HRange );
                index = (DialValue[0]+1.0)*(YRange-48);
                NewScrlY = (index>>1)+16;
                ReDrawFlag |= RFRotateX;
            }
        } else if( status & Button3Mask )
            if( dx ) /* Z Rotation Horizontal */
            {   WrapDial( 2, (Real)dx/WRange );
                ReDrawFlag |= RFRotateZ;
            }
    } else /* MMInsight */
        switch( status & (Button1Mask|Button2Mask|Button3Mask) )
        {   case( Button1Mask ):
                    if( dx ) /* Rotate Y Horizontal */
                    {   WrapDial( 1, (Real)dx/WRange );
                        index = (DialValue[1]+1.0)*(XRange-48);
                        NewScrlX = (index>>1)+16;
                        ReDrawFlag |= RFRotateY;
                    }

                    if( dy ) /* Rotate X Vertical */
                    {   WrapDial( 0, (Real)dy/HRange );
                        index = (DialValue[0]+1.0)*(YRange-48);
                        NewScrlY = (index>>1)+16;
                        ReDrawFlag |= RFRotateX;
                    }
                    break;

            case( Button2Mask ):
                    if( dx ) /* Translate X Horizontal */
                    {   ClampDial( 4, (Real)dx/XRange );
                        ReDrawFlag |= RFTransX;
                    }

                    if( dy ) /* Translate Y Vertical */
                    {   ClampDial( 5, (Real)dy/YRange );
                        ReDrawFlag |= RFTransY;
                    }
                    break;

            case( Button1Mask|Button2Mask ):
                    ClampDial( 3, (Real)dx/WRange - (Real)dy/HRange );
                    ReDrawFlag |= RFZoom;
                    break;

            case( Button1Mask|Button3Mask ):
                    WrapDial( 2, (Real)dx/WRange - (Real)dy/HRange );
                    ReDrawFlag |= RFRotateZ;
                    break;

            case( Button1Mask|Button2Mask|Button3Mask ):
                    ClampDial( 7, (Real)dx/XRange - (Real)dy/YRange );
                    ReDrawFlag |= RFSlab;
                    break;
        }
}


static void DoneEvents()
{
    register Real temp;
    register int index;

    if( HeldButton == YScrlDial )
    {   index = NewScrlY+HeldStep;
        if( YScrlDial<3 )
        {   if( index<16 )             
            {   index += YRange-48;
            } else if( index>YRange-32 ) 
                index -= YRange-48;
            NewScrlY = index;
        } else NewScrlY = CropRange(index,16,YRange-32);
    }

    if( NewScrlY != ScrlY )
    {   XClearArea(dpy,YScrlWin,0,ScrlY,16,16,False);
        XCopyArea(dpy,Scrl,YScrlWin,gcon,0,0,16,16,0,NewScrlY);

        temp = ((Real)(NewScrlY-16))/(YRange-48);
        DialValue[YScrlDial] = 2.0*temp - 1.0;
        ReDrawFlag |= (1<<YScrlDial);
        ScrlY = NewScrlY;
    }

    if( HeldButton == XScrlDial )
    {   index = NewScrlX+HeldStep;
        if( XScrlDial<3 )
        {   if( index<16 ) 
            {   index += XRange-48;
            } else if( index>XRange-32 ) 
                index -= XRange-48;
            NewScrlX = index;
        } else NewScrlX = CropRange(index,16,XRange-32);
    }

    if( NewScrlX != ScrlX )
    {   XClearArea(dpy,XScrlWin,ScrlX,0,16,16,False);
        XCopyArea(dpy,Scrl,XScrlWin,gcon,0,0,16,16,NewScrlX,0);

        temp = ((Real)(NewScrlX-16))/(XRange-48);
        DialValue[XScrlDial] = 2.0*temp - 1.0;
        ReDrawFlag |= (1<<XScrlDial);
        ScrlX = NewScrlX;
    }
    /* XSync(dpy,False); */
    XFlush(dpy);
}


static int ProcessEvent( event )
    XEvent *event;
{
    register int result;
    register int index;

    result = 0;
    switch( event->type )
    {   case(ButtonPress):
            {   XButtonPressedEvent *ptr;

                HeldButton = -1;
                ptr = (XButtonPressedEvent*)event;

                if( ptr->window==CanvWin )
                {   InitX = PointX = ptr->x;
                    InitY = PointY = ptr->y;
                } else if( ptr->window==MenuWin )
                {   if( !DisableMenu )
                        if( HandleMenuClick(ptr->x) )
                            result = HandleMenuLoop();
                } else if( ptr->window==XScrlWin )
                {   ReDrawFlag |= RFRotateY;
                    if( ptr->x<16 )
                    {   HeldButton = XScrlDial;
                        HeldStep = -XScrlSkip;
                    } else if( ptr->x>=XRange-16 )
                    {   HeldButton = XScrlDial;
                        HeldStep = XScrlSkip;
                    } else
                    {   index = ptr->x-8;
                        if( XScrlDial<3 )
                        {   if( index>XRange-32 ) index -= XRange-48;
                            else if( index<16 ) index += XRange-48;
                            NewScrlX = index;
                        } else NewScrlX = CropRange(index,16,XRange-32);
                    }

                } else if( ptr->window==YScrlWin )
                {   ReDrawFlag |= RFRotateX;
                    if( ptr->y<16 )
                    {   HeldButton = YScrlDial;
                        HeldStep = -YScrlSkip;
                    } else if( ptr->y>=YRange-16 )
                    {   HeldButton = YScrlDial;
                        HeldStep = YScrlSkip;
                    } else
                    {   index = ptr->y-8;
                        if( YScrlDial<3 )
                        {   if( index>YRange-32 ) index -= YRange-48;
                            else if( index<16 ) index += YRange-48;
                            NewScrlY = index;
                        } else NewScrlY = CropRange(index,16,YRange-32);
                    }

                } 
            } break;

        case(MotionNotify):
            {   XMotionEvent *ptr;
                int dx, dy;

                ptr = (XMotionEvent*)event;
                if( ptr->window==CanvWin )
                {   if( !IsClose(ptr->x,InitX) || !IsClose(ptr->y,InitY) )
                    {   dx = ptr->x-PointX;  dy = ptr->y-PointY;
                        MouseMove( ptr->state, dx, dy );

                        PointX = ptr->x;
                        PointY = ptr->y;
                    }
                } else if( HeldButton == -1 )
                {   if( ptr->window==XScrlWin )
                    {   index = ptr->x-8;
                        NewScrlX = CropRange(index,16,XRange-32);
                    } else /* if( ptr->window==YScrlWin ) */
                    {   index = ptr->y-8;
                        NewScrlY = CropRange(index,16,YRange-32);
                    }
                }
            } break;
             
        case(ButtonRelease):
            {   XButtonReleasedEvent *ptr;

                if( HeldButton != -1 )
                {   /* Three button emulation fix! */
                    DoneEvents();  HeldButton = -1;
                }

                ptr = (XButtonReleasedEvent*)event;
                if( ptr->window==CanvWin )
                {   PointX = ptr->x;  PointY = ptr->y;
                    if( IsClose(PointX,InitX) && IsClose(PointY,InitY) )
                        if( ptr->state & (ShiftMask|ControlMask) )
                        {      ReDrawFlag |= RFPoint1;
                        } else ReDrawFlag |= RFPoint2;
                }
            } break;

        case(KeyPress):
            {   XKeyPressedEvent *ptr;
                static KeySym symbol;
                static char keychar;

                keychar = '\0';
                ptr = (XKeyPressedEvent*)event;
                index = XLookupString(ptr,&keychar,1,&symbol,NULL);
                switch( symbol )
                {   case(XK_Begin):
                    case(XK_Home):  ProcessCharacter(0x01);  break;
                    case(XK_Right): ProcessCharacter(0x06);  break;
                    case(XK_Left):  ProcessCharacter(0x02);  break;
                    case(XK_End):   ProcessCharacter(0x05);  break;
                    case(XK_Up):
                    case(XK_Prior): ProcessCharacter(0x10);  break;
                    case(XK_Down):
                    case(XK_Next):  ProcessCharacter(0x0e);  break;

                    case(XK_F10):   if( !DisableMenu )
                                    {   SelectMenu(0);
                                        result = HandleMenuLoop();
                                    }
                                    break;

                    default:        if( index == 1 )
                                        if( !(ptr->state&Mod1Mask) )
                                        {   if( ProcessCharacter(keychar) )
                                            {   if( ProcessCommand() )
                                                    RasMolExit();

                                                if( !CommandActive )
                                                    ResetCommandLine(0);
                                            }
                                        } else if( !DisableMenu )
                                            if( HandleMenuKey(keychar) )
                                                result = HandleMenuLoop();
                }
            } break;


        case(Expose):
            {   XExposeEvent *ptr;

                ptr = (XExposeEvent*)event;
                if( ptr->window==CanvWin )
                {   if( image ) {
#ifdef MITSHM
                        if( SharedMemFlag )
                        {   XShmPutImage( dpy, CanvWin, gcon, image,
                                          ptr->x, ptr->y, ptr->x, ptr->y,
                                          ptr->width, ptr->height, False);
                            XSync(dpy,False);
                        } else
#endif 
                        XPutImage( dpy, CanvWin, gcon, image,
                                   ptr->x, ptr->y, ptr->x, ptr->y,
                                   ptr->width, ptr->height );
                    } else XClearWindow( dpy, CanvWin );
		    
                } else if( ptr->window==MainWin )
                {   DrawMainWin();
                    DrawMenuBar();
                } else if( ptr->window==XScrlWin )
                {   DrawXScroll();
                } else if( ptr->window==YScrlWin )
                    DrawYScroll();
                XFlush(dpy);
            } break;

        case(EnterNotify):
            {   XCrossingEvent *ptr;

                ptr = (XCrossingEvent*)event;
                if( ptr->detail != NotifyInferior )
                {   if( LocalMap )
                        XInstallColormap(dpy,lmap);
                }
#ifdef DIALBOX
                if( UseDials )
                    GetDialState();
#endif
            }
            break;

        case(LeaveNotify):
            if( LocalMap )
            {   XCrossingEvent *ptr;

                ptr = (XCrossingEvent*)event;
                if( ptr->detail != NotifyInferior )
                    XUninstallColormap(dpy,lmap);
            }
            break;

        case(ConfigureNotify):
            {   XConfigureEvent *ptr;
                register int wide,high;

                ptr = (XConfigureEvent*)event;
                high = CropRange(ptr->height,MinHeight,MaxHeight);
                wide = CropRange(ptr->width, MinWidth, MaxWidth );

                if( (wide!=MainWide) || (high!=MainHigh) )
                    ReSizeWindow(wide,high);
            } break;

        case(ClientMessage):
            {   XClientMessageEvent *ptr;

                ptr = (XClientMessageEvent*)event;
                if( (ptr->message_type==ProtoXAtom) && 
                    (ptr->data.l[0]==DelWinXAtom) )
                    RasMolExit();
            } break;

        case(PropertyNotify):
            {   XPropertyEvent *ptr;

                ptr = (XPropertyEvent*)event;
                if( (ptr->atom==CommAtom) &&
                    (ptr->state==PropertyNewValue) )
                    HandleIPCCommand();
            } break;

        case(MapNotify):
            DrawXScroll();
            DrawYScroll();
            DrawMainWin();
            DrawMenuBar();
            break;

        default:  
#ifdef DIALBOX
            if( event->type == DialEvent )
                HandleDialEvent( event );
#endif
            break;
    }
    return( result );
}


/*************************/
/* Modal Dialog Handling */
/*************************/


static int HandleMenuLoop()
{
    register unsigned int mask;
    register int result;
    register int done;
    auto XEvent event;

    /* Passive Pointer Grab */
    mask = ButtonPressMask | ButtonReleaseMask | ButtonMotionMask;
    XGrabPointer(dpy,MenuWin,False,mask,
                 GrabModeAsync,GrabModeAsync,
                 None,None,CurrentTime);

    HeldButton = -1;
    MenuFocus = True;
    DrawMenuBar();

    result = 0;
    done = False;
    while( !done )
    {   XNextEvent( dpy, &event );
        switch( event.type )
        {   case(Expose):
                {   XExposeEvent *ptr;

                    ptr = (XExposeEvent*)&event;
                    if( ptr->window==PopUpWin )
                    {   DrawPopUpMenu();
                    } else ProcessEvent(&event);
                } break;

            case(ButtonPress): 
                {   XButtonPressedEvent *ptr;

                    ptr = (XButtonPressedEvent*)&event;
                    /* All Events Relative to MenuWin */
                    if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                    {   HandleMenuClick(ptr->x);
                    } else if( PopUpFlag )
                    {   HandleItemClick(ptr->x,ptr->y);
                    } else done = True;
                } break;

            case(MotionNotify):
                    if( ItemFocus )
                    {   XMotionEvent *ptr;

                        ptr = (XMotionEvent*)&event;
                        /* All Events Relative to MenuWin */
                        if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                        {   HandleMenuClick( ptr->x );
                        } else if( PopUpFlag )
                            HandleItemMove(ptr->x,ptr->y);
                    } break;

            case(ButtonRelease):
                    {   XButtonReleasedEvent *ptr;

                        ptr = (XButtonReleasedEvent*)&event;
                        /* All Events Relative to MenuWin */
                        if( (ptr->y>=0) && (ptr->y<=FontHigh+5) )
                        {   if( HandleMenuClick( ptr->x ) )
                            {   SelectFirstItem(MenuBarSelect);
                                DrawPopUpMenu();
                            } else done = True;
                        } else if( PopUpFlag )
                        {   if( ItemFocus )
                                HandleItemClick(ptr->x,ptr->y);
                            if( ItemFlag )
                                result = (MenuBarSelect<<8) +
                                         MenuItemSelect+1;
                            done = True;
                        } else done = False;
                        ItemFocus = False;
                    }
                    break;
     
            case(KeyPress):
                if( !ItemFocus )
                {   XKeyPressedEvent *ptr;
                    static KeySym symbol;
                    static char keychar;
                    register int index;

                    ptr = (XKeyPressedEvent*)&event;
                    index = XLookupString(ptr,&keychar,1,&symbol,NULL);
                    switch( symbol )
                    {   case(XK_Right): index = MenuBarSelect+1;
                                        if( index != MenuBarCount )
                                        {   SelectMenu( index );
                                        } else SelectMenu( 0 );
                                        break;

                        case(XK_Left):  if( MenuBarSelect )
                                        {   SelectMenu( MenuBarSelect-1 );
                                        } else SelectMenu( MenuBarCount-1 );
                                        break;

                        case(XK_Up):    if( !PopUpFlag )
                                        {   PopUpFlag = True;
                                            SelectMenu(MenuBarSelect);
                                        } else SelectPrevItem();
                                        break;

                        case(XK_Down):  if( !PopUpFlag )
                                        {   PopUpFlag = True;
                                            SelectMenu(MenuBarSelect);
                                        } else SelectNextItem();
                                        break;

                        case(XK_KP_Enter):
                        case(XK_Linefeed):
                        case(XK_Return):   if( PopUpFlag && ItemFlag )
                                               result = (MenuBarSelect<<8) +
                                                        MenuItemSelect+1;
                                           done = True;
                                           break;

                        default:    if( (index==1) && (keychar>=' ') ) 
                                    {   if( !(ptr->state&Mod1Mask) )
                                        {   if( PopUpFlag )
                                            {   result = HandleItemKey(keychar);
                                                if( result ) done = True;
                                            } else HandleMenuKey(keychar);
                                        } else HandleMenuKey(keychar);
                                    }
                    }
                } break;


            case(ConfigureNotify):  /* done = True; */
            default:                ProcessEvent(&event);
        }
    }

    /* Passive Grab Release */
    XUngrabPointer(dpy,CurrentTime);


    XUnmapWindow(dpy,PopUpWin);
    PopUpFlag = False;
    MenuFocus = False;
    DrawMenuBar();
    return( result );
}


int FetchEvent( wait )
    int wait;
{
    register int result;
    auto XEvent event;


    NewScrlX = ScrlX;
    NewScrlY = ScrlY;

    if( HeldButton != -1 ) wait = False;
    while( XPending(dpy) || (wait && !ReDrawFlag) )
    {   XNextEvent( dpy, &event );
        result = ProcessEvent(&event);
        if( result ) return( result );
    }
    DoneEvents();
    return( 0 );
}


int LookUpColour( name, red, grn, blu )
    char *name; int *red, *grn, *blu;
{
    static XColor exact, close;
    register Colormap map;

    map = (LocalMap)? lmap : cmap;
    if( XLookupColor(dpy,map,name,&exact,&close) )
    {   *red = exact.red>>8;
        *grn = exact.green>>8;
        *blu = exact.blue>>8;
        return(True);
    } else 
        return(False);
}


void BeginWait()
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,hglass);
        XDefineCursor(dpy,MainWin,hglass);
        XFlush(dpy);
    }
}


void EndWait()
{
    if( UseHourGlass )
    {   XDefineCursor(dpy,CanvWin,cross);
        XDefineCursor(dpy,MainWin,arrow);
        XFlush(dpy);
    }
}


void CloseDisplay()
{
#ifdef DIALBOX
    register int num;
#endif

    /* FatalXError! */
    if( !dpy ) return;

    if( image ) 
    {
#ifdef MITSHM
        if( SharedMemFlag )
        {   XShmDetach( dpy, &xshminfo );
            image->data = (char*)NULL;
            shmdt( xshminfo.shmaddr );
        }
#endif
        XDestroyImage( image );
    }

    if( *TkInterp )
    {   XGrabServer( dpy );
        DeRegisterInterpName(TkInterp);
        XUngrabServer( dpy );
    }

#ifdef DIALBOX
    if( UseDials )
    {   if( UseDialLEDs )
            for( num=0; num<8; num++ )
                SetDialLabel(num,"");
        XCloseDevice(dpy,Dials);
    }
#endif
    XCloseDisplay( dpy );
}
