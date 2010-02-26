/* graphics.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#ifdef APPLEMAC      
#define DefaultWide  400
#define DefaultHigh  400
#endif

#ifdef IBMPC
#define DefaultWide  480
#define DefaultHigh  480
#endif

#ifndef DefaultWide
#define DefaultWide  576
#define DefaultHigh  576
#endif

#ifdef EIGHTBIT
#define LutSize  256
#else
#define LutSize  1024
#endif

#define RFRotateX  0x0001
#define RFRotateY  0x0002
#define RFRotateZ  0x0004
#define RFZoom     0x0008
#define RFTransX   0x0010
#define RFTransY   0x0020
#define RFTransZ   0x0040
#define RFSlab     0x0080
#define RFReSize   0x0100
#define RFColour   0x0200
#define RFRefresh  0x0400
#define RFPoint1   0x1000
#define RFPoint2   0x2000

#define RFTrans    0x0070
#define RFRotate   0x0007
#define RFApply    0x017F
#define RFDials    0x00FF
#define RFMagnify  0x0108
#define RFInitial  0x01FF
#define RFPoint    0x3000

#define MMRasMol   0x00
#define MMInsight  0x01
#define MMQuanta   0x02

#define ButMax   8

#ifdef GRAPHICS
double DialValue[8];
int WinHigh, WinWide;
int PointX, PointY;
int XRange, WRange;
int YRange, HRange;
int UseHourGlass;
int DisableMenu;
int ReDrawFlag;
int MouseMode;
int Range;

Pixel __huge *FBuffer;
short __huge *DBuffer;

Pixel Lut[LutSize];
Byte RLut[LutSize];
Byte GLut[LutSize];
Byte BLut[LutSize];
Byte ULut[LutSize];

#ifdef IBMPC
LOGPALETTE __far *Palette;
HPALETTE ColourMap;
HGLOBAL FBufHandle;
HGLOBAL DBufHandle;
HBITMAP PixMap;
HWND CanvWin;
#endif /* IBMPC */

#ifdef APPLEMAC
ControlHandle HScroll;
ControlHandle VScroll;
CursHandle CanvCursor;
CursHandle CmndCursor;
CursHandle WaitCursor;
WindowPtr CanvWin;
WindowPtr CmndWin;
THPrint PrintHand;
Handle FBufHandle;
Handle DBufHandle;
#endif /* APPLEMAC */

#else /* GRAPHICS */
extern double DialValue[8];
extern int WinHigh, WinWide;
extern int PointX, PointY;
extern int XRange, WRange;
extern int YRange, HRange;
extern int UseHourGlass;
extern int DisableMenu;
extern int ReDrawFlag;
extern int MouseMode;
extern int Range;

extern Pixel __huge *FBuffer;
extern short __huge *DBuffer;

extern Pixel Lut[LutSize];
extern Byte RLut[LutSize];
extern Byte GLut[LutSize];
extern Byte BLut[LutSize];
extern Byte ULut[LutSize];

#ifdef IBMPC
extern LOGPALETTE __far *Palette;
extern HPALETTE ColourMap;
extern HGLOBAL FBufHandle;
extern HGLOBAL DBufHandle;
extern HBITMAP PixMap;
extern HWND CanvWin;
#endif /* IBMPC */

#ifdef APPLEMAC
extern ControlHandle HScroll;
extern ControlHandle VScroll;
extern CursHandle CanvCursor;
extern CursHandle CmndCursor;
extern CursHandle WaitCursor;
extern WindowPtr CanvWin;
extern WindowPtr CmndWin;
extern THPrint PrintHand;
extern Handle FBufHandle;
extern Handle DBufHandle;
#endif /* APPLEMAC */

#ifdef FUNCPROTO
int CreateImage();
void TransferImage();
int ClipboardImage();
void ClearImage();
int PrintImage();

void AllocateColourMap();
void UpdateScrollBars();
int LookUpColour( char*, int*, int*, int* );
void SetMouseMode( int );
void EnableMenus( int );
void CloseDisplay();
void BeginWait();
void EndWait();

#ifdef IBMPC
int OpenDisplay( HANDLE, int );
#else
int OpenDisplay( int, int );
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
int FetchEvent( int );
#endif

#else /* non-ANSI C compiler */
int CreateImage();
void TransferImage();
int ClipboardImage();
void ClearImage();
int PrintImage();

int OpenDisplay();
void AllocateColourMap();
void UpdateScrollBars();
int LookUpColour();
void SetMouseMode();
void EnableMenus();
void CloseDisplay();
void BeginWait();
void EndWait();

#if !defined(IBMPC) && !defined(APPLEMAC)
int FetchEvent();
#endif

#endif
#endif /* GRAPHICS */
