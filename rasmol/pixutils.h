/* pixutils.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

typedef struct {
        Pixel __huge *fbuf;
        short __huge *dbuf;
        int xmax, ymax;
        int yskip;
    } ViewStruct;

#define ZValid(z)     ((!UseSlabPlane) || ((z)<SlabValue))
#define XValid(x)     (((x)>=0)&&((x)<View.xmax))
#define YValid(y)     (((y)>=0)&&((y)<View.ymax))

#ifdef PIXUTILS
ViewStruct View;
int SplineCount;
int FontSize;

#else
extern ViewStruct View;
extern int SplineCount;
extern int FontSize;

#ifdef FUNCPROTO
void PlotDeepPoint( int, int, int, int );
void ClipDeepPoint( int, int, int, int );
void DrawTwinLine( int, int, int, int, int, int, int, int );
void ClipTwinLine( int, int, int, int, int, int, int, int );
void DrawTwinVector( int, int, int, int, int, int, int, int );
void ClipTwinVector( int, int, int, int, int, int, int, int );
void ClipDashVector( int, int, int, int, int, int, int, int );

void DrawCylinder( int, int, int, int, int, int, int, int, int );
void ClipCylinder( int, int, int, int, int, int, int, int, int );
void DrawSphere( int, int, int, int, int );
void ClipSphere( int, int, int, int, int );

void SetFontSize( int );
void DisplayString( int, int, int, char*, int );
void InitialisePixUtils();

#else /* non-ANSI C compiler */
void PlotDeepPoint();
void ClipDeepPoint();
void DrawTwinLine();
void ClipTwinLine();
void DrawTwinVector();
void ClipTwinVector();
void ClipDashVector();

void DrawCylinder();
void ClipCylinder();
void DrawSphere();
void ClipSphere();

void SetFontSize();
void DisplayString();
void InitialisePixUtils();

#endif
#endif
