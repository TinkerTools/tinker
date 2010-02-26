/* transfor.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#define MaskColourFlag  0x01
#define MaskRadiusFlag  0x02
#define MaskBothFlag    0x03

#define RibColInside    0x01
#define RibColOutside   0x02
#define RibColBoth      0x03

#ifdef EIGHTBIT
#define DefaultAmbient    0.6
#define ColourDepth       16
#define ColourMask        15
#ifdef APPLEMAC
#define LastShade         14
#else
#define LastShade         15
#endif
#else
#define DefaultAmbient    0.05
#define ColourDepth       32
#define ColourMask        31
#define LastShade         31
#endif

typedef struct { 
        Long refcount;
        unsigned char r;
        unsigned char g;
        unsigned char b;
    } ShadeDesc;

#ifdef IBMPC
#define BackCol    0
#define BoxCol     1
#define LabelCol   2
#define FirstCol   3
#endif

#ifdef APPLEMAC
#define BackCol    1
#define BoxCol     2
#define LabelCol   3
#define FirstCol   4
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
#define BackCol    5
#define BoxCol     6
#define LabelCol   7
#define FirstCol   8
#endif

#define Colour2Shade(x)  ((int)((x)-FirstCol)/ColourDepth)
#define Shade2Colour(x)  ((x)*ColourDepth+FirstCol)

#ifdef TRANSFORM
ShadeDesc Shade[LastShade];
Real RotX[3],RotY[3],RotZ[3];
Real MatX[3],MatY[3],MatZ[3];
Real InvX[3],InvY[3],InvZ[3];
Long OrigCX,OrigCY,OrigCZ;
Long CenX, CenY, CenZ;

int FakeSpecular,SpecPower;
int BackR,BackG,BackB;
int LabR,LabG,LabB;
int BoxR,BoxG,BoxB;
int UseLabelCol;
int UseBackFade;
Real Ambient;

Real Scale,MaxZoom;
Real DScale,IScale;
Long SideLen,Offset;
Card WorldRadius,WorldSize;
int XOffset,YOffset,ZOffset;
int UseScreenClip;
int ZoomRange;

int Hydrogens;
int DrawAtoms,MaxAtomRadius;
int DrawBonds,MaxBondRadius;
int ZoneBoth;

#else
extern ShadeDesc Shade[LastShade];
extern Real RotX[3],RotY[3],RotZ[3];
extern Real MatX[3],MatY[3],MatZ[3];
extern Real InvX[3],InvY[3],InvZ[3];
extern Long OrigCX, OrigCY, OrigCZ;
extern Long CenX, CenY, CenZ;

extern int FakeSpecular,SpecPower;
extern int BackR,BackG,BackB;
extern int LabR,LabG,LabB;
extern int BoxR,BoxG,BoxB;
extern int UseLabelCol;
extern int UseBackFade;
extern Real Ambient;

extern Real Scale,MaxZoom;
extern Real DScale,IScale;
extern Long SideLen,Offset;
extern Card WorldRadius,WorldSize;
extern int XOffset,YOffset,ZOffset;
extern int UseScreenClip;
extern int ZoomRange;

extern int Hydrogens;
extern int DrawAtoms,MaxAtomRadius;
extern int DrawBonds,MaxBondRadius;
extern int ZoneBoth;

#ifdef FUNCPROTO
void SetRadiusValue( int );
void SetVanWaalRadius();
void DisableSpacefill();
void EnableWireframe( int, int );
void DisableWireframe();

void SelectZoneExpr( Expr* );
void RestrictZoneExpr( Expr* );
void RestrictZone( int );
void SelectZone( int );

int IsCPKColour( Atom __far * );
int IsVDWRadius( Atom __far * );

void DefineColourMap();
void ResetColourMap();

void ColourBondNone();
void ColourMonitNone();
void ColourBondAttrib( int, int, int );
void ColourMonitAttrib( int, int, int );
void MonoColourAttrib( int, int, int );
void CPKColourAttrib();
void ShapelyColourAttrib();

void DefaultRepresentation();

void DetermineClipping();
void InitialiseTransform();
void InitialTransform();
void PrepareTransform();
void ReviseInvMatrix();
void ApplyTransform();
void ResetTransform();

#else /* non-ANSI C compiler */
void SetRadiusValue();
void SetVanWaalRadius();
void DisableSpacefill();
void EnableWireframe();
void DisableWireframe();

void SelectZoneExpr();
void RestrictZoneExpr();
void RestrictZone();
void SelectZone();

int IsCPKColour();
int IsVDWRadius();

void DefineColourMap();
void ResetColourMap();

void ColourBondNone();
void ColourMonitNone();
void ColourBondAttrib();
void ColourMonitAttrib();
void MonoColourAttrib();
void CPKColourAttrib();
void ShapelyColourAttrib();

void DefaultRepresentation();

void DetermineClipping();
void InitialiseTransform();
void InitialTransform();
void PrepareTransform();
void ReviseInvMatrix();
void ApplyTransform();
void ResetTransform();

#endif
#endif
