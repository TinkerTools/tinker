/* render.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

/* These values set the sizes of the sphere rendering
 * tables. The first value, maxrad, is the maximum
 * sphere radius and the second value is the table
 * size = (maxrad*(maxrad+1))/2 + 1
 */
/* #define MAXRAD    120   256   */
/* #define MAXTABLE  7261  32897 */
#define MAXRAD    255
#define MAXTABLE  32641

#define SlabReject       0x00
#define SlabHalf         0x01
#define SlabHollow       0x02
#define SlabFinal        0x03
#define SlabClose        0x04
#define SlabSection      0x05

#define PickNone         0x00
#define PickIdent        0x01
#define PickDist         0x02
#define PickAngle        0x03
#define PickTorsn        0x04
#define PickLabel        0x05
#define PickMonit        0x06
#define PickCentr        0x07

#define ViewLeft         0
#define ViewRight        1

#define ColBits          24

#define VOXORDER       21
#define VOXORDER2      (VOXORDER*VOXORDER)
#define VOXSIZE        (VOXORDER2*VOXORDER)

typedef struct _Item {
        struct _Item __far *list;
        Atom  __far *data;
    } Item;

#ifdef RENDER
int UseDepthCue;
int UseStereo,StereoView;
int UseShadow,DisplayMode;
int UseClipping,UseSlabPlane;
int SlabMode,SlabValue;
int SlabInten,SliceValue;
int ImageRadius,ImageSize;

double StereoAngle;

int DrawBoundBox,DrawAxes;
int DrawUnitCell;

Real IVoxRatio;
int VoxelsClean;
int BucketFlag;
int FBClear;

Card __far *ColConst;
#if defined(IBMPC) || defined(APPLEMAC)
void __far * __far *HashTable;
Byte __far * __far *LookUp;
Byte __far *Array;

#else /* UNIX */
void *HashTable[VOXSIZE];
Byte *LookUp[MAXRAD];
Byte Array[MAXTABLE];
#endif

#else
extern int UseDepthCue;
extern int UseStereo,StereoView;
extern int UseShadow, DisplayMode;
extern int UseClipping,UseSlabPlane;
extern int SlabMode,SlabValue;
extern int SlabInten,SliceValue;
extern int ImageRadius,ImageSize;

extern double StereoAngle;

extern int DrawBoundBox,DrawAxes;
extern int DrawUnitCell;

extern Real IVoxRatio;
extern int VoxelsClean;
extern int BucketFlag;
extern int FBClear;

extern Card __far *ColConst;
#if defined(IBMPC) || defined(APPLEMAC)
extern void __far * __far *HashTable;
extern Byte __far * __far *LookUp;
extern Byte __far *Array;

#else /* UNIX */
extern void *HashTable[VOXSIZE];
extern Byte *LookUp[MAXRAD];
extern Byte Array[MAXTABLE];
#endif

#ifdef FUNCPROTO
void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ShadowTransform();

void ResetVoxelData();
void CreateVoxelData( int );

void DrawFrame();
void ResetRenderer();
void InitialiseRenderer();
void SetStereoMode( int );
void SetPickMode( int );
void PickAtom( int, int, int );
unsigned int isqrt( Card );

#else /* non-ANSI C compiler */
void ClearBuffers();
void ReSizeScreen();
void ReAllocBuffers();
void ShadowTransform();

void ResetVoxelData();
void CreateVoxelData();

void DrawFrame();
void ResetRenderer();
void InitialiseRenderer();
void SetStereoMode();
void SetPickMode();
void PickAtom();
unsigned int isqrt();

#endif
#endif
