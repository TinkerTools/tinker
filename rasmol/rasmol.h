/* rasmol.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

/*===========================*/
/*  User Definable Options!  */
/*===========================*/

/* #define IBMPC        */
/* #define APPLEMAC     */
/* #define DIALBOX      */
#define TERMIOS
#define SOCKETS
#define MITSHM

/* #define MMIOLIB      */
/* #define OMLCOMMS     */
/* #define CEXIOLIB     */

#if !defined(EIGHTBIT) && !defined(THIRTYTWOBIT) && !defined(SIXTEENBIT)
/* #define THIRTYTWOBIT */
/* #define SIXTEENBIT   */
/* #define EIGHTBIT     */
#endif

/*========================*/
/*  Default User Options! */
/*========================*/

#ifdef IBMPC
#undef THIRTYTWOBIT
#undef SIXTEENBIT
#endif

#if !defined(EIGHTBIT) && !defined(THIRTYTWOBIT) && !defined(SIXTEENBIT)
#define EIGHTBIT
#endif

#ifndef RASMOLDIR
#ifdef IBMPC
#define RASMOLDIR  "C:\\RASWIN\\"
#endif

#if !defined(IBMPC) && !defined(APPLEMAC)
#define RASMOLDIR "/usr/local/lib/rasmol/"
#endif
#endif

#if !defined(INVERT) && defined(IBMPC)
#define INVERT
#endif

#ifndef FUNCPROTO
#if defined(__STDC__) || defined(IBMPC) || defined(APPLEMAC)
#define FUNCPROTO
#endif
#endif

/*==============================*/
/*  Application-wide Constants  */
/*==============================*/

#ifndef True
#define True  1
#define False 0
#endif

#ifndef PI   /* Avoid Linux Warnings! */
#define PI   3.14159265358979323846
#endif

typedef double Real;
#ifndef APPLEMAC
typedef unsigned char Byte;
#endif

#ifdef __STDC__
typedef signed char Char;
#else
typedef char Char;
#endif

#ifdef _LONGLONG
typedef unsigned int Card;
typedef int Long;
#else
typedef unsigned long Card;
typedef long Long;
#endif

#ifdef EIGHTBIT
typedef unsigned char Pixel;
#else
#ifdef THIRTYTWOBIT
typedef Long Pixel;
#else
typedef short Pixel;
#endif
#endif

#define Rad2Deg      (180.0/PI)
#define Deg2Rad      (PI/180.0)
#define AbsFun(a)    (((a)<0)? -(a) : (a))
#define MinFun(a,b)  (((a)<(b))? (a) : (b) )
#define MaxFun(a,b)  (((a)>(b))? (a) : (b) )

#if defined(__STDC__) || defined(IBMPC)  || defined(APPLEMAC) || defined(__sgi)
#define ToUpper(x)   (toupper((x)))
#else
#define ToUpper(x)   (islower((x))?toupper((x)):(x))
#endif

#if !defined(IBMPC) || defined(_WIN32)
#ifdef APPLEMAC
#define _fmalloc   NewPtrSys
#define _ffree(x)  DisposePtr((Ptr)(x))
#else
#define _fmalloc   malloc
#define _ffree     free
#endif
#define _fstrnicmp strnicmp
#define _fstrcmp   strcmp
#define _fmemset   memset
#define __huge
#define __far
#endif

#define ItemCount       8
#define AdvPickAtom     0
#define AdvPickNumber   1
#define AdvSelectCount  2
#define AdvName         3
#define AdvIdent        4
#define AdvClass        5
#define AdvImage        6
#define AdvPickCoord    7

#ifndef RASMOL
#ifdef FUNCPROTO
void WriteChar( int );
void WriteString( char* );
void RasMolFatalExit( char* );
void AdviseUpdate( int );
void RefreshScreen();
void RasMolExit();

#else /* non-ANSI C compiler */
void WriteChar();
void WriteString();
void RasMolFatalExit();
void RefreshScreen();
void AdviseUpdate();
void RasMolExit();

#endif
#endif
