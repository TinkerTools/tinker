/* rasmol.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

/* Includes a patch for Linux due
 * to Reece Hart of the Ponder Lab
 */

#ifndef sun386
#include <stdlib.h>
#endif
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include <signal.h>

#define RASMOL
#include "rasmol.h"
#include "graphics.h"
#include "molecule.h"
#include "infile.h"
#include "abstree.h"
#include "transfor.h"
#include "command.h"
#include "render.h"
#include "repres.h"
#include "pixutils.h"
#include "outfile.h"

#ifdef TERMIOS
#include <sys/types.h>
#include <sys/time.h>

#ifdef esv
#include <sysv/sys/termio.h>
#else
#ifdef __FreeBSD__
#include <sys/ioctl.h>
#include <sys/termios.h>
#define TCSETAW TIOCSETAW
#define TCGETA  TIOCGETA
#else
#ifdef _CONVEX_SOURCE
#include <sys/ioctl.h>
#include "/usr/sys/base/h/ioctl.h"
#define TCSETAW TIOCSETAW
#define TCGETA  TIOCGETA
#else
#if defined(linux)
#include <termios.h>
#include <unistd.h>
#else
#include <sys/termio.h>
#endif /* linux */
#endif /* _CONVEX_SOURCE */
#endif /* __FreeBSD__ */
#endif /* esv */

#ifdef esv
#include <sysv/unistd.h>
#else
#include <unistd.h>
#endif

#if defined(_SEQUENT_) || defined(_AIX)
#include <sys/select.h>
#endif

#ifdef __sgi
/* Avoid 'bzero' Compiler Warnings! */
#include <bstring.h>
#endif
#endif /* TERMIOS */

#ifdef SOCKETS
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#endif

#define TwoPi           6.28318531
#define IsIdentChar(x)  ((isalnum(x))||((x)=='_')||((x)=='$'))

#ifdef TERMIOS
#if defined(linux) || defined(__FreeBSD__)
static struct termios OrigTerm;
static struct termios IntrTerm;
#else
static struct termio OrigTerm;
static struct termio IntrTerm;
#endif

#if defined(linux)
static fd_set OrigWaitSet;
static fd_set WaitSet;
#else
static struct fd_set OrigWaitSet;
static struct fd_set WaitSet;
#endif
static struct timeval TimeOut;
static int WaitWidth;
static int FileNo;

#ifdef SOCKETS
/* Supported Protocols */
#define ProtoRasMol   0x01

#define AMNone         0x00
#define AMPickAtom     0x01
#define AMPickNumber   0x02
#define AMSelectCount  0x04
#define AMMolName      0x08
#define AMPickCoord    0x10

typedef struct {
    int bitmask;
    char *upper;
    char *name;
    } AdviseType;

static AdviseType AdviseMap[ItemCount] = {
    { AMPickAtom,    "PICK",    "Pick"    },  /* AdvPickAtom    */
    { AMPickNumber,  "PICKNO",  "PickNo"  },  /* AdvPickNumber  */
    { AMSelectCount, "COUNT",   "Count"   },  /* AdvSelectCount */
    { AMMolName,     "NAME",    "Name"    },  /* AdvName        */
    { AMNone,        "IDENT",   "Ident"   },  /* AdvIdent       */
    { AMNone,        "CLASS",   "Class"   },  /* AdvClass       */
    { AMNone,        "IMAGE",   "Image"   },  /* AdvImage       */
    { AMPickCoord,   "PICKXYZ", "PickXYZ" }   /* AdvPickCoord   */
        };

typedef struct {
        int protocol;
        int socket;
        int advise;
    } IPCConv;

#define MaxConvNum     8
static IPCConv ConvData[MaxConvNum];

static int ServerPort;
static int UseSockets;
static int SocketNo;
#endif  /* SOCKETS */
#endif  /* TERMIOS */

static int InitialWide;
static int InitialHigh;

static char *FileNamePtr;
static char *ScriptNamePtr;
static int LabelOptFlag;
static int FileFormat;
static int ProfCount;
static int LexState;

/* Function Prototype */
#ifdef FUNCPROTO
int HandleEvents( int );
#else
int HandleEvents();
#endif
int ProcessCommand();
void RasMolExit();

/* Either stdout or stderr */
#define OutFp stdout


void WriteChar( ch )
    char ch;
{   putc(ch,OutFp);
}


void WriteString( ptr )
    char *ptr;
{   fputs(ptr,OutFp);
}


static void ResetTerminal()
{
#ifdef SOCKETS
    register int i;
#endif

#ifdef TERMIOS
    if( isatty(FileNo) )                    /* drain terminal */
#if defined(linux)
      tcdrain ( FileNo );
#else
      ioctl(FileNo, TCSETAW, &OrigTerm);
#endif
#endif

#ifdef SOCKETS
    if( UseSockets )
    {   close(SocketNo);

        for( i=0; i<MaxConvNum; i++ )
            if( ConvData[i].protocol )
            {   close(ConvData[i].socket);
                ConvData[i].protocol = 0;
            }
    }
#endif
}


void RasMolExit()
{
    WriteChar('\n');
    if( CommandActive )
        WriteChar('\n');

    if( Interactive )
        CloseDisplay();
    ResetTerminal();
    exit(0);
}


void RasMolFatalExit( msg )
    char *msg;
{
    WriteChar('\n');
    WriteString(msg);
    WriteChar('\n');
    WriteChar(0x07);

    if( Interactive )
        CloseDisplay();
    ResetTerminal();
    exit(1);
}


#ifdef SOCKETS
static int OpenSocket()
{
    struct sockaddr_in addr;
    auto int length;
    register int i;

    UseSockets = False;
    SocketNo = socket(AF_INET, SOCK_STREAM, 0);
    if( SocketNo < 0 ) return( False );

    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(ServerPort);

    if( bind(SocketNo, (struct sockaddr*)&addr, sizeof(addr)) )
    {   close(SocketNo);
        return( False );
    }

    if( !ServerPort )
    {   length = sizeof(addr);
        if( !getsockname(SocketNo, (struct sockaddr*)&addr, &length) )
        {   ServerPort = ntohs(addr.sin_port);
            fprintf(stderr,"RasMol Server TCP/IP Port: %d\n",ServerPort);
        }
    }

    UseSockets = True;
    for( i=0; i<MaxConvNum; i++ )
        ConvData[i].protocol = 0;

    listen( SocketNo, 5 );
    return( True );
}


static int OpenConnection( socket )
    int socket;
{
    register int i;

    for( i=0; i<MaxConvNum; i++ )
        if( !ConvData[i].protocol )
        {   FD_SET(socket,&OrigWaitSet);
            if( socket >= WaitWidth )
                WaitWidth = socket+1;

            ConvData[i].protocol = ProtoRasMol;
            ConvData[i].socket = socket;
            ConvData[i].advise = AMNone;
            return( True );
        }

    close( socket );
    return( False );
}


static void CloseConnection( conv )
    int conv;
{
    FD_CLR(ConvData[conv].socket,&OrigWaitSet);
    close( ConvData[conv].socket );
    ConvData[conv].protocol = 0;
}


static int IsAdviseRequest( ptr, conv )
    char *ptr;  int conv;
{
    static char item[34];
    register char *dst;
    register char *src;
    register int i,ch;

    if( !strncmp(ptr,"Advise:",7) )
    {   src = ptr+7;

        while( True )
        {   ch = *src++;
            if( isspace(ch) )
                continue;

            if( isalpha(ch) )
            {   dst = item;
                *dst++ = ToUpper(ch);
                while( IsIdentChar(*src) )
                {   if( dst < item+32 )
                    {   ch = *src++;
                        *dst++ = ToUpper(ch);
                    } else src++;
                }
                *dst = '\0';

                for( i=0; i<ItemCount; i++ )
                    if( !strcmp(item,AdviseMap[i].upper) )
                    {   ConvData[conv].advise |= AdviseMap[i].bitmask;
                        break;
                    }

                /* Warning: Unknown Advise Item! */
            } else if( ch != ',' )
                break;
        }
        return( True );
    }
    return( False );
}


static void HandleSocketData( conv )
    int conv;
{
    register char *src,*dst;
    register int result;
    register int ch,len;
    char buffer[4097];

    len = read( ConvData[conv].socket, buffer, 4096 );
    if( len > 0 )
    {   buffer[len] = '\0';
        src = dst = buffer;
        while( (ch = *src++) )
            if( (ch>=' ') && (ch<='~') )
                *dst++ = ch;
        *dst = '\0';

        if( !IsAdviseRequest(buffer,conv) )
        {   result = ExecuteIPCCommand(buffer);
            if( result == IPC_Exit )
            {   CloseConnection( conv );
            } else if( result == IPC_Quit )
                RasMolExit();
        }
    } else CloseConnection( conv );
}
#endif /* SOCKETS */


static void InitTerminal(sockets)
    int sockets;
{
#ifdef TERMIOS
    register int i;
#endif

#ifdef SIGTTIN
    signal(SIGTTIN,SIG_IGN);
#endif
#ifdef SIGTTOU
    signal(SIGTTOU,SIG_IGN);
#endif

#ifdef TERMIOS
    FileNo = fileno(stdin);
    FD_ZERO(&OrigWaitSet);
    FD_SET(FileNo,&OrigWaitSet);
    WaitWidth = FileNo+1;

#ifdef SOCKETS
    OpenSocket();
    if( UseSockets )
    {   FD_SET(SocketNo,&OrigWaitSet);
        if( SocketNo >= WaitWidth )
            WaitWidth = SocketNo+1;
    }
#endif

    for( i=0; i<32; i++ )
        if( sockets & (1<<i) )
        {   if( i >= WaitWidth )
                WaitWidth = i+1;
            FD_SET(i,&OrigWaitSet);
        }


    if( isatty(FileNo) )
      {   
#if defined(linux)
      tcgetattr(FileNo, &OrigTerm);
#else
      ioctl(FileNo, TCGETA, &OrigTerm);
#endif
      IntrTerm = OrigTerm;

      IntrTerm.c_iflag |= IGNBRK|IGNPAR;
      IntrTerm.c_iflag &= ~(BRKINT|PARMRK|INPCK|IXON|IXOFF);
      IntrTerm.c_lflag &= ~(ICANON|ECHO|ECHOE|ECHOK|ECHONL|NOFLSH);
      IntrTerm.c_lflag |= ISIG;
      
      IntrTerm.c_cc[VMIN] = 1;
      IntrTerm.c_cc[VTIME] = 0;

#ifdef VSUSP /* Disable ^Z */
      IntrTerm.c_cc[VSUSP] = 0;
#endif

#if defined(linux)
      tcdrain ( FileNo );
#else
      ioctl(FileNo, TCSETAW, &IntrTerm);
#endif
    }
#endif /* TERMIOS */

    setbuf(stdin,(char*)NULL);
}


static int FetchCharacter()
{
#ifdef TERMIOS
    register int status;
#ifdef SOCKETS
    register int i;
#endif

#ifdef SOCKETS
    if( Interactive || UseSockets )
#else
    if( Interactive )
#endif
        do {
            if( !CommandActive )
                ResetCommandLine(0);

            if( Interactive )
            {   HandleEvents(False);

                /* avoid slow response time */
                if( !CommandActive )
                    ResetCommandLine(0);
            }

            TimeOut.tv_sec = 1;
            TimeOut.tv_usec = 0;
            WaitSet = OrigWaitSet;
#ifdef __hpux
            status = select( WaitWidth, (int*)&WaitSet, (int*)NULL, 
                                        (int*)NULL, &TimeOut );
#elif defined(linux)
            status = select( WaitWidth, (fd_set*)&WaitSet, (fd_set*)NULL, 
                                        (fd_set*)NULL, &TimeOut );
#else
            status = select( WaitWidth, &WaitSet, (struct fd_set*)NULL, 
                                        (struct fd_set*)NULL, &TimeOut );
#endif

#ifdef SOCKETS
            if( UseSockets )
            {   if( FD_ISSET(SocketNo,&WaitSet) )
                {   OpenConnection( accept(SocketNo,0,0) );
                } else for( i=0; i<MaxConvNum; i++ )
                    if( ConvData[i].protocol )
                        if( FD_ISSET(ConvData[i].socket,&WaitSet) )
                            HandleSocketData( i );
            }
#endif
        } while( (status<1) || !FD_ISSET(FileNo,&WaitSet) );
#endif /* TERMIOS */

    if( !CommandActive )
        ResetCommandLine(0);

    return( getc(stdin) );
}


static int ReadCharacter()
{
    register int tmp;
    register int ch;

    if( LexState )
    {   ch = LexState;
        LexState = 0;
        return( ch );
    }

    ch = FetchCharacter();
    if( ch!=0x1b ) return( ch );

    ch = FetchCharacter();
    if( (ch!='[') && (ch!='O') ) 
        return( ch );

    switch( tmp=FetchCharacter() )
    {   case('A'): return( 0x10 );
        case('B'): return( 0x0e );
        case('C'): return( 0x06 );
        case('D'): return( 0x02 );
    }
    LexState = tmp;
    return(ch);
}


void RasMolSignalExit( i )
    int i;
{
    RasMolFatalExit("*** Quit ***");
}


static void LoadInitFile()
{
    register char *src,*dst;
    register FILE *initrc;
    register char *fname;
    char fnamebuf[128];

    fname = ".rasmolrc";

    initrc = fopen(fname,"r");
    if( !initrc && (src=(char*)getenv("HOME")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '/';

        src = fname; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( !initrc && (src=(char*)getenv("RASMOLPATH")) )
    {   dst = fnamebuf; 
        while( *src )
            *dst++ = *src++;
        *dst++ = '/';

        src = "rasmolrc"; fname = fnamebuf;
        while( (*dst++ = *src++) );
        initrc = fopen(fname,"r");
    }

    if( initrc )
        LoadScriptFile(initrc,fname);
}


static void HandleMenu( hand )
     int hand;
{
    register int menu;
    register int item;
    register int mask;

    menu = hand>>8;
    item = hand&0xff;
    switch( menu )
    {   case(0):  /* File Menu */
                  switch( item )
                  {   case(1):  /* Open */
                                if( !Database )
                                    ResetCommandLine(2);
                                break;

                      case(2):  /* Save As */
                                if( Database )
                                    ResetCommandLine(4);
                                break;

                      case(3):  /* Close */
                                ZapDatabase();
                                break;

                      case(5):  /* Exit */
                                RasMolExit();
                                break;
                  } 
                  break;

        case(1):  /* Display Menu */
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

        case(2):  /* Colours Menu */
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

        case(3):  /* Option Menu */
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

        case(4):  /* Export Menu */
                  ResetCommandLine(3);
                  StateOption = item;
                  break;

        case(5):  /* Help Menu */
                  break;
    }
}


void RefreshScreen()
{
    if( !UseSlabPlane )
    {   ReDrawFlag &= ~(RFTransZ|RFSlab|RFPoint);
    } else ReDrawFlag &= ~(RFTransZ|RFPoint);

    if( ReDrawFlag )
    {   if( ReDrawFlag & RFReSize )
            ReSizeScreen();

        if( ReDrawFlag & RFColour )
        {   if( Interactive ) 
                ClearImage();
            DefineColourMap();
        }

        if( Database )
        {   if( Interactive )
                BeginWait();
            if( ReDrawFlag & RFApply ) 
                ApplyTransform();
            DrawFrame();
            if( Interactive )
            {   TransferImage();
                EndWait();
            }
        } else if( Interactive )
        {   ClearBuffers();
            TransferImage();
        }
        ReDrawFlag = 0;
    }
}


#ifdef SOCKETS
static char AdviseBuffer[256];
static int AdviseLen;

static void PrepareAdviseItem( item )
    int item;
{
    register char *src, *dst;
    register int i,flag;

    dst = AdviseBuffer;
    src = AdviseMap[item].name;
    while( *src ) *dst++ = *src++;
    *dst++ = ':';  *dst++ = ' ';

    switch( item )
    {   case(AdvPickAtom):
                  if( QAtom )
                  {   src = Residue[QGroup->refno];
                      flag = False;
                      for( i=0; i<3; i++ )
                          if( (src[i]!=' ') && !isalpha(src[i]) )
                              flag = True;

                      if( flag ) *dst++ = '[';
                      for( i=0; i<3; i++ )
                          if( src[i]!=' ' )
                              *dst++ = src[i];
                      if( flag ) *dst++ = ']';
                      sprintf(dst,"%d",QGroup->serno);
                      for( dst=AdviseBuffer; *dst; dst++ );
                      if( QChain->ident!=' ' )
                      {   if( isdigit(QChain->ident) ) *dst++ = ':';
                          *dst++ = QChain->ident;
                      }
                      *dst++ = '.';

                      src = ElemDesc[QAtom->refno];
                      for( i=0; i<4; i++ )
                          if( src[i]!=' ' ) 
                              *dst++ = src[i];
                  }
                  *dst++ = '\n';
                  *dst = '\0';
                  break;

        case(AdvPickNumber):
                  if( !QAtom )
                  {   *dst++ = '\n'; *dst = '\0';
                  } else sprintf(dst,"%ld\n",(long)QAtom->serno);
                  break;

        case(AdvSelectCount):
                  sprintf(dst,"%ld\n",(long)SelectCount);
                  break;
                 
        case(AdvName):
                  src = Info.moleculename;
                  while( *src ) *dst++ = *src++;
                  *dst++ = '\n'; *dst = '\0';
                  break;

        case(AdvPickCoord):
                  if( !QAtom )
                  {   *dst++ = '\n'; *dst = '\0';
                  } else sprintf( dst, "%ld\t%ld\t%ld\n",
                             (long)QAtom->xorg, 
                             (long)QAtom->yorg, 
                             (long)QAtom->zorg);
                  break;

        default:  *dst++ = '\n';
                  *dst = '\0';
                  break;
    }
    AdviseLen = strlen(AdviseBuffer)+1;
}
#endif


void AdviseUpdate( item )
    int item;
{
#ifdef SOCKETS
    register int mask;
    register int i;

    AdviseLen = 0;
    if( UseSockets && (mask=AdviseMap[item].bitmask) )
    {   for( i=0; i<MaxConvNum; i++ )
            if( ConvData[i].protocol && (ConvData[i].advise&mask) )
            {   if( !AdviseLen ) PrepareAdviseItem(item);
                write(ConvData[i].socket,AdviseBuffer,AdviseLen);
            }
    }
#endif
}


int ProcessCommand()
{
    switch(CurState)
    {   case(1):  /* RasMol Prompt */
                  return( ExecuteCommand() );

        case(2):  /* XYZ Filename */
                  if( *CurLine && FetchFile(FormatXYZ,False,CurLine) )
                      DefaultRepresentation();
                  ResetCommandLine(1);
                  break;

        case(3):  /* Export Image Filename */
                  if( *CurLine ) switch( StateOption )
                  {   case(1):   WriteGIFFile(CurLine);            break;
                      case(2):   WriteEPSFFile(CurLine,True,True); break;
                      case(3):   WriteIRISFile(CurLine);           break;
                      case(4):   WritePICTFile(CurLine);           break;
                  }
                  ResetCommandLine(1);
                  break;

        case(4):  /* Save Molecule Filename */
                  if( *CurLine )
                      SaveXYZMolecule(CurLine);
                  ResetCommandLine(1);
                  break;
    }
    return( False );
}


int HandleEvents( wait )
    int wait;
{
    register int result;

    result = FetchEvent( wait );
    while( ReDrawFlag || result )
    {   if( !result )
        {   if( ReDrawFlag&RFPoint )
            {   if( Database )
                {   if( ReDrawFlag & RFPoint1 )
                    {      PickAtom(True,PointX,PointY);
                    } else PickAtom(False,PointX,PointY);
#ifdef SOCKETS
                    AdviseUpdate(AdvPickNumber);
                    AdviseUpdate(AdvPickAtom);
                    AdviseUpdate(AdvPickCoord);
#endif
                }
                ReDrawFlag &= ~RFPoint;
            }

            if( ReDrawFlag )
                RefreshScreen();
        } else if( !IsPaused )
            HandleMenu( result );
        result = FetchEvent( False );
    }
    return( True );
}


static void ProfileExecution()
{
#ifdef TIME
    register long start,stop;
#else
    static struct timeval start;
    static struct timeval stop;
    register double secs;
#endif
    register Real delta;
    register int i;

    delta = TwoPi/ProfCount;

    printf("Profiling Execution!\n");

#ifdef TIME
    start = time((time_t *)NULL);
#else
    gettimeofday(&start,(struct timezone *)NULL);
#endif

    for( i=0; i<ProfCount; i++ )
    {   DrawFrame();
        if( Interactive )
            TransferImage();
        ReDrawFlag |= RFRotateY;
        DialValue[1] += delta;
        /* UpdateScrollBars(); */
        ApplyTransform();
    }

#ifdef TIME
    stop = time((time_t *)NULL);
    fprintf(stderr,"Execution of %d frames\n",ProfCount);
    fprintf(stderr,"Duration = %ld seconds\n",stop-start);
#else
    gettimeofday(&stop,(struct timezone *)NULL);
    secs = (stop.tv_sec - start.tv_sec) + (double)
           (stop.tv_usec - start.tv_usec)/1000000.0;
    fprintf(stderr,"Execution of %d frames\n",ProfCount);
    fprintf(stderr,"Duration = %g seconds\n",secs);
#endif
}


static void InitDefaultValues()
{
    Interactive = True;

    LabelOptFlag = False;

    FileNamePtr = NULL;
    ScriptNamePtr = NULL;
    InitialWide = DefaultWide;
    InitialHigh = DefaultHigh;
    ServerPort = 21069;
    ProfCount = 0;

    FileFormat = FormatXYZ;

    CalcBondsFlag = True;
}


static void DisplayUsage()
{
    fputs("usage: rasmol [-nodisplay] [-script scriptfile] ",OutFp);
    fputs("[[-format] file]\n    formats: -xyz \n\n",OutFp);
    exit(1);
}


#define FORMATOPTMAX   1
static struct {
        char *ident;
        int format;
    } FormatOpt[FORMATOPTMAX] = { 
            { "xyz",        FormatXYZ      }
                                };


static void ProcessOptions(argc,argv)
    int argc;  char *argv[];
{
    register char *ptr;
    register int i,j;

    for( i=1; i<argc; i++ )
    {   ptr = argv[i];
        if( (*ptr=='-') && ptr[1] )
        {   ptr++;

            if( !strcmp(ptr,"nodisplay") )
            {   Interactive = False;
            } else if( !strcmp(ptr,"prof") ||
                       !strcmp(ptr,"profile") )
            {   ProfCount = 200;
            } else if( !strcmp(ptr,"noconnect") )
            {   CalcBondsFlag = False;
            } else if( !strcmp(ptr,"connect") )
            {   CalcBondsFlag = True;

            } else if( !strcmp(ptr,"script") )
            {   if( i == argc-1 ) DisplayUsage();
                ScriptNamePtr = argv[++i];
            } else if( !strcmp(ptr,"width") ||
                       !strcmp(ptr,"wide") )
            {   if( i == argc-1 ) DisplayUsage();
                InitialWide = atoi(argv[++i]);
                if( InitialWide < 48 )
                    InitialWide = 48;
            } else if( !strcmp(ptr,"height") ||
                       !strcmp(ptr,"high") )
            {   if( i == argc-1 ) DisplayUsage();
                InitialHigh = atoi(argv[++i]);
                if( InitialHigh < 48 )
                    InitialHigh = 48;

#ifdef SOCKETS
            } else if( !strcmp(ptr,"port") )
            {   if( i == argc-1 ) DisplayUsage();
                ServerPort = atoi(argv[++i]);
#endif

            } else  /* File Formats! */
            {   for( j=0; j<FORMATOPTMAX; j++ )
                    if( !strcmp(ptr,FormatOpt[j].ident) )
                    {   FileFormat = FormatOpt[j].format;
                        break;
                    }

                if( j==FORMATOPTMAX )
                    DisplayUsage();
            }
        } else
            if( !FileNamePtr )
            {   FileNamePtr = ptr;
            } else DisplayUsage();
    }
}


int main(argc,argv)
int argc; char *argv[];
{
    register FILE *fp;
    register int temp;
    register int done;
    register char ch;

    InitDefaultValues();
    ProcessOptions(argc,argv);
    ReDrawFlag = 0;
    
    temp = Interactive;
    setbuf(OutFp,(char *)NULL);
    Interactive = OpenDisplay(InitialWide,InitialHigh);
    InitTerminal(Interactive);

    signal(SIGINT,RasMolSignalExit);
    signal(SIGPIPE,SIG_IGN);

    WriteString("TINKER RasMol Molecule Viewer\n");
    WriteString("Version 4.0, October 2002\n");
    WriteString("Based on RasMol 2.6 by Roger Sayle\n\n");

#ifdef EIGHTBIT
    WriteString("[8bit version]\n\n");
#else
#ifdef SIXTEENBIT
    WriteString("[16bit version]\n\n");
#else
    WriteString("[32bit version]\n\n");
#endif
#endif

    if( !Interactive )
    {   if( temp )
        {   WriteString("No suitable display detected!\n");
        } else WriteString("Display window disabled!\n");
    }

    InitialiseCommand();
    InitialiseTransform();
    InitialiseDatabase();
    InitialiseRenderer();
    InitialisePixUtils();
    InitialiseAbstree();
    InitialiseOutFile();
    InitialiseRepres();

    if( ProfCount )
    {   if( FileNamePtr )
        {   strcpy(DataFileName,FileNamePtr);

            if( strcmp(FileNamePtr,"-") )
            {   done = FetchFile(FileFormat,True,FileNamePtr);
            } else done = ProcessFile(FileFormat,True,stdin);
            if( !done )
                RasMolFatalExit("Profile Error: Unable to read data file!");
        } else RasMolFatalExit("Profile Error: No molecule filename!");

        ReDrawFlag |= RFRefresh | RFColour;

        /* SetVanWaalRadius(); */
        /* CPKColourAttrib(); */

        FakeSpecular = True;

        /* Default is WireFrame Representation
        EnableWireframe(WireFlag,0); */

        /* Default is Sticks Representation */
        if( MainAtomCount<256 )
        {   EnableWireframe(CylinderFlag,40);
        } else EnableWireframe(CylinderFlag,80);

        /* Default is Ball & Stick Representation
        SetRadiusValue(120);
        EnableWireframe(CylinderFlag,40); */

        RefreshScreen();

        /* Avoid Pending Events */
        if( Interactive ) 
        {   FetchEvent(False);
            if( ReDrawFlag )
                RefreshScreen();
        }
        ProfileExecution();
        RasMolExit();
    }

    if( FileNamePtr )
    {   strcpy(DataFileName,FileNamePtr);

        if( !strcmp(FileNamePtr,"-") )
        {   done = ProcessFile(FileFormat,True,stdin);
        } else done = FetchFile(FileFormat,True,FileNamePtr);

        if( done )
        {   DefaultRepresentation();

    if( MainAtomCount<256 )
    {   EnableWireframe(CylinderFlag,40);
    } else EnableWireframe(CylinderFlag,80);

            RefreshScreen();
        }
    }

    LexState = 0;
    ResetCommandLine(1);

    LoadInitFile();
    if( ScriptNamePtr )
    {   if( !(fp=fopen(ScriptNamePtr,"r")) )
        {   fprintf(OutFp,"Error: File '%s' not found!\n",ScriptNamePtr);
        } else LoadScriptFile(fp,ScriptNamePtr);
    }

    if( FileNamePtr && !strcmp(FileNamePtr,"-") )
    {   /* Finished Processing after stdin? */

#ifdef TERMIOS
        if( Interactive )
        {   FD_CLR(FileNo,&OrigWaitSet);
            fclose(stdin);
        } else RasMolExit();
#else
        if( !Interactive )
            RasMolExit();
#endif
    }

#ifdef TERMIOS
    while( True )
    {   ch = ReadCharacter();
        if( ProcessCharacter(ch) )
            if( ProcessCommand() )
                break;
        RefreshScreen();
    }
#else /* !TERMIOS */
    if( Interactive )
    {   /* X Windows Only! */
        while( HandleEvents(True) )
            if( !CommandActive )
                ResetCommandLine(0);
    } else while( True )
    {   /* Command Line Only! */
        ch = ReadCharacter();
        if( ProcessCharacter(ch) )
            if( ProcessCommand() )
                break;
        RefreshScreen();
    }
#endif
    RasMolExit();
    return( 0 );
}
