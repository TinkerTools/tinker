/* repres.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

typedef struct _Monitor {
        struct _Monitor *next;
        Atom __far *src;
        Atom __far *dst;
        unsigned short dist;
        short col;
    } Monitor;

typedef struct _Label {
        struct _Label *next;
        Long  refcount;
        char *label;
    } Label;

#ifdef REPRES
Monitor *MonitList;
Label *LabelList;

int ProbeRadius;
int DrawLabels;
int DrawMonitDistance;

#else
extern Monitor *MonitList;
extern Label *LabelList;

extern int ProbeRadius;
extern int DrawLabels;
extern int DrawMonitDistance;

#ifdef FUNCPROTO
int DeleteLabels();
void DeleteLabel( Label* );
Label *CreateLabel( char*, int );
void DefineLabels( char* );
void DefaultLabels( int );
void DisplayLabels();

void DeleteMonitors();
void AddMonitors( Atom __far*, Atom __far* );
void CreateMonitor( Long, Long );
void DisplayMonitors();

void DeleteSurface();
void CalculateSurface( int );
void DisplaySurface();

void ResetRepres();
void InitialiseRepres();

#else /* non-ANSI C compiler */
int DeleteLabels();
void DeleteLabel();
Label *CreateLabel();
void DefineLabels();
void DefaultLabels();
void DisplayLabels();

void DeleteMonitors();
void AddMonitors();
void CreateMonitor();
void DisplayMonitors();

void DeleteSurface();
void CalculateSurface();
void DisplaySurface();

void ResetRepres();
void InitialiseRepres();

#endif
#endif
