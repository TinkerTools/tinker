/* outfile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#ifdef OUTFILE
int UseTransparent;
int UseOutLine;

#else
extern int UseTransparent;
extern int UseOutLine;

#ifdef FUNCPROTO
int WriteVectPSFile( char* );
int WriteEPSFFile( char*, int, int );
int WriteRastFile( char*, int );
int WritePICTFile( char* );
int WriteIRISFile( char* );
int WritePPMFile( char*, int );
int WriteGIFFile( char* );
int WriteBMPFile( char* );
void InitialiseOutFile();

#else /* non-ANSI C compiler */
int WriteVectPSFile();
int WriteEPSFFile();
int WriteRastFile();
int WritePICTFile();
int WriteIRISFile();
int WritePPMFile();
int WriteGIFFile();
int WriteBMPFile();
void InitialiseOutFile();

#endif
#endif
