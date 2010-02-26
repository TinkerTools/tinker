# Microsoft Visual C++ generated build script - Do not modify

PROJ = RASWIN
DEBUG = 0
PROGTYPE = 0
CALLER = 
ARGS = 
DLLS = 
D_RCDEFINES = /d_DEBUG 
R_RCDEFINES = /dNDEBUG 
ORIGIN = MSVC
ORIGIN_VER = 1.00
PROJPATH = C:\RASWIN\SOURCE\
USEMFC = 0
CC = cl
CPP = cl
CXX = cl
CCREATEPCHFLAG = 
CPPCREATEPCHFLAG = 
CUSEPCHFLAG = 
CPPUSEPCHFLAG = 
FIRSTC = ABSTREE.C   
FIRSTCPP =             
RC = rc
CFLAGS_D_WEXE = /nologo /G2 /W3 /Zi /AM /Od /D "_DEBUG" /FR /GA /Fd"RASWIN.PDB"
CFLAGS_R_WEXE = /nologo /f /G2 /FPi87 /Zp4 /Gy /W3 /Gf /vmg /vms /vd0 /AM /O2 /Oa /Oe /Og /Oi /Ol /Ot /Ox /Oz /Ob2 /OV9 /D "NDEBUG" /FR /GA 
LFLAGS_D_WEXE = /NOLOGO /NOD /PACKC:61440 /STACK:10240 /ALIGN:16 /ONERROR:NOEXE /CO  
LFLAGS_R_WEXE = /NOLOGO /NOD /PACKC:61440 /STACK:10240 /ALIGN:16 /ONERROR:NOEXE  
LIBS_D_WEXE = libw mlibcew commdlg.lib shell.lib  oldnames
LIBS_R_WEXE = libw mlibcew commdlg.lib shell.lib   oldnames
RCFLAGS = /nologo 
RESFLAGS = /nologo /k
RUNFLAGS = 
DEFFILE = RASWIN.DEF
OBJS_EXT = 
LIBS_EXT = 
!if "$(DEBUG)" == "1"
CFLAGS = $(CFLAGS_D_WEXE)
LFLAGS = $(LFLAGS_D_WEXE)
LIBS = $(LIBS_D_WEXE)
MAPFILE = nul
RCDEFINES = $(D_RCDEFINES)
!else
CFLAGS = $(CFLAGS_R_WEXE)
LFLAGS = $(LFLAGS_R_WEXE)
LIBS = $(LIBS_R_WEXE)
MAPFILE = nul
RCDEFINES = $(R_RCDEFINES)
!endif
!if [if exist MSVC.BND del MSVC.BND]
!endif
SBRS = ABSTREE.SBR \
		MOLECULE.SBR \
		COMMAND.SBR \
		MSWIN31.SBR \
		OUTFILE.SBR \
		PIXUTILS.SBR \
		RASWIN.SBR \
		SCRIPT.SBR \
		RENDER.SBR \
		TRANSFOR.SBR


ABSTREE_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\abstree.h


MOLECULE_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\render.h


COMMAND_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\command.h \
	c:\raswin\source\tokens.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\render.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\pixutils.h \
	c:\raswin\source\outfile.h \
	c:\raswin\source\script.h


MSWIN31_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\render.h


OUTFILE_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\outfile.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\render.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\pixutils.h \
	c:\raswin\source\script.h


PIXUTILS_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\pixutils.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\render.h \
	c:\raswin\source\font.h


RASWIN_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\raswin.idm \
	c:\raswin\source\molecule.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\pixutils.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\render.h \
	c:\raswin\source\outfile.h


SCRIPT_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\script.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\render.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\pixutils.h


RENDER_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\graphics.h \
	c:\raswin\source\render.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\pixutils.h


TRANSFOR_DEP = c:\raswin\source\rasmol.h \
	c:\raswin\source\transfor.h \
	c:\raswin\source\molecule.h \
	c:\raswin\source\command.h \
	c:\raswin\source\abstree.h \
	c:\raswin\source\render.h \
	c:\raswin\source\graphics.h


RASWIN_RCDEP = c:\raswin\source\raswin.idm \
	c:\raswin\source\raswin.cur \
	c:\raswin\source\raswin.ico


all:	$(PROJ).EXE $(PROJ).BSC

ABSTREE.OBJ:	ABSTREE.C $(ABSTREE_DEP)
	$(CC) $(CFLAGS) $(CCREATEPCHFLAG) /c ABSTREE.C

MOLECULE.OBJ:	MOLECULE.C $(MOLECULE_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c MOLECULE.C

COMMAND.OBJ:	COMMAND.C $(COMMAND_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c COMMAND.C

MSWIN31.OBJ:	MSWIN31.C $(MSWIN31_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c MSWIN31.C

OUTFILE.OBJ:	OUTFILE.C $(OUTFILE_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c OUTFILE.C

PIXUTILS.OBJ:	PIXUTILS.C $(PIXUTILS_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c PIXUTILS.C

RASWIN.OBJ:	RASWIN.C $(RASWIN_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c RASWIN.C

SCRIPT.OBJ:	SCRIPT.C $(SCRIPT_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c SCRIPT.C

RENDER.OBJ:	RENDER.C $(RENDER_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c RENDER.C

TRANSFOR.OBJ:	TRANSFOR.C $(TRANSFOR_DEP)
	$(CC) $(CFLAGS) $(CUSEPCHFLAG) /c TRANSFOR.C

RASWIN.RES:	RASWIN.RC $(RASWIN_RCDEP)
	$(RC) $(RCFLAGS) $(RCDEFINES) -r RASWIN.RC


$(PROJ).EXE::	RASWIN.RES

$(PROJ).EXE::	ABSTREE.OBJ MOLECULE.OBJ COMMAND.OBJ MSWIN31.OBJ OUTFILE.OBJ PIXUTILS.OBJ \
	RASWIN.OBJ SCRIPT.OBJ RENDER.OBJ TRANSFOR.OBJ $(OBJS_EXT) $(DEFFILE)
	echo >NUL @<<$(PROJ).CRF
ABSTREE.OBJ +
MOLECULE.OBJ +
COMMAND.OBJ +
MSWIN31.OBJ +
OUTFILE.OBJ +
PIXUTILS.OBJ +
RASWIN.OBJ +
SCRIPT.OBJ +
RENDER.OBJ +
TRANSFOR.OBJ +
$(OBJS_EXT)
$(PROJ).EXE
$(MAPFILE)
c:\msvc\lib\+
$(LIBS)
$(DEFFILE);
<<
	link $(LFLAGS) @$(PROJ).CRF
	$(RC) $(RESFLAGS) RASWIN.RES $@
	@copy $(PROJ).CRF MSVC.BND

$(PROJ).EXE::	RASWIN.RES
	if not exist MSVC.BND 	$(RC) $(RESFLAGS) RASWIN.RES $@

run: $(PROJ).EXE
	$(PROJ) $(RUNFLAGS)


$(PROJ).BSC: $(SBRS)
	bscmake @<<
/o$@ $(SBRS)
<<
