# Microsoft Visual C++ Generated NMAKE File, Format Version 2.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

!IF "$(CFG)" == ""
CFG=Win32 Debug
!MESSAGE No configuration specified.  Defaulting to Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "Win32 Release" && "$(CFG)" != "Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "raswin32.mak" CFG="Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

################################################################################
# Begin Project
# PROP Target_Last_Scanned "Win32 Debug"
MTL=MkTypLib.exe
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Win32 Release"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "WinRel"
# PROP BASE Intermediate_Dir "WinRel"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "WinRel"
# PROP Intermediate_Dir "WinRel"
OUTDIR=.\WinRel
INTDIR=.\WinRel

ALL : $(OUTDIR)/raswin32.exe $(OUTDIR)/raswin32.bsc

$(OUTDIR) : 
    if not exist $(OUTDIR)/nul mkdir $(OUTDIR)

# ADD BASE MTL /nologo /D "NDEBUG" /win32
# ADD MTL /nologo /D "NDEBUG" /win32
MTL_PROJ=/nologo /D "NDEBUG" /win32 
# ADD BASE CPP /nologo /Zp4 /MT /W3 /vmg /vms /vd0 /GX /YX /Ox /Ot /Oa /Og /Oi /Ob2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR /c
# ADD CPP /nologo /Zp4 /ML /W3 /vmg /vms /vd0 /GX /YX /Ox /Ot /Oa /Og /Oi /Ob2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR /c
CPP_PROJ=/nologo /Zp4 /ML /W3 /vmg /vms /vd0 /GX /YX /Ox /Ot /Oa /Og /Oi /Ob2\
 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /FR$(INTDIR)/\
 /Fp$(OUTDIR)/"raswin32.pch" /Fo$(INTDIR)/ /c 
CPP_OBJS=.\WinRel/
# ADD BASE RSC /l 0x809 /d "NDEBUG"
# ADD RSC /l 0x809 /d "NDEBUG"
RSC_PROJ=/l 0x809 /fo$(INTDIR)/"RASWIN.res" /d "NDEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o$(OUTDIR)/"raswin32.bsc" 
BSC32_SBRS= \
	$(INTDIR)/ABSTREE.SBR \
	$(INTDIR)/MOLECULE.SBR \
	$(INTDIR)/COMMAND.SBR \
	$(INTDIR)/MSWIN31.SBR \
	$(INTDIR)/OUTFILE.SBR \
	$(INTDIR)/PIXUTILS.SBR \
	$(INTDIR)/RASWIN.SBR \
	$(INTDIR)/SCRIPT.SBR \
	$(INTDIR)/RENDER.SBR \
	$(INTDIR)/TRANSFOR.SBR

$(OUTDIR)/raswin32.bsc : $(OUTDIR)  $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 oldnames.lib /NOLOGO /STACK:0x10240 /SUBSYSTEM:windows /MACHINE:IX86
# ADD LINK32 oldnames.lib user32.lib gdi32.lib comdlg32.lib shell32.lib /NOLOGO /STACK:0x10240 /SUBSYSTEM:windows /MACHINE:IX86
LINK32_FLAGS=oldnames.lib user32.lib gdi32.lib comdlg32.lib shell32.lib /NOLOGO\
 /STACK:0x10240 /SUBSYSTEM:windows /INCREMENTAL:no /PDB:$(OUTDIR)/"raswin32.pdb"\
 /MACHINE:IX86 /DEF:".\RASWIN32.DEF" /OUT:$(OUTDIR)/"raswin32.exe" 
DEF_FILE=.\RASWIN32.DEF
LINK32_OBJS= \
	$(INTDIR)/ABSTREE.OBJ \
	$(INTDIR)/MOLECULE.OBJ \
	$(INTDIR)/COMMAND.OBJ \
	$(INTDIR)/MSWIN31.OBJ \
	$(INTDIR)/OUTFILE.OBJ \
	$(INTDIR)/PIXUTILS.OBJ \
	$(INTDIR)/RASWIN.OBJ \
	$(INTDIR)/SCRIPT.OBJ \
	$(INTDIR)/RENDER.OBJ \
	$(INTDIR)/TRANSFOR.OBJ \
	$(INTDIR)/RASWIN.res

$(OUTDIR)/raswin32.exe : $(OUTDIR)  $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "Win32 Debug"

# PROP BASE Use_MFC 1
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "WinDebug"
# PROP BASE Intermediate_Dir "WinDebug"
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "WinDebug"
# PROP Intermediate_Dir "WinDebug"
OUTDIR=.\WinDebug
INTDIR=.\WinDebug

ALL : $(OUTDIR)/raswin32.exe $(OUTDIR)/raswin32.bsc

$(OUTDIR) : 
    if not exist $(OUTDIR)/nul mkdir $(OUTDIR)

# ADD BASE MTL /nologo /D "_DEBUG" /win32
# ADD MTL /nologo /D "_DEBUG" /win32
MTL_PROJ=/nologo /D "_DEBUG" /win32 
# ADD BASE CPP /nologo /MT /W3 /GX /Zi /YX /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /Fd"RASWIN32.PDB" /c
# ADD CPP /nologo /ML /W3 /GX /Zi /YX /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /FR /Fd"RASWIN32.PDB" /c
CPP_PROJ=/nologo /ML /W3 /GX /Zi /YX /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /D "_MBCS" /FR$(INTDIR)/ /Fp$(OUTDIR)/"raswin32.pch" /Fo$(INTDIR)/\
 /Fd"RASWIN32.PDB" /c 
CPP_OBJS=.\WinDebug/
# ADD BASE RSC /l 0x809 /d "_DEBUG"
# ADD RSC /l 0x809 /d "_DEBUG"
RSC_PROJ=/l 0x809 /fo$(INTDIR)/"RASWIN.res" /d "_DEBUG" 
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o$(OUTDIR)/"raswin32.bsc" 
BSC32_SBRS= \
	$(INTDIR)/ABSTREE.SBR \
	$(INTDIR)/MOLECULE.SBR \
	$(INTDIR)/COMMAND.SBR \
	$(INTDIR)/MSWIN31.SBR \
	$(INTDIR)/OUTFILE.SBR \
	$(INTDIR)/PIXUTILS.SBR \
	$(INTDIR)/RASWIN.SBR \
	$(INTDIR)/SCRIPT.SBR \
	$(INTDIR)/RENDER.SBR \
	$(INTDIR)/TRANSFOR.SBR

$(OUTDIR)/raswin32.bsc : $(OUTDIR)  $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
# ADD BASE LINK32 oldnames.lib /NOLOGO /STACK:0x10240 /SUBSYSTEM:windows /DEBUG /MACHINE:IX86
# ADD LINK32 oldnames.lib user32.lib gdi32.lib comdlg32.lib shell32.lib /NOLOGO /STACK:0x10240 /SUBSYSTEM:windows /DEBUG /MACHINE:IX86
LINK32_FLAGS=oldnames.lib user32.lib gdi32.lib comdlg32.lib shell32.lib /NOLOGO\
 /STACK:0x10240 /SUBSYSTEM:windows /INCREMENTAL:yes\
 /PDB:$(OUTDIR)/"raswin32.pdb" /DEBUG /MACHINE:IX86 /DEF:".\RASWIN32.DEF"\
 /OUT:$(OUTDIR)/"raswin32.exe" 
DEF_FILE=.\RASWIN32.DEF
LINK32_OBJS= \
	$(INTDIR)/ABSTREE.OBJ \
	$(INTDIR)/MOLECULE.OBJ \
	$(INTDIR)/COMMAND.OBJ \
	$(INTDIR)/MSWIN31.OBJ \
	$(INTDIR)/OUTFILE.OBJ \
	$(INTDIR)/PIXUTILS.OBJ \
	$(INTDIR)/RASWIN.OBJ \
	$(INTDIR)/SCRIPT.OBJ \
	$(INTDIR)/RENDER.OBJ \
	$(INTDIR)/TRANSFOR.OBJ \
	$(INTDIR)/RASWIN.res

$(OUTDIR)/raswin32.exe : $(OUTDIR)  $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Group "Source Files"

################################################################################
# Begin Source File

SOURCE=.\ABSTREE.C
DEP_ABSTR=\
	.\RASMOL.H\
	.\MOLECULE.H\
	.\ABSTREE.H

$(INTDIR)/ABSTREE.OBJ :  $(SOURCE)  $(DEP_ABSTR) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\MOLECULE.C
DEP_MOLEC=\
	.\RASMOL.H\
	.\MOLECULE.H\
	.\COMMAND.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\RENDER.H

$(INTDIR)/MOLECULE.OBJ :  $(SOURCE)  $(DEP_MOLEC) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\COMMAND.C
DEP_COMMA=\
	.\RASMOL.H\
	.\COMMAND.H\
	.\TOKENS.H\
	.\MOLECULE.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\RENDER.H\
	.\GRAPHICS.H\
	.\PIXUTILS.H\
	.\OUTFILE.H\
	.\SCRIPT.H

$(INTDIR)/COMMAND.OBJ :  $(SOURCE)  $(DEP_COMMA) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\MSWIN31.C
DEP_MSWIN=\
	.\RASMOL.H\
	.\GRAPHICS.H

$(INTDIR)/MSWIN31.OBJ :  $(SOURCE)  $(DEP_MSWIN) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\OUTFILE.C
DEP_OUTFI=\
	.\RASMOL.H\
	.\OUTFILE.H\
	.\MOLECULE.H\
	.\COMMAND.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\RENDER.H\
	.\GRAPHICS.H\
	.\PIXUTILS.H\
	.\SCRIPT.H

$(INTDIR)/OUTFILE.OBJ :  $(SOURCE)  $(DEP_OUTFI) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\PIXUTILS.C
DEP_PIXUT=\
	.\RASMOL.H\
	.\PIXUTILS.H\
	.\GRAPHICS.H\
	.\MOLECULE.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\RENDER.H\
	.\FONT.H

$(INTDIR)/PIXUTILS.OBJ :  $(SOURCE)  $(DEP_PIXUT) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\RASWIN.C
DEP_RASWI=\
	.\RASMOL.H\
	.\RASWIN.IDM\
	.\MOLECULE.H\
	.\ABSTREE.H\
	.\GRAPHICS.H\
	.\PIXUTILS.H\
	.\TRANSFOR.H\
	.\COMMAND.H\
	.\RENDER.H\
	.\OUTFILE.H

$(INTDIR)/RASWIN.OBJ :  $(SOURCE)  $(DEP_RASWI) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\SCRIPT.C
DEP_SCRIP=\
	.\RASMOL.H\
	.\SCRIPT.H\
	.\MOLECULE.H\
	.\COMMAND.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\RENDER.H\
	.\GRAPHICS.H\
	.\PIXUTILS.H

$(INTDIR)/SCRIPT.OBJ :  $(SOURCE)  $(DEP_SCRIP) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\RENDER.C
DEP_RENDE=\
	.\RASMOL.H\
	.\MOLECULE.H\
	.\GRAPHICS.H\
	.\RENDER.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\COMMAND.H\
	.\PIXUTILS.H

$(INTDIR)/RENDER.OBJ :  $(SOURCE)  $(DEP_RENDE) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\TRANSFOR.C
DEP_TRANS=\
	.\RASMOL.H\
	.\MOLECULE.H\
	.\ABSTREE.H\
	.\TRANSFOR.H\
	.\COMMAND.H\
	.\RENDER.H\
	.\GRAPHICS.H

$(INTDIR)/TRANSFOR.OBJ :  $(SOURCE)  $(DEP_TRANS) $(INTDIR)

# End Source File
################################################################################
# Begin Source File

SOURCE=.\RASWIN.RC
DEP_RASWIN=\
	.\RASWIN.CUR\
	.\RASWIN.ICO\
	.\RASWIN.IDM

$(INTDIR)/RASWIN.res :  $(SOURCE)  $(DEP_RASWIN) $(INTDIR)
   $(RSC) $(RSC_PROJ)  $(SOURCE) 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\RASWIN32.DEF
# End Source File
# End Group
# End Project
################################################################################
