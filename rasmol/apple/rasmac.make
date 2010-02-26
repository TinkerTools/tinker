#   File:       RasMol.symantec.make
#   Target:     RasMol.symantec
#   Sources:
#               abstree.c applemac.c command.c molecule.c outfile.c pixutils.c  
#		rasmac.c render.c script.c transfor.c 
#   Created:    Tue Oct 10 15:08:59 1995


PPC_C = SMrC
PPC_OPTIONS = 
OPTIONS = -ansi relaxed -align power -opt all
C_OPTIONS = {OPTIONS} 



COBJECTS = abstree.c.o applemac.c.o command.c.o infile.c.o ¶
           molecule.c.o outfile.c.o pixutils.c.o rasmac.c.o ¶
		   render.c.o repres.c.o script.c.o transfor.c.o 

RasMol.symantec Ä RasMol.symantec.make RasMac.rsrc {COBJECTS} 
	PPCLink   -warn  {COBJECTS} ¶
		"{SharedLibraries}"InterfaceLib ¶
		"{SharedLibraries}"MathLib  ¶
		"{SharedLibraries}"StdCLib  ¶
		"{PPCLibraries}"StdCRuntime.o  ¶
		"{PPCLibraries}"PPCCRuntime.o ¶
		-o RasMol.symantec
	Echo "include ¶"RasMac.rsrc¶" ;" ¶
		| Rez -a  -o RasMol.symantec
	SetFile -c RSML RasMol.symantec


abstree.c.o Ä abstree.c
	{PPC_C} abstree.c {PPC_OPTIONS} {C_OPTIONS}

applemac.c.o Ä applemac.c
	{PPC_C} applemac.c {PPC_OPTIONS} {C_OPTIONS}

command.c.o Ä command.c
	{PPC_C} command.c {PPC_OPTIONS} {C_OPTIONS}

infile.c.o Ä infile.c
	{PPC_C} infile.c {PPC_OPTIONS} {C_OPTIONS}

molecule.c.o Ä molecule.c
	{PPC_C} molecule.c {PPC_OPTIONS} {C_OPTIONS}

outfile.c.o Ä  outfile.c
	{PPC_C} outfile.c {PPC_OPTIONS} {C_OPTIONS}

pixutils.c.o Ä  pixutils.c
	{PPC_C} pixutils.c {PPC_OPTIONS} {C_OPTIONS}

rasmac.c.o Ä  rasmac.c
	{PPC_C} rasmac.c {PPC_OPTIONS} {C_OPTIONS}

render.c.o Ä  render.c
	{PPC_C} render.c {PPC_OPTIONS} {C_OPTIONS}

repres.c.o Ä  repres.c
	{PPC_C} repres.c {PPC_OPTIONS} {C_OPTIONS}

script.c.o Ä  script.c
	{PPC_C} script.c {PPC_OPTIONS} {C_OPTIONS}

transfor.c.o Ä  transfor.c
	{PPC_C} transfor.c {PPC_OPTIONS} {C_OPTIONS}

