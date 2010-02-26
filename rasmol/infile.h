/* infile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

#ifdef INFILE

#else

#ifdef FUNCPROTO
int LoadPDBMolecule( FILE*, int );
int LoadMacroModelMolecule( FILE* );
int LoadAlchemyMolecule( FILE* );
int LoadCharmmMolecule( FILE* );
int LoadBiosymMolecule( FILE* );
int LoadMOPACMolecule( FILE* );
int LoadSHELXMolecule( FILE* );
int LoadMol2Molecule( FILE* );
int LoadFDATMolecule( FILE* );
int LoadMDLMolecule( FILE* );
int LoadXYZMolecule( FILE* );
int LoadCEXMolecule( FILE* );

int SaveAlchemyMolecule( char* );
int SavePDBMolecule( char* );
int SaveMDLMolecule( char* );
int SaveXYZMolecule( char* );
int SaveCIFMolecule( char* );
int SaveCEXMolecule( char* );

#else /* non-ANSI C compiler */
int LoadPDBMolecule();
int LoadMacroModelMolecule();
int LoadAlchemyMolecule();
int LoadCharmmMolecule();
int LoadMOPACMolecule();
int LoadMol2Molecule();
int LoadXYZMolecule();
int LoadMDLMolecule();

int SaveAlchemyMolecule();
int SaveMol2Molecule();
int SavePDBMolecule();
int SaveMDLMolecule();
int SaveXYZMolecule();
int SaveCIFMolecule();
int SaveCEXMolecule();

#endif
#endif
