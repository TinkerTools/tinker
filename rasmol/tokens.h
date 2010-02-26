/* tokens.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 * TINKER Viewer Version
 */

/* Lexeme Tokens */
#define IdentTok       256
#define NumberTok      257
#define FloatTok       258
#define StringTok      259

/* Command Tokens */
#define AdviseTok      260
#define BackboneTok    261
#define CartoonTok     262
#define CentreTok      263
#define ClipboardTok   264
#define ColourTok      265
#define ConnectTok     266
#define DashTok        267
#define DefineTok      268
#define DisplayTok     269
#define EchoTok        270
#define ExitTok        271
#define HelpTok        272
#define LabelTok       273
#define LoadTok        274
#define MonitorTok     275
#define PrintTok       276
#define QuitTok        277
#define RefreshTok     278
#define RenumTok       279
#define ResetTok       280
#define ResizeTok      281
#define RestrictTok    282
#define RotateTok      283
#define SaveTok        284
#define ScriptTok      285
#define SelectTok      286
#define SetTok         287
#define ShowTok        288
#define SlabTok        289
#define SourceTok      290
#define SpacefillTok   291
#define StructureTok   292
#define SymmetryTok    293
#define TraceTok       294
#define TranslateTok   295
#define WaitTok        296
#define WireframeTok   297
#define WriteTok       298
#define ZapTok         299
#define ZoomTok        300

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=310) && ((x)<=348))
#define PredTokOrd(x)  ((x)-310)
#define PredTokChr(x)  ((x)+310)

#define AlphaTok       310
#define AminoTok       311
#define ATTok          312
#define BondedTok      313
#define CGTok          314
#define CystineTok     315
#define DNATok         316
#define HelixTok       317
#define HeteroTok      318
#define HydrogenTok    319
#define IonTok         320
#define LigandTok      321
#define MainChainTok   322
#define NucleicTok     323
#define ProteinTok     324
#define PurineTok      325
#define PyrimidineTok  326
#define RNATok         327
#define SelectedTok    328
#define SheetTok       329
#define SidechainTok   330
#define SolventTok     331
#define TurnTok        332
#define WaterTok       333

#define AcidicTok      334
#define AcyclicTok     335
#define AliphaticTok   336
#define AromaticTok    337
#define BasicTok       338
#define BuriedTok      339
#define ChargedTok     340
#define CyclicTok      341
#define HydrophobicTok 342
#define LargeTok       343
#define MediumTok      344
#define NeutralTok     345
#define PolarTok       346
#define SmallTok       347
#define SurfaceTok     348

/* Property Tokens */
#define IsPropTok(x)   (((x)>=350) && ((x)<=355))
#define TemperatureTok 350
#define RadiusTok      351
#define AtomNoTok      352
#define ElemNoTok      353
#define ModelTok       354
#define ResNoTok       355

/* File Format Tokens */
/* Warning! Tokens are related to Format values */
#define IsMoleculeFormat(x)  (((x)>=360) && ((x)<=375))

#define PDBTok         360
#define MacroModelTok  361
#define GaussianTok    362
#define AlchemyTok     363
#define NMRPDBTok      364
#define CharmmTok      365
#define BiosymTok      366
#define MOPACTok       367
#define SHELXTok       368
#define Mol2Tok        369
#define FDATTok        370
#define MMDBTok        371
#define MDLTok         372
#define XYZTok         373
#define CIFTok         374
#define CEXTok         375

/* Raster Tokens */
#define IsImageFormat(x) (((x)>=366) && ((x)<=389))
#define GIFTok         376
#define PPMTok         377
#define SUNTok         378
#define SUNRLETok      379
#define EPSFTok        380
#define PICTTok        381
#define IRISTok        382
#define BMPTok         383
#define MonoPSTok      384
#define VectPSTok      385
#define KinemageTok    386
#define MolScriptTok   387
#define POVRayTok      388
#define VRMLTok        389

/* Feature Tokens */
#define AtomTok        390
#define BondTok        391
#define DotsTok        392
#define HBondTok       393
#define RibbonTok      394
#define SSBondTok      395
#define Ribbon1Tok     396
#define Ribbon2Tok     397

/* Expression Tokens */
#define TrueTok        400
#define FalseTok       401
#define AllTok         402
#define NoneTok        403
#define AndTok         404
#define OrTok          405
#define NotTok         406
#define WithinTok      407
#define XorTok         408

/* Colour Tokens */
#define BlueTok        410
#define BlueTintTok    411
#define BlackTok       412
#define BrownTok       413
#define CyanTok        414
#define GoldTok        415
#define GrayTok        416
#define GreenTok       417
#define GreenblueTok   418
#define GreenTintTok   419
#define HotPinkTok     420
#define MagentaTok     421
#define OrangeTok      422
#define PinkTok        423
#define PinkTintTok    424
#define PurpleTok      425
#define RedTok         426
#define RedorangeTok   427
#define SeaTok         428
#define SkyTok         429
#define VioletTok      430
#define WhiteTok       431
#define YellowTok      432
#define YellowTintTok  433

#define CPKTok         434
#define ShapelyTok     435
#define ResidueTok     436
#define UserTok        437
#define GroupTok       438
#define ChainTok       439
#define TypeTok        440
#define PotentialTok   441
#define ChargeTok      442

/* Variable Tokens */
#define AmbientTok     450
#define AxesTok        451
#define BackFadeTok    452
#define BackgroundTok  453
#define BondModeTok    454
#define BoundBoxTok    455
#define DepthCueTok    456
#define FontSizeTok    457
#define HourGlassTok   458
#define MenusTok       459
#define MouseTok       460
#define PickingTok     461
#define ShadowTok      462
#define SlabModeTok    463
#define SpecularTok    464
#define SpecPowerTok   465
#define StrandsTok     466
#define TransparentTok 467
#define UnitCellTok    468

/* SlabMode Tokens */
#define RejectTok      470
#define HalfTok        471
#define HollowTok      472
#define SolidTok       473
#define SectionTok     474

/* MouseMode Tokens */
#define RasMolTok      475
#define InsightTok     476
#define QuantaTok      477
#define SybylTok       478

/* Information Tokens */
#define InfoTok        480
#define SequenceTok    481
#define VersionTok     482

/* Display Mode Tokens */
#define NormalTok      485
#define StereoTok      486
#define MonoTok        487
#define HardwareTok    488

/* Axis Tokens */
#define XTok           490
#define YTok           491
#define ZTok           492

/* Picking Tokens */
#define IdentifyTok    495
#define DistanceTok    496
#define AngleTok       497
#define TorsionTok     498

/* Misc Tokens */
#define InLineTok      500
#define VDWTok         501

typedef struct {
                char *ident;
                int token;
               } KeywordEntry;

#define MAXKEYLEN 11
static int KeyLen[MAXKEYLEN+1] = {
        0, 3, 8, 30, 64, 101, 152, 185, 214, 230, 234, 241 };

static KeywordEntry Keyword[] = {
            { "X",  XTok },
            { "Y",  YTok },
            { "Z",  ZTok },  /* 3 */

            { "AT", ATTok   },
            { "CG", CGTok   },
            { "ON", TrueTok },
            { "OR", OrTok   },
            { "PS", EPSFTok },  /* 8 */

            { "ALL", AllTok   },
            { "AND", AndTok   },
            { "BMP", BMPTok   },
            { "CEX", CEXTok   },
            { "CIF", CIFTok   },
            { "CPK", CPKTok   },
            { "DNA", DNATok   },
            { "GIF", GIFTok   },
            { "ION", IonTok   },
            { "MDL", MDLTok   },
            { "NOT", NotTok   },
            { "OFF", FalseTok },
            { "PDB", PDBTok   },
            { "PPM", PPMTok   },
            { "RED", RedTok   },
            { "RGB", IRISTok  },
            { "RNA", RNATok   },
            { "SET", SetTok   },
            { "SUN", SUNTok   },
            { "VDW", VDWTok   },
            { "XYZ", XYZTok   },
            { "ZAP", ZapTok   }, /* 30 */

            { "ATOM", AtomTok },
            { "AXES", AxesTok },
            { "AXIS", AxesTok },
            { "BLUE", BlueTok },
            { "BOND", BondTok },
            { "CYAN", CyanTok },
            { "DASH", DashTok },
            { "DOTS", DotsTok },
            { "ECHO", EchoTok },
            { "EPSF", EPSFTok },
            { "EXIT", ExitTok },
            { "FDAT", FDATTok },
            { "HALF", HalfTok },
            { "HELP", HelpTok },
            { "INFO", InfoTok },
            { "IONS", IonTok  },
            { "IRIS", IRISTok },
            { "LOAD", LoadTok },
            { "MMDB", MMDBTok },
            { "MOL2", Mol2Tok },
            { "MONO", MonoTok },
            { "NONE", NoneTok },
            { "PICT", PICTTok },
            { "QUIT", QuitTok },
            { "SAVE", SaveTok },
            { "SHOW", ShowTok },
            { "SLAB", SlabTok },
            { "TRUE", TrueTok },
            { "TURN", TurnTok },
            { "TYPE", TypeTok },
            { "USER", UserTok },
            { "VRML", VRMLTok },
            { "WAIT", WaitTok },
            { "ZOOM", ZoomTok }, /* 64 */

            { "ALPHA", AlphaTok    },
            { "AMINO", AminoTok    },
            { "ANGLE", AngleTok    },
            { "ATOMS", AtomTok     },
            { "BASIC", BasicTok    },
            { "BLACK", BlackTok    },
            { "BONDS", BondTok     },
            { "CHAIN", ChainTok    },
            { "COLOR", ColourTok   },
            { "FALSE", FalseTok    },
            { "GREEN", GreenTok    },
            { "GROUP", GroupTok    },
            { "HBOND", HBondTok    },
            { "HELIX", HelixTok    },
            { "IDENT", IdentifyTok },
            { "LABEL", LabelTok    },
            { "LARGE", LargeTok    },
            { "MENUS", MenusTok    },
            { "MODEL", ModelTok    },
            { "MOPAC", MOPACTok    },
            { "MOUSE", MouseTok    },
            { "PAUSE", WaitTok     },
            { "POLAR", PolarTok    },
            { "PRINT", PrintTok    },
            { "RENUM", RenumTok    },
            { "RESET", ResetTok    },
            { "RESNO", ResNoTok    },
            { "SHEET", SheetTok    },
            { "SHELX", SHELXTok    },
            { "SMALL", SmallTok    },
            { "SOLID", SolidTok    },
            { "SYBYL", SybylTok    },
            { "TRACE", TraceTok    },
            { "TURNS", TurnTok     },
            { "WATER", WaterTok    },
            { "WHITE", WhiteTok    },
            { "WRITE", WriteTok    },  /* 101 */

            { "ACIDIC", AcidicTok },
            { "ANGLES", AngleTok  },
            { "ATOMNO", AtomNoTok },
            { "BIOSYM", BiosymTok },
            { "BONDED", BondedTok },
            { "BURIED", BuriedTok },
            { "CENTER", CentreTok },
            { "CENTRE", CentreTok },
            { "CHARGE", ChargeTok },
            { "CHARMM", CharmmTok },
            { "COLORS", ColourTok },
            { "COLOUR", ColourTok },
            { "CYCLIC", CyclicTok },
            { "DASHES", DashTok   },
            { "DEFINE", DefineTok },
            { "ELEMNO", ElemNoTok },
            { "HBONDS", HBondTok  },
            { "HETERO", HeteroTok },
            { "HOLLOW", HollowTok },
            { "INLINE", InLineTok },
            { "LABELS", LabelTok  },
            { "LIGAND", LigandTok },
            { "MEDIUM", MediumTok },
            { "MONOPS", MonoPSTok },
            { "NMRPDB", NMRPDBTok },
            { "NORMAL", NormalTok },
            { "ORANGE", OrangeTok },
            { "POVRAY", POVRayTok },
            { "PURINE", PurineTok },
            { "PURPLE", PurpleTok },
            { "QUANTA", QuantaTok },
            { "RADIUS", RadiusTok },
            { "RASMOL", RasMolTok },
            { "RASWIN", RasMolTok },
            { "REJECT", RejectTok },
            { "RESIZE", ResizeTok },
            { "RIBBON", RibbonTok },
            { "ROTATE", RotateTok },
            { "SCRIPT", ScriptTok },
            { "SELECT", SelectTok },
            { "SHADOW", ShadowTok },
            { "SHEETS", SheetTok  },
            { "SOURCE", SourceTok },
            { "SSBOND", SSBondTok },
            { "STEREO", StereoTok },
            { "SUNRLE", SUNRLETok },
            { "VECTPS", VectPSTok },
            { "VIOLET", VioletTok },
            { "WATERS", WaterTok  },
            { "WITHIN", WithinTok },
            { "YELLOW", YellowTok },  /* 152 */

            { "ACYCLIC", AcyclicTok },
            { "ALCHEMY", AlchemyTok },
            { "AMBIENT", AmbientTok },
            { "CARTOON", CartoonTok },
            { "CHARGED", ChargedTok },
            { "CHARGES", ChargeTok  },
            { "COLOURS", ColourTok  },
            { "CONNECT", ConnectTok },
            { "CYSTINE", CystineTok },
            { "DISPLAY", DisplayTok },
            { "HELICES", HelixTok   },
            { "INSIGHT", InsightTok },
            { "LIGANDS", LigandTok  },
            { "MAGENTA", MagentaTok },
            { "MONITOR", MonitorTok },
            { "NEUTRAL", NeutralTok },
            { "NUCLEIC", NucleicTok },
            { "PICKING", PickingTok },
            { "PROTEIN", ProteinTok },
            { "PURINES", PurineTok  },
            { "REFRESH", RefreshTok },
            { "RESIDUE", ResidueTok },
            { "RIBBON1", Ribbon1Tok },
            { "RIBBON2", Ribbon2Tok },
            { "RIBBONS", RibbonTok  },
            { "SECTION", SectionTok },
            { "SHADOWS", ShadowTok  },
            { "SHAPELY", ShapelyTok },
            { "SOLVENT", SolventTok },
            { "SSBONDS", SSBondTok  },
            { "STRANDS", StrandsTok },
            { "SURFACE", SurfaceTok },  
            { "TORSION", TorsionTok }, /* 185 */

            { "AROMATIC", AromaticTok },
            { "BACKBONE", BackboneTok },
            { "BACKFADE", BackFadeTok },
            { "BONDMODE", BondModeTok },
            { "BOUNDBOX", BoundBoxTok },
            { "CARTOONS", CartoonTok  },
            { "DEPTHCUE", DepthCueTok },
            { "DISTANCE", DistanceTok },
            { "FONTSIZE", FontSizeTok },
            { "GAUSSIAN", GaussianTok },
            { "HARDWARE", HardwareTok },
            { "HYDROGEN", HydrogenTok },
            { "IDENTIFY", IdentifyTok },
            { "KINEMAGE", KinemageTok },
            { "MONITORS", MonitorTok  },
            { "NEGATIVE", AcidicTok   },
            { "POSITIVE", BasicTok    },
            { "RENUMBER", RenumTok    },
            { "RESTRICT", RestrictTok },
            { "RIBBONS1", Ribbon1Tok  },
            { "RIBBONS2", Ribbon2Tok  },
            { "SELECTED", SelectedTok },
            { "SEQUENCE", SequenceTok },
            { "SLABMODE", SlabModeTok },
            { "SOLVENTS", SolventTok  },
            { "SPECULAR", SpecularTok }, 
            { "SYMMETRY", SymmetryTok },
            { "TORSIONS", TorsionTok  },
            { "UNITCELL", UnitCellTok },  /* 214 */

            { "ALIPHATIC", AliphaticTok },
            { "CLIPBOARD", ClipboardTok },
            { "DISTANCES", DistanceTok  },
            { "GREENBLUE", GreenblueTok },
            { "HOURGLASS", HourGlassTok },
            { "MAINCHAIN", MainChainTok },
            { "MOLSCRIPT", MolScriptTok },
            { "MOUSEMODE", MouseTok     },
            { "POTENTIAL", PotentialTok },
            { "REDORANGE", RedorangeTok },
            { "SIDECHAIN", SidechainTok },
            { "SPACEFILL", SpacefillTok },
            { "SPECPOWER", SpecPowerTok },
            { "STRUCTURE", StructureTok },
            { "TRANSLATE", TranslateTok },
            { "WIREFRAME", WireframeTok },  /* 230 */

            { "BACKGROUND", BackgroundTok },
            { "MACROMODEL", MacroModelTok },
            { "MONOCHROME", MonoTok       },
            { "PYRIMIDINE", PyrimidineTok },  /* 234 */

            { "BOUNDINGBOX", BoundBoxTok    },
            { "HYDROPHOBIC", HydrophobicTok },
            { "INFORMATION", InfoTok        },
            { "PYRIMIDINES", PyrimidineTok, },
            { "TEMPERATURE", TemperatureTok },
            { "TRANSPARENT", TransparentTok }  /* 241 */
                };
