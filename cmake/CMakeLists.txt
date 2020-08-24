cmake_minimum_required(VERSION 3.9)

project(tinker VERSION 8.8.1 LANGUAGES Fortran)

include(GNUInstallDirs)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
    set(
        CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
    )
endif()

# Determine main Tinker directory
find_path(
    TINKER_DIR
    "source/dynamic.f"
    HINTS ../
    DOC "Main Tinker directory"
    REQUIRED
)

find_package(OpenMP REQUIRED)

set(BUILD_FFTW OFF CACHE BOOL "Build included fftw package")

if(BUILD_FFTW)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")
    set(ENABLE_OPENMP ON CACHE BOOL "Use OpenMP for multithreading")
    set(ENABLE_THREADS ON CACHE BOOL "Use pthread for multithreading")
    add_subdirectory(${TINKER_DIR}/fftw ${CMAKE_BINARY_DIR}/fftw)
    set(FFT_LIB fftw3)
else()
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        # Intel compilers have an FFT library built-in
        set(FFT_LIB "")
    else()
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
        set(FFT_LIB PkgConfig::FFTW3)
    endif()
endif()

if(APPLE)
    set(ARCH_FLAGS "-mtune=native")
else()
    set(ARCH_FLAGS "-march=native")
endif()

if(WIN32)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS_RELEASE "/O3 /QxHost /Qip- /Qprec-div- /w")
        set(CMAKE_Fortran_FLAGS_DEBUG "/debug:full")
    endif()
elseif(UNIX) # includes cygwin
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast ${ARCH_FLAGS}")
        set(CMAKE_Fortran_FLAGS_DEBUG "-Og -fbacktrace -fcheck=bounds -Wall")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost -no-ipo -no-prec-div")
        set(CMAKE_Fortran_FLAGS_DEBUG "-g -warn all -check all")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        set(CMAKE_Fortran_FLAGS_RELEASE "-fast -O3")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Mdclchk -Mbounds -Mchkptr -Mchkstk")
    endif()
else()
    message(WARNING "OS-specific flags were not set properly!")
endif()

set(_FILES
    action active align alterchg analysis analyz angang angbnd angles angpot
    angtor argue ascii atmlst atomid atoms attach baoab basefile bath beeman
    bicubic bitor bitors bndpot bndstr bonds born bound bounds boxes bussi
    calendar cell center cflux charge chgpen chgpot chgtrn chkpole chkring
    chkxyz cholesky chrono chunks clock cluster column command connect connolly
    control couple cspline ctrpot cutoffs damping deflate delete deriv dcflux
    diagq diffeq dipole disgeo disp dma domega dsppot eangang eangang1 eangang2
    eangang3 eangle eangle1 eangle2 eangle3 eangtor eangtor1 eangtor2 eangtor3
    ebond ebond1 ebond2 ebond3 ebuck ebuck1 ebuck2 ebuck3 echarge echarge1
    echarge2 echarge3 echgdpl echgdpl1 echgdpl2 echgdpl3 echgtrn echgtrn1
    echgtrn2 echgtrn3 edipole edipole1 edipole2 edipole3 edisp edisp1 edisp2
    edisp3 egauss egauss1 egauss2 egauss3 egeom egeom1 egeom2 egeom3 ehal ehal1
    ehal2 ehal3 eimprop eimprop1 eimprop2 eimprop3 eimptor eimptor1  eimptor2
    eimptor3 elj elj1 elj2 elj3 embed emetal emetal1 emetal2 emetal3 emm3hb
    emm3hb1 emm3hb2 emm3hb3 empole empole1 empole2 empole3 energi energy eopbend
    eopbend1 eopbend2 eopbend3 eopdist eopdist1 eopdist2 eopdist3 epitors
    epitors1 epitors2 epitors3 epolar epolar1 epolar2 epolar3 erepel erepel1
    erepel2 erepel3 erf erxnfld erxnfld1 erxnfld2 erxnfld3 esolv esolv1 esolv2
    esolv3 estrbnd estrbnd1 estrbnd2 estrbnd3 estrtor estrtor1 estrtor2 estrtor3
    etors etors1 etors2 etors3 etortor etortor1 etortor2 etortor3 eurey eurey1
    eurey2 eurey3 evcorr ewald extra extra1 extra2 extra3 faces fatal fft fft3d
    fftpack field fields files final flatten fracs freeunit freeze geometry
    getarc getint getkey getmol getmol2 getnumb getpdb getprm getref getstring
    gettext getword getxyz ghmcstep gkstuf gradient gradrgd gradrot group groups
    grpline gyrate hescut hessian hessn hessrgd hessrot hpmf hybrid ielscf image
    impose improp imptor induce inertia inform initatom initial initprm initres
    initrot insert inter invbeta invert iounit jacobi kanang kangang kangle
    kangs kangtor kantor katom katoms kbond kbonds kcflux kcharge kchgflx
    kchgtrn kchrge kcpen kctrn kdipol kdipole kdisp kdsp kewald kextra keys
    kgeom khbond kimprop kimptor kinetic kiprop kitors kmetal kmpole kmulti
    kopbend kopbnd kopdist kopdst korbit korbs kpitor kpitors kpolar kpolr
    krepel krepl ksolut ksolv kstbnd kstrbnd kstrtor ksttor ktors ktorsn ktortor
    ktrtor kurey kurybr kvdw kvdwpr kvdws lattice lbfgs light lights limits
    linmin makeint makeref makexyz math maxwell mdinit mdrest mdsave mdstat
    mdstuf mechanic merck merge minima molcul moldyn molecule moment moments
    mplpot mpole mrecip mutant mutate nblist neigh nextarg nexttext nonpol nose
    nspline nucleo number numeral numgrad ocvm omega opbend opdist openend
    openmp optinit optsave orbital orbits orient orthog output overlap params
    paths pbstuf pdb phipsi picalc piorbs pistuf pitors pme pmestuf pmpb polar
    polgrp polopt polpcg polpot poltcg polymer potent potfit pressure prmkey
    promo prtdyn prterr prtint prtmol2 prtpdb prtprm prtseq prtxyz ptable qmstuf
    qrfact quatfit random rattle readdyn readgau readgdma readint readmol
    readmol2 readpdb readprm readseq readxyz refer repel replica reppot resdue
    respa restrn rgddyn rgdstep rigid ring rings rmsfit rotbnd rotlist rotpole
    rxnfld rxnpot scales sdstep search sequen server shakeup shunt sigmoid
    simplex sizes sktstuf socket solpot solute sort square stodyn strbnd strtor
    suffix surface surfatom switch syntrn tarray tcgstuf temper titles tncg
    torphase torpot torque tors torsions tortor tree trimtext unitcell units
    uprior urey urypot usage valfit vdw vdwpot verlet version vibs virial volume
    warp xtals xyzatm zatom zclose zcoord
)
set(_LIB "")
foreach(_FILE IN LISTS _FILES)
    list(APPEND _LIB ${TINKER_DIR}/source/${_FILE}.f)
endforeach()
add_library(tinker STATIC ${_LIB})
install(
    TARGETS tinker
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set(_BINS
    alchemy analyze anneal archive bar correlate crystal diffuse
    distgeom document dynamic gda intedit intxyz minimize minirot
    minrigid mol2xyz molxyz monte newton newtrot nucleic optimize
    optirot optrigid path pdbxyz polarize poledit potential
    prmedit protein pss pssrigid pssrot radial saddle scan sniffer
    spacefill spectrum superpose testgrad testhess testpair
    testpol testrot testvir timer timerot torsfit valence vibbig
    vibrate vibrot xtalfit xtalmin xyzedit xyzint xyzmol2 xyzpdb
)
foreach(_BIN IN LISTS _BINS)
    add_executable(${_BIN} ${TINKER_DIR}/source/${_BIN}.f)
    target_link_libraries(${_BIN} tinker ${FFT_LIB} OpenMP::OpenMP_Fortran)
    install(
        TARGETS ${_BIN}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endforeach()

