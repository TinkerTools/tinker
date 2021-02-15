enable_language(CUDA)
enable_language(CXX)

# Find OpenMM Libraries

set(OPENMM_DIR /usr/local/openmm CACHE PATH "Directory of OpenMM libraries.")

if(DEFINED ENV{OPENMM_PLUGIN_DIR})
    get_filename_component(OPENMM_DIR "$ENV{OPENMM_PLUGIN_DIR}/../../" ABSOLUTE)
else()
    message(WARNING "Ensure that the environment variable OPENMM_PLUGIN_DIR is set before running Tinker-OpenMM.")
endif()

find_library(OPENMM_LIBRARY NAMES libOpenMM.so OpenMM HINTS ${OPENMM_DIR}/lib REQUIRED)
find_library(OPENMM_AMOEBA_LIBRARY NAMES libOpenMMAmoeba.so OpenMMAmoeba HINTS ${OPENMM_DIR}/lib REQUIRED)
include_directories(${OPENMM_DIR}/include)

# Build Tinker's OpenMM library

set(OPENMM_FILES
    gpu_cards.cu
    kopenmm.f
    ommdata.f
    ommstuf.cpp
    openmm.f
)

set(_LIB "")
foreach(_FILE IN LISTS OPENMM_FILES)
    list(APPEND _LIB ${TINKER_DIR}/openmm/${_FILE})
endforeach()

add_library(tinker_omm STATIC ${_LIB})
target_link_libraries(
    tinker_omm
    tinker
    OpenMP::OpenMP_Fortran
    ${FFTW_LIB} ${FFTW_THREADED_LIB}
    ${CUDA_LIBRARIES}
    ${OPENMM_LIBRARY} ${OPENMM_AMOEBA_LIBRARY}
    cuda
    nvidia-ml
)

install(
    TARGETS tinker_omm
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Build Tinker's OpenMM executables

set(OPENMM_BINS
    analyze_omm
    bar_omm
    dynamic_omm
)

foreach(_BIN IN LISTS OPENMM_BINS)
    add_executable(${_BIN} ${TINKER_DIR}/openmm/${_BIN}.f)
    set_property(TARGET ${_BIN} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(
        ${_BIN}
        tinker tinker_omm
        OpenMP::OpenMP_Fortran
        ${FFTW_LIB} ${FFTW_THREADED_LIB}
        ${CUDA_LIBRARIES}
        ${OPENMM_LIBRARY} ${OPENMM_AMOEBA_LIBRARY}
        cuda
        nvidia-ml
    )
    install(
        TARGETS ${_BIN}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endforeach()
