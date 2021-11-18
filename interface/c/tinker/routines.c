#include "macro.hh"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef TINKER_GFORTRAN
void _gfortran_set_args(int, char**);
void tinkerFortranRuntimeBegin(int argc, char** argv) { _gfortran_set_args(argc, argv); }
void tinkerFortranRuntimeEnd() {}
#endif

#ifdef TINKER_IFORT
void for_rtl_init_(int*, char**);
int for_rtl_finish_();
extern int for__l_argc;
extern char** for__a_argv;
void tinkerFortranRuntimeBegin(int argc, char** argv) { for__l_argc = argc; for__a_argv = argv; for_rtl_init_(&argc, argv); }
void tinkerFortranRuntimeEnd() { for_rtl_finish_(); }
#endif

#ifdef __cplusplus
}
#endif
