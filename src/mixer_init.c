#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void VEstep_(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void init_ermg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void main_ermg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void main_ermgo(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void vertex_coord_C(void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"VEstep_",        (DL_FUNC) &VEstep_,        13},
    {"init_ermg",      (DL_FUNC) &init_ermg,      12},
    {"main_ermg",      (DL_FUNC) &main_ermg,      20},
    {"main_ermgo",     (DL_FUNC) &main_ermgo,     14},
    {"vertex_coord_C", (DL_FUNC) &vertex_coord_C,  9},
    {NULL, NULL, 0}
};

void R_init_mixer(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
