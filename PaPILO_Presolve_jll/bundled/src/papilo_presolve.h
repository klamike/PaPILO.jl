/*
 * Derived from PaPILO (https://github.com/scipopt/papilo), Apache-2.0.
 */

#ifndef PAPILO_PRESOLVE_H
#define PAPILO_PRESOLVE_H

#include "papilo_presolve_export.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Papilo_Presolve PAPILO_PRESOLVE;
typedef struct Papilo_PresolveResult PAPILO_PRESOLVE_RESULT;

typedef enum {
    PAPILO_STATUS_UNCHANGED = 0,
    PAPILO_STATUS_REDUCED = 1,
    PAPILO_STATUS_INFEASIBLE = 2,
    PAPILO_STATUS_UNBOUNDED = 3,
    PAPILO_STATUS_UNBND_OR_INFEAS = 4
} PAPILO_STATUS;

PAPILO_PRESOLVE_EXPORT PAPILO_PRESOLVE* papilo_create(void);
PAPILO_PRESOLVE_EXPORT const char* papilo_last_error(void);
PAPILO_PRESOLVE_EXPORT void papilo_free(PAPILO_PRESOLVE* p);
PAPILO_PRESOLVE_EXPORT void papilo_set_threads(PAPILO_PRESOLVE* p, int n);
PAPILO_PRESOLVE_EXPORT void papilo_set_verbosity(PAPILO_PRESOLVE* p, int level);

PAPILO_PRESOLVE_EXPORT PAPILO_PRESOLVE_RESULT* papilo_apply(
    PAPILO_PRESOLVE* p, int ncols, int nrows,
    const double* obj, double obj_offset,
    const double* col_lb, const double* col_ub,
    const unsigned char* col_lb_inf, const unsigned char* col_ub_inf,
    const unsigned char* col_integral,
    const double* row_lhs, const double* row_rhs,
    const unsigned char* row_lhs_inf, const unsigned char* row_rhs_inf,
    int nnz, const int* row_start, const int* col_idx, const double* vals);

PAPILO_PRESOLVE_EXPORT PAPILO_STATUS papilo_status(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT int papilo_num_cols(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT int papilo_num_rows(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT int papilo_nnz(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT int papilo_orig_cols(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT int papilo_orig_rows(const PAPILO_PRESOLVE_RESULT* r);
PAPILO_PRESOLVE_EXPORT void papilo_get_obj(const PAPILO_PRESOLVE_RESULT* r, double* coeffs, double* offset);
PAPILO_PRESOLVE_EXPORT void papilo_get_col_bounds(const PAPILO_PRESOLVE_RESULT* r, double* lb, double* ub, unsigned char* lb_inf, unsigned char* ub_inf, unsigned char* integral);
PAPILO_PRESOLVE_EXPORT void papilo_get_row_bounds(const PAPILO_PRESOLVE_RESULT* r, double* lhs, double* rhs, unsigned char* lhs_inf, unsigned char* rhs_inf);
PAPILO_PRESOLVE_EXPORT void papilo_get_matrix(const PAPILO_PRESOLVE_RESULT* r, int* row_start, int* col_idx, double* vals);
PAPILO_PRESOLVE_EXPORT void papilo_get_col_map(const PAPILO_PRESOLVE_RESULT* r, int* map);
PAPILO_PRESOLVE_EXPORT void papilo_get_row_map(const PAPILO_PRESOLVE_RESULT* r, int* map);
PAPILO_PRESOLVE_EXPORT int papilo_postsolve(const PAPILO_PRESOLVE_RESULT* r, const double* reduced, double* original);
PAPILO_PRESOLVE_EXPORT void papilo_map_primal(const PAPILO_PRESOLVE_RESULT* r, const double* x_orig, double* x_red);
PAPILO_PRESOLVE_EXPORT void papilo_result_free(PAPILO_PRESOLVE_RESULT* r);

#ifdef __cplusplus
}
#endif

#endif
