/*
 * Derived from PaPILO (https://github.com/scipopt/papilo), Apache-2.0.
 */

#include "papilo_presolve.h"
#include "papilo/core/Presolve.hpp"
#include "papilo/core/ProblemBuilder.hpp"
#include "papilo/core/postsolve/Postsolve.hpp"
#include "papilo/core/postsolve/ReductionType.hpp"
#include <exception>
#include <string>
#include <vector>

using namespace papilo;

struct Papilo_Presolve {
   Presolve<double> presolve;
   Papilo_Presolve() { presolve.addDefaultPresolvers(); }
};

struct Papilo_PresolveResult {
   PresolveStatus status;
   Problem<double> problem;
   PostsolveStorage<double> postsolve;
   int origCols, origRows;
};

static thread_local std::string papilo_error;

static PAPILO_STATUS papilo_status_from_cpp(PresolveStatus s) {
   switch (s) {
      case PresolveStatus::kUnchanged: return PAPILO_STATUS_UNCHANGED;
      case PresolveStatus::kReduced: return PAPILO_STATUS_REDUCED;
      case PresolveStatus::kInfeasible: return PAPILO_STATUS_INFEASIBLE;
      case PresolveStatus::kUnbounded: return PAPILO_STATUS_UNBOUNDED;
      case PresolveStatus::kUnbndOrInfeas: return PAPILO_STATUS_UNBND_OR_INFEAS;
   }
   return PAPILO_STATUS_UNBND_OR_INFEAS;
}

static bool papilo_ptr_ok(const void* p, int n) { return n == 0 || p != nullptr; }
static void papilo_clear_error() { papilo_error.clear(); }
static void papilo_set_error(const char* msg) { papilo_error = msg ? msg : "papilo error"; }
static void papilo_set_error(const std::exception& e) { papilo_error = e.what(); }

template <typename Ret, typename Fn>
static Ret papilo_try(const char* fallback_error, Ret fail_value, const Fn& fn) {
   try { return fn(); }
   catch (const std::exception& e) { papilo_set_error(e); return fail_value; }
   catch (...) { papilo_set_error(fallback_error); return fail_value; }
}

static PAPILO_PRESOLVE_RESULT* papilo_apply_impl(
    PAPILO_PRESOLVE* p, int ncols, int nrows,
    const double* obj, double obj_offset,
    const double* col_lb, const double* col_ub,
    const unsigned char* col_lb_inf, const unsigned char* col_ub_inf,
    const unsigned char* col_integral,
    const double* row_lhs, const double* row_rhs,
    const unsigned char* row_lhs_inf, const unsigned char* row_rhs_inf,
    int nnz, const int* row_start, const int* col_idx, const double* vals)
{
   ProblemBuilder<double> b;
   b.reserve(nnz, nrows, ncols);
   b.setNumCols(ncols);
   b.setNumRows(nrows);
   b.setObjOffset(obj_offset);

   for (int j = 0; j < ncols; ++j) {
      b.setObj(j, obj[j]);
      b.setColLb(j, col_lb[j]);
      b.setColUb(j, col_ub[j]);
      b.setColLbInf(j, col_lb_inf[j]);
      b.setColUbInf(j, col_ub_inf[j]);
      b.setColIntegral(j, col_integral[j]);
   }

   for (int i = 0; i < nrows; ++i) {
      if (row_start[i] > row_start[i + 1] || row_start[i] < 0 || row_start[i + 1] > nnz) {
         papilo_set_error("invalid row_start");
         return nullptr;
      }
      b.setRowLhs(i, row_lhs[i]);
      b.setRowRhs(i, row_rhs[i]);
      b.setRowLhsInf(i, row_lhs_inf[i]);
      b.setRowRhsInf(i, row_rhs_inf[i]);
      for (int k = row_start[i]; k < row_start[i + 1]; ++k) {
         if (col_idx[k] < 0 || col_idx[k] >= ncols) {
            papilo_set_error("invalid col_idx");
            return nullptr;
         }
         if (vals[k] != 0.0) b.addEntry(i, col_idx[k], vals[k]);
      }
   }

   Problem<double> problem = b.build();
   int origCols = problem.getNCols(), origRows = problem.getNRows();
   PresolveResult<double> res = p->presolve.apply(problem, false);

   auto* r = new Papilo_PresolveResult();
   r->status = res.status;
   r->problem = std::move(problem);
   r->postsolve = std::move(res.postsolve);
   r->origCols = origCols;
   r->origRows = origRows;
   return r;
}

static int papilo_postsolve_impl(const PAPILO_PRESOLVE_RESULT* r, const double* reduced, double* original) {
   int ncols = static_cast<int>(r->postsolve.origcol_mapping.size());
   Solution<double> redSol(SolutionType::kPrimal);
   redSol.primal.assign(reduced, reduced + ncols);

   Solution<double> origSol(SolutionType::kPrimal);
   Num<double> num;
   Message msg;
   msg.setVerbosityLevel(VerbosityLevel::kQuiet);

   Postsolve<double> ps(msg, num);
   // we pass isoptimal = false since CoolPDLP has looser tolerances than papilo.
   ps.undo(redSol, origSol, r->postsolve, false);
   for (size_t j = 0; j < origSol.primal.size(); ++j) original[j] = origSol.primal[j];
   return 0;
}

extern "C" {

PAPILO_PRESOLVE* papilo_create() {
   papilo_clear_error();
   return papilo_try("papilo_create failed", static_cast<PAPILO_PRESOLVE*>(nullptr),
                     []() { return new Papilo_Presolve(); });
}
const char* papilo_last_error() { return papilo_error.c_str(); }
void papilo_free(PAPILO_PRESOLVE* p) { delete p; }
void papilo_set_threads(PAPILO_PRESOLVE* p, int n) {
   papilo_clear_error();
   if (!p) { papilo_set_error("papilo_set_threads: null presolve"); return; }
   p->presolve.getPresolveOptions().threads = n;
}
// TODO: time limit seems to be not working
void papilo_set_verbosity(PAPILO_PRESOLVE* p, int level) {
   papilo_clear_error();
   if (!p) { papilo_set_error("papilo_set_verbosity: null presolve"); return; }
   VerbosityLevel vl = level == 0 ? VerbosityLevel::kQuiet :
                       level == 1 ? VerbosityLevel::kError :
                       level == 2 ? VerbosityLevel::kWarning :
                       level == 3 ? VerbosityLevel::kInfo : VerbosityLevel::kDetailed;
   p->presolve.setVerbosityLevel(vl);
}

PAPILO_PRESOLVE_RESULT* papilo_apply(
    PAPILO_PRESOLVE* p, int ncols, int nrows,
    const double* obj, double obj_offset,
    const double* col_lb, const double* col_ub,
    const unsigned char* col_lb_inf, const unsigned char* col_ub_inf,
    const unsigned char* col_integral,
    const double* row_lhs, const double* row_rhs,
    const unsigned char* row_lhs_inf, const unsigned char* row_rhs_inf,
    int nnz, const int* row_start, const int* col_idx, const double* vals)
{
   papilo_clear_error();
   if (!p || ncols < 0 || nrows < 0 || nnz < 0) {
      papilo_set_error("papilo_apply: invalid dimensions"); return nullptr;
   }
   if (!papilo_ptr_ok(obj, ncols) || !papilo_ptr_ok(col_lb, ncols) || !papilo_ptr_ok(col_ub, ncols) ||
       !papilo_ptr_ok(col_lb_inf, ncols) || !papilo_ptr_ok(col_ub_inf, ncols) || !papilo_ptr_ok(col_integral, ncols) ||
       !papilo_ptr_ok(row_lhs, nrows) || !papilo_ptr_ok(row_rhs, nrows) || !papilo_ptr_ok(row_lhs_inf, nrows) || !papilo_ptr_ok(row_rhs_inf, nrows) ||
       !papilo_ptr_ok(row_start, nrows + 1) || !papilo_ptr_ok(col_idx, nnz) || !papilo_ptr_ok(vals, nnz)) {
      papilo_set_error("papilo_apply: null input pointer"); return nullptr;
   }
   if (row_start[0] != 0 || row_start[nrows] != nnz) {
      papilo_set_error("papilo_apply: invalid row_start bounds"); return nullptr;
   }
   return papilo_try("papilo_apply failed", static_cast<PAPILO_PRESOLVE_RESULT*>(nullptr), [&]() {
      return papilo_apply_impl(p, ncols, nrows, obj, obj_offset, col_lb, col_ub, col_lb_inf, col_ub_inf, col_integral, row_lhs, row_rhs, row_lhs_inf, row_rhs_inf, nnz, row_start, col_idx, vals);
   });
}

PAPILO_STATUS papilo_status(const PAPILO_PRESOLVE_RESULT* r) { return r ? papilo_status_from_cpp(r->status) : PAPILO_STATUS_UNBND_OR_INFEAS; }
int papilo_num_cols(const PAPILO_PRESOLVE_RESULT* r) { return r ? r->problem.getNCols() : -1; }
int papilo_num_rows(const PAPILO_PRESOLVE_RESULT* r) { return r ? r->problem.getNRows() : -1; }
int papilo_nnz(const PAPILO_PRESOLVE_RESULT* r) { return r ? r->problem.getConstraintMatrix().getNnz() : -1; }
int papilo_orig_cols(const PAPILO_PRESOLVE_RESULT* r) { return r ? r->origCols : -1; }
int papilo_orig_rows(const PAPILO_PRESOLVE_RESULT* r) { return r ? r->origRows : -1; }
void papilo_get_obj(const PAPILO_PRESOLVE_RESULT* r, double* coeffs, double* offset) {
   if (!r) return;
   const auto& obj = r->problem.getObjective();
   for (int j = 0; j < r->problem.getNCols(); ++j) coeffs[j] = obj.coefficients[j];
   *offset = obj.offset;
}
void papilo_get_col_bounds(const PAPILO_PRESOLVE_RESULT* r, double* lb, double* ub,
                           unsigned char* lb_inf, unsigned char* ub_inf, unsigned char* integral) {
   if (!r) return;
   const auto& d = r->problem.getVariableDomains();
   const auto& f = r->problem.getColFlags();
   for (int j = 0; j < r->problem.getNCols(); ++j) {
      lb[j] = d.lower_bounds[j];
      ub[j] = d.upper_bounds[j];
      lb_inf[j] = f[j].test(ColFlag::kLbInf);
      ub_inf[j] = f[j].test(ColFlag::kUbInf);
      integral[j] = f[j].test(ColFlag::kIntegral);
   }
}
void papilo_get_row_bounds(const PAPILO_PRESOLVE_RESULT* r, double* lhs, double* rhs,
                           unsigned char* lhs_inf, unsigned char* rhs_inf) {
   if (!r) return;
   const auto& m = r->problem.getConstraintMatrix();
   const auto& f = m.getRowFlags();
   for (int i = 0; i < r->problem.getNRows(); ++i) {
      lhs[i] = m.getLeftHandSides()[i];
      rhs[i] = m.getRightHandSides()[i];
      lhs_inf[i] = f[i].test(RowFlag::kLhsInf);
      rhs_inf[i] = f[i].test(RowFlag::kRhsInf);
   }
}
void papilo_get_matrix(const PAPILO_PRESOLVE_RESULT* r, int* row_start, int* col_idx, double* vals) {
   if (!r) return;
   const auto& m = r->problem.getConstraintMatrix();
   int k = 0;
   for (int i = 0; i < r->problem.getNRows(); ++i) {
      row_start[i] = k;
      auto row = m.getRowCoefficients(i);
      for (int j = 0; j < row.getLength(); ++j) {
         col_idx[k] = row.getIndices()[j];
         vals[k] = row.getValues()[j];
         ++k;
      }
   }
   row_start[r->problem.getNRows()] = k;
}
void papilo_get_col_map(const PAPILO_PRESOLVE_RESULT* r, int* map) {
   if (!r) return;
   const auto& m = r->postsolve.origcol_mapping;
   for (size_t j = 0; j < m.size(); ++j) map[j] = m[j];
}
void papilo_get_row_map(const PAPILO_PRESOLVE_RESULT* r, int* map) {
   if (!r) return;
   const auto& m = r->postsolve.origrow_mapping;
   for (size_t i = 0; i < m.size(); ++i) map[i] = m[i];
}
int papilo_postsolve(const PAPILO_PRESOLVE_RESULT* r, const double* reduced, double* original) {
   papilo_clear_error();
   if (!r || !reduced || !original) { papilo_set_error("papilo_postsolve: null input"); return 1; }
   return papilo_try("papilo_postsolve failed", 1, [&]() { return papilo_postsolve_impl(r, reduced, original); });
}

void papilo_map_primal(const PAPILO_PRESOLVE_RESULT* r, const double* x_orig, double* x_red) {
   if (!r || !x_orig || !x_red) return;

   const auto& cm     = r->postsolve.origcol_mapping;   // reduced → original
   const int   n_red  = static_cast<int>(cm.size());
   const int   n_orig = r->origCols;

   // Step 1: build inverse map: original → reduced (-1 if not present)
   std::vector<int> orig_to_red(n_orig, -1);
   for (int j = 0; j < n_red; ++j) orig_to_red[cm[j]] = j;

   // Step 2: identity gather
   for (int j = 0; j < n_red; ++j) x_red[j] = x_orig[cm[j]];

   // Step 3: parallel-column fixup.
   // For kParallelCol: y = col2 + col2scale * col1  (col2 survives, col1 eliminated)
   //   indices[s+0] = orig_col1 (eliminated)
   //   indices[s+2] = orig_col2 (surviving)
   //   values[s+4]  = col2scale
   // Step 2 already placed x_orig[orig_col2] in x_red[j_red].
   // We add col2scale * x_orig[orig_col1] to complete the inverse.
   const auto& types   = r->postsolve.types;
   const auto& indices = r->postsolve.indices;
   const auto& values  = r->postsolve.values;
   const auto& start   = r->postsolve.start;

   for (int t = 0; t < static_cast<int>(types.size()); ++t) {
      if (types[t] != ReductionType::kParallelCol) continue;
      const int s        = start[t];
      const int orig_col1 = indices[s];      // eliminated
      const int orig_col2 = indices[s + 2];  // surviving
      const double col2scale = values[s + 4];

      const int j_red = orig_to_red[orig_col2];
      if (j_red >= 0)
         x_red[j_red] += col2scale * x_orig[orig_col1];
   }
}

void papilo_result_free(PAPILO_PRESOLVE_RESULT* r) { delete r; }

} // extern "C"
