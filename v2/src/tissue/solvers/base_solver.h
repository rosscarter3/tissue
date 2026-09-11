//
// Solver base: owns the working state matrices, the print schedule, and the
// shared per-step bookkeeping. Concrete solvers implement simulate().
// Solver parameter files are backwards compatible with the legacy format
// (first token = solver id, then the id-specific fields; '#' comments).
//
#ifndef TISSUE2_SOLVERS_BASE_SOLVER_H
#define TISSUE2_SOLVERS_BASE_SOLVER_H

#include <iostream>
#include <memory>
#include <string>

#include "tissue/core/matrix.h"
#include "tissue/core/tissue.h"

namespace tissue {

class BaseSolver {
public:
  explicit BaseSolver(Tissue *T) : T_(T) {}
  virtual ~BaseSolver() = default;

  // Reads the solver file (comment-filtered), dispatches on the solver id
  // (RK5Adaptive, RK4, Euler, HeunIto) and returns the configured solver.
  static std::unique_ptr<BaseSolver> getSolver(Tissue *T,
                                               const std::string &file);

  // Copies the tissue state mirror into the solver's working matrices.
  void getInit();
  // Copies the working matrices back into the tissue mirror; numCellVariable
  // = size_t(-1) copies whole cell rows (adapting their length).
  void setTissueVariables(size_t numCellVariable = static_cast<size_t>(-1));

  virtual void simulate() = 0;

  // Per-print-point output, dispatched on printFlag (see solver files).
  void print(std::ostream &os = std::cout);
  // Final state in init format (used by -init_output).
  void printInit(std::ostream &os) const;

  double startTime() const { return startTime_; }
  double endTime() const { return endTime_; }

protected:
  // Shared end-of-step sequence: direction update, reaction updates,
  // compartment changes, connectivity check.
  void postStep(double h);
  // Print-schedule setup shared by all solvers (legacy semantics: numPrint<=0
  // none, 1 final only, 2 first+last, else evenly spaced accumulating).
  void initPrintSchedule(double tiny);

  Tissue *T_;
  Matrix cellData_, wallData_, vertexData_;
  Matrix cellDerivs_, wallDerivs_, vertexDerivs_;
  double t_ = 0.0;
  double startTime_ = 0.0, endTime_ = 0.0;
  int printFlag_ = 0;
  int numPrint_ = 0;
  unsigned int numOk_ = 0, numBad_ = 0;

  // print schedule state
  bool doPrint_ = true;
  double printTime_ = 0.0;
  double printDeltaTime_ = 0.0;

private:
  // print() bookkeeping (legacy function statics)
  int tCount_ = 0;
  int nOld_ = 0;
  unsigned int okOld_ = 0, badOld_ = 0;
  double tOld_ = 0.0;
};

} // namespace tissue

#endif
