//
// Concrete solvers. Parameter-file layouts are identical to the legacy ones
// (whitespace tokens after the solver id):
//   RK5Adaptive: startTime endTime printFlag numPrint maxStep errorTolerance
//   RK4/Euler:   startTime endTime printFlag numPrint step
//   HeunIto:     startTime endTime printFlag numPrint step volume volumeFlag
//
// Numerical notes vs legacy (deliberate fixes, see v2/README.md):
//  - RK5 Cash-Karp stage 3 writes the wall derivatives into ak3 (legacy wrote
//    them into ak2, corrupting wall-variable integration).
//  - RK4 accumulates k3 into the vertex weighted sum (legacy dropped it).
//  - The redundant duplicate derivative evaluation per step in Euler/RK4 is
//    removed (identical results for side-effect-free reactions).
//
#ifndef TISSUE2_SOLVERS_SOLVERS_H
#define TISSUE2_SOLVERS_SOLVERS_H

#include <istream>

#include <vector>

#include "tissue/solvers/base_solver.h"

namespace tissue {

class Euler : public BaseSolver {
public:
  Euler(Tissue *T, std::istream &in);
  void simulate() override;

private:
  double h_ = 0.0;
};

// Growth stepping with mechanical equilibrium solved (FIRE) rather than
// integrated; see quasi_static.cpp for why.
class QuasiStatic : public BaseSolver {
public:
  QuasiStatic(Tissue *T, std::istream &in);
  void simulate() override;

private:
  size_t relax();
  size_t relaxBB();   // Barzilai-Borwein alternative; see quasi_static.cpp
  double maxForce() const;
  void calibrateStep();
  static void addScaled(Matrix &y, const Matrix &k, double h);
  static void combineHeun(Matrix &y, const Matrix &y0, const Matrix &k1,
                          const Matrix &k2, double h);
  void restorePositional(Matrix &y, const Matrix &y0) const;
  void forceOnly();
  void applyPrescribed(double h);

  double hGrowth_ = 0.0;
  double forceTol_ = 1e-3;
  int maxRelax_ = 2000;
  double dt0_ = 0.0, dtMax_ = 0.0, scaleHint_ = 1.0;
  double firstProbe_ = 0.0;  // single-probe estimate, for the diagnostic
  unsigned int relaxNotConverged_ = 0;
  int relaxMethod_ = 0;   // 0 = FIRE, 1 = Barzilai-Borwein
  Matrix velocity_;
  Matrix vertexVel_;              // prescribed (non-relaxable) velocity
  bool hasPrescribed_ = false;
  Matrix cellStart_, wallStart_, cellK1_, wallK1_;
  std::vector<size_t> posIndex_;   // flat cellData indices that are positions
  std::vector<double> posVel_;     // their FIRE velocities
};

class RK4 : public BaseSolver {
public:
  RK4(Tissue *T, std::istream &in);
  void simulate() override;

private:
  void rk4Step();
  double h_ = 0.0;
  Matrix ytCell_, ytWall_, ytVertex_;
  Matrix dytCell_, dytWall_, dytVertex_;
  Matrix dymCell_, dymWall_, dymVertex_;
};

class RK5Adaptive : public BaseSolver {
public:
  RK5Adaptive(Tissue *T, std::istream &in);
  void simulate() override;

private:
  void rkqs(double hTry, double &hDid, double &hNext);
  void rkck(double h);
  void allocateScratch();

  double h1_ = 0.0;  // maximal (and initial) step
  double eps_ = 0.0; // error tolerance

  Matrix yScalC_, yScalW_, yScalV_;
  Matrix yTempC_, yTempW_, yTempV_;
  Matrix yErrC_, yErrW_, yErrV_;
  Matrix ak2C_, ak2W_, ak2V_;
  Matrix ak3C_, ak3W_, ak3V_;
  Matrix ak4C_, ak4W_, ak4V_;
  Matrix ak5C_, ak5W_, ak5V_;
  Matrix ak6C_, ak6W_, ak6V_;
  Matrix yTrialC_, yTrialW_, yTrialV_;
};

class HeunIto : public BaseSolver {
public:
  HeunIto(Tissue *T, std::istream &in);
  void simulate() override;

private:
  void heunStep();
  double h_ = 0.0;
  double vol_ = 1.0;
  int volFlag_ = 0;
  Matrix sdydtC_, sdydtW_, sdydtV_;
  Matrix stC_, stW_, stV_;
  Matrix randC_, randW_, randV_;
  Matrix y1C_, y1W_, y1V_;
  Matrix dydt2C_, dydt2W_, dydt2V_;
};

} // namespace tissue

#endif
