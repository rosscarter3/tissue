#include <cmath>
#include <cstdlib>
#include <iostream>

#include "tissue/solvers/detail.h"
#include "tissue/solvers/solvers.h"

namespace tissue {

// --- RK4 (fixed step) ---------------------------------------------------------

RK4::RK4(Tissue *T, std::istream &in) : BaseSolver(T) {
  getInit();
  in >> startTime_ >> endTime_;
  t_ = startTime_;
  in >> printFlag_ >> numPrint_;
  in >> h_;
}

void RK4::simulate() {
  if (h_ <= 0.0 || endTime_ - startTime_ <= 0.0) {
    std::cerr << "RK4::simulate() Wrong time borders or time step for "
                 "simulation. No simulation performed." << std::endl;
    return;
  }
  std::cerr << "Simulating using fourth-order Runge-Kutta" << std::endl;
  T_->initiateReactions(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  cellDerivs_.reshapeLike(cellData_);
  wallDerivs_.reshapeLike(wallData_);
  vertexDerivs_.reshapeLike(vertexData_);
  T_->initiateDirection(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  initPrintSchedule(1e-10);
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  while (t_ < endTime_) {
    if (doPrint_ && t_ >= printTime_) {
      printTime_ += printDeltaTime_;
      print();
    }
    rk4Step();
    numOk_++;
    postStep(h_);
    if ((t_ + h_) == t_) {
      std::cerr << "RK4::simulate() Step size too small." << std::endl;
      std::exit(-1);
    }
    t_ += h_;
  }
  if (doPrint_) {
    T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
               vertexDerivs_);
    print();
  }
  std::cerr << "Simulation done." << std::endl;
}

void RK4::rk4Step() {
  const double hh = 0.5 * h_;
  const double h6 = h_ / 6.0;
  ytCell_.reshapeLike(cellData_);
  ytWall_.reshapeLike(wallData_);
  ytVertex_.reshapeLike(vertexData_);
  dytCell_.reshapeLike(cellData_);
  dytWall_.reshapeLike(wallData_);
  dytVertex_.reshapeLike(vertexData_);
  dymCell_.reshapeLike(cellData_);
  dymWall_.reshapeLike(wallData_);
  dymVertex_.reshapeLike(vertexData_);

  // k1
  T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
             vertexDerivs_);
  detail::assignAxpy(ytCell_, cellData_, cellDerivs_, hh);
  detail::assignAxpy(ytWall_, wallData_, wallDerivs_, hh);
  detail::assignAxpy(ytVertex_, vertexData_, vertexDerivs_, hh);
  // k2
  T_->derivs(ytCell_, ytWall_, ytVertex_, dytCell_, dytWall_, dytVertex_);
  detail::assignAxpy(ytCell_, cellData_, dytCell_, hh);
  detail::assignAxpy(ytWall_, wallData_, dytWall_, hh);
  detail::assignAxpy(ytVertex_, vertexData_, dytVertex_, hh);
  // k3
  T_->derivs(ytCell_, ytWall_, ytVertex_, dymCell_, dymWall_, dymVertex_);
  detail::assignAxpy(ytCell_, cellData_, dymCell_, h_);
  detail::assignAxpy(ytWall_, wallData_, dymWall_, h_);
  detail::assignAxpy(ytVertex_, vertexData_, dymVertex_, h_);
  // dym := k2 + k3 (legacy dropped k3 for vertices; fixed here)
  detail::addScaled(dymCell_, dytCell_, 1.0);
  detail::addScaled(dymWall_, dytWall_, 1.0);
  detail::addScaled(dymVertex_, dytVertex_, 1.0);
  // k4
  T_->derivs(ytCell_, ytWall_, ytVertex_, dytCell_, dytWall_, dytVertex_);
  // y += h/6 * (k1 + k4 + 2*(k2 + k3))
  detail::addScaled(cellData_, cellDerivs_, h6);
  detail::addScaled(cellData_, dytCell_, h6);
  detail::addScaled(cellData_, dymCell_, 2.0 * h6);
  detail::addScaled(wallData_, wallDerivs_, h6);
  detail::addScaled(wallData_, dytWall_, h6);
  detail::addScaled(wallData_, dymWall_, 2.0 * h6);
  detail::addScaled(vertexData_, vertexDerivs_, h6);
  detail::addScaled(vertexData_, dytVertex_, h6);
  detail::addScaled(vertexData_, dymVertex_, 2.0 * h6);
}

// --- RK5Adaptive (Cash-Karp with step control) -----------------------------------

RK5Adaptive::RK5Adaptive(Tissue *T, std::istream &in) : BaseSolver(T) {
  getInit();
  in >> startTime_ >> endTime_;
  t_ = startTime_;
  in >> printFlag_ >> numPrint_;
  in >> h1_ >> eps_;
}

void RK5Adaptive::allocateScratch() {
  for (Matrix *m : {&yScalC_, &yTempC_, &yErrC_, &ak2C_, &ak3C_, &ak4C_,
                    &ak5C_, &ak6C_, &yTrialC_})
    m->reshapeLike(cellData_);
  for (Matrix *m : {&yScalW_, &yTempW_, &yErrW_, &ak2W_, &ak3W_, &ak4W_,
                    &ak5W_, &ak6W_, &yTrialW_})
    m->reshapeLike(wallData_);
  for (Matrix *m : {&yScalV_, &yTempV_, &yErrV_, &ak2V_, &ak3V_, &ak4V_,
                    &ak5V_, &ak6V_, &yTrialV_})
    m->reshapeLike(vertexData_);
}

void RK5Adaptive::simulate() {
  if (h1_ <= 0.0 || endTime_ - startTime_ <= 0.0) {
    std::cerr << "RK5Adaptive::simulate() Wrong time borders or time step for "
                 "simulation. No simulation performed." << std::endl;
    std::exit(-1);
  }
  const double tiny = 1e-9 * eps_; // legacy value (NR uses 1e-30)
  T_->initiateReactions(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  cellDerivs_.reshapeLike(cellData_);
  wallDerivs_.reshapeLike(wallData_);
  vertexDerivs_.reshapeLike(vertexData_);
  T_->initiateDirection(cellData_, wallData_, vertexData_, cellDerivs_,
                        wallDerivs_, vertexDerivs_);
  initPrintSchedule(tiny);
  allocateScratch();
  double h = h1_;
  t_ = startTime_;
  numOk_ = numBad_ = 0;
  for (;;) {
    // k1 for this step attempt.
    T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
               vertexDerivs_);
    // Error scaling with the pre-clamp step (legacy order).
    auto fillScal = [&](Matrix &scal, const Matrix &y, const Matrix &dydx) {
      auto sf = scal.flat();
      auto yf = y.flat();
      auto df = dydx.flat();
      parallelFor(
          sf.size(),
          [&](size_t b, size_t e) {
            for (size_t k = b; k < e; ++k)
              sf[k] = std::fabs(yf[k]) + std::fabs(df[k] * h) + tiny;
          },
          detail::kGrain);
    };
    fillScal(yScalC_, cellData_, cellDerivs_);
    fillScal(yScalW_, wallData_, wallDerivs_);
    fillScal(yScalV_, vertexData_, vertexDerivs_);

    if (doPrint_ && t_ >= printTime_) {
      printTime_ += printDeltaTime_;
      print();
    }
    double tMin = endTime_ < printTime_ ? endTime_ : printTime_;
    if (t_ + h > tMin)
      h = tMin - t_;

    double hDid = 0.0, hNext = 0.0;
    rkqs(h, hDid, hNext);
    if (hDid == h)
      ++numOk_;
    else
      ++numBad_;

    postStep(h);
    if (!cellData_.sameShape(yScalC_) || !wallData_.sameShape(yScalW_) ||
        !vertexData_.sameShape(yScalV_))
      allocateScratch();

    if (t_ >= endTime_) {
      if (doPrint_) {
        T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
                   vertexDerivs_);
        print();
      }
      std::cerr << "Simulation done." << std::endl;
      return;
    }
    h = hNext;
    if (h > h1_)
      h = h1_;
  }
}

void RK5Adaptive::rkqs(double hTry, double &hDid, double &hNext) {
  constexpr double SAFETY = 0.9;
  constexpr double PGROW = -0.2;
  constexpr double PSHRNK = -0.25;
  constexpr double ERRCON = 1.89e-4;

  double h = hTry;
  double errMax;
  for (;;) {
    rkck(h);
    errMax = detail::maxErrRatio(yErrC_, yScalC_);
    errMax = std::max(errMax, detail::maxErrRatio(yErrW_, yScalW_));
    errMax = std::max(errMax, detail::maxErrRatio(yErrV_, yScalV_));
    errMax /= eps_;
    if (errMax <= 1.0)
      break;
    double hTemp = SAFETY * h * std::pow(errMax, PSHRNK);
    if (h >= 0.0)
      h = hTemp > 0.1 * h ? hTemp : 0.1 * h;
    else
      h = hTemp > 0.1 * h ? 0.1 * h : hTemp;
    if (t_ + h == t_) {
      std::cerr << "Warning stepsize underflow in RK5Adaptive::rkqs"
                << std::endl;
      std::exit(-1);
    }
  }
  hNext = errMax > ERRCON ? SAFETY * h * std::pow(errMax, PGROW) : 5.0 * h;
  t_ += (hDid = h);
  cellData_.copyFrom(yTempC_);
  wallData_.copyFrom(yTempW_);
  vertexData_.copyFrom(yTempV_);
}

void RK5Adaptive::rkck(double h) {
  // Cash-Karp coefficients (Numerical Recipes).
  constexpr double b21 = 0.2;
  constexpr double b31 = 3.0 / 40.0, b32 = 9.0 / 40.0;
  constexpr double b41 = 0.3, b42 = -0.9, b43 = 1.2;
  constexpr double b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
                   b54 = 35.0 / 27.0;
  constexpr double b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
                   b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0,
                   b65 = 253.0 / 4096.0;
  constexpr double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0,
                   c6 = 512.0 / 1771.0;
  constexpr double dc5 = -277.00 / 14336.0;
  constexpr double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
                   dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;

  // Per-group linear combination out = y + h * sum(a_i * k_i), fused into a
  // single parallel pass (the hot memory-bound loop of the integrator).
  auto combine = [&](Matrix &out, const Matrix &y,
                     std::initializer_list<std::pair<double, const Matrix *>>
                         terms) {
    auto of = out.flat();
    auto yf = y.flat();
    double a[5];
    const double *k[5];
    size_t nt = 0;
    for (auto &[c, m] : terms) {
      a[nt] = c * h;
      k[nt] = m->flat().data();
      ++nt;
    }
    parallelFor(
        of.size(),
        [&](size_t b, size_t e) {
          switch (nt) {
          case 1:
            for (size_t i = b; i < e; ++i)
              of[i] = yf[i] + a[0] * k[0][i];
            break;
          case 2:
            for (size_t i = b; i < e; ++i)
              of[i] = yf[i] + a[0] * k[0][i] + a[1] * k[1][i];
            break;
          case 3:
            for (size_t i = b; i < e; ++i)
              of[i] = yf[i] + a[0] * k[0][i] + a[1] * k[1][i] + a[2] * k[2][i];
            break;
          case 4:
            for (size_t i = b; i < e; ++i)
              of[i] = yf[i] + a[0] * k[0][i] + a[1] * k[1][i] +
                      a[2] * k[2][i] + a[3] * k[3][i];
            break;
          default:
            for (size_t i = b; i < e; ++i)
              of[i] = yf[i] + a[0] * k[0][i] + a[1] * k[1][i] +
                      a[2] * k[2][i] + a[3] * k[3][i] + a[4] * k[4][i];
          }
        },
        detail::kGrain);
  };

  const Matrix &k1C = cellDerivs_, &k1W = wallDerivs_, &k1V = vertexDerivs_;

  combine(yTrialC_, cellData_, {{b21, &k1C}});
  combine(yTrialW_, wallData_, {{b21, &k1W}});
  combine(yTrialV_, vertexData_, {{b21, &k1V}});
  T_->derivs(yTrialC_, yTrialW_, yTrialV_, ak2C_, ak2W_, ak2V_);

  combine(yTrialC_, cellData_, {{b31, &k1C}, {b32, &ak2C_}});
  combine(yTrialW_, wallData_, {{b31, &k1W}, {b32, &ak2W_}});
  combine(yTrialV_, vertexData_, {{b31, &k1V}, {b32, &ak2V_}});
  // Legacy wrote the wall output of this stage into ak2W (clobbering stage 2
  // and leaving ak3W stale); fixed to ak3W here.
  T_->derivs(yTrialC_, yTrialW_, yTrialV_, ak3C_, ak3W_, ak3V_);

  combine(yTrialC_, cellData_, {{b41, &k1C}, {b42, &ak2C_}, {b43, &ak3C_}});
  combine(yTrialW_, wallData_, {{b41, &k1W}, {b42, &ak2W_}, {b43, &ak3W_}});
  combine(yTrialV_, vertexData_, {{b41, &k1V}, {b42, &ak2V_}, {b43, &ak3V_}});
  T_->derivs(yTrialC_, yTrialW_, yTrialV_, ak4C_, ak4W_, ak4V_);

  combine(yTrialC_, cellData_,
          {{b51, &k1C}, {b52, &ak2C_}, {b53, &ak3C_}, {b54, &ak4C_}});
  combine(yTrialW_, wallData_,
          {{b51, &k1W}, {b52, &ak2W_}, {b53, &ak3W_}, {b54, &ak4W_}});
  combine(yTrialV_, vertexData_,
          {{b51, &k1V}, {b52, &ak2V_}, {b53, &ak3V_}, {b54, &ak4V_}});
  T_->derivs(yTrialC_, yTrialW_, yTrialV_, ak5C_, ak5W_, ak5V_);

  combine(yTrialC_, cellData_,
          {{b61, &k1C}, {b62, &ak2C_}, {b63, &ak3C_}, {b64, &ak4C_}, {b65, &ak5C_}});
  combine(yTrialW_, wallData_,
          {{b61, &k1W}, {b62, &ak2W_}, {b63, &ak3W_}, {b64, &ak4W_}, {b65, &ak5W_}});
  combine(yTrialV_, vertexData_,
          {{b61, &k1V}, {b62, &ak2V_}, {b63, &ak3V_}, {b64, &ak4V_}, {b65, &ak5V_}});
  T_->derivs(yTrialC_, yTrialW_, yTrialV_, ak6C_, ak6W_, ak6V_);

  combine(yTempC_, cellData_,
          {{c1, &k1C}, {c3, &ak3C_}, {c4, &ak4C_}, {c6, &ak6C_}});
  combine(yTempW_, wallData_,
          {{c1, &k1W}, {c3, &ak3W_}, {c4, &ak4W_}, {c6, &ak6W_}});
  combine(yTempV_, vertexData_,
          {{c1, &k1V}, {c3, &ak3V_}, {c4, &ak4V_}, {c6, &ak6V_}});

  auto errCombine = [&](Matrix &err,
                        std::initializer_list<std::pair<double, const Matrix *>>
                            terms) {
    auto ef = err.flat();
    double a[5];
    const double *k[5];
    size_t nt = 0;
    for (auto &[c, m] : terms) {
      a[nt] = c * h;
      k[nt] = m->flat().data();
      ++nt;
    }
    parallelFor(
        ef.size(),
        [&](size_t b, size_t e) {
          for (size_t i = b; i < e; ++i)
            ef[i] = a[0] * k[0][i] + a[1] * k[1][i] + a[2] * k[2][i] +
                    a[3] * k[3][i] + a[4] * k[4][i];
        },
        detail::kGrain);
  };
  errCombine(yErrC_,
             {{dc1, &k1C}, {dc3, &ak3C_}, {dc4, &ak4C_}, {dc5, &ak5C_}, {dc6, &ak6C_}});
  errCombine(yErrW_,
             {{dc1, &k1W}, {dc3, &ak3W_}, {dc4, &ak4W_}, {dc5, &ak5W_}, {dc6, &ak6W_}});
  errCombine(yErrV_,
             {{dc1, &k1V}, {dc3, &ak3V_}, {dc4, &ak4V_}, {dc5, &ak5V_}, {dc6, &ak6V_}});
}

} // namespace tissue
