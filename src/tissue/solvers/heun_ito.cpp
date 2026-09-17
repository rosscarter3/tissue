//
// Heun predictor-corrector with Itô noise (Carrillo et al. 2003 PRE):
// noise amplitude per element is sqrt(|flux sum| * h / volume), the same
// Gaussian increment is used in predictor and corrector, and cell/wall
// variables have an absorbing barrier at 0.
//
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "tissue/core/random.h"
#include "tissue/solvers/detail.h"
#include "tissue/solvers/solvers.h"

namespace tissue {

HeunIto::HeunIto(Tissue *T, std::istream &in) : BaseSolver(T) {
  getInit();
  in >> startTime_ >> endTime_;
  t_ = startTime_;
  in >> printFlag_ >> numPrint_;
  in >> h_ >> vol_ >> volFlag_;
}

void HeunIto::simulate() {
  if (h_ <= 0.0 || endTime_ - startTime_ <= 0.0) {
    std::cerr << "HeunIto::simulate() Wrong time borders or time step for "
                 "simulation. No simulation performed." << std::endl;
    return;
  }
  std::cerr << "Simulating using an HeunIto solver" << std::endl;
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
    heunStep();
    numOk_++;
    postStep(h_);
    if ((t_ + h_) == t_) {
      std::cerr << "HeunIto::simulate() Step size too small." << std::endl;
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

void HeunIto::heunStep() {
  const double hh = 0.5 * h_;
  for (Matrix *m : {&sdydtC_, &stC_, &randC_, &y1C_, &dydt2C_})
    m->reshapeLike(cellData_);
  for (Matrix *m : {&sdydtW_, &stW_, &randW_, &y1W_, &dydt2W_})
    m->reshapeLike(wallData_);
  for (Matrix *m : {&sdydtV_, &stV_, &randV_, &y1V_, &dydt2V_})
    m->reshapeLike(vertexData_);

  T_->derivsWithAbs(cellData_, wallData_, vertexData_, cellDerivs_,
                    wallDerivs_, vertexDerivs_, sdydtC_, sdydtW_, sdydtV_);

  // Predictor with noise. RNG draw order (cells, walls, vertices; row-major)
  // matches legacy for reproducibility, so this stays serial.
  for (size_t i = 0; i < cellData_.rows(); ++i) {
    double volume = 1.0;
    if (volFlag_ == 1)
      volume = T_->cellVolume(i, vertexData_);
    else if (volFlag_ > 1) {
      std::cerr << "HeunIto::heunStep() volFlag " << volFlag_
                << " not supported (0 or 1)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (size_t j = 0; j < cellData_.rowSize(i); ++j) {
      randC_[i][j] = random::Grand();
      stC_[i][j] = std::sqrt(sdydtC_[i][j] * h_ / (vol_ * volume));
      y1C_[i][j] = cellData_[i][j] + h_ * cellDerivs_[i][j] +
                   stC_[i][j] * randC_[i][j];
      if (y1C_[i][j] < 0.0)
        y1C_[i][j] = 0.0;
    }
  }
  for (size_t i = 0; i < wallData_.rows(); ++i)
    for (size_t j = 0; j < wallData_.rowSize(i); ++j) {
      randW_[i][j] = random::Grand();
      stW_[i][j] = std::sqrt(sdydtW_[i][j] * h_ / vol_);
      y1W_[i][j] = wallData_[i][j] + h_ * wallDerivs_[i][j] +
                   stW_[i][j] * randW_[i][j];
      if (y1W_[i][j] < 0.0)
        y1W_[i][j] = 0.0;
    }
  for (size_t i = 0; i < vertexData_.rows(); ++i)
    for (size_t j = 0; j < vertexData_.rowSize(i); ++j) {
      randV_[i][j] = random::Grand();
      stV_[i][j] = std::sqrt(sdydtV_[i][j] * h_ / vol_);
      y1V_[i][j] = vertexData_[i][j] + h_ * vertexDerivs_[i][j] +
                   stV_[i][j] * randV_[i][j];
    }

  // Corrector reusing the same noise increment.
  T_->derivs(y1C_, y1W_, y1V_, dydt2C_, dydt2W_, dydt2V_);
  auto correct = [&](Matrix &y, const Matrix &k1, const Matrix &k2,
                     const Matrix &st, const Matrix &rnd, bool barrier) {
    auto yf = y.flat();
    auto k1f = k1.flat();
    auto k2f = k2.flat();
    auto stf = st.flat();
    auto rf = rnd.flat();
    for (size_t k = 0; k < yf.size(); ++k) {
      yf[k] = yf[k] + hh * (k1f[k] + k2f[k]) + stf[k] * rf[k];
      if (barrier && yf[k] < 0.0)
        yf[k] = 0.0;
    }
  };
  correct(cellData_, cellDerivs_, dydt2C_, stC_, randC_, true);
  correct(wallData_, wallDerivs_, dydt2W_, stW_, randW_, true);
  correct(vertexData_, vertexDerivs_, dydt2V_, stV_, randV_, false);
}

} // namespace tissue
