#include <cstdlib>
#include <iostream>

#include "tissue/solvers/detail.h"
#include "tissue/solvers/solvers.h"

namespace tissue {

Euler::Euler(Tissue *T, std::istream &in) : BaseSolver(T) {
  getInit();
  in >> startTime_ >> endTime_;
  t_ = startTime_;
  in >> printFlag_ >> numPrint_;
  in >> h_;
}

void Euler::simulate() {
  if (h_ <= 0.0 || endTime_ - startTime_ <= 0.0) {
    std::cerr << "Euler::simulate() Wrong time borders or time step for "
                 "simulation. No simulation performed." << std::endl;
    return;
  }
  std::cerr << "Simulating using explicit Euler." << std::endl;
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
    T_->derivs(cellData_, wallData_, vertexData_, cellDerivs_, wallDerivs_,
               vertexDerivs_);
    detail::addScaled(cellData_, cellDerivs_, h_);
    detail::addScaled(wallData_, wallDerivs_, h_);
    detail::addScaled(vertexData_, vertexDerivs_, h_);
    numOk_++;
    postStep(h_);
    cellDerivs_.reshapeLike(cellData_);
    wallDerivs_.reshapeLike(wallData_);
    vertexDerivs_.reshapeLike(vertexData_);
    if ((t_ + h_) == t_) {
      std::cerr << "Euler::simulate() Step size too small." << std::endl;
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

} // namespace tissue
