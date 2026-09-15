#include "tissue/solvers/base_solver.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sys/resource.h>

#include "tissue/io/config.h"
#include "tissue/io/file_utils.h"
#include "tissue/io/vtk.h"
#include "tissue/solvers/solvers.h"

namespace tissue {

namespace {
// User CPU seconds since the previous call (legacy myTimes::getDiffTime).
double cpuDiffTime() {
  static double previous = 0.0;
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  double now = static_cast<double>(usage.ru_utime.tv_sec) +
               1e-6 * static_cast<double>(usage.ru_utime.tv_usec);
  double diff = now - previous;
  previous = now;
  return diff;
}
} // namespace

std::unique_ptr<BaseSolver> BaseSolver::getSolver(Tissue *T,
                                                  const std::string &file) {
  auto in = io::openCommentFiltered(file);
  if (!in) {
    std::cerr << "BaseSolver::getSolver() - Cannot open file " << file
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string idValue;
  *in >> idValue;
  if (idValue == "RK5Adaptive")
    return std::make_unique<RK5Adaptive>(T, *in);
  if (idValue == "QuasiStatic")
    return std::make_unique<QuasiStatic>(T, *in);
  if (idValue == "RK4")
    return std::make_unique<RK4>(T, *in);
  if (idValue == "Euler")
    return std::make_unique<Euler>(T, *in);
  if (idValue == "HeunIto")
    return std::make_unique<HeunIto>(T, *in);
  std::cerr << "BaseSolver::getSolver() - Unknown solver: " << idValue
            << " (known: RK5Adaptive, QuasiStatic, RK4, Euler, HeunIto)" 
            << std::endl;
  std::exit(EXIT_FAILURE);
}

void BaseSolver::getInit() {
  cellData_ = T_->cellState();
  wallData_ = T_->wallState();
  vertexData_ = T_->vertexState();
  cellDerivs_.reshapeLike(cellData_);
  wallDerivs_.reshapeLike(wallData_);
  vertexDerivs_.reshapeLike(vertexData_);
}

void BaseSolver::setTissueVariables(size_t numCellVariable) {
  if (cellData_.rows() != T_->numCell() || wallData_.rows() != T_->numWall() ||
      vertexData_.rows() != T_->numVertex()) {
    std::cerr << "BaseSolver::setTissueVariables wrong size of data:"
              << std::endl
              << "Cells: " << T_->numCell() << " vs " << cellData_.rows()
              << ", Walls: " << T_->numWall() << " vs " << wallData_.rows()
              << ", Vertices: " << T_->numVertex() << " vs "
              << vertexData_.rows() << std::endl;
    std::exit(-1);
  }
  Matrix &cellState = T_->cellState();
  if (numCellVariable == static_cast<size_t>(-1)) {
    cellState = cellData_;
  } else {
    for (size_t i = 0; i < cellData_.rows(); ++i)
      for (size_t j = 0; j < numCellVariable; ++j)
        cellState[i][j] = cellData_[i][j];
  }
  T_->wallState() = wallData_;
  T_->vertexState() = vertexData_;
}

void BaseSolver::postStep(double h) {
  T_->updateDirection(h, cellData_, wallData_, vertexData_, cellDerivs_,
                      wallDerivs_, vertexDerivs_);
  T_->updateReactions(cellData_, wallData_, vertexData_, h);
  T_->checkCompartmentChange(cellData_, wallData_, vertexData_, cellDerivs_,
                             wallDerivs_, vertexDerivs_);
  T_->checkConnectivity(1);
}

void BaseSolver::initPrintSchedule(double tiny) {
  printTime_ = endTime_ + tiny;
  printDeltaTime_ = endTime_ + 2.0 * tiny;
  doPrint_ = true;
  if (numPrint_ <= 0)
    doPrint_ = false;
  else if (numPrint_ == 1) {
    // only the final print
  } else if (numPrint_ == 2) {
    printTime_ = startTime_ - tiny;
  } else {
    printTime_ = startTime_ - tiny;
    printDeltaTime_ = (endTime_ - startTime_) / static_cast<double>(numPrint_ - 1);
  }
}

void BaseSolver::print(std::ostream &os) {
  double time = cpuDiffTime();
  std::cerr << tCount_ << " " << t_ << " " << cellData_.rows() << " "
            << wallData_.rows() << " " << vertexData_.rows() << " " << numOk_
            << " " << numBad_ << "  " << t_ - tOld_ << " "
            << static_cast<int>(cellData_.rows()) - nOld_ << " "
            << numOk_ - okOld_ << " " << numBad_ - badOld_ << " " << time
            << std::endl;
  tOld_ = t_;
  nOld_ = static_cast<int>(cellData_.rows());
  okOld_ = numOk_;
  badOld_ = numBad_;

  std::string vtkOutputFolder = config::getValue("vtk_output", 0);
  if (vtkOutputFolder.empty())
    vtkOutputFolder = "vtk";

  const size_t dimension = vertexData_.cols();

  auto printVertexBlock = [&]() {
    os << vertexData_.rows() << " " << dimension << std::endl;
    for (size_t i = 0; i < vertexData_.rows(); ++i) {
      for (size_t d = 0; d < dimension; ++d)
        os << vertexData_[i][d] << " ";
      os << std::endl;
    }
  };

  if (printFlag_ == 0) {
    if (tCount_ == 0)
      os << numPrint_ << "\n";
    if (vertexData_.rows() == 0) {
      os << "0 0" << std::endl << "0 0" << std::endl;
      return; // legacy: tCount not incremented
    }
    printVertexBlock();
    const size_t numCellVar = T_->numCellVariable();
    os << cellData_.rows() << " " << numCellVar + 2 << std::endl;
    for (size_t i = 0; i < cellData_.rows(); ++i) {
      const CellTopo &cell = T_->cell(i);
      os << cell.numVertex() << " ";
      for (size_t v : cell.vertices)
        os << v << " ";
      for (size_t j = 0; j < numCellVar; ++j)
        os << cellData_[i][j] << " ";
      os << T_->cellVolume(i, vertexData_) << " " << cell.numWall()
         << std::endl;
    }
    os << wallData_.rows() << " " << (wallData_.cols() - 1) + 5 << std::endl;
    for (size_t i = 0; i < wallData_.rows(); ++i) {
      double distance = T_->wallLengthFromVertices(i, vertexData_);
      os << T_->wall(i).vertex1 << " " << T_->wall(i).vertex2 << " ";
      for (size_t j = 0; j < wallData_.rowSize(i); ++j)
        os << wallData_[i][j] << " ";
      os << i << " " << distance << " " << distance - wallData_[i][0] << " "
         << (distance - wallData_[i][0]) / wallData_[i][0] << std::endl;
    }
    os << std::endl;
  } else if (printFlag_ == 1 || printFlag_ == 2) {
    std::filesystem::create_directories(vtkOutputFolder);
    std::string pvdFile = vtkOutputFolder + "/tissue.pvd";
    std::string cellFile = vtkOutputFolder + "/VTK_cells.vtu";
    std::string wallFile = vtkOutputFolder + "/VTK_walls.vtu";
    setTissueVariables(T_->numCellVariable());
    if (tCount_ == 0)
      io::writeFullPvd(pvdFile, cellFile, wallFile,
                       static_cast<size_t>(numPrint_));
    io::writeVtu(*T_, cellFile, wallFile, static_cast<size_t>(tCount_),
                 printFlag_ == 1);
  } else if (printFlag_ == 3) {
    if (tCount_ == 0)
      os << numPrint_ << "\n";
    printVertexBlock();
    os << cellData_.rows() << " " << cellData_.cols() + 3 << std::endl;
    for (size_t i = 0; i < cellData_.rows(); ++i) {
      const CellTopo &cell = T_->cell(i);
      os << cell.numVertex() << " ";
      for (size_t v : cell.vertices)
        os << v << " ";
      for (size_t j = 0; j < cellData_.rowSize(i); ++j)
        os << cellData_[i][j] << " ";
      os << i << " " << T_->cellVolume(i, vertexData_) << " " << cell.numWall()
         << std::endl;
    }
  } else if (printFlag_ == 4) {
    if (tCount_ == 0)
      os << numPrint_ << "\n";
    printVertexBlock();
    os << wallData_.rows() << " " << (wallData_.cols() - 1) + 5 << std::endl;
    for (size_t i = 0; i < wallData_.rows(); ++i) {
      double distance = T_->wallLengthFromVertices(i, vertexData_);
      os << "2 " << T_->wall(i).vertex1 << " " << T_->wall(i).vertex2 << " ";
      for (size_t j = 0; j < wallData_.rowSize(i); ++j)
        os << wallData_[i][j] << " ";
      os << i << " " << distance << " " << distance - wallData_[i][0] << " "
         << (distance - wallData_[i][0]) / wallData_[i][0] << std::endl;
    }
  } else if (printFlag_ == 5) {
    for (size_t i = 0; i < cellData_.rows(); ++i) {
      os << "0 " << i << " " << t_ << " ";
      for (size_t j = 0; j < cellData_.rowSize(i); ++j)
        os << cellData_[i][j] << " ";
      os << i << " " << T_->cellVolume(i, vertexData_) << " "
         << T_->cell(i).numWall() << std::endl;
    }
    for (size_t i = 0; i < wallData_.rows(); ++i) {
      os << "1 " << i << " " << t_ << " ";
      for (size_t j = 0; j < wallData_.rowSize(i); ++j)
        os << wallData_[i][j] << " ";
      double distance = T_->wallLengthFromVertices(i, vertexData_);
      os << i << " " << distance << " " << distance - wallData_[i][0]
         << std::endl;
    }
    os << std::endl;
  } else if (printFlag_ == 77) {
    for (size_t i = 0; i < cellData_.rows(); ++i) {
      Vec3 o = T_->cellPosition(i, vertexData_);
      os << i << " " << T_->cellVolume(i, vertexData_) << " "
         << T_->cell(i).numVertex() << " " << o[0] << " " << o[1] << " ";
      for (size_t v : T_->cell(i).vertices)
        os << vertexData_[v][0] << " " << vertexData_[v][1] << " ";
      os << std::endl;
    }
  } else if (printFlag_ == 107) {
    printInit(os);
  } else {
    std::cerr << "BaseSolver::print() Wrong printFlag value (" << printFlag_
              << " not implemented in tissue v2)." << std::endl;
  }
  tCount_++;
}

void BaseSolver::printInit(std::ostream &os) const {
  T_->printInit(cellData_, wallData_, vertexData_, os);
}

void BaseSolver::printInitCenterTri(std::ostream &os) const {
  T_->printInitCenterTri(cellData_, wallData_, vertexData_, os);
}

} // namespace tissue
