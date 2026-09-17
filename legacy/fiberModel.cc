//
// Filename     : fiberModel.cc
// Description  : Classes describing updates of mechanical properties, e.g.
// Young's modulii given input
//                from elsewhere
// Author(s)    : Henrik Jonsson (henrik.jonsson@slcu.cam.ac.uk)
// Created      : January 2020
// Revision     : $Id:$
//

#include "baseReaction.h"
#include "fiberModel.h"
#include "tissue.h"

namespace FiberModel {

double isotropicFunction(double youngMatrix, double youngFiber) {
  return youngMatrix +
         0.5 * youngFiber; // youngL = youngT = youngMatrix + 0.5*youngFiber
}
double linearFunction(double input, double youngMatrix, double youngFiber) {
  return youngMatrix + 0.5 * (1.0 + input) * youngFiber;
}
double hillFunction(double input, double youngMatrix, double youngFiber,
                    double Kh, double Nh) {
  return youngMatrix +
         0.5 *
             (1 + (std::pow(input, Nh) /
                   (std::pow((1. - input), Nh) * std::pow(Kh, Nh) +
                    std::pow(input, Nh)))) *
             youngFiber;
}
double transverseFromLongitudinal(double youngL, double youngMatrix,
                                  double youngFiber) {
  return 2 * youngMatrix + youngFiber - youngL;
}

General::General(std::vector<double> &paraValue,
                 std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() != 8) {
    std::cerr
        << "FiberModel::General::General() Uses eight parameters:" << std::endl
        << "k_rate, equilibrium_threshold (only update if vertex moving "
           "slower), "
        << std::endl
        << "linear-hill flag (0=linear,1=Hill,2=Hill_direct)," << std::endl
        << "K_hill, n_hill (only used if Hill versions selected)," << std::endl
        << "young_matrix, young_fiber, (material limits, should be same as for "
           "the mechanical (TRBS) model used)"
        << std::endl
        << "initialization flag (0=no initiation/1=isotropic/2=anisotropic"
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if (indValue.size() != 3 || indValue[0].size() != 1 ||
      (indValue[1].size() != 1 && indValue[1].size() != 2) ||
      indValue[2].size() != 1) {
    std::cerr << "FiberModel::General::General() " << std::endl
              << "First level gives stress/strain anisotropy index to read."
              << std::endl
              << "Second level gives Young_Longitudinal index to update "
                 "(p_2=1,2) or Y_L and F_L (p_2=2) "
              << "where the first one is read and the second updated."
              << std::endl
              << "Third level gives index for velocity to be compared with p_1."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::General");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "velocitythreshold";
  tmp[2] = "linear-hill-flag";
  tmp[3] = "k_hill";
  tmp[4] = "n_hill";
  tmp[5] = "Y_matrix";
  tmp[6] = "Y_fiber";
  tmp[7] = "init_flag";
  // if (parameter(7)==3)
  //   tmp[8] = "Poisson";

  setParameterId(tmp);
}

void General::initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                       DataMatrix &vertexData, DataMatrix &cellDerivs,
                       DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double Kh = parameter(3);
  double Nh = parameter(4);
  double youngMatrix = parameter(5);
  double youngFiber = parameter(6);

  if (parameter(7) == 1) { // initialize with iso material
    for (size_t cellIndex = 0; cellIndex < numCell;
         ++cellIndex) { // initiating with 0 anisotropy and isotropic material
      cellData[cellIndex][AnisoIndex] = 0;
      cellData[cellIndex][YoungLIndex] =
          youngMatrix +
          0.5 * youngFiber; // youngL = youngMatrix + 0.5*youngFiber;
    }
  }
  if (parameter(7) == 2) {
    for (size_t cellIndex = 0; cellIndex < numCell;
         ++cellIndex) { // initiating with current anisotropy
      double anisotropy = cellData[cellIndex][AnisoIndex];
      if (parameter(2) == 0) // linear
        cellData[cellIndex][YoungLIndex] =
            youngMatrix + 0.5 * (1 + anisotropy) * youngFiber;

      if (parameter(2) == 1 || parameter(2) == 2) // Hill
        cellData[cellIndex][YoungLIndex] =
            youngMatrix +
            0.5 *
                (1 + (std::pow(anisotropy, Nh) /
                      (std::pow((1 - anisotropy), Nh) * std::pow(Kh, Nh) +
                       std::pow(anisotropy, Nh)))) *
                youngFiber;
      // std::cerr<< cellData[cellIndex][variableIndex(1,0)] << std::endl;
    }
  }
}

void General::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                     DataMatrix &vertexData, DataMatrix &cellDerivs,
                     DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {} //

void General::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                     DataMatrix &vertexData, double h) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  size_t velocityIndex = variableIndex(2, 0);
  double Kh = parameter(3);
  double Nh = parameter(4);
  double youngMatrix = parameter(5);
  double youngFiber = parameter(6);
  // double poisson=parameter(7);

  if (parameter(0) == 0.0)
    return;

  for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
    double anisotropy = cellData[cellIndex][AnisoIndex];

    // --- Linear feedback rule ---
    if (parameter(2) == 0 // linear
        && cellData[cellIndex][velocityIndex] < parameter(1)
        //&& cellData[cellIndex][YoungLIndex] < youngMatrix+0.5*(1+anisotropy)*
        //youngFiber
        && cellData[cellIndex][YoungLIndex] <
               youngMatrix + youngFiber // not to exceed maximum when using
                                        // absolute stress anisotropy
    ) {                                 // linear Fiber model
      double deltaYoungL = youngMatrix + 0.5 * (1 + anisotropy) * youngFiber -
                           cellData[cellIndex][YoungLIndex];
      cellData[cellIndex][YoungLIndex] += parameter(0) * h * deltaYoungL;
      // cellData[cellIndex][YoungLIndex] = youngMatrix+0.5*(1+anisotropy)*
      // youngFiber;
    }

    if (parameter(2) == 1 // Hill
        && cellData[cellIndex][velocityIndex] < parameter(1)
        // && cellData[cellIndex][YoungLIndex] <
        // youngMatrix+
        // 0.5*(1+(std::pow(anisotropy,Nh)
        //         /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))*
        //         youngFiber
        && cellData[cellIndex][YoungLIndex] <
               youngMatrix + youngFiber // not to exceed maximum when using
                                        // absolute stress anisotropy
    ) {                                 // Hill-like Fiber model

      double deltaYoungL =
          youngMatrix +
          0.5 *
              (1 + (std::pow(anisotropy, Nh) /
                    (std::pow((1 - anisotropy), Nh) * std::pow(Kh, Nh) +
                     std::pow(anisotropy, Nh)))) *
              youngFiber -
          cellData[cellIndex][YoungLIndex];
      cellData[cellIndex][YoungLIndex] += parameter(0) * h * deltaYoungL;
      // cellData[cellIndex][YoungLIndex] = youngMatrix+
      //   0.5*(1+(std::pow(anisotropy,Nh)/(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)+std::pow(anisotropy,Nh))))*
      //   youngFiber;
    }
    if (parameter(2) == 2) { // Hill-like Fiber model instantenous update for
                             // fiber deposition model
      size_t fiberLIndex = variableIndex(1, 1);
      // youngFiber= cellData[cellIndex][YoungLIndex];

      cellData[cellIndex][YoungLIndex] =
          0.5 *
          (1 + (std::pow(anisotropy, Nh) /
                (std::pow((1 - anisotropy), Nh) * std::pow(Kh, Nh) +
                 std::pow(anisotropy, Nh)))) *
          youngFiber;
    }
  }
}

Linear::Linear(std::vector<double> &paraValue,
               std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() < 2 || paraValue.size() > 4) {
    std::cerr
        << "FiberModel::Linear::Linear() Uses at least 2 and up to 4 "
           "parameters:"
        << std::endl
        << "young_matrix, young_fiber, (material limits, should be same as for "
           "the mechanical "
        << "(TRBS) model used)" << std::endl
        << "[initialization_flag] (0 - no initiation; 1 - isotropic; 2 - "
           "anisotropic"
        << std::endl
        << "[equilibrium_threshold] (only update if vertex moving slower), "
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() != 2 && indValue.size() != 3) ||
      indValue[0].size() != 1 ||
      (indValue[1].size() != 1 && indValue[1].size() != 2) ||
      (indValue.size() == 3 && indValue[2].size() != 1)) {
    std::cerr << "FiberModel::Linear::Linear() Error in indices given:"
              << std::endl
              << "First level gives stress/strain anisotropy index to read."
              << std::endl
              << "Second level gives Young_Longitudinal index to update, or "
                 "Y_L and Y_T indices."
              << std::endl
              << "Third optional level gives index for velocity to be compared "
                 "with p_3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() == 3 && paraValue.size() < 4) ||
      (indValue.size() < 3 && paraValue.size() == 4)) {
    std::cerr << "FiberModel::Linear::Linear() 4th optional [eq_threshold] "
                 "parameter must be comparable "
              << "with variable index variable (optional) given at level 3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::Linear");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "Y_matrix";
  tmp[1] = "Y_fiber";
  if (numParameter() > 2)
    tmp[2] = "init_flag";
  if (numParameter() > 3)
    tmp[3] = "velocity_threshold";

  setParameterId(tmp);
}

void Linear::initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                      DataMatrix &vertexData, DataMatrix &cellDerivs,
                      DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {

  if (numParameter() < 3 || parameter(2) == 0)
    return;

  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(0);
  double youngFiber = parameter(1);

  if (parameter(2) == 1) { // initialize with iso material
    // initiating with 0 anisotropy and isotropic material
    double youngIso = FiberModel::isotropicFunction(youngMatrix, youngFiber);
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      cellData[cellIndex][AnisoIndex] = 0.0;
      cellData[cellIndex][YoungLIndex] = youngIso;
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] = youngIso;
      }
    }
  } else if (parameter(2) == 2) { // initiating with current anisotropy
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      double anisotropy = cellData[cellIndex][AnisoIndex];
      cellData[cellIndex][YoungLIndex] =
          FiberModel::linearFunction(anisotropy, youngMatrix, youngFiber);
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

void Linear::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                    DataMatrix &vertexData, DataMatrix &cellDerivs,
                    DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
} // do nothing

void Linear::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                    DataMatrix &vertexData, double h) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(0);
  double youngFiber = parameter(1);

  for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
    double anisotropy = cellData[cellIndex][AnisoIndex];

    // only update cell if v<v_th or v_th not given
    if (numParameter() < 4 || numVariableIndexLevel() < 3 ||
        cellData[cellIndex][variableIndex(2, 0)] < parameter(4)) {

      double youngL =
          FiberModel::linearFunction(anisotropy, youngMatrix, youngFiber);

      if (youngL > youngMatrix + youngFiber) // not to exceed maximum when using
                                             // absolute stress anisotropy
        youngL = youngMatrix + youngFiber;

      cellData[cellIndex][YoungLIndex] = youngL;
      if (numVariableIndex(1) > 1) { // Also save Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(youngL, youngMatrix,
                                                   youngFiber);
      }
    }
  }
}

LinearEquilibrium::LinearEquilibrium(
    std::vector<double> &paraValue,
    std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() < 3 || paraValue.size() > 5) {
    std::cerr
        << "FiberModel::LinearEquilibrium::LinearEquilibrium() Uses at least 3 "
           "and up to 5 parameters:"
        << std::endl
        << "k_rate, (rate of material update (towards target)" << std::endl
        << "young_matrix, young_fiber, (material limits, should be same as for "
           "the mechanical "
        << "(TRBS) model used)" << std::endl
        << "[initialization_flag] (0 - no initiation; 1 - isotropic; 2 - "
           "anisotropic"
        << std::endl
        << "[equilibrium_threshold] (only update if vertex moving slower), "
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() != 2 && indValue.size() != 3) ||
      indValue[0].size() != 1 ||
      (indValue[1].size() != 1 && indValue[1].size() != 2) ||
      (indValue.size() == 3 && indValue[2].size() != 1)) {
    std::cerr << "FiberModel::LinearEquilibrium::LinearEquilibrium() Error in "
                 "indices given:"
              << std::endl
              << "First level gives stress/strain anisotropy index to read."
              << std::endl
              << "Second level gives Young_Longitudinal index to update, or "
                 "Y_L and Y_T indices."
              << std::endl
              << "Third optional level gives index for velocity to be compared "
                 "with p_4."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() == 3 && paraValue.size() < 5) ||
      (indValue.size() < 3 && paraValue.size() == 5)) {
    std::cerr << "FiberModel::LinearEquilibrium::LinearEquilibrium() "
              << "5th optional [eq_threshold] parameter must be comparable "
              << "with variable index variable (optional) given at level 3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::LinearEquilibrium");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "Y_matrix";
  tmp[2] = "Y_fiber";
  if (numParameter() > 3)
    tmp[3] = "init_flag";
  if (numParameter() > 4)
    tmp[4] = "velocity_threshold";

  setParameterId(tmp);
}

void LinearEquilibrium::initiate(Tissue &T, DataMatrix &cellData,
                                 DataMatrix &wallData, DataMatrix &vertexData,
                                 DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                                 DataMatrix &vertexDerivs) {

  if (numParameter() < 4 || parameter(3) == 0)
    return;

  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(1);
  double youngFiber = parameter(2);

  if (parameter(3) == 1) { // initialize with iso material
    // initiating with 0 anisotropy and isotropic material
    double youngIso = FiberModel::isotropicFunction(youngMatrix, youngFiber);
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      cellData[cellIndex][AnisoIndex] = 0.0;
      cellData[cellIndex][YoungLIndex] = youngIso;
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] = youngIso;
      }
    }
  } else if (parameter(3) == 2) { // initiating with current anisotropy

    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      double anisotropy = cellData[cellIndex][AnisoIndex];
      cellData[cellIndex][YoungLIndex] =
          FiberModel::linearFunction(anisotropy, youngMatrix, youngFiber);
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

void LinearEquilibrium::derivs(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                               DataMatrix &vertexDerivs) {} // do nothing

void LinearEquilibrium::update(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               double h) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double kRate = parameter(0);
  double youngMatrix = parameter(1);
  double youngFiber = parameter(2);

  for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
    double anisotropy = cellData[cellIndex][AnisoIndex];

    // only update cell if v<v_th or v_th not given
    if (numParameter() < 5 || numVariableIndexLevel() < 3 ||
        cellData[cellIndex][variableIndex(2, 0)] < parameter(5)) {

      double youngLTarget =
          FiberModel::linearFunction(anisotropy, youngMatrix, youngFiber);

      if (youngLTarget >
          youngMatrix + youngFiber) // not to exceed maximum when using absolute
                                    // stress anisotropy
        youngLTarget = youngMatrix + youngFiber;

      double update =
          kRate * h * (youngLTarget - cellData[cellIndex][YoungLIndex]);
      cellData[cellIndex][YoungLIndex] += update;

      // if maximal value is exceded, the rate parameter is too large
      if (cellData[cellIndex][YoungLIndex] > youngMatrix + youngFiber) {
        std::cerr << "FiberModel::LinearEquilibrium::update ERROR: "
                  << "Update rate too large - Youngs modulus exceded maximal "
                     "value in cell "
                  << cellIndex << "(" << cellData[cellIndex][YoungLIndex]
                  << " > " << youngMatrix << " + " << youngFiber << "."
                  << std::endl;
      }
      if (numVariableIndex(1) > 1) { // Also save Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

Hill::Hill(std::vector<double> &paraValue,
           std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() < 4 || paraValue.size() > 6) {
    std::cerr
        << "FiberModel::Hill::Hill() Uses at least 4 and up to 6 parameters:"
        << std::endl
        << "young_matrix, young_fiber, (material limits, should be same as for "
           "the mechanical "
        << "(TRBS) model used." << std::endl
        << "K_hill, n_hill, (parameters for feedback function, should be same "
           "as for the mechanical "
        << "(TRBS) model used." << std::endl
        << "[initialization_flag] (0 - no initiation; 1 - isotropic; 2 - "
           "anisotropic"
        << std::endl
        << "[equilibrium_threshold] (only update if vertex moving slower), "
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() != 2 && indValue.size() != 3) ||
      indValue[0].size() != 1 ||
      (indValue[1].size() != 1 && indValue[1].size() != 2) ||
      (indValue.size() == 3 && indValue[2].size() != 1)) {
    std::cerr << "FiberModel::Hill::Hill() Error in indices given:" << std::endl
              << "First level gives stress/strain anisotropy index to read."
              << std::endl
              << "Second level gives Young_Longitudinal index to update, or "
                 "Y_L and Y_T indices."
              << std::endl
              << "Third optional level gives index for velocity to be compared "
                 "with p_5."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() == 3 && paraValue.size() < 6) ||
      (indValue.size() < 3 && paraValue.size() == 6)) {
    std::cerr << "FiberModel::Hill::Hill() 6th optional [eq_threshold] "
                 "parameter must be comparable "
              << "with variable index variable (optional) given at level 3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::Hill");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "Y_matrix";
  tmp[1] = "Y_fiber";
  tmp[2] = "K_hill";
  tmp[3] = "n_hill";
  if (numParameter() > 4)
    tmp[4] = "init_flag";
  if (numParameter() > 5)
    tmp[5] = "velocity_threshold";

  setParameterId(tmp);
}

void Hill::initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                    DataMatrix &vertexData, DataMatrix &cellDerivs,
                    DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {

  if (numParameter() < 5 || parameter(4) == 0)
    return;

  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(0);
  double youngFiber = parameter(1);
  double KHill = parameter(2);
  double nHill = parameter(3);

  if (parameter(4) == 1) { // initialize with iso material
    // initiating with 0 anisotropy and isotropic material
    double youngIso = FiberModel::isotropicFunction(youngMatrix, youngFiber);
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      cellData[cellIndex][AnisoIndex] = 0.0;
      cellData[cellIndex][YoungLIndex] = youngIso;
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] = youngIso;
      }
    }
  } else if (parameter(4) == 2) { // initiating with current anisotropy
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      double anisotropy = cellData[cellIndex][AnisoIndex];
      cellData[cellIndex][YoungLIndex] = FiberModel::hillFunction(
          anisotropy, youngMatrix, youngFiber, KHill, nHill);
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

void Hill::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, DataMatrix &cellDerivs,
                  DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
} // do nothing

void Hill::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                  DataMatrix &vertexData, double h) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(0);
  double youngFiber = parameter(1);
  double KHill = parameter(2);
  double nHill = parameter(3);

  for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
    double anisotropy = cellData[cellIndex][AnisoIndex];

    // only update cell if v<v_th or v_th not given
    if (numParameter() < 6 || numVariableIndexLevel() < 3 ||
        cellData[cellIndex][variableIndex(2, 0)] < parameter(6)) {

      double youngL = FiberModel::hillFunction(anisotropy, youngMatrix,
                                               youngFiber, KHill, nHill);

      if (youngL > youngMatrix + youngFiber) // not to exceed maximum when using
                                             // absolute stress anisotropy
        youngL = youngMatrix + youngFiber;

      cellData[cellIndex][YoungLIndex] = youngL;
      if (numVariableIndex(1) > 1) { // Also save Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(youngL, youngMatrix,
                                                   youngFiber);
      }
    }
  }
}

HillEquilibrium::HillEquilibrium(std::vector<double> &paraValue,
                                 std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() < 5 || paraValue.size() > 7) {
    std::cerr
        << "FiberModel::HillEquilibrium::HillEquilibrium() Uses at least 5 and "
           "up to 7 parameters:"
        << std::endl
        << "k_rate, (rate of update towards target" << std::endl
        << "young_matrix, young_fiber, (material limits, should be same as for "
           "the mechanical "
        << "(TRBS) model used." << std::endl
        << "K_hill, n_hill, (parameters for feedback function, should be same "
           "as for the mechanical "
        << "(TRBS) model used." << std::endl
        << "[initialization_flag] (0 - no initiation; 1 - isotropic; 2 - "
           "anisotropic"
        << std::endl
        << "[equilibrium_threshold] (only update if vertex moving slower), "
        << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() != 2 && indValue.size() != 3) ||
      indValue[0].size() != 1 ||
      (indValue[1].size() != 1 && indValue[1].size() != 2) ||
      (indValue.size() == 3 && indValue[2].size() != 1)) {
    std::cerr << "FiberModel::HillEquilibrium::HillEquilibrium() Error in "
                 "indices given:"
              << std::endl
              << "First level gives stress/strain anisotropy index to read."
              << std::endl
              << "Second level gives Young_Longitudinal index to update, or "
                 "Y_L and Y_T indices."
              << std::endl
              << "Third optional level gives index for velocity to be compared "
                 "with p_6."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  if ((indValue.size() == 3 && paraValue.size() < 7) ||
      (indValue.size() < 3 && paraValue.size() == 7)) {
    std::cerr << "FiberModel::HillEquilibrium::HillEquilibrium()"
              << "7th optional [eq_threshold] parameter must be comparable "
              << "with variable index variable (optional) given at level 3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::HillEquilibrium");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "Y_matrix";
  tmp[2] = "Y_fiber";
  tmp[3] = "K_hill";
  tmp[4] = "n_hill";
  if (numParameter() > 5)
    tmp[5] = "init_flag";
  if (numParameter() > 6)
    tmp[6] = "velocity_threshold";

  setParameterId(tmp);
}

void HillEquilibrium::initiate(Tissue &T, DataMatrix &cellData,
                               DataMatrix &wallData, DataMatrix &vertexData,
                               DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                               DataMatrix &vertexDerivs) {

  if (numParameter() < 6 || parameter(5) == 0)
    return;

  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double youngMatrix = parameter(1);
  double youngFiber = parameter(2);
  double KHill = parameter(3);
  double nHill = parameter(4);

  if (parameter(5) == 1) { // initialize with iso material
    // initiating with 0 anisotropy and isotropic material
    double youngIso = FiberModel::isotropicFunction(youngMatrix, youngFiber);
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      cellData[cellIndex][AnisoIndex] = 0.0;
      cellData[cellIndex][YoungLIndex] = youngIso;
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] = youngIso;
      }
    }
  } else if (parameter(5) == 2) { // initiating with current anisotropy
    for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {
      double anisotropy = cellData[cellIndex][AnisoIndex];
      cellData[cellIndex][YoungLIndex] = FiberModel::hillFunction(
          anisotropy, youngMatrix, youngFiber, KHill, nHill);
      if (numVariableIndex(1) > 1) { // Also initiate Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

void HillEquilibrium::derivs(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             DataMatrix &cellDerivs, DataMatrix &wallDerivs,
                             DataMatrix &vertexDerivs) {} // do nothing

void HillEquilibrium::update(Tissue &T, DataMatrix &cellData,
                             DataMatrix &wallData, DataMatrix &vertexData,
                             double h) {
  size_t numCell = cellData.size();
  size_t AnisoIndex = variableIndex(0, 0);
  size_t YoungLIndex = variableIndex(1, 0);
  double kRate = parameter(0);
  double youngMatrix = parameter(1);
  double youngFiber = parameter(2);
  double KHill = parameter(3);
  double nHill = parameter(4);

  for (size_t cellIndex = 0; cellIndex < numCell; ++cellIndex) {

    // only update cell if v<v_th or v_th not given
    if (numParameter() < 7 || numVariableIndexLevel() < 3 ||
        cellData[cellIndex][variableIndex(2, 0)] < parameter(7)) {

      double anisotropy = cellData[cellIndex][AnisoIndex];
      double youngLTarget = FiberModel::hillFunction(anisotropy, youngMatrix,
                                                     youngFiber, KHill, nHill);

      if (youngLTarget >
          youngMatrix + youngFiber) // not to exceed maximum when using absolute
                                    // stress anisotropy
        youngLTarget = youngMatrix + youngFiber;

      double update =
          kRate * h * (youngLTarget - cellData[cellIndex][YoungLIndex]);
      cellData[cellIndex][YoungLIndex] += update;

      // if maximal value is exceded, the rate parameter is too large
      if (cellData[cellIndex][YoungLIndex] > youngMatrix + youngFiber) {
        std::cerr << "FiberModel::HillEquilibrium::update ERROR: "
                  << "Update rate too large - Youngs modulus exceded maximal "
                     "value in cell "
                  << cellIndex << "(" << cellData[cellIndex][YoungLIndex]
                  << " > " << youngMatrix << " + " << youngFiber << "."
                  << std::endl;
      }
      if (numVariableIndex(1) > 1) { // Also save Y_T
        cellData[cellIndex][variableIndex(1, 1)] =
            FiberModel::transverseFromLongitudinal(
                cellData[cellIndex][YoungLIndex], youngMatrix, youngFiber);
      }
    }
  }
}

Deposition::Deposition(std::vector<double> &paraValue,
                       std::vector<std::vector<size_t> > &indValue) {
  if (paraValue.size() != 7) {
    std::cerr << "FiberModel::Deposition::Deposition() "
              << "Uses seven parameters: k_rate, initial randomness(percent),"
              << "velocity threshold, initiation flag(0/1), k_hill and n_hill"
              << "fiber strength" << std::endl;
    exit(EXIT_FAILURE);
  }

  if (indValue.size() != 1 || indValue[0].size() != 6) {
    std::cerr << "FiberModel::Deposition::Deposition() "
              << "First index: anisotropy."
              << "second index:  misses stress."
              << "third index: cell area."
              << "fourth index: fiber."
              << "fifth index: longitudinal young fiber"
              << "fourth index: velocity" << std::endl;
    exit(EXIT_FAILURE);
  }

  setId("FiberModel::Deposition");
  setParameter(paraValue);
  setVariableIndex(indValue);

  std::vector<std::string> tmp(numParameter());
  tmp[0] = "k_rate";
  tmp[1] = "randomness";
  tmp[2] = "velocity_threshold";
  tmp[3] = "init_flag";
  tmp[4] = "k_Hill";
  tmp[5] = "n_Hill";
  tmp[6] = "fiber";

  setParameterId(tmp);
}

void Deposition::initiate(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                          DataMatrix &vertexData, DataMatrix &cellDerivs,
                          DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {
  size_t numCells = cellData.size();
  size_t YoungFiberIndex = variableIndex(0, 3);
  size_t YoungFiberLIndex = variableIndex(0, 4);

  double randCoef = parameter(1);
  double youngFiber = parameter(6);
  size_t initFlag = parameter(3);

  if (initFlag == 1) // initialize with uniform material
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      cellData[cellIndex][YoungFiberIndex] = youngFiber;
      cellData[cellIndex][YoungFiberLIndex] = youngFiber / 2;
    }
  if (initFlag == 2 ||
      initFlag == 3) { // initialize with random material stiffness
    srand((unsigned)time(NULL));
    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      cellData[cellIndex][YoungFiberIndex] =
          youngFiber *
          (1 + randCoef * 2 * (0.5 - ((double)rand() / (RAND_MAX))));
      cellData[cellIndex][YoungFiberLIndex] =
          cellData[cellIndex][YoungFiberIndex] / 2;

      std::cerr << cellData[cellIndex][YoungFiberIndex] << "  "
                << cellData[cellIndex][YoungFiberLIndex] << std::endl;
    }
  }
  return;
}

void Deposition::derivs(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                        DataMatrix &vertexData, DataMatrix &cellDerivs,
                        DataMatrix &wallDerivs, DataMatrix &vertexDerivs) {}

void Deposition::update(Tissue &T, DataMatrix &cellData, DataMatrix &wallData,
                        DataMatrix &vertexData, double h) {
  size_t numCells = cellData.size();

  size_t AnisoIndex = variableIndex(0, 0);
  size_t MstressIndex = variableIndex(0, 1);
  size_t areaIndex = variableIndex(0, 2);
  size_t YoungFiberIndex = variableIndex(0, 3);
  size_t YoungFiberLIndex = variableIndex(0, 4);
  size_t velocityIndex = variableIndex(0, 5);

  double k_rate = parameter(0);
  // double randCoef=parameter(1);
  double vThresh = parameter(2);
  double Kh = parameter(4);
  double Nh = parameter(5);
  double YoungFiber = parameter(6);
  if (parameter(0) == 0.0)
    return;

  double smax = 14; // 28;
  double smin = 8;  //;

  bool equil = true;
  for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex)
    if (cellData[cellIndex][velocityIndex] > vThresh)
      equil = false;

  if (equil) {
    double totalStressArea = 0;
    double totalArea = 0;

    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      if (cellData[cellIndex][MstressIndex] < smax &&
          cellData[cellIndex][MstressIndex] > smin) {
        totalStressArea += (cellData[cellIndex][MstressIndex] - smin) *
                           cellData[cellIndex][areaIndex];
      } else if (cellData[cellIndex][MstressIndex] > smax) {
        totalStressArea += (smax - smin) * cellData[cellIndex][areaIndex];
      }
      totalArea += cellData[cellIndex][areaIndex];
    }

    for (size_t cellIndex = 0; cellIndex < numCells; ++cellIndex) {
      double targetFiber = 0;
      if (cellData[cellIndex][MstressIndex] < smax &&
          cellData[cellIndex][MstressIndex] > smin) {
        targetFiber = (YoungFiber * totalArea *
                       (cellData[cellIndex][MstressIndex] - smin)) /
                      totalStressArea;
      } else if (cellData[cellIndex][MstressIndex] > smax) {
        targetFiber =
            (YoungFiber * totalArea * (smax - smin)) / totalStressArea;
      } else if (cellData[cellIndex][MstressIndex] < smin) {
        targetFiber = 0;
      }
      cellData[cellIndex][YoungFiberIndex] +=
          k_rate * (targetFiber - cellData[cellIndex][YoungFiberIndex]);

      if (parameter(3) == 3) {
        double anisotropy = cellData[cellIndex][AnisoIndex];
        cellData[cellIndex][YoungFiberLIndex] +=
            k_rate *
            ((0.5 *
              (1 + (std::pow(anisotropy, Nh) /
                    (std::pow((1 - anisotropy), Nh) * std::pow(Kh, Nh) +
                     std::pow(anisotropy, Nh)))) *
              cellData[cellIndex][YoungFiberIndex]) -
             cellData[cellIndex][YoungFiberLIndex]);
      } else
        cellData[cellIndex][YoungFiberLIndex] =
            cellData[cellIndex][YoungFiberIndex] / 2;
    }
    // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
    //   totalStressArea+=cellData[cellIndex][MstressIndex]
    //     *cellData[cellIndex][areaIndex];
    //   totalArea+=cellData[cellIndex][areaIndex];
    // }
    // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
    //   cellData[cellIndex][YoungFiberIndex]+=k_rate*
    //     ((YoungFiber*totalArea*cellData[cellIndex][MstressIndex])
    //      /totalStressArea
    //      -cellData[cellIndex][YoungFiberIndex]);
    //   cellData[cellIndex][YoungFiberLIndex]=cellData[cellIndex][YoungFiberIndex]/2;
    //   // double anisotropy= cellData[cellIndex][AnisoIndex];
    //   // cellData[cellIndex][YoungFiberLIndex]+=k_rate*
    //   //   ((0.5*(1+(std::pow(anisotropy,Nh)
    //   //           /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
    //   //             +std::pow(anisotropy,Nh))))*
    //   cellData[cellIndex][YoungFiberIndex])
    //   //   -cellData[cellIndex][YoungFiberLIndex]);
    // }
    // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
    //   totalStressArea+=cellData[cellIndex][MstressIndex]*cellData[cellIndex][areaIndex];
    //   totalArea+=cellData[cellIndex][areaIndex];
    // }
    // for (size_t cellIndex=0; cellIndex<numCells; ++cellIndex){
    //   cellData[cellIndex][YoungFiberIndex]+=k_rate*
    //     ((YoungFiber*totalArea*cellData[cellIndex][MstressIndex])/totalStressArea
    //      -cellData[cellIndex][YoungFiberIndex]);
    //   cellData[cellIndex][YoungFiberLIndex]=cellData[cellIndex][YoungFiberIndex]/2;
    //   // double anisotropy= cellData[cellIndex][AnisoIndex];
    //   // cellData[cellIndex][YoungFiberLIndex]+=k_rate*
    //   //   ((0.5*(1+(std::pow(anisotropy,Nh)
    //   //           /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
    //   //             +std::pow(anisotropy,Nh))))*
    //   cellData[cellIndex][YoungFiberIndex])
    //   //   -cellData[cellIndex][YoungFiberLIndex]);
    // }
  }
  // double tmp;
  // //std::cerr<<"total    "<<totalStressArea<<std::endl;

  // if(totalStressArea !=0){
  //   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex){
  //     tmp=h*k_rate*cellData[cellIndex][YoungFiberIndex];
  //     cellData[cellIndex][YoungFiberIndex]-=tmp;
  //     totalFiber+=tmp*cellData[cellIndex][areaIndex];
  //   }

  //   for (size_t cellIndex=0; cellIndex<numCell; ++cellIndex) {
  //     double anisotropy= cellData[cellIndex][AnisoIndex];
  //     cellData[cellIndex][YoungFiberIndex]+=(totalFiber*cellData[cellIndex][MstressIndex])/totalStressArea;
  //     cellData[cellIndex][YoungFiberLIndex] =
  //       0.5*(1+(std::pow(anisotropy,Nh)
  //               /(std::pow((1-anisotropy),Nh)*std::pow(Kh,Nh)
  //                 +std::pow(anisotropy,Nh))))*
  //                 cellData[cellIndex][YoungFiberIndex];
  //   }
  // }
}
} // namespace FiberModel
