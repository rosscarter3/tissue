//
// Tissue v2 simulator: command-line compatible with the legacy simulator.
//   simulator modelFile initFile simulatorParaFile [flags]
// Flags: -init_output file, -init_output_format format, -verbose flag,
//        -vtk_output dir, -centerTri_init, -debug_output file (accepted, not
//        yet used), -help.
//
#include <csignal>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>

#include "tissue/core/tissue.h"
#include "tissue/io/config.h"
#include "tissue/solvers/base_solver.h"

namespace {

tissue::BaseSolver *signalSolver = nullptr;

void solverSignalHandler(int signal) {
  std::cerr << "simulator: received signal " << signal << std::endl;
  if (signalSolver) {
    std::string fileName = tissue::config::getValue("init_output", 0);
    if (!fileName.empty()) {
      std::ofstream out(fileName);
      if (out) {
        signalSolver->printInit(out);
        std::cerr << "simulator: state written to " << fileName << std::endl;
      }
    }
  }
  std::exit(EXIT_FAILURE);
}

} // namespace

int main(int argc, char *argv[]) {
  using namespace tissue;

  config::registerOption("init_output", 1);
  config::registerOption("init_output_format", 1);
  config::registerOption("help", 0);
  config::registerOption("centerTri_init", 0);
  config::registerOption("verbose", 1);
  config::registerOption("debug_output", 1);
  config::registerOption("vtk_output", 1);

  int verboseFlag = 1;
  std::string configFile;
  if (const char *home = std::getenv("HOME"))
    configFile = std::string(home) + "/.tissue";
  config::init(argc, argv, configFile);

  std::string verboseString = config::getValue("verbose", 0);
  if (!verboseString.empty())
    verboseFlag = std::atoi(verboseString.c_str());

  if (config::getBooleanValue("help")) {
    std::cerr << std::endl
              << "Usage: " << argv[0] << " modelFile initFile simulatorParaFile."
              << std::endl
              << std::endl;
    std::cerr << "Possible additional flags are:" << std::endl;
    std::cerr << "-centerTri_init - Init file is assumed to have cell "
                 "variables storing a central triangulation." << std::endl;
    std::cerr << "-init_output file - Set filename for output of final state "
                 "in init file format." << std::endl;
    std::cerr << "-init_output_format format - Sets format for output of "
                 "final state (tissue supported in v2)." << std::endl;
    std::cerr << "-verbose flag - Set flag for verbose (1) or silent (0) "
                 "output mode to stderr." << std::endl;
    std::cerr << "-vtk_output dir - Directory name for vtk-format output. "
                 "Default is \"vtk\"." << std::endl;
    std::cerr << "-help - Shows this message." << std::endl;
    std::exit(EXIT_FAILURE);
  } else if (config::argc() != 4) {
    std::cerr << "Type '" << argv[0] << " -help' for usage." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string modelFile = config::argv(1);
  std::string initFile = config::argv(2);
  std::string simPara = config::argv(3);

  try {
    Tissue T;
    if (verboseFlag)
      std::cerr << "Reading model file " << modelFile << std::endl;
    T.readModel(modelFile, verboseFlag);
    if (verboseFlag)
      std::cerr << "Reading init file " << initFile << std::endl;
    if (!config::getBooleanValue("centerTri_init")) {
      T.readInit(initFile, verboseFlag);
    } else {
      std::cerr << "Assuming init file format with central triangulation "
                   "stored in cell variables." << std::endl;
      T.readInitCenterTri(initFile, verboseFlag);
    }

    if (verboseFlag)
      std::cerr << "Generating solver from file " << simPara << std::endl;
    std::unique_ptr<BaseSolver> solver = BaseSolver::getSolver(&T, simPara);
    signalSolver = solver.get();
    std::signal(SIGINT, solverSignalHandler);
    std::signal(SIGHUP, solverSignalHandler);
    std::signal(SIGTERM, solverSignalHandler);

    if (verboseFlag)
      std::cerr << "Initiating solver from tissue." << std::endl;
    solver->getInit();
    std::cerr << "Start simulation." << std::endl;
    solver->simulate();

    std::string fileName = config::getValue("init_output", 0);
    if (!fileName.empty()) {
      std::ofstream out(fileName);
      if (!out) {
        std::cerr << "Warning: main() - Cannot open file for init output ("
                  << fileName << ")" << std::endl;
      } else {
        std::cerr << "Setting tissue variables from simulator data."
                  << std::endl;
        solver->setTissueVariables(T.numCellVariable());
        std::string initFormat = config::getValue("init_output_format", 0);
        if (initFormat.empty() || initFormat == "tissue") {
          std::cerr << "Printing init in file " << fileName
                    << " using tissue format." << std::endl;
          solver->printInit(out);
        } else {
          std::cerr << "Warning: main() - init output format '" << initFormat
                    << "' is not ported to tissue v2 yet. No init file "
                       "written." << std::endl;
        }
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "simulator: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
