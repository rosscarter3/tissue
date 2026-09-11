# Regression tests: run every tutorial model with the v2 simulator.
# Numerical comparison against the legacy simulator is done manually (see
# v2/README.md); these tests assert the runs complete successfully.

set(TUTORIALS ${CMAKE_CURRENT_SOURCE_DIR}/../../examples/tutorials)

function(tissue_tutorial_test name model init solver)
  add_test(NAME ${name}
    COMMAND simulator ${TUTORIALS}/${model} ${TUTORIALS}/${init}
            ${TUTORIALS}/${solver})
  set_tests_properties(${name} PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

tissue_tutorial_test(tutorial_mass_action
  mass_action/massAction.model mass_action/twoSquare.init
  mass_action/solver.rk5)
tissue_tutorial_test(tutorial_growth_division_random
  growth_and_division/random.model growth_and_division/oneSquare.init
  growth_and_division/solver.rk5)
tissue_tutorial_test(tutorial_mesh_growth
  mesh/growth.model mesh/twoSquare.init mesh/solver.rk5)
tissue_tutorial_test(tutorial_mesh_growth_ct
  mesh/growthCT.model mesh/twoSquare3D.init mesh/solverCT.rk5)
tissue_tutorial_test(tutorial_mesh_growth_sv
  mesh/growthSV.model mesh/twoSquareSV.init mesh/solver.rk5)
