# Regression tests: run every tutorial model with the v2 simulator.
# Numerical comparison against the legacy simulator is done manually (see
# README.md); these tests assert the runs complete successfully.

set(TUTORIALS ${CMAKE_CURRENT_SOURCE_DIR}/../examples/tutorials)

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
# oneSquare never reaches random.model's volume threshold of 40, so this pair
# runs the whole 2000 time units and ends with the one cell it started with:
# it covers growth and the mechanics, but not a single division. Kept as that
# smoke test, with the division itself covered separately below.
tissue_tutorial_test(tutorial_growth_division_random
  growth_and_division/random.model growth_and_division/oneSquare.init
  growth_and_division/solver.rk5)

# Division actually happening. The meristem starts with 267 cells large enough
# to divide, so this reaches 268 within two time units - enough to exercise the
# topology surgery, and fast because it stops right after.
function(tissue_division_test name model init)
  add_test(NAME ${name}
    COMMAND simulator ${TUTORIALS}/${model} ${TUTORIALS}/${init}
            ${CMAKE_CURRENT_SOURCE_DIR}/division.rk5)
  set_tests_properties(${name} PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()

tissue_division_test(division_random_direction
  growth_and_division/random.model growth_and_division/meristem.init)
tissue_tutorial_test(tutorial_mesh_growth
  mesh/growth.model mesh/twoSquare.init mesh/solver.rk5)
tissue_tutorial_test(tutorial_mesh_growth_ct
  mesh/growthCT.model mesh/twoSquare3D.init mesh/solverCT.rk5)
tissue_tutorial_test(tutorial_mesh_growth_sv
  mesh/growthSV.model mesh/twoSquareSV.init mesh/solver.rk5)
