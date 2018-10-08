list(APPEND PDEKIT_Solver_HEADERS
 SolverBase.hpp
 ElementalMatrixOperator.hpp
 SolverSetupAlgorithm.hpp
 ExplicitSolver.hpp
 ImplicitSolver.hpp
 FEMetric.hpp
 FEMetricData.hpp
 InitialCondition.hpp
 NumFlux.hpp
 NumFluxAUSM.hpp
 NumFluxLaxFriedrichs.hpp
 Postprocessing.hpp
 SolverIO.hpp
 art_visc/ArtificialViscosity.hpp
)

list(APPEND PDEKIT_Solver_SOURCES
 SolverBase.cpp
 ExplicitSolver.cpp
 ImplicitSolver.cpp
 SolverIO.cpp
)



