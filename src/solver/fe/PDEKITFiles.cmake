list(APPEND PDEKIT_Solver_FE_HEADERS
  HelmholtzSolverCG.hpp
  HelmholtzSolverHDG.hpp
  HelmholtzSolverCGHDG.hpp
  InteriorSolverCGHDG.hpp
  DGSolver.hpp
  DGHyperbolicSolver.hpp
  DGEllipticSolver.hpp
  DGMethodConstData.hpp
  DGMethodScratchData.hpp
  DGTimeUpdate.hpp
  assembly/DGExplicitCellWorker.hpp
  assembly/DGExplicitFacetWorker.hpp
  assembly/HDGCellWorker.hpp
  assembly/HDGTraceWorker.hpp
  cell_terms/DGCellTerm.hpp
  facet_terms/DGFacetTerm.hpp
)

list(APPEND PDEKIT_Solver_FE_SOURCES
  DGTimeUpdate.cpp
)

