list(APPEND PDEKIT_Mesh_Adaptation_HEADERS
  adaptation/CellAdaptOpBase.hpp
  adaptation/CellAdaptOpFactory.hpp
  adaptation/CellAdaptOp.hpp
  adaptation/CellAdaptOpTag.hpp
  adaptation/CellAdaptOpTriag.hpp
  adaptation/CellAdaptOpQuad.hpp
  adaptation/DoNothingCellAdaptOp.hpp
  adaptation/GeometryAdapter.hpp
  adaptation/LocalInterpolator.hpp
  adaptation/MeshAdaptSequence.hpp
  adaptation/MeshAdaptStep.hpp
)

list(APPEND PDEKIT_Mesh_Adaptation_SOURCES
  adaptation/CellAdaptOpBase.cpp
  adaptation/CellAdaptOpFactory.cpp
  adaptation/CellAdaptOp.cpp
  adaptation/CellAdaptOpTag.cpp
  adaptation/CellAdaptOpTriag.cpp
  adaptation/CellAdaptOpQuad.cpp
  adaptation/DoNothingCellAdaptOp.cpp
  adaptation/LocalInterpolator.cpp
  adaptation/MeshAdaptStep.cpp
)
