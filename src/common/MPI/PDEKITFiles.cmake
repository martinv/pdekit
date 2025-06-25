#CmakeLists.txt in InterpolationPointSet dir

list(APPEND PDEKIT_MPI_HEADERS
  MPI/CommMap.hpp
  MPI/MPITypes.hpp
  MPI/MPIEnv.hpp
  MPI/MPIRecv.hpp
  MPI/MPISend.hpp
)

list(APPEND PDEKIT_MPI_SOURCES
  MPI/CommMap.cpp
  MPI/MPIEnv.cpp
  MPI/MPITypes.cpp
  MPI/MPIRecv.cpp
  MPI/MPISend.cpp
)
