list(APPEND PDEKIT_Mesh_Shape_Function_HEADERS
  shape_function/ModalExpansion.hpp
  shape_function/DubinerExpansionTriag.hpp
  shape_function/DubinerExpansionTetra.hpp
  shape_function/CarnevaliExpansionLine.hpp
  shape_function/CarnevaliExpansionTriag.hpp
  shape_function/C0ExpansionTriag.hpp
  shape_function/C0ExpansionTetra.hpp
  shape_function/ModalExpansionHexa.hpp
  shape_function/ModalExpansionLine.hpp
  shape_function/ModalExpansionQuad.hpp
  shape_function/ModalBasisTag.hpp
  shape_function/ModalBasisFactory.hpp
  shape_function/ShapeFunction.hpp
  shape_function/SFTag.hpp
)

list(APPEND PDEKIT_Mesh_Shape_Function_SOURCES
  shape_function/ModalExpansion.cpp
  shape_function/DubinerExpansionTriag.cpp
  shape_function/DubinerExpansionTetra.cpp
  shape_function/CarnevaliExpansionLine.cpp
  shape_function/CarnevaliExpansionTriag.cpp
  shape_function/ModalExpansionHexa.cpp
  shape_function/ModalExpansionLine.cpp
  shape_function/ModalExpansionQuad.cpp
  shape_function/ModalBasisTag.cpp
  shape_function/ModalBasisFactory.cpp
  shape_function/ShapeFunction.cpp
  shape_function/SFTag.cpp
)
