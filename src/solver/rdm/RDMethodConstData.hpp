#ifndef PDEKIT_RD_Method_Const_Data_hpp
#define PDEKIT_RD_Method_Const_Data_hpp

namespace pdekit
{

namespace solver
{

namespace rdm
{

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
class RDMethodConstData
{
  public:
  /// Default constructor
  RDMethodConstData();

  /// Copy constructor
  RDMethodConstData(const RDMethodConstData &other_data);

  /// Destructor
  ~RDMethodConstData();

  /// Assignment operator
  RDMethodConstData &operator=(const RDMethodConstData &other_data);

  /// Metric data (geometric factors - jacobians, normals etc.) for element
  /// geometry
  CellGeoMetType CGM;

  /// Cell solution metric - values of u and derivatives o u in quadrature
  /// points
  CellSolMetType CSM;

  /// Flux metric data in one element/on one face
  CellFluxMetType CFM;

  /// Metric data for source term
  CellSolMetType CSrcM;

  private:
};

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::RDMethodConstData()
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::RDMethodConstData(
    const RDMethodConstData &other_data)
    : CGM(other_data.CGM), CSM(other_data.CSM), CFM(other_data.CFM), CSrcM(other_data.CSrcM)
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::~RDMethodConstData()
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>
    &RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::operator=(
        const RDMethodConstData &other_data)
{
  CGM   = other_data.CGM;
  CSM   = other_data.CSM;
  CFM   = other_data.CFM;
  CSrcM = other_data.CSrcM;
  return *this;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
