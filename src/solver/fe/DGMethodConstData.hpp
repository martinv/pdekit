#ifndef PDEKIT_DG_Method_Const_Data_hpp
#define PDEKIT_DG_Method_Const_Data_hpp

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
class DGMethodConstData
{
  public:
  /// Default constructor
  DGMethodConstData();

  /// Copy constructor
  DGMethodConstData(const DGMethodConstData &other_data);

  /// Destructor
  ~DGMethodConstData();

  /// Assignment operator
  DGMethodConstData &operator=(const DGMethodConstData &other_data);

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
DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::DGMethodConstData()
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::DGMethodConstData(
    const DGMethodConstData &other_data)
    : CGM(other_data.CGM), CSM(other_data.CSM), CFM(other_data.CFM), CSrcM(other_data.CSrcM)
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::~DGMethodConstData()
{
}

// ----------------------------------------------------------------------------

template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>
    &DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType>::operator=(
        const DGMethodConstData &other_data)
{
  CGM   = other_data.CGM;
  CSM   = other_data.CSM;
  CFM   = other_data.CFM;
  CSrcM = other_data.CSrcM;
  return *this;
}

// ----------------------------------------------------------------------------

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
