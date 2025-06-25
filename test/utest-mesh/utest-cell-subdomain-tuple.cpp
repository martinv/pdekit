/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_subdomain_tuple_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers

#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/local_topology/TraceEntityTuple.hpp"
#include "mesh/local_topology/TraceTupleFactory.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(line_facet_subdomain_tuple_utest)
{
  TraceTupleFactory::instance_type &inc_factory = TraceTupleFactory::instance();

  /*
  std::shared_ptr<internal::CellSubdomainTupleBase> inc_base =
      inc_factory.create(std::pair<Uint, Uint>(NO_SPLIT, NO_SPLIT));

  CellSubdomainTuple fpattern;

  fpattern.change_type(LINE_TO_TWO_SEGMENTS, NO_SPLIT);

  std::cout << fpattern.get().adapt_operations() << std::endl;

  std::cout << "Nb. elem on left facet side = " <<
  fpattern.get().nb_elem_on_facet_side(LEFT)
            << std::endl;
  std::cout << "Nb. elem on right facet side = " <<
  fpattern.get().nb_elem_on_facet_side(RIGHT)
            << std::endl;

  std::cout << fpattern.get() << std::endl;
  */
}

// ----------------------------------------------------------------------------
