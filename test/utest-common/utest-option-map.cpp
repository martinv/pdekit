/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE options_test
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iostream>

#include "common/OptionMap.hpp"
#include "common/PDEKit.hpp"

using namespace pdekit;
using namespace pdekit::common;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(option_map_insert_modify_option)
{
  common::OptionMap opt_map;
  opt_map.create<Uint>("age");
  opt_map.create<Uint>("age2", 35);

  const Uint age2_a = opt_map.get<Uint>("age2");

  opt_map.set<Uint>("age2", 36);
  const Uint age2_b = opt_map.get<Uint>("age2");

  BOOST_CHECK_EQUAL(age2_a, 35);
  BOOST_CHECK_EQUAL(age2_b, 36);

  opt_map.create<Real>("length");
  opt_map.set<Real>("length", 1.6);

  const Float lf = 0.8;
  opt_map.create<Float>("length");
  opt_map.set<Float>("length", lf);

  const Real length_real   = opt_map.get<Real>("length");
  const Float length_float = opt_map.get<Float>("length");

  BOOST_CHECK_CLOSE(length_real, 1.6, 1.e-10);
  BOOST_CHECK_CLOSE(length_float, 0.8, 1.e-5);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(option_map_initialize_from_file)
{
  std::ofstream outfile;
  outfile.open("options.dat");

  /*
  outfile << "Param1 Real 0.1" << std::endl;
  outfile << "Param2 Float 0.2" << std::endl;
  outfile << "Param3 Real 0.3" << std::endl;
  outfile << "Param4 Int 1" << std::endl;
  outfile << "Param5 Uint 2" << std::endl;
  outfile << "Param6 String Milie" << std::endl;
  outfile << "Param7 StringTuple2 Matyas Vymazal" << std::endl;
  outfile << "Param8 StringTuple3 2000 Space Odyssey" << std::endl;
  */

  outfile << "Param1 Real" << std::endl;
  outfile << "Param2 Float" << std::endl;
  outfile << "Param3 Real" << std::endl;
  outfile << "Param4 Int" << std::endl;
  outfile << "Param5 Uint" << std::endl;
  outfile << "Param6 String" << std::endl;
  outfile << "Param7 StringTuple2" << std::endl;
  outfile << "Param8 StringTuple3" << std::endl;

  outfile.close();

  OptionMap opt_map;
  opt_map.initialize_from_file("options.dat");
  opt_map.print();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(option_map_configure_from_file)
{
  std::ofstream outfile;
  outfile.open("options_init.dat");

  outfile << "Param1 Real" << std::endl;
  outfile << "Param2 Float" << std::endl;
  outfile << "Param3 Real" << std::endl;
  outfile << "Param4 Int" << std::endl;
  outfile << "Param5 Uint" << std::endl;
  outfile << "Param6 String" << std::endl;
  outfile << "Param7 StringTuple2" << std::endl;
  outfile << "Param8 StringTuple3" << std::endl;

  outfile.close();

  outfile.open("option_values.dat");

  outfile << "Param1 Real 0.1" << std::endl;
  outfile << "Param2 Float 0.2" << std::endl;
  outfile << "Param3 Real 0.3" << std::endl;
  outfile << "Param4 Int 1" << std::endl;
  outfile << "Param5 Uint 2" << std::endl;
  outfile << "Param6 String Milie" << std::endl;
  outfile << "Param7 StringTuple2 Matyas Vymazal" << std::endl;
  outfile << "Param8 StringTuple3 2000 Space Odyssey" << std::endl;

  OptionMap opt_map;
  opt_map.initialize_from_file("options_init.dat");
  // opt_map.print();
  opt_map.configure_from_file("option_values.dat");
  // opt_map.print();
}

// ----------------------------------------------------------------------------
