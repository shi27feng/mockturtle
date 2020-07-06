/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "experiments.hpp"

#include <mockturtle/networks/aig.hpp>
#include <mockturtle/algorithms/multithreaded_rewriting.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <lorina/aiger.hpp>
#include <fmt/format.h>

#include <string>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint64_t, uint64_t, double, bool> exp( "multithreaded_rewriting", "benchmark", "size_before", "size_after", "runtime", "equivalent" );
  for ( auto const& benchmark : epfl_benchmarks() )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[e] could not parse benchmark {}\n", benchmark );
      continue;
    }

    uint64_t const size_before = aig.num_gates();

    multithreaded_cut_enumeration_params ps;
    multithreaded_cut_enumeration_stats st;
    multithreaded_cut_enumeration( aig, ps, &st );
    aig = cleanup_dangling( aig );

    uint64_t const size_after = aig.num_gates();

    auto const cec = abc_cec( aig, benchmark );
    exp( benchmark, size_before, size_after, to_seconds( st.time_total ), cec );
  }

  exp.save();
  exp.table();
}
