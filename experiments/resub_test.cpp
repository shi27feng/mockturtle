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
#include <mockturtle/io/verilog_reader.hpp>
#include <lorina/aiger.hpp>
#include <fmt/format.h>

#include <string>

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  std::string verilog_file =
    "module top( x0 , x1 , x2 , y0 , y1 );\n"
    "  input x0 , x1 , x2 ;\n"
    "  output y0 , y1 ;\n"
    "   wire n4 , n5 ;\n"
    "   assign n4 = x0 & x1 ;\n"
    "   assign y0 = x2 ;\n"
    "   assign y1 = n4 ;\n"
    "endmodule\n";

  std::stringstream ss;
  ss << verilog_file;
  
  aig_network aig;
  if ( lorina::read_verilog( ss, verilog_reader( aig ) ) != lorina::return_code::success )
  {
    fmt::print( "[e] could not parse benchmark\n" );
    std::abort();
  }
  
  index_list indices;
  encode( indices, aig );
  
  aig_network window_aig;
  decode( window_aig, indices );
  write_verilog( window_aig, "win.v" );
  
  /* optimize index list using the resubstitution algorithm */
  auto const raw_array = indices.raw_data();
  abcresub::Abc_ResubPrepareManager( 1 );
  int num_resubs;
  int *new_indices_raw;
  uint64_t new_entries = abcresub::Abc_ResubComputeWindow( raw_array.first, raw_array.second, 1000, -1, 0, 0, 0, 0, &new_indices_raw, &num_resubs );
  fmt::print( "Performed resub {} times.  Reduced {} nodes.\n", num_resubs, new_entries > 0 ? raw_array.second - new_entries : 0 );
  abcresub::Abc_ResubPrepareManager( 0 );
  if ( new_entries == 0 )
  {
    return -1; /* next */
  }

  index_list new_indices( new_indices_raw, 2*new_entries );
  
  /* convert to mini AIG */
  aig_network window_aig_new;
  decode( window_aig_new, new_indices );
  write_verilog( window_aig_new, "win_opt.v" );

  /* TODO: verify that windows are equivalent */
  system( "abc -c \"cec -n win.v win_opt.v\"" );

  return 0;
}
