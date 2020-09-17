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

namespace experiments
{

struct abc_script_stats
{
  double time{0.0};
  uint64_t size{0};
  uint64_t depth{0};
};

abc_script_stats abc_script( std::string const& benchmark, std::string const& script )
{
  std::string command = fmt::format( "abc -q \"{}; {}; print_stats;\"", benchmark_path( benchmark ), script );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }
  
  std::stringstream ss( result );
  std::string line;

  std::smatch sm;

  abc_script_stats st;
  while ( std::getline( ss, line, '\n' ) )
  {
    if ( std::regex_match( line, sm, std::regex( R"(TOTAL\s+=\s+([0-9\.]+).*))" ) ) )
    {
      // std::cout << "MATCHED: " << std::stod( sm[1] ) << std::endl;
      st.time = std::stod( sm[1] );
    }

    if ( std::regex_search( line, sm, std::regex( R"(and\s*=\s*([0-9]+))" ) ) )
    {
      // std::cout << "MATCHED: " << std::stoul( sm[1] ) << ' ' << std::endl;
      st.size = std::stoul( sm[1] );
    }

    if ( std::regex_search( line, sm, std::regex( R"(lev\s*=\s*([0-9]+))" ) ) )
    {
      // std::cout << "MATCHED: " << std::stoul( sm[1] ) << ' ' << std::endl;
      st.depth = std::stoul( sm[1] );
    }
  }

  return st;
}

}

int main()
{
  using namespace experiments;
  using namespace mockturtle;

  experiment<std::string, uint64_t, uint64_t, uint64_t, uint64_t, double, uint64_t, uint64_t, double, uint64_t, uint64_t, double, uint64_t, uint64_t, double, bool>
    exp( "multithreaded_rewriting", "benchmark",
         "size", "depth",
         "rw.size",  "rw.depth",  "rw.time",
         "irw.size", "irw.depth", "irw.time",
         "drw.size", "drw.depth", "drw.time",
         "rrw.size", "rrw.depth", "rrw.time",
         "cec" );

  for ( auto const& benchmark : epfl_benchmarks( experiments::div ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[e] could not parse benchmark {}\n", benchmark );
      continue;
    }
    depth_view depth_aig0{aig};
    
    uint64_t const size_before = aig.num_gates();

    auto const st_rw = experiments::abc_script( benchmark, "rewrite; rewrite; rewrite -v" );
    auto const st_irw = experiments::abc_script( benchmark, "irw; irw; irw -v" );
    auto const st_drw = experiments::abc_script( benchmark, "drw; drw; drw -v" );

    multithreaded_cut_enumeration_params ps;
    multithreaded_cut_enumeration_stats st;
    multithreaded_cut_enumeration( aig, ps, &st );
    aig = cleanup_dangling( aig );
    multithreaded_cut_enumeration( aig, ps, &st );
    aig = cleanup_dangling( aig );
    multithreaded_cut_enumeration( aig, ps, &st );
    aig = cleanup_dangling( aig );    
    multithreaded_cut_enumeration( aig, ps, &st );
    aig = cleanup_dangling( aig );

    depth_view depth_aig1{aig};
    
    uint64_t const size_after = aig.num_gates();

    auto const cec = benchmark == "hyp" ? true : abc_cec( aig, benchmark );
    exp( benchmark, size_before, depth_aig0.depth(),
         st_rw.size, st_rw.depth, st_rw.time,
         st_irw.size, st_irw.depth, st_irw.time,
         st_drw.size, st_drw.depth, st_drw.time,
         size_after, depth_aig1.depth(), mockturtle::to_seconds( st.time_total ), cec );

    exp.save();
    exp.table();
  }

  exp.save();
  exp.table();

  return 0;
}
