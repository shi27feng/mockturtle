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
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/multithreaded_rewriting.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/node_resynthesis/xag_npn.hpp>
#include <mockturtle/algorithms/node_resynthesis/shannon.hpp>
#include <mockturtle/algorithms/node_resynthesis/dsd.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <lorina/aiger.hpp>

#include <fmt/format.h>

#include <string>

template<typename Plot, typename Ntk>
void run_rewrite( Plot& p, Ntk const& ntk )
{
  mockturtle::write_verilog( ntk, "temp.v" );

  uint64_t times{0};
  p( times++, ntk.num_gates() );

  experiments::abc_script( "temp.v", "strash" );
      
  while ( true )
  {
    fmt::print( "ABC rewrite {}: nodes={}\n", p.back().first, p.back().second );
    auto const st = experiments::abc_script( "temp.aig", "rewrite -v; cec -n temp.aig" );

    if ( p.back().second == st.size )
    {
      break;
    }

    p( times++, st.size );
  }
}

template<typename Plot, typename Ntk>
void run_irw( Plot& p, Ntk const& ntk )
{
  mockturtle::write_verilog( ntk, "temp.v" );

  uint64_t times{0};
  p( times++, ntk.num_gates() );

  experiments::abc_script( "temp.v", "strash" );
      
  while ( true )
  {
    fmt::print( "ABC irw {}: nodes={}\n", p.back().first, p.back().second );
    auto const st = experiments::abc_script( "temp.aig", "irw -v; cec -n temp.aig" );

    if ( p.back().second == st.size )
    {
      break;
    }

    p( times++, st.size );
  }
}

template<typename Plot, typename Ntk>
void run_drw( Plot& p, Ntk const& ntk )
{
  mockturtle::write_verilog( ntk, "temp.v" );

  uint64_t times{0};
  p( times++, ntk.num_gates() );

  experiments::abc_script( "temp.v", "strash" );
      
  while ( true )
  {
    fmt::print( "ABC drw {}: nodes={}\n", p.back().first, p.back().second );
    auto const st = experiments::abc_script( "temp.aig", "drw -v; cec -n temp.aig" );

    if ( p.back().second == st.size )
    {
      break;
    }

    p( times++, st.size );
  }
}

template<typename Plot, typename Ntk>
void run_rrw( Plot& p, Ntk const& ntk )
{
  using namespace mockturtle;
  multithreaded_cut_enumeration_params ps;
  multithreaded_cut_enumeration_stats st;

  uint64_t times{0};
  p( times++, ntk.num_gates() );

  /* make a copy of ntk */
  auto copy_ntk = cleanup_dangling( ntk );
  
  while ( true )
  { 
    fmt::print( "rrw {}: nodes={}\n", p.back().first, p.back().second );

    if ( has_cycle( copy_ntk ) )
    {
      std::cout << "cycle check: failed" << std::endl;
      std::abort();
    }
    multithreaded_cut_enumeration( copy_ntk, ps, &st );
    if ( has_cycle( copy_ntk ) )
    {
      std::cout << "cycle check: failed" << std::endl;
      std::abort();
    }
    copy_ntk = cleanup_dangling( copy_ntk );

    /* check equivalence */
    if ( !experiments::abc_cec2( copy_ntk, ntk ) )
    {
      std::cout << "cec: failed" << std::endl;
      std::abort();
    }
    
    uint64_t const size_after = copy_ntk.num_gates();    
    if ( p.back().second == size_after )
    {
      break;
    }
    
    p( times++, size_after );
  }

  mockturtle::write_verilog( copy_ntk, "final.v" );
}

template<typename Plot, typename Ntk>
void run_cut_rewrite4( Plot& p, Ntk const& ntk )
{
  using namespace mockturtle;
  xag_npn_resynthesis<aig_network> resyn;
  cut_rewriting_params ps;
  ps.cut_enumeration_ps.cut_size = 4;

  uint64_t times{0};
  p( times++, ntk.num_gates() );

  /* make a copy of ntk */
  auto copy_ntk = cleanup_dangling( ntk );
  
  while ( true )
  { 
    fmt::print( "mockturtle cut-rewrite {}: nodes={}\n", p.back().first, p.back().second );

    copy_ntk = cut_rewriting( copy_ntk, resyn, ps );    
    copy_ntk = cleanup_dangling( copy_ntk );

    /* check equivalence */
    if ( !experiments::abc_cec2( copy_ntk, ntk ) )
    {
      std::cout << "cec: failed" << std::endl;
      std::abort();
    }
    
    uint64_t const size_after = copy_ntk.num_gates();    
    if ( p.back().second == size_after )
    {
      break;
    }
    
    p( times++, size_after );
  }
}

template<typename Plot, typename Ntk>
void run_cut_rewrite6_shannon( Plot& p, Ntk const& ntk )
{
  using namespace mockturtle;

  xag_npn_resynthesis<aig_network> resyn;
  cut_rewriting_params ps;
  ps.cut_enumeration_ps.cut_size = 4;

  shannon_resynthesis<aig_network,xag_npn_resynthesis<aig_network>> resyn2( 4, &resyn );
  
  uint64_t times{0};
  p( times++, ntk.num_gates() );

  /* make a copy of ntk */
  auto copy_ntk = cleanup_dangling( ntk );
  
  while ( true )
  { 
    fmt::print( "mockturtle cut-rewrite {}: nodes={}\n", p.back().first, p.back().second );

    copy_ntk = cut_rewriting( copy_ntk, resyn2, ps );    
    copy_ntk = cleanup_dangling( copy_ntk );

    /* check equivalence */
    if ( !experiments::abc_cec2( copy_ntk, ntk ) )
    {
      std::cout << "cec: failed" << std::endl;
      std::abort();
    }
    
    uint64_t const size_after = copy_ntk.num_gates();    
    if ( p.back().second == size_after )
    {
      break;
    }
    
    p( times++, size_after );
  }
}

template<typename Plot, typename Ntk>
void run_cut_rewrite6_dsd( Plot& p, Ntk const& ntk )
{
  using namespace mockturtle;

  xag_npn_resynthesis<aig_network> resyn;
  cut_rewriting_params ps;
  ps.cut_enumeration_ps.cut_size = 4;

  dsd_resynthesis_params dsd_ps;
  dsd_ps.prime_input_limit = 4;
  dsd_resynthesis<aig_network,xag_npn_resynthesis<aig_network>> resyn2( resyn, dsd_ps );
  
  uint64_t times{0};
  p( times++, ntk.num_gates() );

  /* make a copy of ntk */
  auto copy_ntk = cleanup_dangling( ntk );
  
  while ( true )
  { 
    fmt::print( "mockturtle cut-rewrite {}: nodes={}\n", p.back().first, p.back().second );

    copy_ntk = cut_rewriting( copy_ntk, resyn2, ps );    
    copy_ntk = cleanup_dangling( copy_ntk );

    /* check equivalence */
    if ( !experiments::abc_cec2( copy_ntk, ntk ) )
    {
      std::cout << "cec: failed" << std::endl;
      std::abort();
    }
    
    uint64_t const size_after = copy_ntk.num_gates();    
    if ( p.back().second == size_after )
    {
      break;
    }
    
    p( times++, size_after );
  }
}


// adder.png
// bar.png
// div.png

int main()
{
  using namespace experiments;
  using namespace mockturtle;
  namespace plt = matplotlibcpp;

  // & ~experiments::hyp
  for ( auto const& benchmark : epfl_benchmarks( experiments::div ) )
  {
    fmt::print( "[i] processing {}\n", benchmark );

    aig_network aig;
    if ( lorina::read_aiger( benchmark_path( benchmark ), aiger_reader( aig ) ) != lorina::return_code::success )
    {
      fmt::print( "[e] could not parse benchmark {}\n", benchmark );
      continue;
    }

    experiments::plot<uint64_t, uint64_t> plot_rrw( "rrw" );
    run_rrw( plot_rrw, aig );
    
    experiments::plot<uint64_t, uint64_t> plot_rewrite( "rw" );
    run_rewrite( plot_rewrite, aig );

    experiments::plot<uint64_t, uint64_t> plot_irw( "irw" );
    run_irw( plot_irw, aig );

    experiments::plot<uint64_t, uint64_t> plot_drw( "drw" );
    run_drw( plot_drw, aig );

    experiments::plot<uint64_t, uint64_t> plot_cut_rewrite4( "crw4" );
    run_cut_rewrite4( plot_cut_rewrite4, aig );

    // experiments::plot<uint64_t, uint64_t> plot_cut_rewrite6_shannon( "crw6s" );
    // run_cut_rewrite6_shannon( plot_cut_rewrite6_shannon, aig );

    experiments::plot<uint64_t, uint64_t> plot_cut_rewrite6_dsd( "crw6d" );
    run_cut_rewrite6_dsd( plot_cut_rewrite6_dsd, aig );
    
    fmt::print( "[i] write graph {}.png\n", benchmark );    
    plt::xkcd();
    plot_rewrite.show();
    plot_irw.show();
    plot_drw.show();
    plot_cut_rewrite4.show();
    // plot_cut_rewrite6_shannon.show();
    plot_cut_rewrite6_dsd.show();
    plot_rrw.show();
    write_plot( fmt::format( "{}.png", benchmark ), benchmark );
  }
  
  return 0;
}
