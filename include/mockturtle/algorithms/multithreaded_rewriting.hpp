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

/*!
  \file multithreaded_rewritin.hpp
  \brief Boolean rewriting using multiple threads

  \author Heinz Riener
*/

#pragma once

#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"
#include "../io/write_verilog.hpp"
#include "simulation.hpp"

#include "rewriting/index_list.hpp"
#include "rewriting/thread_pool.hpp"
#include "rewriting/window_view.hpp"
#include "rewriting/debug.hpp"

#include "../utils/abc_resub.hpp"

#include <kitty/kitty.hpp>

#include <random>
#include <limits>

namespace mockturtle
{

/*! \brief Parameters for multithreaded_cut_enumeration.
 *
 */
struct multithreaded_cut_enumeration_params
{
  /*! \brief Number of threads for cut enumeration. */
  uint64_t num_threads{1u};

  /*! \brief Be verbose. */
  bool verbose{true};

  /*! \brief Be very verbose. */
  bool very_verbose{false};

  /* \brief Enable additional checks to verify that the results are correct. */
  bool verify{false};
};

/*! \brief Statistics for multithreaded_cut_enumeration.
 *
 */
struct multithreaded_cut_enumeration_stats
{
  /*! \brief Total time. */
  stopwatch<>::duration time_total{0};

  /*! \brief Time for creating windows. */
  stopwatch<>::duration time_create_window{0};
  stopwatch<>::duration time_initialize_window{0};
  stopwatch<>::duration time_window_inputs{0};
  stopwatch<>::duration time_add_side_inputs{0};
  stopwatch<>::duration time_explore{0};
  stopwatch<>::duration time_expand_inputs{0};
  stopwatch<>::duration time_grow{0};
  stopwatch<>::duration time_try_adding_node{0};

  /*! \brief Total number of created windows. */
  uint64_t total_num_windows{0};

  /*! \brief Total number of candidate windows. */
  uint64_t total_num_candidates{0};

  /*! \brief Total number of leaves. */
  uint64_t total_num_leaves{0};

  /*! \brief Total number of nodes. */
  uint64_t total_num_nodes{0};

  uint64_t count_success_resub{0};
  uint64_t count_failed_resub{0};

  /*! \brief Prints report. */
  void report() const
  {
    fmt::print( "[i] total time                = {:>7.2f} secs\n", to_seconds( time_total ) );
    fmt::print( "[i] time create_window        = {:>7.2f} secs\n", to_seconds( time_create_window ) );
    fmt::print( "[i]   time initialize_window  = {:>7.2f} secs\n", to_seconds( time_initialize_window ) );
    fmt::print( "[i]     time explore          = {:>7.2f} secs\n", to_seconds( time_explore ) );
    fmt::print( "[i]   time window_inputs      = {:>7.2f} secs\n", to_seconds( time_window_inputs ) );
    fmt::print( "[i]   time grow               = {:>7.2f} secs\n", to_seconds( time_grow ) );
    fmt::print( "[i]     time add_side_inputs  = {:>7.2f} secs\n", to_seconds( time_add_side_inputs ) );
    fmt::print( "[i]     time expand_inputs    = {:>7.2f} secs\n", to_seconds( time_expand_inputs ) );

    fmt::print( "[i] successful resubs  = {:>5d}\n", count_success_resub );
    fmt::print( "[i] failed resubs      = {:>5d}\n", count_failed_resub );

    fmt::print( "[i] number of leaves = {:>5d} (avg. leaves per window {:>5.2f})\n",
                total_num_windows, double(total_num_leaves) / total_num_windows );
    fmt::print( "[i] number of nodes = {:>5d} (avg. nodes per window {:>5.2f})\n",
                total_num_nodes, double(total_num_nodes) / total_num_windows );
    fmt::print( "[i] number of candidates = {:>5d}\n",
                total_num_candidates );
    fmt::print( "[i] number of windows = {:>5d} (avg. success rate {:>5.2f})\n",
                total_num_windows, double(total_num_windows) / total_num_candidates );
  }
};

namespace detail
{

template<class Ntk>
class window_manager
{
public:
  /*! \brief Maximum number of inputs. */
  static constexpr uint64_t const max_inputs = 6u;

  /*! \brief Maximum number of levels. */
  static constexpr uint64_t const max_levels = 5u;

public:
  using node   = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit window_manager( Ntk const& ntk, multithreaded_cut_enumeration_stats& st )
    : ntk( ntk )
    , st( st )
    , levels( ntk.depth() + 1u )
    , paths( ntk.size() )
    , size ( ntk.size() )
  {
  }

  std::optional<std::pair<std::vector<node>,std::vector<node>>> create_window( node const& pivot )
  {
    stopwatch t( st.time_create_window );

    /* resize paths */
    levels.resize( ntk.depth() + 1 );
    // std::cout << "set up levels with " << ntk.depth() << std::endl;
    paths.resize( ntk.size() );

    std::optional<std::vector<node>> window_nodes = initialize_window( pivot );
    if ( !window_nodes )
    {
      return std::nullopt;
    }

    std::vector<node> window_leaves = create_window_inputs( *window_nodes );

    /* consider a window, which initially has a larger input space */
    if ( window_leaves.size() <= max_inputs + 3 )
    {
      grow( *window_nodes, window_leaves );
    }

    if ( window_leaves.size() <= max_inputs )
    {
      std::sort( std::begin( window_leaves ), std::end( window_leaves ) );

      /* ensure that no window node is dead */
      for ( const auto& n : *window_nodes )
      {
        if ( ntk.is_dead( n ) || n >= size )
        {
          return std::nullopt;
        }

        bool dead_window = false;
        ntk.foreach_fanin( n, [&]( signal const& fi ){
          if ( ntk.is_dead( ntk.get_node( fi ) ) )
          {
            dead_window = true;
            return false;
          }
          return true;
        });

        if ( dead_window )
        {
          return std::nullopt;
        }
      }

      /* topologically sort the nodes */
      auto const sorted_nodes = topo_sort( *window_nodes, window_leaves );
      // std::sort( std::begin( *window_nodes ), std::end( *window_nodes ) );

      /* ensure that no window leave is dead */
      for ( const auto& n : window_leaves )
      {
        if ( ntk.is_dead( n ) )
          return std::nullopt;

        bool dead_window = false;
        ntk.foreach_fanin( n, [&]( signal const& fi ){
          if ( ntk.is_dead( ntk.get_node( fi ) ) )
          {
            dead_window = true;
            return false;
          }
          return true;
        });

        if ( dead_window )
        {
          return std::nullopt;
        }
      }

      return std::make_pair( sorted_nodes, window_leaves );
    }
    else
    {
      return std::nullopt;
    }
  }

  std::vector<node> topo_sort( std::vector<node> const& nodes, std::vector<node> const& leaves )
  {
    std::vector<node> topo_order;
    topo_order.reserve( nodes.size() );

    ntk.incr_trav_id();
    ntk.set_visited( ntk.get_node( ntk.get_constant( false ) ), ntk.trav_id() );
    topo_order.emplace_back( ntk.get_node( ntk.get_constant( false ) ) );

    for ( const auto& l : leaves )
    {
      ntk.set_visited( l, ntk.trav_id() );
      topo_order.emplace_back( l );
    }

    for ( const auto& n : nodes )
    {
      topo_sort_recur( topo_order, n );
    }
    return topo_order;
  }

  void topo_sort_recur( std::vector<node>& topo_order, node const& n )
  {
    if ( ntk.is_dead( n ) || ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.foreach_fanin( n, [&]( signal const& s ){
      if ( ntk.visited( ntk.get_node( s ) ) != ntk.trav_id() )
      {
        topo_sort_recur( topo_order, ntk.get_node( s ) );
      }
    });

    /* all have been visited */
    ntk.set_visited( n, ntk.trav_id() );
    topo_order.emplace_back( n );
  }

private:
  /* collect nodes on the path recursively from the meeting point to the root node, excluding the meeting point */
  void gather( node const& n, std::vector<node>& visited )
  {
    if ( n == 0u || ntk.is_dead( n ) )
      return;
    visited.emplace_back( n );

    auto const prev = paths[n];
    if ( prev == 0 )
      return;

    assert( ntk.visited( prev ) == ntk.visited( n ) );
    gather( prev, visited );
  }

  /* explore the frontier of nodes in breath-first traversal */
  bool explore( std::vector<node>& visited, uint64_t start, node& pi_meet, node& pi_node )
  {
    stopwatch t( st.time_explore );

    uint64_t const end = visited.size();
    pi_meet = ntk.get_node( ntk.get_constant( false ) );
    pi_node = ntk.get_node( ntk.get_constant( false ) );

    bool result = false;
    for ( uint64_t i = start; i < end; ++i )
    {
      node const& n = visited.at( i );
      if ( ntk.is_constant( n ) || ntk.is_ci( n ) || ntk.is_dead( n ) )
        continue;

      ntk.foreach_fanin( n, [&]( signal const& fi ){
        if ( ntk.is_dead( ntk.get_node( fi ) ) )
          return true;

        /* if the node was visited on the paths to both fanins, collect it */
        if ( ntk.visited( n ) >= ntk.trav_id() - 1u &&
             ntk.visited( ntk.get_node( fi ) ) >= ntk.trav_id() - 1u &&
             ntk.visited( n ) != ntk.visited( ntk.get_node( fi ) ) )
        {
          pi_meet = ntk.get_node( fi );
          pi_node = n;
          result = true;
          return false; /* terminate and return TRUE */
        }

        /* if the node was visited already on this path, skip it */
        if ( ntk.visited( ntk.get_node( fi ) ) >= ntk.trav_id() - 1u )
        {
          assert( ntk.visited( n ) == ntk.visited( ntk.get_node( fi ) ) );
          return true; /* next */
        }

        /* label the node as visited */
        ntk.set_visited( ntk.get_node( fi ), ntk.visited( n ) );
        assert( paths.size() > ntk.get_node( fi ) );
        paths[ntk.get_node( fi )] = n;

        visited.push_back( ntk.get_node( fi ) );
        return true; /* next */
      });

      if ( result )
      {
        return true;
      }
    }
    return false;
  }

  std::optional<std::vector<node>> initialize_window( node const& pivot )
  {
    assert( !ntk.is_constant( pivot ) && !ntk.is_ci( pivot ) && !ntk.is_dead( pivot ) );

    stopwatch t( st.time_initialize_window );

    std::vector<node> visited;
    visited.reserve( 150u ); /* reserve more memory to avoid memory issues */

    uint64_t start{0};

    /* start paths for both fanins of the pivot node */
    ntk.foreach_fanin( pivot, [&]( signal const& fi ){
      if ( ntk.is_dead( ntk.get_node( fi ) ) )
        return true;

      ntk.incr_trav_id();
      visited.push_back( ntk.get_node( fi ) );
      assert( paths.size() > ntk.get_node( fi ) );
      paths[ntk.get_node( fi )] = 0;
      ntk.set_visited( ntk.get_node( fi ), ntk.trav_id() );
      return true;
    });

    /* perform several iterations of breath-first search */
    uint64_t i = 0u;
    for ( ; i < max_levels; ++i )
    {
      uint64_t next = visited.size();

      node meet, n;
      if ( explore( visited, start, meet, n ) )
      {
        /* found the shared path */
        // FIXME: skipping some error checks */

        /* collect the initial window */
        visited.clear();
        gather( paths[meet], visited );
        gather( n, visited );
        assert( !ntk.is_dead( pivot ) );
        visited.push_back( pivot );
        break;
      }

      start = next;
    }

    /* if no meeting point is found, make sure to return NULL */
    if ( i == max_levels )
    {
      return std::nullopt;
    }
    else
    {
      return visited;
    }
  }

  std::vector<node> create_window_inputs( std::vector<node>& window_nodes )
  {
    stopwatch t( st.time_window_inputs );

    std::vector<node> inputs;
    inputs.reserve( 10 );

    /* mark all window nodes as visited */
    ntk.incr_trav_id();
    for ( auto const& n : window_nodes )
    {
      ntk.set_visited( n, ntk.trav_id() );
    }

    /* collet fanins of these nodes as inputs */
    for ( auto const& n : window_nodes )
    {
      assert( !ntk.is_dead( n ) );
      assert( !ntk.is_constant( n ) && !ntk.is_ci( n ) && "node is not a gate" );

      ntk.foreach_fanin( n, [&]( signal const& fi ){
        if ( ntk.is_dead( ntk.get_node( fi ) ) )
          return true;

        node const& fanin_node = ntk.get_node( fi );
        if ( ntk.visited( fanin_node ) == ntk.trav_id() )
          return true;

        if ( std::find( std::begin( inputs ), std::end( inputs ), fanin_node ) == std::end( inputs ) )
        {
          inputs.emplace_back( fanin_node );
        }
        return true;
      });
    }

    /* mark inputs and add them to the window nodes */
    for ( const auto& i : inputs )
    {
      ntk.set_visited( i, ntk.trav_id() );
      window_nodes.emplace_back( i );
    }

    return inputs;
  }

  void add_side_inputs( std::vector<node>& window_nodes )
  {
    stopwatch t( st.time_add_side_inputs );

    /* ensure that levels is empty */
    assert( levels.size() > 0u );

    /* window nodes are labeled with the current traversal id */
    for ( node const& n : window_nodes )
    {
      if ( ntk.is_dead( n ) )
        continue;

      assert( ntk.visited( n ) == ntk.trav_id() );
      assert( !ntk.is_dead( n ) );
      assert( ntk.level( n ) >= 0 );
      assert( ntk.level( n ) < levels.size() );
      levels[ntk.level( n )].push_back( n );
    }

    /* iterate through all objects and explore their fanouts */
    for ( std::vector<node>& level : levels )
    {
      for ( node const& n : level )
      {
        assert( !ntk.is_dead( n ) );
        ntk.foreach_fanout( n, [&]( node const& fo, uint64_t index ){
          /* explore first 5 fanouts of the node */
          if ( index == 5u )
            return false;

          if ( ntk.is_dead( fo ) )
            return true;

          /* ensure that fo is an internal node */
          if ( ntk.is_constant( fo ) || ntk.is_ci( fo ) )
            return true;

          /* ensure that fo is in the window */
          if ( ntk.visited( fo ) == ntk.trav_id() )
            return true;

          /* ensure that fanins are in the window */
          bool fanins_are_in_the_window = true;
          ntk.foreach_fanin( fo, [&]( signal const& fi ){
            if ( ntk.is_dead( ntk.get_node( fi ) ) )
              return true;

            if ( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() )
            {
              fanins_are_in_the_window = false;
              return false;
            }
            return true;
          });
          if ( !fanins_are_in_the_window )
            return true;

          /* add fanout to the window and to the levelized structure */
          ntk.set_visited( fo, ntk.trav_id() );
          levels[ntk.level( fo )].emplace_back( fo );
          assert( !ntk.is_dead( fo ) );
          window_nodes.emplace_back( fo );
          return true;
        });
      }
    }

    /* iterate through the nodes in the levelized structure */
    for ( std::vector<node>& level : levels )
    {
      level.clear();
    }
  }

  void expand_inputs( std::vector<node>& window_nodes, std::vector<node>& inputs )
  {
    stopwatch t( st.time_expand_inputs );

    bool changed = true;
    while ( changed )
    {
      changed = false;

      auto it = std::begin( inputs );
      while ( it != std::end( inputs ) )
      {
        node const input = *it;

        /* skip all inputs that are not gates */
        if ( ntk.is_constant( input ) || ntk.is_ci( input ) || ntk.is_dead( input ) )
        {
          ++it;
          continue; /* next */
        }

        /* skip if none of the fanins has been marked */
        bool no_fanin_marked = true;
        ntk.foreach_fanin( input, [&]( signal const& fi ){
          if ( ntk.is_dead( ntk.get_node( fi ) ) )
            return true;

          if ( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() )
          {
            no_fanin_marked = false;
            return false;
          }
          return true;
        });
        if ( no_fanin_marked )
        {
          ++it;
          continue; /* next */
        }

        /* remove the current input from the vector */
        it = inputs.erase( it );

        /* ensure that the input is within the window nodes */
        assert( std::find( std::begin( window_nodes ), std::end( window_nodes ), input ) != std::end( window_nodes ) );

        ntk.foreach_fanin( input, [&]( signal const& fi ){
          if ( ntk.is_dead( ntk.get_node( fi ) ) )
            return true;

          if ( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() )
          {
            return true; /* next */
          }

            assert( std::find( std::begin( inputs ), std::end( inputs ), ntk.get_node( fi ) ) == std::end( inputs ) );
            inputs.emplace_back( ntk.get_node( fi ) );

            try_adding_node( { ntk.get_node( fi ) }, window_nodes );
            assert( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() || ntk.is_dead( ntk.get_node( fi ) ) );
            return true;
          });

        changed = true;
      }
    }
  }

  std::optional<node> select_one_input( std::vector<node>& inputs )
  {
    std::optional<node> best_input{std::nullopt};
    std::optional<uint64_t> best_weight{std::nullopt};

    for ( node const& i : inputs )
    {
      if ( ntk.is_constant( i ) || ntk.is_ci( i ) || ntk.is_dead( i ) )
        continue;

      std::vector<node> fis;
      ntk.foreach_fanin( i, [&]( signal const& fi ){
        if ( ntk.is_dead( ntk.get_node( fi ) ) )
          return true;

        assert( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() );
        fis.emplace_back( ntk.get_node( fi ) );
        return true;
      });

      std::vector<node> t;
      uint64_t const weight = try_adding_node( fis, t, false );
      if ( !best_weight || *best_weight < weight )
      {
        best_weight = weight;
        best_input = i;
      }
    }
    return best_input;
  }

  void grow( std::vector<node>& window_nodes, std::vector<node>& inputs )
  {
    stopwatch t( st.time_grow );

    add_side_inputs( window_nodes );
    expand_inputs( window_nodes, inputs );

    std::optional<node> n;
    while ( inputs.size() < max_inputs && ( n = select_one_input( inputs ) ) )
    {
      if ( ntk.is_dead( *n ) )
        continue;

      std::vector<node> fanin_nodes;
      ntk.foreach_fanin( *n, [&]( signal const& fi ){
        if ( ntk.is_dead( ntk.get_node( fi ) ) )
          return;

        assert( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() );
        assert( !ntk.is_dead( ntk.get_node( fi ) ) );
        fanin_nodes.emplace_back( ntk.get_node( fi ) );
      });

      try_adding_node( fanin_nodes, window_nodes );
      inputs.erase( std::remove( std::begin( inputs ), std::end( inputs ), *n ), std::end( inputs ) );

      ntk.foreach_fanin( *n, [&]( signal const& fi ){
        if ( ntk.is_dead( ntk.get_node( fi ) ) )
          return;

        assert( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() );
        assert( !ntk.is_dead( ntk.get_node( fi ) ) );
        inputs.emplace_back( ntk.get_node( fi ) );
      });

      expand_inputs( window_nodes, inputs );
    }
  }

  uint64_t try_adding_node( std::vector<node> const& pivots, std::vector<node>& nodes, bool collect_nodes = true )
  {
    stopwatch t( st.time_try_adding_node );

    /* ensure that levels is empty */
    assert( levels.size() > 0u );

    uint64_t count{0};

    /* add the pivots to the window and to the levelized structure */
    for ( node const& pivot : pivots )
    {
      if ( ntk.is_dead( pivot ) )
        continue;

      assert( ntk.visited( pivot ) != ntk.trav_id() && "pivot must not be part of the window" );
      ntk.set_visited( pivot, ntk.trav_id() );
      assert( !ntk.is_dead( pivot ) );
      assert( levels.size() > ntk.level( pivot ) );
      levels[ntk.level( pivot )].emplace_back( pivot );
    }

    /* iterate through all objects and explore their fanouts */
    for ( std::vector<node>& level : levels )
    {
      for ( node const& n : level )
      {
        ntk.foreach_fanout( n, [&]( node const& fo, uint64_t index ){
            /* explore first 5 fanouts of the node */
            if ( index == 5u )
              return false;

            /* ensure that fo is an internal node */
            if ( ntk.is_constant( fo ) || ntk.is_ci( fo ) )
              return true;

            /* ensure that fo is in the window */
            if ( ntk.visited( fo ) == ntk.trav_id() )
              return true;

            /* ensure that fanins are in the window */
            bool fanins_are_in_the_window = true;
            ntk.foreach_fanin( fo, [&]( signal const& fi ){
                if ( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() )
                {
                  fanins_are_in_the_window = false;
                  return false;
                }
                return true;
              });
            if ( !fanins_are_in_the_window )
              return true;

            /* add fanout to the window and to the levelized structure */
            ntk.set_visited( fo, ntk.trav_id() );
            // if ( ntk.level( fo ) >= levels.size() )
            // {
            //   std::cout << fo << ' ' << ntk.level( fo ) << ' ' << levels.size() << std::endl;
            //   std::cout << ntk.depth() << std::endl;
            // }
            assert( ntk.level( fo ) < levels.size() );
            levels[ntk.level( fo )].emplace_back( fo );
            ++count;
            return true;
          });
      }
    }

    /* iterate through the nodes in the levelized structure */
    for ( std::vector<node>& level : levels )
    {
      for ( node& l : level )
      {
        if ( !collect_nodes ) /* it was a test run - unmark the node */
        {
          ntk.set_visited( l, ntk.visited( l ) - 1u ); /* unmark node */
        }
        else /* it was a real run - permanently add to the node to the window */
        {
          if ( !ntk.is_dead( l ) )
          {
            nodes.emplace_back( l );
          }
        }
      }
      level.clear();
    }

    return count;
  }

private:
  Ntk const& ntk;
  multithreaded_cut_enumeration_stats& st;

  std::vector<std::vector<node>> levels;
  std::vector<node> paths;
  uint64_t size{0};
};

} /* namespace detail */

template<typename Ntk>
std::vector<typename Ntk::signal> find_roots( Ntk const& ntk, std::vector<typename Ntk::node> const& nodes, std::vector<typename Ntk::node> const& leaves )
{
  using signal = typename Ntk::signal;

  std::vector<int32_t> refs( ntk.size() );

  /* reference nodes */
  for ( const auto& n : nodes )
  {
    ntk.foreach_fanin( n, [&]( signal const& fi ){
      ++refs[ntk.get_node( fi )];
    });
  }

  refs[0] = ntk.fanout_size( 0 );
  for ( const auto& l : leaves )
  {
    refs[l] = ntk.fanout_size( l );
  }

  std::vector<signal> roots;
  for ( const auto& n : nodes )
  {
    if ( refs[n] != int32_t( ntk.fanout_size( n ) ) )
      roots.emplace_back( ntk.make_signal( n ) );
  }

  return roots;
}

namespace detail
{

template<class Ntk>
class multithreaded_cut_enumeration_impl
{
public:
  using node   = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit multithreaded_cut_enumeration_impl( Ntk& ntk, multithreaded_cut_enumeration_params const& ps, multithreaded_cut_enumeration_stats& st )
    : ntk( ntk )
    , ps( ps )
    , st( st )
  {
    auto const update_level_of_new_node = [&]( const auto& n ) {
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_existing_node = [&]( node const& n, const auto& old_children ) {
      (void)old_children;
      ntk.resize_levels();
      update_node_level( n );
    };

    auto const update_level_of_deleted_node = [&]( const auto& n ) {
      ntk.set_level( n, -1 );
      // std::cout << "set_level " << n << ' ' << -1 << std::endl;
    };

    ntk._events->on_add.emplace_back( update_level_of_new_node );
    ntk._events->on_modified.emplace_back( update_level_of_existing_node );
    ntk._events->on_delete.emplace_back( update_level_of_deleted_node );
  }

  ~multithreaded_cut_enumeration_impl()
  {
  }

  void kresub_aig_test()
  {
    static uint64_t counter = 0ul;

    window_manager windows( ntk, st );
    ntk.foreach_gate( [&]( node const& n ){
      auto const result = windows.create_window( n );
      if ( !result )
      {
        return true;
      }

      window_view2 window( ntk, result->second, find_roots( ntk, result->first, result->second ), result->first );
      index_list indices;
      encode( indices, window );

      index_list new_indices;
      {
        int *indices_raw = ABC_CALLOC( int, 2*indices.num_entries()+1 );
        uint64_t pos = 0;
        indices.foreach_entry( [&]( uint32_t i, uint32_t j ){
          indices_raw[pos] = i;
          indices_raw[pos+1] = j;
          pos+=2;
        });

        abcresub::Abc_ResubPrepareManager( 1 );
        int *new_indices_raw = nullptr;
        int num_resubs = 0;
        uint64_t new_entries = abcresub::Abc_ResubComputeWindow( indices_raw, indices.num_entries(), 1000, -1, 0, 0, 0, 0, &new_indices_raw, &num_resubs );
        abcresub::Abc_ResubPrepareManager( 0 );

        fmt::print( "Performed resub {} times.  Reduced {} nodes.\n", num_resubs, new_entries > 0 ? indices.num_entries() - new_entries : 0 );          

        if ( new_entries > 0 )
        {
          for ( uint64_t i = 0; i < 2*new_entries; i+=2 )
          {
            new_indices.push2( new_indices_raw[i], new_indices_raw[i+1] );
          }
        }
        else if ( new_entries == 0 )
        {
          aig_network window_aig;
          decode( window_aig, indices );
          write_verilog( window_aig, fmt::format( "win{}.v", counter++ ) );
        }

        if ( indices_raw )
        {
          ABC_FREE( indices_raw );
        }
        if ( new_indices_raw )
        {
          ABC_FREE( new_indices_raw );
        }
      }

      return true;
    });
  }

  void enumerate_windows_test()
  {
    window_manager windows( ntk, st );

    // progress_bar pbar{ntk.size(), "resub |{0}| node = {1:>4}   cand = {2:>4}   est. gain = {3:>5}", ps.progress};

    auto const size = ntk.num_gates();
    ntk.foreach_gate( [&]( node const& n, auto i ) {
      if ( i >= size )
      {
        return false; /* terminate */
      }

      if ( ntk.is_dead( n ) )
      {
        return true; /* next */
      }

      ++st.total_num_candidates;

      auto const result = windows.create_window( n );
      if ( !result )
      {
        return true;
      }

#if 0
      std::cout << "leaves = ";
      for ( const auto& node : result->second )
      {
        std::cout << node << ' ';
      }
      std::cout << std::endl;

      std::cout << "nodes = ";
      for ( const auto& node : result->first )
      {
        assert( !ntk.is_dead( node ) );
        std::cout << node << ' ';
      }
      std::cout << std::endl;
#endif

      /* TOOD: ensure that the constant false is included in the window */
      /* TODO: ensure that the nodes are topologically sorted */

      std::vector<signal> const roots = find_roots( ntk, result->first, result->second );
#if 0
      std::cout << "root = ";
      for ( const auto& signal : roots )
      {
        std::cout << ntk.get_node( signal ) << ':' << ntk.is_complemented( signal ) << ' ';
      }
      std::cout << std::endl;
#endif

      /* make a view on the window */
      window_view2 win( ntk, result->second, roots, result->first );
      assert( win.size() == result->first.size() );
      assert( win.num_pis() == result->second.size() );

#if 0
      std::cout << win.size() << std::endl;
      win.foreach_node( [&]( auto const& n ){
        std::cout << n << ' ';
      });
      std::cout << std::endl;

      win.foreach_gate( [&]( auto const& n ){
        std::cout << n << ' ';
      });
      std::cout << std::endl;

      win.foreach_pi( [&]( auto const& n ){
        std::cout << n << ' ';
      });
      std::cout << std::endl;

      win.foreach_po( [&]( auto const& n ){
        std::cout << ntk.get_node( n ) << ' ';
      });
      std::cout << std::endl;
#endif

      assert( check_size( win ) );
      // assert( !has_dangling_roots( win ) ); /* fanouts are not update */
      // write_verilog( win, std::cout );

      ++st.total_num_windows;
      st.total_num_nodes += result->first.size();
      st.total_num_leaves += result->second.size();

      /* convert window into index list */
      index_list indices;
      encode( indices, win );
      // indices.print();
      // indices.print_raw();

      /* use ABC to optimize the window */
      index_list new_indices;
      {
        int *indices_raw = ABC_CALLOC( int, 2*indices.num_entries()+1 );
        uint64_t pos = 0;
        indices.foreach_entry( [&]( uint32_t i, uint32_t j ){
          indices_raw[pos] = i;
          indices_raw[pos+1] = j;
          pos+=2;
        });

        abcresub::Abc_ResubPrepareManager( 1 );
        int *new_indices_raw = nullptr;
        int num_resubs = 0;
        uint64_t new_entries = abcresub::Abc_ResubComputeWindow( indices_raw, indices.num_entries(), 1000, -1, 0, 0, 0, 0, &new_indices_raw, &num_resubs );
        abcresub::Abc_ResubPrepareManager( 0 );

        if ( new_entries > 0 )
        {
          // fmt::print( "Performed resub {} times.  Reduced {} nodes.\n", num_resubs, new_entries > 0 ? indices.num_entries() - new_entries : 0 );
          for ( uint64_t i = 0; i < 2*new_entries; i+=2 )
          {
            new_indices.push2( new_indices_raw[i], new_indices_raw[i+1] );
          }
        }

        if ( indices_raw )
        {
          ABC_FREE( indices_raw );
        }
        if ( new_indices_raw )
        {
          ABC_FREE( new_indices_raw );
        }

        if ( new_entries == 0 )
        {
          return true; /* next window */
        }
      }

      if ( ps.verify )
      {
        check_window_equivalence( indices, new_indices );
      }

#if 0
      /* old code: a much shorter integration of ABC */
      /* optimize index list using the resubstitution algorithm */
      auto const raw_array = indices.raw_data();

      int num_resubs;
      int *new_indices_raw;
      uint64_t new_entries = abcresub::Abc_ResubComputeWindow( raw_array.first, raw_array.second, 1000, -1, 0, 0, 0, 0, &new_indices_raw, &num_resubs );
      // fmt::print( "Performed resub {} times.  Reduced {} nodes.\n", num_resubs, new_entries > 0 ? raw_array.second - new_entries : 0 );
      if ( new_entries == 0 )
      {
        return true; /* next */
      }

      index_list new_indices( new_indices_raw, 2*new_entries );
      ABC_FREE( new_indices_raw );
#endif

      assert( win.num_pos() == new_indices.num_pos() );

      /* substitute optimized window into the large network */
      std::vector<signal> inputs;
      win.foreach_pi( [&]( node const& n ){
          inputs.push_back( win.make_signal( n ) );
        });

      /* collect outputs and ensure that they are regular */
      std::vector<node> outputs;
      win.foreach_po( [&]( signal const& s ){
          if ( win.is_complemented( s ) )
          {
            std::cout << "[e] window outputs are not normalized" << std::endl;
            std::abort();
          }
          else
          {
            outputs.push_back( win.get_node( s ) );
          }
        });

      uint64_t counter = 0;
      insert( ntk, std::begin( inputs ), std::end( inputs ), new_indices,
              [&]( signal const& s ){
                auto const output = outputs.at( counter++ );

                /* only substitute if the output is different */
                if ( output == ntk.get_node( s ) && !ntk.is_complemented( s ) )
                {
                  return true;
                }

                // fmt::print( "substitute node {} with signal {}{}\n",
                //             output, ntk.is_complemented( s ) ? "~" : "", ntk.get_node( s ) );
                ntk.substitute_node( output, s );
                return true;
              });

      if ( ps.verify && has_cycle( ntk ) )
      {
        std::cout << "[e] cycle detected in DAG" << std::endl;
        std::abort();
      }

      return true;
    });
  }

  void check_window_equivalence( index_list const& a, index_list const& b)
  {
    aig_network window_aig;
    decode( window_aig, a );
    write_verilog( window_aig, "win.v" );

    aig_network window_aig_new;
    decode( window_aig_new, b );
    write_verilog( window_aig_new, "win_opt.v" );

    system( "abc -q \"cec -n win.v win_opt.v\" >> verify.log" );
  }

  void multithreaded_test()
  {
#if 0
    thread_pool threads{ps.num_threads};
    uint64_t const size = ntk.size();
    uint64_t num_processed_nodes = 0u;

    std::random_device rand;
    std::mt19937 gen( rand() ); /* seed random generator */
    std::uniform_int_distribution<> distribution( 0, size - 1u ); /* define the range */

    window_manager windows( ntk, st );
    while ( num_processed_nodes < size )
    {
      auto const index = distribution( gen );
      auto const pivot = ntk.index_to_node( index );
      assert( pivot < ntk.size() );

      /* ensure that each node is evaluated only once */
      if ( mark0( pivot ) )
      {
        continue;
      }
      set_mark0( pivot );

      /* create window */
      std::vector<std::vector<node>> levels( ntk.depth() + 1u );
      std::vector<node> paths( ntk.size() );

      std::vector<node> nodes;
      std::vector<node> leaves;
      if ( windows.create_window( pivot, levels, paths, nodes, leaves ) )
      {
        std::cout << "success" << std::endl;
      }
      else
      {
        // std::cout << "failed" << std::endl;
      }

      threads.enqueue( [=]{
          /* evaluate all cuts of the current node concurrently */
          // if ( ps.very_verbose )
          // {
          //   fmt::print( "[i] evaluate cut for node at index {} ({}/{})\n", index, num_processed_nodes, ntk.size() );
          // }
        } );

      ++num_processed_nodes;
    }
#endif
  }

  void run()
  {
    stopwatch t( st.time_total );
    enumerate_windows_test();
  }

  bool mark0( node const& n )
  {
    return ntk.value( n ) & 1;
  }

  void set_mark0( node const& n )
  {
    ntk.set_value( n, ntk.value( n ) | 1 );
  }

  void clear_mark0( node const& n )
  {
    ntk.set_value( n, ntk.value( n ) & 0 );
  }

  /* maybe should move to depth_view */
  void update_node_level( node const& n, bool top_most = true )
  {
    uint32_t curr_level = ntk.level( n );

    uint32_t max_level = 0;
    ntk.foreach_fanin( n, [&]( const auto& f ) {
      auto const p = ntk.get_node( f );
      auto const fanin_level = ntk.level( p );
      if ( fanin_level > max_level )
      {
        max_level = fanin_level;
      }
    } );
    ++max_level;

    if ( curr_level != max_level )
    {
      ntk.set_level( n, max_level );
      // std::cout << "set_level " << n << ' ' << max_level << std::endl;

      if ( max_level > ntk.depth() )
      {
        ntk.set_depth( max_level + 1 ); /* FIXME */
      }

      /* update only one more level */
      if ( top_most )
      {
        ntk.foreach_fanout( n, [&]( const auto& p ) {
          update_node_level( p, false );
        } );
      }
    }
  }

private:
  Ntk& ntk;
  multithreaded_cut_enumeration_params const& ps;
  multithreaded_cut_enumeration_stats& st;
}; /* multithreaded_cut_enumeration_impl */

} /* namespace detail */

template<class Ntk>
void multithreaded_cut_enumeration( Ntk const& ntk, multithreaded_cut_enumeration_params const& ps = {}, multithreaded_cut_enumeration_stats* pst = nullptr )
{
  depth_view<Ntk> ntk2{ntk};
  fanout_view<decltype( ntk2 )> ntk3{ntk2};

  multithreaded_cut_enumeration_stats st;
  detail::multithreaded_cut_enumeration_impl cut_enum( ntk3, ps, st );
  cut_enum.run();

  if ( pst )
  {
    *pst = st;
  }
}

template<class Ntk>
void kresub_aig_test( Ntk const& ntk )
{
  depth_view<Ntk> ntk2{ntk};
  fanout_view<decltype( ntk2 )> ntk3{ntk2};

  multithreaded_cut_enumeration_params ps;
  multithreaded_cut_enumeration_stats st;
  detail::multithreaded_cut_enumeration_impl cut_enum( ntk3, ps, st );
  cut_enum.kresub_aig_test();
}

} /* namespace mockturtle */
