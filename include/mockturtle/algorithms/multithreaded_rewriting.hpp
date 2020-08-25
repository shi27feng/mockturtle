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

#include <condition_variable>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>
#include <random>
#include <limits>

namespace mockturtle
{

namespace detail
{

class thread_pool
{
public:
  using task_type = std::function<void()>;

public:
  explicit thread_pool( uint64_t num_threads )
    : num_threads( num_threads )
  {
    start();
  }

  ~thread_pool()
  {
    stop();
  }

  template<class T>
  auto enqueue( T task ) -> std::future<decltype( task() )>
  {
    auto wrapper = std::make_shared<std::packaged_task<decltype( task() )()>>( std::move( task ) );
    {
      std::unique_lock<std::mutex> lock{event_mutex};
      tasks.emplace( [=]{
          (*wrapper)();
        } );
    }

    event.notify_one();
    return wrapper->get_future();
  }

private:
  void start()
  {
    for ( auto i = 0ul; i < num_threads; ++i )
    {
      threads.emplace_back( [=]{
          while ( true )
          {
            task_type task;

            {
              std::unique_lock<std::mutex> lock{event_mutex};
              event.wait( lock, [=]{ return stopping || !tasks.empty(); } );

              if ( stopping && tasks.empty() )
              {
                break;
              }

              /* take first task from the queue */
              task = std::move( tasks.front() );
              tasks.pop();
            }

            /* execute task */
            task();
          }
        } );
    }
  }

  void stop() noexcept
  {
    {
      std::unique_lock<std::mutex> lock{event_mutex};
      stopping = true;
    }
    event.notify_all();

    for ( auto &thread : threads )
    {
      thread.join();
    }
  }

private:
  uint64_t num_threads;
  std::vector<std::thread> threads;
  std::queue<task_type> tasks;

  std::condition_variable event;
  std::mutex event_mutex;
  bool stopping = false;
};

} /* namespace detail */

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
  using node   = node<Ntk>;
  using signal = signal<Ntk>;

public:
  explicit window_manager( Ntk const& ntk, multithreaded_cut_enumeration_stats& st )
    : ntk( ntk )
    , st( st )
    , levels( ntk.depth() + 1u )
    , paths( ntk.size() )
  {
  }

  std::optional<std::pair<std::vector<node>,std::vector<node>>> create_window( node const& pivot )
  {
    stopwatch t( st.time_create_window );

    std::optional<std::vector<node>> window_nodes = initialize_window( pivot );
    if ( !window_nodes )
    {
      return std::nullopt;
    }

    std::vector<node> window_leaves = create_window_inputs( *window_nodes );

    /* consider a window, which initially has a larger input space */
    if ( window_leaves.size() <= max_inputs + 2 )
    {
      grow( *window_nodes, window_leaves );
    }

    if ( window_leaves.size() <= max_inputs )
    {
      std::sort( std::begin( *window_nodes ), std::end( *window_nodes ) );
      std::sort( std::begin( window_leaves ), std::end( window_leaves ) );
      return std::make_pair( *window_nodes, window_leaves );
    }
    else
    {
      return std::nullopt;
    }
  }

private:
  /* collect nodes on the path recursively from the meeting point to the root node, excluding the meeting point */
  void gather( node const& n, std::vector<node>& visited )
  {
    if ( n == 0u )
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
      if ( ntk.is_constant( n ) || ntk.is_ci( n ) )
        continue;

      ntk.foreach_fanin( n, [&]( signal const& fi ){
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
          paths[ntk.get_node( fi )] = n;
          visited.emplace_back( ntk.get_node( fi ) );
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
    stopwatch t( st.time_initialize_window );

    std::vector<node> visited( 100u );
    assert( !ntk.is_constant( pivot ) && !ntk.is_ci( pivot ) );

    uint64_t start{0};

    /* start paths for both fanins of the pivot node */
    ntk.foreach_fanin( pivot, [&]( signal const& fi ){
        ntk.incr_trav_id();
        visited.emplace_back( ntk.get_node( fi ) );
        paths[ntk.get_node( fi )] = 0;
        ntk.set_visited( ntk.get_node( fi ), ntk.trav_id() );
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
        visited.emplace_back( pivot );
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
      assert( !ntk.is_constant( n ) && !ntk.is_ci( n ) && "node is not a gate" );

      ntk.foreach_fanin( n, [&]( signal const& fi ){
          node const& fanin_node = ntk.get_node( fi );
          if ( ntk.visited( fanin_node ) == ntk.trav_id() )
          {
            return true;
          }

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
      assert( ntk.visited( n ) == ntk.trav_id() );
      levels[ntk.level( n )].emplace_back( n );
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
            levels[ntk.level( fo )].emplace_back( fo );
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
        if ( ntk.is_constant( input ) || ntk.is_ci( input ) )
        {
          ++it;
          continue; /* next */
        }

        /* skip if none of the fanins has been marked */
        bool no_fanin_marked = true;
        ntk.foreach_fanin( input, [&]( signal const& fi ){
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
            if ( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() )
            {
              return true; /* next */
            }

            assert( std::find( std::begin( inputs ), std::end( inputs ), ntk.get_node( fi ) ) == std::end( inputs ) );
            inputs.emplace_back( ntk.get_node( fi ) );

            try_adding_node( { ntk.get_node( fi ) }, window_nodes );
            assert( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() );
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
      if ( ntk.is_constant( i ) || ntk.is_ci( i ) )
        continue;

      std::vector<node> fis;
      ntk.foreach_fanin( i, [&]( signal const& fi ){
          assert( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() );
          fis.emplace_back( ntk.get_node( fi ) );
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
      std::vector<node> fanin_nodes;
      ntk.foreach_fanin( *n, [&]( signal const& fi ){
          assert( ntk.visited( ntk.get_node( fi ) ) != ntk.trav_id() );
          fanin_nodes.emplace_back( ntk.get_node( fi ) );
        });

      try_adding_node( fanin_nodes, window_nodes );
      inputs.erase( std::remove( std::begin( inputs ), std::end( inputs ), *n ), std::end( inputs ) );

      ntk.foreach_fanin( *n, [&]( signal const& fi ){
          assert( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() );
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
      assert( ntk.visited( pivot ) != ntk.trav_id() && "pivot must not be part of the window" );

      ntk.set_visited( pivot, ntk.trav_id() );
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
          nodes.emplace_back( l );
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
};

} /* namespace detail */

namespace detail
{

template<class Ntk>
class multithreaded_cut_enumeration_impl
{
public:
  using node = node<Ntk>;
  using signal = signal<Ntk>;

public:
  explicit multithreaded_cut_enumeration_impl( Ntk const& ntk, multithreaded_cut_enumeration_params const& ps, multithreaded_cut_enumeration_stats& st )
    : ntk( ntk )
    , ps( ps )
    , st( st )
  {
  }

  void enumerate_windows_test()
  {
    window_manager windows( ntk, st );

    ntk.foreach_gate( [&]( node const& n ){
        auto const result = windows.create_window( n );
        if ( result )
        {
          ++st.total_num_windows;
          st.total_num_nodes += result->first.size();
          st.total_num_leaves += result->second.size();
        }
        ++st.total_num_candidates;
        return true;
      });
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
          if ( ps.very_verbose )
          {
            fmt::print( "[i] evaluate cut for node at index {} ({}/{})\n", index, num_processed_nodes, ntk.size() );
          }
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

private:
  Ntk const& ntk;
  multithreaded_cut_enumeration_params const& ps;
  multithreaded_cut_enumeration_stats& st;
}; /* multithreaded_cut_enumeration_impl */

} /* namespace detail */

template<class Ntk>
void multithreaded_cut_enumeration( Ntk const& ntk, multithreaded_cut_enumeration_params const& ps = {}, multithreaded_cut_enumeration_stats *pst = nullptr )
{
  depth_view<Ntk> ntk2{ntk};
  fanout_view<decltype( ntk2 )> ntk3{ntk2};

  multithreaded_cut_enumeration_stats st;
  detail::multithreaded_cut_enumeration_impl cut_enum( ntk3, ps, st );
  cut_enum.run();
  if ( ps.verbose )
  {
    st.report();
  }
  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
