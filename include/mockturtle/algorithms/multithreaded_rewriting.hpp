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
  /*! \brief Maximum number of leaves for a cut. */
  uint32_t cut_size{6u};

  /*! \brief Maximum number of cuts for a node. */
  uint32_t cut_limit{16u};

  /*! \brief Number of threads for cut enumeration. */
  uint32_t num_threads{32u};

  /*! \brief Be verbose. */
  bool verbose{true};

  /*! \brief Be very verbose. */
  bool very_verbose{true};
};

/*! \brief Statistics for multithreaded_cut_enumeration.
 *
 */
struct multithreaded_cut_enumeration_stats
{
  /*! \brief Total time. */
  stopwatch<>::duration time_total{0};

  /*! \brief Prints report. */
  void report() const
  {
    std::cout << fmt::format( "[i] total time       = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{


template<class Ntk>
class window_manager
{
public:
  static constexpr uint64_t const max_leaves = 6u;
  static constexpr uint64_t const max_fanouts = 1000u;

public:
  using node = node<Ntk>;
  using signal = signal<Ntk>;

public:
  explicit window_manager( Ntk const& ntk )
    : ntk( ntk )
    , levels( ntk.depth() + 1u )
  {
  }

  void reconvergence_driven_cut( node const& pivot )
  {
    /* clean up */
    leaves.clear();
    nodes.clear();

    nodes.emplace_back( pivot );
    ntk.set_visited( pivot, ntk.trav_id() );

    ntk.foreach_fanin( pivot, [&]( signal const& fi ){
        node const& n = ntk.get_node( fi );
        if ( n == 0 )
          return;

        leaves.emplace_back( n );
        nodes.emplace_back( n );
        ntk.set_visited( n, ntk.trav_id() );
        return;
      } );

    if ( leaves.size() > max_leaves )
    {
      /* special case: cut already overflows at the current node bc the cut size limit is very low */
      leaves.clear();
    }
    else
    {
      /* compute the cut */
      while ( construct_cut() );
      assert( leaves.size() <= max_leaves );
    }
  }

  bool construct_cut()
  {
    uint64_t best_cost{std::numeric_limits<uint64_t>::max()};
    std::optional<node> best_fanin;
    uint64_t best_position;

    uint64_t position = 0;
    for ( const auto& l : leaves )
    {
      uint64_t const current_cost = cost( l );
      if ( best_cost > current_cost ||
           ( best_cost == current_cost && best_fanin && ntk.level( l ) > ntk.level( *best_fanin ) ) )
      {
        best_cost = current_cost;
        best_fanin = std::make_optional( l );
        best_position = position;
      }

      if ( best_cost == 0u )
        break;

      ++position;
    }

    if ( !best_fanin )
    {
      return false;
    }

    if ( leaves.size() - 1 + best_cost > max_leaves )
    {
      return false;
    }

    /* remove the best node from the array */
    leaves.erase( std::begin( leaves ) + best_position );

    /* add the fanins of best to leaves and nodes */
    ntk.foreach_fanin( *best_fanin, [&]( signal const& fi ){
        node const& n = ntk.get_node( fi );
        if ( n != 0 && ( ntk.visited( n ) != ntk.trav_id() ) )
        {
          ntk.set_visited( n, ntk.trav_id() );
          nodes.emplace_back( n );
          leaves.emplace_back( n );
        }
      });

    assert( leaves.size() <=  max_leaves );
    return true;
  }

  uint64_t cost( node const &node )
  {
    /* make sure the node is in the construction zone */
    assert( ntk.visited( node ) == ntk.trav_id() );

    /* cannot expand over the CI node */
    if ( ntk.is_constant( node ) || ntk.is_ci( node ) )
    {
      return std::numeric_limits<uint64_t>::max();
    }

    /* get the cost of the cone */
    uint64_t cost = 0u;
    ntk.foreach_fanin( node, [&]( const auto& fi ){
        cost += ( ntk.visited( ntk.get_node( fi ) ) == ntk.trav_id() ) ? 0 : 1;
      });

    /* always accept if the number of leaves does not increase */
    if ( cost < ntk.fanin_size( node ) )
    {
      return std::numeric_limits<uint64_t>::max();
    }

    /* skip nodes with many fanouts */
    if ( ntk.fanout_size( node ) > max_fanouts )
    {
      return std::numeric_limits<uint64_t>::max();
    }

    /* return the number of nodes that will be on the leaves if this node is removed */
    return cost;
  }

  std::vector<node> create_window( node const& pivot )
  {
    ntk.incr_trav_id();

    reconvergence_driven_cut( pivot );

    std::cout << "[i] leaves = ";
    for ( const auto& l : leaves )
    {
      std::cout << l << ' ';
    }
    std::cout << std::endl;

    std::cout << "[i] nodes = ";
    for ( const auto& n : nodes )
    {
      std::cout << n << ' ';
    }
    std::cout << std::endl;


    return nodes;
  }


private:
  Ntk const& ntk;

  /* levelized structure to temporarily store nodes */
  std::vector<std::vector<node>> levels;

  /* number of nodes stored in levels */
  uint64_t num_nodes{0u};

  std::vector<node> leaves;
  std::vector<node> nodes;
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

  void run()
  {
    stopwatch t( st.time_total );

    thread_pool threads{ps.num_threads};
    uint64_t const size = ntk.size();
    uint64_t num_processed_nodes = 0u;

    std::random_device rand;
    std::mt19937 gen( rand() ); /* seed random generator */
    std::uniform_int_distribution<> distribution( 0, size - 1u ); /* define the range */

    window_manager windows( ntk );

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
      windows.create_window( pivot );

      threads.enqueue( [=]{
          /* evaluate all cuts of the current node concurrently */
          if ( ps.very_verbose )
          {
            fmt::print( "[i] evaluate cut for node at index {} ({}/{})\n", index, num_processed_nodes, ntk.size() );
          }
        } );

      ++num_processed_nodes;
    }
  }

  void create_window( node const& pivot )
  {
    (void)pivot;


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
