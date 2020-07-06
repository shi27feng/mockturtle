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

namespace detail
{

template<class Ntk>
class multithreaded_cut_enumeration_impl
{
public:
  explicit multithreaded_cut_enumeration_impl( Ntk const& ntk, multithreaded_cut_enumeration_params const& ps )
    : ntk( ntk )
    , ps( ps )
  {
  }

  void run()
  {
    ntk.incr_trav_id();
    auto const trav_id = ntk.trav_id();

    thread_pool threads{ps.num_threads};
    uint64_t const size = ntk.size();
    uint64_t num_processed_nodes = 0u;

    std::random_device rand;
    std::mt19937 gen( rand() ); /* seed random generator */
    std::uniform_int_distribution<> distribution( 0, size ); /* define the range */

    while ( num_processed_nodes < size )
    {
      auto const index = distribution( gen );

      auto const node = ntk.index_to_node( index );
      if ( ntk.visited( node ) == trav_id )
        continue;
      ntk.set_visited( node, trav_id );

      threads.enqueue( [=]{
          if ( ps.very_verbose )
          {
            std::cout << fmt::format( "[i] compute cut for node at index {} ({})\n", index, ntk.size() );
          }
        } );

      ++num_processed_nodes;
    }
  }

private:
  Ntk const& ntk;
  multithreaded_cut_enumeration_params const& ps;
}; /* multithreaded_cut_enumeration_impl */

} /* namespace detail */

template<class Ntk>
void multithreaded_cut_enumeration( Ntk const& ntk, multithreaded_cut_enumeration_params const& ps = {} )
{
  detail::multithreaded_cut_enumeration_impl cut_enum( ntk, ps );
  cut_enum.run();
}

} /* namespace mockturtle */
