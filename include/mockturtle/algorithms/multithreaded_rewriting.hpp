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
#include <mutex>
#include <thread>
#include <vector>

namespace mockturtle
{

namespace detail
{

class thread_pool
{
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

private:
  void start()
  {
    for ( auto i = 0ul; i < num_threads; ++i )
    {
      threads.emplace_back( [=]{
          while ( true )
          {
            std::unique_lock<std::mutex> lock{stop_mutex};
            stop_signal.wait( lock, [=]{ return stopping; } );

            if ( stopping )
            {
              break;
            }
          }
        } );
    }
  }

  void stop() noexcept
  {
    {
      std::unique_lock<std::mutex> lock{stop_mutex};
      stopping = true;
    }
    stop_signal.notify_all();

    for ( auto &thread : threads )
    {
      thread.join();
    }
  }

private:
  uint64_t num_threads;
  std::vector<std::thread> threads;

  std::condition_variable stop_signal;
  std::mutex stop_mutex;
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
    thread_pool threads{ps.num_threads};
    ntk.foreach_node( [this]( auto const node ) {
      const auto index = ntk.node_to_index( node );

      if ( ps.very_verbose )
      {
        std::cout << fmt::format( "[i] compute cut for node at index {}\n", index );
      }
      } );
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
