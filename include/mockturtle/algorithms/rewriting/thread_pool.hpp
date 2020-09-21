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
  \file thread_pool.hpp
  \brief Thread pool.

  \author Heinz Riener
*/

#include <condition_variable>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

namespace mockturtle
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

} /* namespace mockturtle */
