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
  \file debug.hpp
  \brief Simple check for examining a logic networks consistency

  \author Heinz Riener
*/

#include <vector>
#include <queue>

/**
 * \brief Returns true iff the network's size and the network's number of nodes are equal.
 */
template<typename Ntk>
bool check_size( Ntk const& ntk )
{
  using node = typename Ntk::node;
  
  uint64_t node_counter = 0u;
  ntk.foreach_node( [&]( node const& ){
    ++node_counter;
  });

  return ( ntk.size() == node_counter );
}

/**
 * \brief Returns the number of dangling roots that are not POs.
 */
template<typename Ntk>
uint64_t count_dangling_roots( Ntk const& ntk )
{
  using node   = typename Ntk::node;
  using signal = typename Ntk::signal;

  std::vector<node> outs;
  ntk.foreach_po( [&]( signal const& o ){
    outs.push_back( ntk.get_node( o ) );
  });
  
  std::vector<int32_t> refs( ntk.size() );
  ntk.foreach_gate( [&]( node const& n ){
    ntk.foreach_fanin( n, [&]( signal const& fi ){
      ++refs[ntk.get_node( fi )];
    });
  });

  uint64_t num_dangling_roots{0};
  ntk.foreach_gate( [&]( node const& n ){
    if ( ( refs[n] != int64_t( ntk.fanout_size( n ) ) ) &&
         ( std::find( std::begin( outs ), std::end( outs ), n ) == std::end( outs ) ) )
    {
      ++num_dangling_roots;
    }
  });

  return num_dangling_roots;
}

/**
 * \brief Returns true iff the network has a cycle.
 */
template<typename Ntk>
bool has_dangling_roots( Ntk const& ntk )
{
  return ( count_dangling_roots( ntk ) != 0 );
}

#if 0
template<typename Ntk>
void test2( Ntk const& ntk )
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  std::vector<std::vector<node>> fanouts( ntk.size() );
  ntk.foreach_gate( [&]( node const& n ){
    ntk.foreach_fanin( n, [&]( signal const& fi ){
      auto& fos = fanouts[ntk.get_node( fi )];
      if ( std::find( std::begin( fos ), std::end( fos ), n ) == std::end( fos ) )
      {
        fos.push_back( n );
      }
    });
  });

  std::vector<node> outs;
  ntk.foreach_co( [&]( signal const& o ){
    outs.push_back( ntk.get_node( o ) );
  });

  uint64_t errors = 0;
  for ( uint64_t index = 0; index < ntk.size(); ++index )
  {
    if ( fanouts[index].size() != ntk.fanout_size( index ) )
    {
      ++errors;
    }
  }

  std::vector<node> topsort; // L
  std::queue<node> fronteer; // S
  std::vector<int32_t> refs( ntk.size() );

  fronteer.push( 0 );
  ntk.foreach_ci( [&]( node const& pi ){
    fronteer.push( pi );
    refs[pi] = 2;
  });

  while ( !fronteer.empty() )
  {    
    node const& n = fronteer.front();
    fronteer.pop();

    topsort.push_back( n );

    for ( auto const& m : fanouts[n] )
    {
      ++refs[m];
      if ( refs[m] == 2 )
      {
        fronteer.push( m );
      }
    }
  }

  std::cout << topsort.size() << " " << ntk.size() << std::endl;
}
#endif

/**
 * \brief Returns true iff the network DAG has a cycle.
 */
template<typename Ntk>
bool has_cycle( Ntk const& ntk )
{
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

  std::vector<std::vector<node>> fanouts( ntk.size() );
  ntk.foreach_gate( [&]( node const& n ){
    ntk.foreach_fanin( n, [&]( signal const& fi ){
      auto& fos = fanouts[ntk.get_node( fi )];
      if ( std::find( std::begin( fos ), std::end( fos ), n ) == std::end( fos ) )
      {
        fos.push_back( n );
      }
    });
  });

  /* Kahn's algorithm for topological sorting */
  std::vector<int32_t> refs( ntk.size() );

  for ( uint32_t i = 0; i < ntk.size(); ++i )
  {
    for ( const auto& fo : fanouts[i] )
    {
      ++refs[fo];
    }
  }

  std::queue<node> Q;
  for ( uint32_t i = 0; i < ntk.size(); ++i )
  {
    if ( refs[i] == 0 )
    {
      Q.push( i );
    }
  }

  uint64_t counter{0};
  std::vector<node> topsort;
  while ( !Q.empty() )
  {
    node const u = Q.front();
    Q.pop();

    topsort.emplace_back( u );

    for ( const auto& fo : fanouts[u] )
    {
      if ( --refs[fo] == 0 )
      {
        Q.push( fo );
      }
    }

    ++counter;
  }

  return ( counter != ntk.size() );
}
