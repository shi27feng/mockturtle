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
  \file window_view.hpp
  \brief Isolated view on a multi-input multi-output window of a network.

  \author Heinz Riener
*/

#include "../../views/immutable_view.hpp"
#include <vector>

namespace mockturtle
{

/*! \brief Implements an isolated view on a window in a network. */
template<typename Ntk>
class window_view2 : public immutable_view<Ntk>
{
public:
  using storage = typename Ntk::storage;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit window_view2( Ntk const& ntk, std::vector<node> const& leaves, std::vector<signal> const& roots, std::vector<node> const& nodes )
    : immutable_view<Ntk>( ntk )
    , leaves( leaves )
    , roots( roots )
    , nodes( nodes )
  {
    /* ensure that all leaves are also stored in nodes */
    assert( nodes.size() >= leaves.size() );

    /* add constant node at the beginning of nodes if missing */
    if ( this->nodes.begin() != this->nodes.end() && *this->nodes.begin() != 0u )
    {
      this->nodes.insert( std::begin( this->nodes ), 0u );
    }

    /* sort nodes topologically */
    /* ... TODO ... */
    
    compute_node_to_index();
    // compute_roots();
    compute_gates();
  }
  
  inline auto size() const
  {
    return nodes.size();
  }

  inline auto num_pis() const
  {
    return leaves.size();
  }

  inline auto num_pos() const
  {
    return roots.size();
  }

  inline auto num_cis() const
  {
    return leaves.size();
  }
  
  inline auto num_cos() const
  {
    return roots.size();
  }
  
  inline auto num_gates() const
  {
    assert( nodes.size() - leaves.size() - 1u == gates.size() );
    return gates.size();
  }

  inline auto node_to_index( node const& n ) const
  {
    assert( node_to_index_map.find( n ) != std::end( node_to_index_map ) );
    return node_to_index_map.at( n );
  }

  inline auto index_to_node( uint32_t index ) const
  {
    return nodes[index];
  }

  inline bool is_pi( node const& pi ) const
  {
    return std::find( std::begin( leaves ), std::end( leaves ), pi ) != std::end( leaves );
  }

  template<typename Fn>
  void foreach_pi( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( leaves ), std::end( leaves ), fn );
  }

  template<typename Fn>
  void foreach_po( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( roots ), std::end( roots ), fn );
  }

  template<typename Fn>
  void foreach_ci( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( leaves ), std::end( leaves ), fn );
  }

  template<typename Fn>
  void foreach_co( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( roots ), std::end( roots ), fn );
  }
  
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( nodes ), std::end( nodes ), fn );
  }

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    detail::foreach_element( std::begin( gates ), std::end( gates ), fn );
  }

protected:
  void compute_node_to_index()
  {
    uint32_t index = 0u;
    for ( const auto& n : nodes )
    {
      node_to_index_map[n] = index++;
    }    
  }

  void compute_gates()
  {
    for ( const auto& n : nodes )
    {
      if ( n == 0 || immutable_view<Ntk>::is_ci( n ) )
      {
        continue;
      }

      /* if the node is not a leave, then it's a gate */
      if ( std::find( std::begin( leaves ), std::end( leaves ), n ) == std::end( leaves ) )
      {
        gates.push_back( n );
      }
    }
  }

#if 0
  void compute_roots()
  {
    std::vector<int32_t> refs( immutable_view<Ntk>::size() );

    /* reference nodes */
    for ( const auto& n : nodes )
    {
      if ( n == 0 || immutable_view<Ntk>::is_ci( n ) )
      {
        continue;
      }

      immutable_view<Ntk>::foreach_fanin( n, [&]( signal const& fi ){
        ++refs[immutable_view<Ntk>::get_node( fi )];
      });
    }
    
    /* ensure that leaves are not marked as outputs */
    for ( const auto& l : leaves )
    {
      refs[l] = immutable_view<Ntk>::fanout_size( l );
    }

    /* collect roots and gates */
    for ( const auto& n : nodes )
    {
      if ( n == 0 || immutable_view<Ntk>::is_ci( n ) )
      {
        continue;
      }

      /* if node has more fanouts than references, we have found an output of the window */
      if ( int32_t( immutable_view<Ntk>::fanout_size( n ) ) != refs[n] )
      {
        // roots.emplace_back( immutable_view<Ntk>::make_signal( n ) );
        std::cout << "root: " << n << std::endl;
      }
    }
  }
#endif

protected:
  std::vector<node> leaves;
  std::vector<signal> roots;
  std::vector<node> nodes;
  std::vector<node> gates;
  spp::sparse_hash_map<node, uint32_t> node_to_index_map;
}; /* window_view */

} /* mockturtle */
