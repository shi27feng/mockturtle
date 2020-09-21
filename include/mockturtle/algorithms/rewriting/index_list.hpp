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
  \file index_list.hpp
  \brief List of indices to represent small networks.

  \author Heinz Riener
*/

#include <fmt/format.h>
#include <vector>

namespace mockturtle
{

/*
 * A wrapper class to represent small networks as a list of indices.
 *
 * The implementation supports AND and XOR gates and with the
 * datastructures used in ABC.
 */
class index_list
{
public:
  using element_type = uint32_t;

public:
  explicit index_list() {}

#if 0
  explicit index_list( int* lits, uint64_t num_lits )
    : literals( lits, lits + num_lits )
  {
    assert( num_lits % 2 == 0 );

    /* count the number of zero entries at the beginning of the list */
    for ( uint64_t i = 0; i < literals.size(); i+=2 )
    {
      if ( literals.at( i ) != 0 || literals.at( i + 1 ) != 0 )
        return;

      ++num_leading_zero_entries;
    }
  }
#endif

  void push2( element_type lit )
  {
    if ( lit == 0 && num_leading_zero_entries == num_entries() )
    {
      ++num_leading_zero_entries;
    }
    else
    {
      ++num_outputs;
    }
    literals.push_back( lit );
    literals.push_back( lit );
  }

  void push2( element_type lit0, element_type lit1 )
  {
    if ( lit0 == 0 && lit1 == 0 && num_leading_zero_entries == num_entries() )
    {
      ++num_leading_zero_entries;
    }
    else if ( lit0 == lit1 )
    {
      ++num_outputs;
    }
    literals.push_back( lit0 );
    literals.push_back( lit1 );
  }

  uint64_t num_entries() const
  {
    return ( literals.size() >> 1 );
  }

  void print_raw() const
  {
    assert( literals.size() % 2 == 0 );
    auto const raw_array = raw_data();
    std::cout << raw_array.second << std::endl;
    for ( uint64_t i = 0; i < raw_array.second*2; ++i )
    {
      std::cout << raw_array.first[i] << ' ';
    }
    std::cout << std::endl;
  }

  void print() const
  {
    assert( literals.size() % 2 == 0 );
    for ( auto i = 0u; i < literals.size(); i += 2 )
    {
      uint64_t const id  = i / 2;
      uint64_t const lit0 = literals.at( i );
      uint64_t const lit1 = literals.at( i+1 );

      /* constant or pi */
      if ( lit0 == 0 && lit1 == 0  )
      {
        fmt::print( "Obj {:4d} : {}\n", id, id == 0 ? "constant" : "PI" );
      }
      else if ( lit0 == lit1 )
      {
        fmt::print( "Obj {:4d} : PO( {:4d} )\n", id, lit0 );
      }
      /* AND gate */
      else if ( lit0 < lit1 )
      {
        fmt::print( "Obj {:4d} : AND( {:4d}, {:4d} )\n", id, lit0, lit1 );
      }
      else
      {
        fmt::print( "Obj {:4d} : XOR( {:4d}, {:4d} )\n", id, lit0, lit1 );
      }
    }
  }

  template<typename Fn>
  void foreach_entry( Fn&& fn ) const
  {
    assert( literals.size() % 2 == 0 );
    for ( auto i = 0u; i < literals.size(); i += 2 )
    {
      fn( literals.at( i ), literals.at( i+1 ) );
    }
  }

  uint64_t num_pis() const
  {
    /* the first entry is reserved for the constant */
    return ( num_leading_zero_entries - 1 );
  }

  uint64_t num_pos() const
  {
    return num_outputs;
  }
  
private:
  std::pair<int*, uint64_t> raw_data() const
  {
    return std::make_pair<int*, uint64_t>( (int*)&literals[0], literals.size() / 2 );
  }

private:
  std::vector<element_type> literals;
  uint64_t num_leading_zero_entries{0};
  uint64_t num_outputs{0};
}; /* index_list */

/* \brief Generates an index_list from a network type */
template<typename Ntk>
void encode( index_list& indices, Ntk const& ntk )
{
  using node   = typename Ntk::node;
  using signal = typename Ntk::signal;

  /* constant */
  // if ( ntk.has_constant() )
  // {
  indices.push2( 0 );
  // }

  /* inputs */
  for ( uint64_t i = 0; i < ntk.num_pis(); ++i )
  {
    indices.push2( 0 );
  }

  /* gates */
  ntk.foreach_gate( [&]( node const& n ){
      std::vector<signal> fanins;
      ntk.foreach_fanin( n, [&]( signal const& fi ){
          fanins.emplace_back( fi );
        });

      uint64_t lit0 = 2*ntk.node_to_index( ntk.get_node( fanins[0] ) ) + ntk.is_complemented( fanins[0] );
      uint64_t lit1 = 2*ntk.node_to_index( ntk.get_node( fanins[1] ) ) + ntk.is_complemented( fanins[1] );
      if ( ( ntk.is_and( n ) && ( lit0 > lit1 ) ) ||
           ( ntk.is_xor( n ) && ( lit0 < lit1 ) ) )
      {
        std::swap( lit0, lit1 );
      }
      indices.push2( lit0, lit1 );
    });

  /* outputs */
  ntk.foreach_po( [&]( signal const& f ){
      indices.push2( 2*ntk.node_to_index( ntk.get_node( f ) ) + ntk.is_complemented( f ) );
    });

  assert( indices.num_entries() == /* ntk.has_constant() */1u + ntk.num_pis() + ntk.num_gates() + ntk.num_pos() );
}

/* \brief Generates a network from an index_list */
template<typename Ntk>
void decode( Ntk& ntk, index_list const& indices )
{
  using signal = typename Ntk::signal;

  std::vector<signal> signals;
  for ( uint64_t i = 0; i < indices.num_pis(); ++i )
  {
    signals.push_back( ntk.create_pi() );
  }

  insert( ntk, std::begin( signals ), std::end( signals ), indices,
          [&]( signal const& s ){ ntk.create_po( s ); });
}

/* \brief Inserts an index_list into an existing network */
template<typename Ntk, typename BeginIter, typename EndIter, typename Fn>
void insert( Ntk& ntk, BeginIter begin, EndIter end, index_list const& indices, Fn&& fn )
{
  using signal = typename Ntk::signal;

  std::vector<signal> signals;
  signals.emplace_back( ntk.get_constant( false ) );
  for ( auto it = begin; it != end; ++it )
  {
    signals.push_back( *it );
  }

  bool skip_leading_zeros = true;
  indices.foreach_entry( [&]( uint32_t lit0, uint32_t lit1 ){
    /* skip constant and inputs at the beginning, but not if an output has been replaced with a constant */
    if ( skip_leading_zeros && ( lit0 == 0 && lit1 == 0 ) )
    {
      return;
    }

    uint64_t const i0 = lit0 >> 1;
    uint64_t const i1 = lit1 >> 1;
    bool const c0 = lit0 % 2;
    bool const c1 = lit1 % 2;
    
    signal const s0 = c0 ? !signals.at( i0 ) : signals.at( i0 );
    signal const s1 = c1 ? !signals.at( i1 ) : signals.at( i1 );

    skip_leading_zeros = false;
    
    /* outputs */
    if ( lit0 == lit1 )
    {
      fn( s0 );
    }
    /* AND gates */
    else if ( i0 < i1 )
    {
      signals.push_back( ntk.create_and( s0, s1 ) );
    }
    /* XOR gates */
    else
    {
      signals.push_back( ntk.create_xor( s0, s1 ) );
    }
  });
}

} /* mockturtle */
