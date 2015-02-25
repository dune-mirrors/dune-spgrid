#ifndef DUNE_SPGRID_PARTITION_HH
#define DUNE_SPGRID_PARTITION_HH

#include <array>
#include <limits>

#include <dune/common/iostream.hh>

#include <dune/grid/spgrid/direction.hh>
#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/mesh.hh>
#include <dune/grid/spgrid/normal.hh>

namespace Dune
{

  // SPBasicPartition
  // ----------------

  template< int dim >
  class SPBasicPartition
  {
    typedef SPBasicPartition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    typedef SPDirection< dimension > Direction;

    SPBasicPartition ( const MultiIndex &begin, const MultiIndex &end ) : bound_{{ begin, end }} {}

    const MultiIndex &begin () const { return bound( 0 ); }
    const MultiIndex &end () const { return bound( 1 ); }

    const MultiIndex &bound ( unsigned int b ) const { assert( b == (b & 1) ); return bound_[ b ]; }

    int bound ( unsigned int b, int i, unsigned int d ) const
    {
      assert( d == (d & 1) );
      return bound( b )[ i ] - (2*b-1) * ((bound( b )[ i ] ^ d) & 1);
    }

    int bound ( const SPNormalId< dimension > &id ) const { return bound( id.face() & 1 )[ id.axis() ]; }

    This intersect ( const This &other ) const
    {
      return This( std::max( begin(), other.begin() ), std::min( end(), other.end() ) );
    }

    bool contains ( const MultiIndex &id ) const;

    bool empty () const;
    bool empty ( Direction dir ) const;

    int volume () const;
    MultiIndex width () const;
    int width ( int i ) const { return std::max( (end()[ i ]+1)/2 - begin()[ i ]/2, 0 ); }

    template< class char_type, class traits >
    void print ( std::basic_ostream< char_type, traits > &out ) const;

    template< class char_type, class traits >
    void print ( std::basic_ostream< char_type, traits > &out, const int i ) const;

  private:
    std::array< MultiIndex, 2 > bound_;
  };



  // SPPartition
  // -----------

  template< int dim >
  class SPPartition
    : public SPBasicPartition< dim >
  {
    typedef SPPartition< dim > This;
    typedef SPBasicPartition< dim > Base;

    typedef unsigned int Flags;

  public:
    static const int dimension = dim;

    typedef typename Base::MultiIndex MultiIndex;
    typedef SPMesh< dimension > Mesh;

    SPPartition ( const Base &base, const unsigned int number );
    SPPartition ( const MultiIndex &begin, const MultiIndex &end,
                  const unsigned int number );
    SPPartition ( const MultiIndex &begin, const MultiIndex &end,
                  const Mesh &globalMesh, const unsigned int number );

    unsigned int number () const;
    const unsigned int &neighbor ( const int face ) const;
    unsigned int &neighbor ( const int face );

    bool hasNeighbor ( const int face ) const;
    Flags boundary () const;
    bool boundary ( const int face ) const;

  private:
    unsigned int number_;
    unsigned int neighbor_[ 2*dimension ];
    Flags boundary_;
  };



  // Implementation of SPBasicPartition
  // ----------------------------------

  template< int dim >
  inline bool SPBasicPartition< dim >::contains ( const MultiIndex &id ) const
  {
    bool contains = true;
    for( int i = 0; i < dimension; ++i )
      contains &= (id[ i ] >= begin()[ i ]) && (id[ i ] <= end()[ i ]);
    return contains;
  }


  template< int dim >
  inline bool SPBasicPartition< dim >::empty () const
  {
    bool empty = false;
    for( int i = 0; i < dimension; ++i )
      empty |= (begin()[ i ] > end()[ i ]);
    return empty;
  }


  template< int dim >
  inline bool SPBasicPartition< dim >::empty ( Direction dir ) const
  {
    bool empty = false;
    for( int i = 0; i < dimension; ++i )
      empty |= (bound( 0, i, dir[ i ] ) > bound( 1, i, dir[ i ] ));
    return empty;
  }


  template< int dim >
  inline int SPBasicPartition< dim >::volume () const
  {
    int volume = 1;
    for( int i = 0; i < dimension; ++i )
      volume *= width( i );
    return volume;
  }


  template< int dim >
  inline typename SPBasicPartition< dim >::MultiIndex
  SPBasicPartition< dim >::width () const
  {
    MultiIndex w;
    for( int i = 0; i < dimension; ++i )
      w[ i ] = width( i );
    return w;
  }


  template< int dim >
  template< class char_type, class traits >
  inline void SPBasicPartition< dim >
    ::print ( std::basic_ostream< char_type, traits > &out ) const
  {
    print( out, 0 );
    for( int i = 1; i < dimension; ++i )
    {
      out << " x ";
      print( out, i );
    }
  }


  template< int dim >
  template< class char_type, class traits >
  inline void SPBasicPartition< dim >
    ::print ( std::basic_ostream< char_type, traits > &out, const int i ) const
  {
    const int b = begin()[ i ];
    const int e= end()[ i ];

    const char left = '[' + (b & 1)*(']'-'[');
    const char right = ']' + (e & 1)*('['-']');

    out << left << ' ' << (b/2) << ", " << ((e+1)/2) << ' ' << right;
  }



  // Implementation of SPPartition
  // -----------------------------

  template< int dim >
  SPPartition< dim >::SPPartition ( const Base &base, const unsigned int number )
  : Base( base ),
    number_( number ),
    boundary_( ((Flags( 1 ) << (2*dimension-1))-1) | (Flags( 1 ) << (2*dimension-1)) )
  {
    for( int face = 0; face < 2*dimension; ++face )
      neighbor_[ face ] = std::numeric_limits< unsigned int >::max();
  }


  template< int dim >
  SPPartition< dim >
   ::SPPartition ( const MultiIndex &begin, const MultiIndex &end,
                   const unsigned int number )
  : Base( begin, end ),
    number_( number ),
    boundary_( ((Flags( 1 ) << (2*dimension-1))-1) | (Flags( 1 ) << (2*dimension-1)) )
  {
    for( int face = 0; face < 2*dimension; ++face )
      neighbor_[ face ] = std::numeric_limits< unsigned int >::max();
  }


  template< int dim >
  SPPartition< dim >
   ::SPPartition ( const MultiIndex &begin, const MultiIndex &end,
                   const Mesh &globalMesh, const unsigned int number )
  : Base( begin, end ),
    number_( number ),
    boundary_( 0 )
  {
    for( int i = 0; i < dimension; ++i )
    {
      neighbor_[ 2*i ] = neighbor_[ 2*i+1 ] = std::numeric_limits< unsigned int >::max();
      boundary_ |= Flags( begin[ i ] == 2*globalMesh.begin()[ i ] ) << (2*i);
      boundary_ |= Flags( end[ i ] == 2*globalMesh.end()[ i ] ) << (2*i+1);
    }
  }


  template< int dim >
  inline unsigned int SPPartition< dim >::number () const
  {
    return number_;
  }


  template< int dim >
  inline const unsigned int &
  SPPartition< dim >::neighbor ( const int face ) const
  {
    assert( (face >= 0) && (face < 2*dimension) );
    return neighbor_[ face ];
  }


  template< int dim >
  inline unsigned int &
  SPPartition< dim >::neighbor ( const int face )
  {
    assert( (face >= 0) && (face < 2*dimension) );
    return neighbor_[ face ];
  }


  template< int dim >
  inline bool SPPartition< dim >::hasNeighbor ( const int face ) const
  {
    return (neighbor( face ) < std::numeric_limits< unsigned int >::max());
  }


  template< int dim >
  inline typename SPPartition< dim >::Flags
  SPPartition< dim >::boundary () const
  {
    return boundary_;
  }


  template< int dim >
  inline bool SPPartition< dim >::boundary ( const int face ) const
  {
    assert( (face >= 0) && (face < 2*dimension) );
    return bool( (boundary_ >> face) & 1 );
  }



  // Auxilliary functions for SPBasicPartition
  // -----------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out,
               const SPBasicPartition< dim > &partition )
  {
    partition.print( out );
    return out;
  }

} // namespace Dune

#endif // #ifndef DUNE_SPGRID_PARTITION_HH
