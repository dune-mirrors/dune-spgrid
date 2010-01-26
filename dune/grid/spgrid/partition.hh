#ifndef DUNE_SPGRID_PARTITION_HH
#define DUNE_SPGRID_PARTITION_HH

#include <dune/common/iostream.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/mesh.hh>

namespace Dune
{

  // SPPartition
  // -----------

  template< int dim >
  class SPPartition
  {
    typedef SPPartition< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;
    typedef SPMesh< dimension > Mesh;

    SPPartition ( const Mesh &mesh )
    : begin_( mesh.begin() ),
      end_( mesh.end() )
    {}

  private:
    SPPartition ( const MultiIndex &begin, const MultiIndex &end )
    : begin_( begin ),
      end_( end )
    {}

  public:
    const MultiIndex &begin () const
    {
      return begin_;
    }

    const MultiIndex &end () const
    {
      return end_;
    }

    bool empty () const;

    This grow ( const int amount, const unsigned int dir = ((1 << dimension)-1) ) const;
    This intersect ( const This &other ) const;

    int volume () const;
    MultiIndex width () const;

  private:
    MultiIndex begin_, end_;
  };



  // Implementation of SPPartition
  // -----------------------------

  template< int dim >
  inline bool SPPartition< dim >::empty () const
  {
    bool empty = false;
    for( int i = 0; i < dimension; ++i )
      empty |= (end[ i ] < begin[ i ]);
    return empty;
  }


  template< int dim >
  inline SPPartition< dim >
  SPPartition< dim >::grow ( const int amount, const unsigned int dir ) const
  {
    MultiIndex b, e;
    for( int i = 0; i < dimension; ++i )
    {
      const int s = ((dir >> i) & 1) * amount;
      b[ i ] = begin()[ i ] - s;
      e[ i ] = end()[ i ] + s;
    }
    return This( b, e );
  }


  template< int dim >
  inline SPPartition< dim >
  SPPartition< dim >::intersect ( const This &other ) const
  {
    MultiIndex b, e;
    for( int i = 0; i < dimension; ++i )
    {
      b[ i ] = std::max( begin()[ i ], other.begin()[ i ] );
      e[ i ] = std::min( end()[ i ], other.end()[ i ] );
    }
    return This( b, e );
  }


  template< int dim >
  inline int SPPartition< dim >::volume () const
  {
    const MultiIndex &w = width();
    int volume = 1;
    for( int i = 0; i < dimension; ++i )
      volume *= w[ i ];
    return volume;
  }


  template< int dim >
  inline typename SPPartition< dim >::MultiIndex SPPartition< dim >::width () const
  {
    MultiIndex width;
    for( int i = 0; i < dimension; ++i )
      width[ i ] = std::max( end()[ i ] - begin()[ i ], 0 );
    return width;
  }



  // Auxilliary functions for SPPartition
  // ------------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPPartition< dim > &partition )
  {
    return out << "[ " << partition.begin() << ", " << partition.end() << " [";
  }

}

#endif // #ifndef DUNE_SPGRID_PARTITION_HH
