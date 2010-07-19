#ifndef DUNE_SPGRID_MESH_HH
#define DUNE_SPGRID_MESH_HH

#include <dune/common/iostream.hh>

#include <dune/grid/spgrid/multiindex.hh>
#include <dune/grid/spgrid/refinement.hh>

namespace Dune
{

  // SPMesh
  // ------

  template< int dim >
  class SPMesh
  {
    typedef SPMesh< dim > This;

  public:
    static const int dimension = dim;

    typedef SPMultiIndex< dimension > MultiIndex;

    explicit SPMesh ( const MultiIndex &width );

    SPMesh ( const MultiIndex &begin, const MultiIndex &end );

    const This &operator+= ( const MultiIndex &shift );
    const This &operator-= ( const MultiIndex &shift );

    const MultiIndex &begin () const;
    const MultiIndex &end () const;

    const MultiIndex &bound ( const int b ) const;

    bool empty () const;

    template< SPRefinementStrategy strategy >
    This refine ( const SPRefinement< dimension, strategy > &refinement ) const;

    This grow ( int size ) const;
    This grow ( const MultiIndex &size ) const;

    This intersect ( const This &other ) const;

    std::pair< This, This > split ( const int dir, const int leftWeight, const int rightWeight ) const;

    int volume () const;

    MultiIndex width () const;
    int width ( const int i ) const;

    static This unitMesh ();

  private:
    MultiIndex bound_[ 2 ];
  };



  // Implementation of SPMesh
  // ------------------------

  template< int dim >
  SPMesh< dim >::SPMesh ( const MultiIndex &width )
  {
    bound_[ 0 ] = MultiIndex::zero();
    bound_[ 1 ] = width;
  }


  template< int dim >
  SPMesh< dim >::SPMesh ( const MultiIndex &begin, const MultiIndex &end )
  {
    bound_[ 0 ] = begin;
    bound_[ 1 ] = end;
  }


  template< int dim >
  inline const typename SPMesh< dim >::This &
  SPMesh< dim >::operator+= ( const MultiIndex &shift )
  {
    for( int b = 0; b < 2; ++b )
      bound_[ b ] += shift;
    return *this;
  }


  template< int dim >
  inline const typename SPMesh< dim >::This &
  SPMesh< dim >::operator-= ( const MultiIndex &shift )
  {
    for( int b = 0; b < 2; ++b )
      bound_[ b ] -= shift;
    return *this;
  }


  template< int dim >
  inline const typename SPMesh< dim >::MultiIndex &
  SPMesh< dim >::begin () const
  {
    return bound( 0 );
  }


  template< int dim >
  inline const typename SPMesh< dim >::MultiIndex &
  SPMesh< dim >::end () const
  {
    return bound( 1 );
  }


  template< int dim >
  inline const typename SPMesh< dim >::MultiIndex &
  SPMesh< dim >::bound ( const int b ) const
  {
    assert( (b == 0) || (b == 1) );
    return bound_[ b ];
  }


  template< int dim >
  inline bool SPMesh< dim >::empty () const
  {
    bool empty = false;
    for( int i = 0; i < dimension; ++i )
      empty |= (end()[ i ] < begin()[ i ]);
    return empty;
  }


  template< int dim >
  template< SPRefinementStrategy strategy >
  inline typename SPMesh< dim >::This
  SPMesh< dim >::refine ( const SPRefinement< dimension, strategy > &refinement ) const
  {
    MultiIndex childBegin, childEnd;
    for( int i = 0; i < dimension; ++i )
    {
      const int factor = refinement.factor( i );
      childBegin[ i ] = factor * begin()[ i ];
      childEnd[ i ] = factor * end()[ i ];
    }
    return This( childBegin, childEnd );
  }


  template< int dim >
  inline typename SPMesh< dim >::This
  SPMesh< dim >::grow ( const int size ) const
  {
    MultiIndex begin, end;
    for( int i = 0; i < dim; ++i )
    {
      begin[ i ] = begin()[ i ] - size;
      end[ i ] = end()[ i ] + size;
    }
    return This( begin, end );
  }


  template< int dim >
  inline typename SPMesh< dim >::This
  SPMesh< dim >::grow ( const MultiIndex &size ) const
  {
    return This( begin() - size, end() + size );
  }


  template< int dim >
  inline typename SPMesh< dim >::This
  SPMesh< dim >::intersect ( const This &other ) const
  {
    return This( std::max( begin(), other.begin() ), std::min( end(), other.end() ) );
  }


  template< int dim >
  inline std::pair< typename SPMesh< dim >::This, typename SPMesh< dim >::This >
  SPMesh< dim >::split ( const int dir, const int leftWeight, const int rightWeight ) const
  {
    const MultiIndex &lbegin = begin();
    const MultiIndex &rend = end();

    assert( (dir >= 0) && (dir < dimension) );
    const int width = (rend[ dir ] - lbegin[ dir ]);
    const int leftWidth = (leftWeight * width) / (leftWeight + rightWeight);

    MultiIndex lend = rend;
    MultiIndex rbegin = lbegin;
    rbegin[ dir ] = lend[ dir ] = lbegin[ dir ] + leftWidth;

    return std::make_pair( This( lbegin, lend ), This( rbegin, rend ) );
  }


  template< int dim >
  inline int SPMesh< dim >::volume () const
  {
    const MultiIndex &w = width();
    int volume = 1;
    for( int i = 0; i < dimension; ++i )
      volume *= w[ i ];
    return volume;
  }


  template< int dim >
  inline typename SPMesh< dim >::MultiIndex SPMesh< dim >::width () const
  {
    MultiIndex w;
    for( int i = 0; i < dimension; ++i )
      w[ i ] = width( i );
    return w;
  }


  template< int dim >
  inline int SPMesh< dim >::width ( const int i ) const
  {
    //return std::max( end()[ i ] - begin()[ i ], 0 );
    return end()[ i ] - begin()[ i ];
  }


  template< int dim >
  inline typename SPMesh< dim >::This
  SPMesh< dim >::unitMesh ()
  {
    MultiIndex w;
    for( int i = 0; i < dimension; ++i )
      w[ i ] = 1;
    return This( w );
  }



  // Auxilliary functions for SPMesh
  // -------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPMesh< dim > &mesh )
  {
    return out << "[ " << mesh.begin() << ", " << mesh.end() << " [";
  }


  template< int dim >
  inline SPMesh< dim >
  operator+ ( const SPMesh< dim > &mesh, const SPMultiIndex< dim > &shift )
  {
    SPMesh< dim > copy( mesh );
    return copy += shift;
  }


  template< int dim >
  inline SPMesh< dim >
  operator- ( const SPMesh< dim > &mesh, const SPMultiIndex< dim > &shift )
  {
    SPMesh< dim > copy( mesh );
    return copy -= shift;
  }

}

#endif // #ifndef DUNE_SPGRID_MESH_HH
