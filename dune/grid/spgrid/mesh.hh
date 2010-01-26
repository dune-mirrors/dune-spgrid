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

    SPMesh ( const MultiIndex &width )
    : begin_( MultiIndex::zero() ),
      end_( width )
    {}

  private:
    SPMesh ( const MultiIndex &begin, const MultiIndex &end )
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

    template< SPRefinementStrategy strategy >
    This refine ( const SPRefinement< dimension, strategy > &refinement ) const;

    std::pair< This, This > split ( const int dir, const int leftWeight, const int rightWeight ) const;

    int volume () const;
    MultiIndex width () const;

  private:
    MultiIndex begin_, end_;
  };



  // Implementation of SPMesh
  // ------------------------

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
      w[ i ] = std::max( end()[ i ] - begin()[ i ], 0 );
    return w;
  }



  // Auxilliary functions for SPMesh
  // -------------------------------

  template< class char_type, class traits, int dim >
  inline std::basic_ostream< char_type, traits > &
  operator<< ( std::basic_ostream< char_type, traits > &out, const SPMesh< dim > &mesh )
  {
    return out << "[ " << mesh.begin() << ", " << mesh.end() << " [";
  }

}

#endif // #ifndef DUNE_SPGRID_MESH_HH
